//
//  MM4RigidBody+Velocity.swift
//
//
//  Created by Philip Turner on 11/20/23.
//

import Numerics

extension MM4RigidBody {
  /// The velocity (in nanometers per picosecond) of each atom.
  public var velocities: [SIMD3<Float>] {
    // _modify not supported b/c it requires very complex caching logic.
    // Workaround: use the exposed setVelocities function.
    _read { fatalError("Not implemented.") }
  }
  
  @_specialize(where T == Double)
  @_specialize(where T == Float)
  public mutating func setVelocities<T: BinaryFloatingPoint>(
    _ buffer: UnsafeBufferPointer<SIMD3<T>>
  ) {
    ensureUniquelyReferenced()
    storage.anchorVelocitiesValid = nil
    storage.velocities = nil
    
    guard buffer.count == storage.atoms.count else {
      fatalError("Velocity buffer was not the correct size.")
    }
    let baseAddress = buffer.baseAddress.unsafelyUnwrapped
    
    for vID in 0..<storage.atoms.vectorCount {
      let (x, y, z) = storage.swizzleToVectorWidth(vID, baseAddress)
      storage.vVelocities[vID &* 3 &+ 0] = x
      storage.vVelocities[vID &* 3 &+ 1] = y
      storage.vVelocities[vID &* 3 &+ 2] = z
    }
  }
  
  @_specialize(where T == Double)
  @_specialize(where T == Float)
  public func getVelocities<T: BinaryFloatingPoint>(
    _ buffer: UnsafeMutableBufferPointer<SIMD3<T>>
  ) {
    guard buffer.count == storage.atoms.count else {
      fatalError("Velocity buffer was not the correct size.")
    }
    let baseAddress = buffer.baseAddress.unsafelyUnwrapped
    
    for vID in 0..<storage.atoms.vectorCount {
      let x = storage.vVelocities[vID &* 3 &+ 0]
      let y = storage.vVelocities[vID &* 3 &+ 1]
      let z = storage.vVelocities[vID &* 3 &+ 2]
      storage.swizzleFromVectorWidth((x, y, z), vID, baseAddress)
    }
  }
}

extension MM4RigidBodyStorage {
  func ensureAnchorVelocitiesValid() {
    guard anchorVelocitiesValid == nil else {
      return
    }
    
    let epsilon: Float = 1e-4
    var firstVelocity: SIMD3<Float>?
    var valid = true
    
    for anchor in anchors {
      let velocity = extractScalar(Int(anchor), vVelocities)
      if let firstVelocity {
        var deviation = velocity - firstVelocity
        var accepted = epsilon * firstVelocity
        deviation.replace(with: -deviation, where: deviation .< 0)
        accepted.replace(with: -accepted, where: accepted .< 0)
        if any(deviation .> accepted) {
          valid = false
        }
      } else {
        firstVelocity = velocity
      }
    }
    guard valid else {
      fatalError("Anchor velocities are invalid.")
    }
    anchorVelocitiesValid = valid
  }
  
  func createLinearVelocity() -> SIMD3<Float> {
    ensureAnchorVelocitiesValid()
    guard anchors.count == 0 else {
      let anchor = anchors.first!
      let velocity = extractScalar(Int(anchor), vVelocities)
      return velocity
    }
    guard atoms.count > 0 else {
      return .zero
    }
    
    var momentum: SIMD3<Double> = .zero
    withMasses(nonAnchorMasses) { vMasses in
      withSegmentedLoop(chunk: 256) {
        var vMomentumX: MM4FloatVector = .zero
        var vMomentumY: MM4FloatVector = .zero
        var vMomentumZ: MM4FloatVector = .zero
        for vID in $0 {
          let x = vVelocities[vID &* 3 &+ 0]
          let y = vVelocities[vID &* 3 &+ 1]
          let z = vVelocities[vID &* 3 &+ 2]
          let mass = vMasses[vID]
          vMomentumX.addProduct(mass, x)
          vMomentumY.addProduct(mass, y)
          vMomentumZ.addProduct(mass, z)
        }
        momentum.x += MM4DoubleVector(vMomentumX).sum()
        momentum.y += MM4DoubleVector(vMomentumY).sum()
        momentum.z += MM4DoubleVector(vMomentumZ).sum()
      }
    }
    return SIMD3<Float>(momentum / nonAnchorMass)
  }
  
  func createAngularVelocity() -> Quaternion<Float> {
    ensureAnchorVelocitiesValid()
    guard anchors.count <= 1, atoms.count > 0 else {
      return .zero
    }
    
    var momentum: SIMD3<Double> = .zero
    ensureCenterOfMassCached()
    guard let centerOfMass else {
      fatalError("This should never happen.")
    }
    withMasses(nonAnchorMasses) { vMasses in
      withSegmentedLoop(chunk: 256) {
        var vMomentumX: MM4FloatVector = .zero
        var vMomentumY: MM4FloatVector = .zero
        var vMomentumZ: MM4FloatVector = .zero
        for vID in $0 {
          let rX = vPositions[vID &* 3 &+ 0] - centerOfMass.x
          let rY = vPositions[vID &* 3 &+ 1] - centerOfMass.y
          let rZ = vPositions[vID &* 3 &+ 2] - centerOfMass.z
          let vX = vVelocities[vID &* 3 &+ 0]
          let vY = vVelocities[vID &* 3 &+ 1]
          let vZ = vVelocities[vID &* 3 &+ 2]
          let mass = vMasses[vID]
          vMomentumX.addProduct(mass, rY * vZ - rZ * vY)
          vMomentumY.addProduct(mass, rZ * vX - rX * vZ)
          vMomentumZ.addProduct(mass, rX * vY - rY * vX)
        }
        momentum.x += MM4DoubleVector(vMomentumX).sum()
        momentum.y += MM4DoubleVector(vMomentumY).sum()
        momentum.z += MM4DoubleVector(vMomentumZ).sum()
      }
    }
    
    ensureMomentOfInertiaCached()
    guard let momentOfInertia else {
      fatalError("This should never happen.")
    }
    let inverse = momentOfInertia.inverse
    let velocityX = inverse.0 * momentum.x
    let velocityY = inverse.1 * momentum.y
    let velocityZ = inverse.2 * momentum.z
    let velocity = velocityX + velocityY + velocityZ
    
    if all(velocity .< .leastNormalMagnitude .&
           velocity .> -.leastNormalMagnitude) {
      return .zero
    } else {
      let lengthSquared = (velocity * velocity).sum()
      if abs(lengthSquared) < .leastNormalMagnitude {
        return .zero
      }
      
      let axis = SIMD3<Float>(velocity / lengthSquared.squareRoot())
      let angle = Float(lengthSquared.squareRoot())
      return Quaternion<Float>(angle: angle, axis: axis)
    }
  }
}

extension MM4RigidBodyStorage {
  func createTotalKineticEnergy() -> Double {
    guard atoms.count > 0 else {
      return 0
    }
    
    // Anchors should never be included in kinetic energy, either free or
    // thermal. For thermal, the kinetic energy from anchor velocities will
    // contribute zero to the total energy. This should let velocity rescaling
    // work properly.
    var kinetic: Double = .zero
    withMasses(nonAnchorMasses) { vMasses in
      withSegmentedLoop(chunk: 256) {
        var vKineticX: MM4FloatVector = .zero
        var vKineticY: MM4FloatVector = .zero
        var vKineticZ: MM4FloatVector = .zero
        for vID in $0 {
          let x = vVelocities[vID &* 3 &+ 0]
          let y = vVelocities[vID &* 3 &+ 1]
          let z = vVelocities[vID &* 3 &+ 2]
          let mass = vMasses[vID]
          vKineticX.addProduct(mass, x * x)
          vKineticY.addProduct(mass, y * y)
          vKineticZ.addProduct(mass, z * z)
        }
        kinetic += MM4DoubleVector(vKineticX).sum()
        kinetic += MM4DoubleVector(vKineticY).sum()
        kinetic += MM4DoubleVector(vKineticZ).sum()
      }
    }
    return kinetic / 2
  }
  
  // What is the return value?
  // Does this modify the state in-place?
  //
  // Likely does not handle linear and angular velocities. For example, it may
  // be used to analyze new thermal velocities in isolation during
  // thermalization.
  func createThermalVelocities() {
    /*
     // Generate the list of Gaussian random numbers.
     OpenMM_SFMT::SFMT sfmt;
     init_gen_rand(randomSeed, sfmt);
     std::vector<double> randoms;
     while (randoms.size() < system.getNumParticles()*3) {
         double x, y, r2;
         do {
             x = 2.0*genrand_real2(sfmt)-1.0;
             y = 2.0*genrand_real2(sfmt)-1.0;
             r2 = x*x + y*y;
         } while (r2 >= 1.0 || r2 == 0.0);
         double multiplier = sqrt((-2.0*std::log(r2))/r2);
         randoms.push_back(x*multiplier);
         randoms.push_back(y*multiplier);
     }

     // Assign the velocities.
     std::vector<Vec3> velocities(system.getNumParticles(), Vec3());
     int nextRandom = 0;
     for (int i = 0; i < system.getNumParticles(); i++) {
         double mass = system.getParticleMass(i);
         if (mass != 0) {
             double velocityScale = sqrt(BOLTZ*temperature/mass);
             velocities[i] = Vec3(randoms[nextRandom++], randoms[nextRandom++], randoms[nextRandom++])*velocityScale;
         }
     }
     return velocities;
     */
  }
}
