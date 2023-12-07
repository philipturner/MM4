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
    _read {
      storage.ensureVelocitiesCached()
      yield storage.velocities!
    }
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
  func createVelocities() -> [SIMD3<Float>] {
    let capacity = atoms.vectorCount * MM4VectorWidth
    let output = [SIMD3<Float>](unsafeUninitializedCapacity: capacity) {
      let baseAddress = $0.baseAddress.unsafelyUnwrapped
      for vID in 0..<atoms.vectorCount {
        let x = vVelocities[vID &* 3 &+ 0]
        let y = vVelocities[vID &* 3 &+ 1]
        let z = vVelocities[vID &* 3 &+ 2]
        swizzleFromVectorWidth((x, y, z), vID, baseAddress)
      }
      $1 = atoms.count
    }
    return output
  }
  
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
  
  // WARNING: When there are anchors, this returns something besides the bulk
  // velocity. It is the linear velocity of non-anchors.
  func createLinearVelocity() -> SIMD3<Float> {
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
  
  // WARNING: This returns a nonzero angular velocity, even when we should store
  // zero (anchors > 1). It is the angular velocity of non-anchors.
  func createAngularVelocity() -> Quaternion<Float> {
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
  // Linear kinetic energy from non-anchor atoms.
  func createLinearKineticEnergy() -> Double {
    ensureLinearVelocityCached()
    guard let linearVelocity else {
      fatalError("This should never happen.")
    }
    
    let v = SIMD3<Double>(linearVelocity)
    return nonAnchorMass * (v * v).sum() / 2
  }
  
  // Angular kinetic energy about an angular mass defined by non-anchor atoms.
  func createAngularKineticEnergy() -> Double {
    ensureMomentOfInertiaCached()
    ensureAngularVelocityCached()
    guard let momentOfInertia,
          let angularVelocity else {
      fatalError("This should never happen.")
    }
    
    let I = momentOfInertia
    let w = SIMD3<Double>(quaternionToVector(angularVelocity))
    let velocityX = I.columns.0 * w.x
    let velocityY = I.columns.1 * w.y
    let velocityZ = I.columns.2 * w.z
    let Iw = velocityX + velocityY + velocityZ
    return (w * Iw).sum() / 2
  }
  
  // Total translational kinetic energy from non-anchor atoms.
  func createTotalKineticEnergy() -> Double {
    // This function should return the same value as separating the bulk and
    // thermal velocities, computing the energies separately, then summing.
    
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
  
  func createThermalVelocities() {
    ensureCenterOfMassCached()
    ensureLinearVelocityCached()
    ensureAngularVelocityCached()
    
    // Express that every bulk quantity, except thermal kinetic energy, is the
    // same as before thermalization.
    guard let centerOfMass,
          let constantLinearVelocity = linearVelocity,
          let constantAngularVelocity = angularVelocity,
          let newThermalKineticEnergy = thermalKineticEnergy else {
      fatalError("This should never happen.")
    }
    
    // Handle the special case where velocity rescaling would cause division by
    // zero.
    if anchors.count == atoms.count {
      guard newThermalKineticEnergy == 0 else {
        fatalError(
          "Cannot create thermal velocities when all atoms are anchors.")
      }
      ensureAnchorVelocitiesValid()
      return
    }
    
    // Reference implementation from OpenMM.
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
    @inline(__always)
    func gaussian(_ seed: MM4UInt32Vector) -> (
      x: MM4FloatVector, y: MM4FloatVector, r2: MM4FloatVector
    ) {
      let seedPair = unsafeBitCast(seed, to: MM4UInt16VectorPair.self)
      let floatPair = MM4FloatVectorPair(seedPair) / Float(UInt16.max)
      let x = 2 * floatPair.evenHalf - 1
      let y = 2 * floatPair.oddHalf - 1
      let r2 = x * x + y * y
      return (x, y, r2)
    }
    
    // First, generate a unitless list of velocities. Pad the list to 64 more
    // than required (21 atoms) to decrease the overhead of repeated
    // reinitialization in the final loop iterations.
    var scalarsRequired = 3 * atoms.vectorCount * MM4VectorWidth
    scalarsRequired = (scalarsRequired + 4 - 1) / 4 * 4
    let scalarsCapacity = 64 + scalarsRequired
    let scalarsPointer: UnsafeMutablePointer<UInt16> =
      .allocate(capacity: scalarsCapacity)
    defer { scalarsPointer.deallocate() }
    
    // Repeatedly compact the list, removing pairs that failed.
    var scalarsFinished = 0
    var generator = SystemRandomNumberGenerator()
    while scalarsFinished < scalarsRequired {
      // Round down to UInt64 alignment.
      scalarsFinished = scalarsFinished / 4 * 4
      
      // Fill up to the capacity, rather than the required amount.
      let quadsToGenerate = (scalarsCapacity - scalarsFinished) / 4
      let quadsPointer: UnsafeMutablePointer<UInt64> = .init( OpaquePointer(scalarsPointer + scalarsFinished / 4))
      for i in 0..<quadsToGenerate {
        quadsPointer[i] = generator.next()
      }
      
      // The first of these pointers acts as a cursor.
      var pairsPointer: UnsafeMutablePointer<UInt32> = .init( OpaquePointer(scalarsPointer + scalarsFinished / 4))
      let pairsVectorPointer: UnsafeMutablePointer<MM4UInt32Vector> = .init( OpaquePointer(scalarsPointer + scalarsFinished / 4))
      
      for vID in 0..<quadsToGenerate * 4 / MM4VectorWidth {
        let seed = pairsVectorPointer[vID]
        let (_, _, r2) = gaussian(seed)
        
        let mask = (r2 .< 1) .& (r2 .!= 0)
        for lane in 0..<MM4VectorWidth {
          if mask[lane] {
            pairsPointer.pointee = seed[lane]
            pairsPointer += 1
          }
        }
      }
      let newScalarsPointer: UnsafeMutablePointer<UInt16> = .init(
        OpaquePointer(pairsPointer))
      scalarsFinished = newScalarsPointer - scalarsPointer
    }
    
    // Generate velocities, zeroing out the ones from anchors.
    //
    // E_particle = E_system / (atoms.count - anchors.count)
    // E_openmm = 3/2 kT
    // kT = 2/3 * (translational kinetic energy)
    //
    // mass = anchor ? INF : mass
    // v_rms = sqrt(2/3 * E_particle / mass)
    // v = (gaussian(0, 1), gaussian(0, 1), gaussian(0, 1)) * v_rms
    let xyPointer: UnsafePointer<MM4UInt32Vector> = .init(
      OpaquePointer(scalarsPointer))
    let zPointer: UnsafePointer<MM4UInt16Vector> = .init(
      OpaquePointer(scalarsPointer + 2 * MM4VectorWidth * atoms.vectorCount))
    let particleCount = Double(atoms.count - anchors.count)
    let particleEnergy = newThermalKineticEnergy / particleCount
    
    // Units: zJ -> kJ/mol
    let particleEnergyTerm = Float(
      Double(2.0 / 3) * MM4KJPerMolPerZJ * particleEnergy)
    
    withMasses(nonAnchorMasses) { vMasses in
      for vID in 0..<atoms.vectorCount {
        var xGaussian: MM4FloatVector
        var yGaussian: MM4FloatVector
        var zGaussian: MM4FloatVector
        do {
          let z = MM4FloatVector(zPointer[vID]) / Float(UInt16.max)
          let zLow = 2 * z.evenHalf - 1
          let zHigh = 2 * z.oddHalf - 1
          let zR2 = zLow * zLow + zHigh * zHigh
          let (x, y, xyR2) = gaussian(xyPointer[vID])
          var (xyLog, zLog) = (xyR2, zR2)
          
          // There is no simple way to access vectorized transcendentals on
          // non-Apple platforms. Ideally, one would copy code from Sleef. We
          // can only keep our fingers crossed that the compiler will
          // "auto-vectorize" this transcendental function, which may have
          // control flow operations that prevent it from actually vectorizing.
          for lane in 0..<MM4VectorWidth {
            xyLog[lane] = Float.log(xyLog[lane])
          }
          for lane in 0..<MM4VectorWidth / 2 {
            zLog[lane] = Float.log(zLog[lane])
          }
          
          let xyMultiplier = (-2 * xyLog / xyR2).squareRoot()
          let zMultiplier = (-2 * zLog / zR2).squareRoot()
          xGaussian = x * xyMultiplier
          yGaussian = y * xyMultiplier
          
          var zBroadcasted: MM4FloatVector = .zero
          zBroadcasted.evenHalf = zMultiplier
          zBroadcasted.oddHalf = zMultiplier
          zGaussian = z * zBroadcasted
        }
        
        var mass = vMasses[vID]
        mass.replace(with: .greatestFiniteMagnitude, where: mass .== 0)
        let velocityScale = (particleEnergyTerm / mass).squareRoot()
        vVelocities[vID &* 3 &+ 0] = xGaussian * velocityScale
        vVelocities[vID &* 3 &+ 1] = yGaussian * velocityScale
        vVelocities[vID &* 3 &+ 2] = zGaussian * velocityScale
      }
    }
    
    // Query the bulk linear and angular momentum.
    let linearDrift = createLinearVelocity()
    let angularDrift = createAngularVelocity()
    let wDrift = quaternionToVector(angularDrift)
    
    // Set momentum to zero and calculate the modified thermal energy.
    var correctedThermalKineticEnergy: Double = .zero
    withMasses(nonAnchorMasses) { vMasses in
      withSegmentedLoop(chunk: 256) {
        var vKineticX: MM4FloatVector = .zero
        var vKineticY: MM4FloatVector = .zero
        var vKineticZ: MM4FloatVector = .zero
        for vID in $0 {
          let rX = vPositions[vID &* 3 &+ 0] - centerOfMass.x
          let rY = vPositions[vID &* 3 &+ 1] - centerOfMass.y
          let rZ = vPositions[vID &* 3 &+ 2] - centerOfMass.z
          var vX = vVelocities[vID &* 3 &+ 0]
          var vY = vVelocities[vID &* 3 &+ 1]
          var vZ = vVelocities[vID &* 3 &+ 2]
          
          // Apply the correction to linear velocity.
          vX -= linearDrift.x
          vY -= linearDrift.y
          vZ -= linearDrift.z
          
          // Apply the correction to angular velocity.
          let w = wDrift
          vX -= w.y * rZ - w.z * rY
          vY -= w.z * rX - w.x * rZ
          vZ -= w.x * rY - w.y * rX
          
          // Mask out the changes to anchor velocities.
          let mass = vMasses[vID]
          vX.replace(with: MM4FloatVector.zero, where: mass .== 0)
          vY.replace(with: MM4FloatVector.zero, where: mass .== 0)
          vZ.replace(with: MM4FloatVector.zero, where: mass .== 0)
          vKineticX.addProduct(mass, vX * vX)
          vKineticY.addProduct(mass, vY * vY)
          vKineticZ.addProduct(mass, vZ * vZ)
          vVelocities[vID &* 3 &+ 0] = vX
          vVelocities[vID &* 3 &+ 1] = vY
          vVelocities[vID &* 3 &+ 2] = vZ
        }
        correctedThermalKineticEnergy += MM4DoubleVector(vKineticX).sum()
        correctedThermalKineticEnergy += MM4DoubleVector(vKineticY).sum()
        correctedThermalKineticEnergy += MM4DoubleVector(vKineticZ).sum()
      }
    }
    
    // Rescale thermal velocities and superimpose over bulk velocities.
    let velocityScale = Float((
      newThermalKineticEnergy / correctedThermalKineticEnergy).squareRoot())
    let w = quaternionToVector(constantAngularVelocity)
    
    withMasses(nonAnchorMasses) { vMasses in
      for vID in 0..<atoms.vectorCount {
        let rX = vPositions[vID &* 3 &+ 0] - centerOfMass.x
        let rY = vPositions[vID &* 3 &+ 1] - centerOfMass.y
        let rZ = vPositions[vID &* 3 &+ 2] - centerOfMass.z
        var vX = vVelocities[vID &* 3 &+ 0]
        var vY = vVelocities[vID &* 3 &+ 1]
        var vZ = vVelocities[vID &* 3 &+ 2]
        
        // Apply the correction to thermal velocity.
        vX *= velocityScale
        vY *= velocityScale
        vZ *= velocityScale
        
        // Apply the bulk angular velocity.
        vX += w.y * rZ - w.z * rY
        vY += w.z * rX - w.x * rZ
        vZ += w.x * rY - w.y * rX
        
        // Mask out the changes to angular velocity for anchors.
        let mass = vMasses[vID]
        vX.replace(with: MM4FloatVector.zero, where: mass .== 0)
        vY.replace(with: MM4FloatVector.zero, where: mass .== 0)
        vZ.replace(with: MM4FloatVector.zero, where: mass .== 0)
        
        // Apply the bulk linear velocity.
        vX += constantLinearVelocity.x
        vY += constantLinearVelocity.y
        vZ += constantLinearVelocity.z
        
        vVelocities[vID &* 3 &+ 0] = vX
        vVelocities[vID &* 3 &+ 1] = vY
        vVelocities[vID &* 3 &+ 2] = vZ
      }
    }
  }
}
