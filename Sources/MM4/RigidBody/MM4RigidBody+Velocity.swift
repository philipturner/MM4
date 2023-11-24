//
//  MM4RigidBody+Velocity.swift
//
//
//  Created by Philip Turner on 11/20/23.
//

import Numerics

extension MM4RigidBody {
  /// The bulk + thermal velocity (in nanometers per picosecond) of each atom.
  ///
  /// Velocities attributed to anchors are ignored. They are replaced with a
  /// value determined by the bulk velocities.
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
    
    guard buffer.count == atoms.count else {
      fatalError("Velocity buffer was not the correct size.")
    }
    let baseAddress = buffer.baseAddress.unsafelyUnwrapped
    
    for vID in 0..<atoms.vectorCount {
      let (x, y, z) = swizzleToVectorWidth(vID, baseAddress: baseAddress)
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
    guard buffer.count == atoms.count else {
      fatalError("Velocity buffer was not the correct size.")
    }
    let baseAddress = buffer.baseAddress.unsafelyUnwrapped
    
    for vID in 0..<atoms.vectorCount {
      let x = storage.vVelocities[vID &* 3 &+ 0]
      let y = storage.vVelocities[vID &* 3 &+ 1]
      let z = storage.vVelocities[vID &* 3 &+ 2]
      swizzleFromVectorWidth((x, y, z), vID, baseAddress: baseAddress)
    }
  }
}

extension MM4RigidBody {
  func ensureAnchorVelocitiesValid() {
    guard storage.anchorVelocitiesValid == nil else {
      return
    }
    
    let epsilon: Float = 1e-4
    var firstVelocity: SIMD3<Float>?
    var valid = true
    
    for anchor in storage.anchors {
      let velocity = extractScalar(Int(anchor), array: storage.vVelocities)
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
    storage.anchorVelocitiesValid = valid
  }
  
  func createLinearVelocity() -> SIMD3<Float> {
    ensureAnchorVelocitiesValid()
    guard anchors.count == 0 else {
      let anchor = storage.anchors.first!
      let velocity = extractScalar(Int(anchor), array: storage.vVelocities)
      return velocity
    }
    guard atoms.count > 0 else {
      return .zero
    }
    
    var momentum: SIMD3<Double> = .zero
    withMasses { vMasses in
      withSegmentedLoop(chunk: 256) {
        var vMomentumX: MM4FloatVector = .zero
        var vMomentumY: MM4FloatVector = .zero
        var vMomentumZ: MM4FloatVector = .zero
        for vID in $0 {
          let x = storage.vVelocities[vID &* 3 &+ 0]
          let y = storage.vVelocities[vID &* 3 &+ 1]
          let z = storage.vVelocities[vID &* 3 &+ 2]
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
    return SIMD3<Float>(momentum / storage.totalMass)
  }
  
  func createAngularVelocity() -> Quaternion<Float> {
    ensureAnchorVelocitiesValid()
    ensureAngularMassCached()
    guard anchors.count <= 1, atoms.count > 0 else {
      return .zero
    }
    
    var momentum: SIMD3<Double> = .zero
    ensureCenterOfMassCached()
    withMasses { vMasses in
      let centerOfMass = storage.centerOfMass!
      withSegmentedLoop(chunk: 256) {
        var vMomentumX: MM4FloatVector = .zero
        var vMomentumY: MM4FloatVector = .zero
        var vMomentumZ: MM4FloatVector = .zero
        for vID in $0 {
          let rX = storage.vPositions[vID &* 3 &+ 0] - centerOfMass.x
          let rY = storage.vPositions[vID &* 3 &+ 1] - centerOfMass.y
          let rZ = storage.vPositions[vID &* 3 &+ 2] - centerOfMass.z
          let vX = storage.vVelocities[vID &* 3 &+ 0]
          let vY = storage.vVelocities[vID &* 3 &+ 1]
          let vZ = storage.vVelocities[vID &* 3 &+ 2]
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
    
    let inverse = storage.angularMass!.inverse
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

extension MM4RigidBody {
  func _createTotalKineticEnergy(
    _ vMasses: UnsafePointer<MM4FloatVector>
  ) -> Double {
    // Mask out the anchors while calculating this.
    var kinetic: Double = .zero
    withSegmentedLoop(chunk: 256) {
      var vKineticX: MM4FloatVector = .zero
      var vKineticY: MM4FloatVector = .zero
      var vKineticZ: MM4FloatVector = .zero
      for vID in $0 {
        let x = storage.vVelocities[vID &* 3 &+ 0]
        let y = storage.vVelocities[vID &* 3 &+ 1]
        let z = storage.vVelocities[vID &* 3 &+ 2]
        let mass = vMasses[vID]
        vKineticX.addProduct(mass, x * x)
        vKineticY.addProduct(mass, y * y)
        vKineticZ.addProduct(mass, z * z)
      }
      kinetic += MM4DoubleVector(vKineticX).sum()
      kinetic += MM4DoubleVector(vKineticY).sum()
      kinetic += MM4DoubleVector(vKineticZ).sum()
    }
    return kinetic / 2
  }
  
  func createTotalKineticEnergy() -> Double {
    guard atoms.count > 0 else {
      return 0
    }
    
    if anchors.count == 0 {
      return withMasses(_createTotalKineticEnergy)
    } else {
      let vMasses: UnsafeMutablePointer<MM4FloatVector> =
        .allocate(capacity: atoms.vectorCount)
      withMasses {
        vMasses.initialize(from: $0, count: atoms.count)
      }
      defer { vMasses.deallocate() }
      
      let masses: UnsafeMutablePointer<Float> = .init(OpaquePointer(vMasses))
      for anchor in self.anchors {
        masses[Int(anchor)] = 0
      }
      return _createTotalKineticEnergy(vMasses)
    }
  }
  
  // What is the return value?
  // Does this modify the state in-place?
  //
  // Likely does not handle linear and angular velocities. For example, it may
  // be used to analyze new thermal velocities in isolation during
  // thermalization.
  func createThermalVelocities() {
    
  }
}
