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
    //
    // Workaround: use the velocities property of the descriptor.
    // Workaround: use the exposed setVelocities function.
    _read {
      storage.ensureVelocitiesCached()
      yield storage.velocities!
    }
  }
  
  @_specialize(where T == Double)
  @_specialize(where T == Float)
  public mutating func setVelocities<T: BinaryFloatingPoint>(
    _ array: Array<SIMD3<T>>
  ) {
    array.withUnsafeBufferPointer {
      setVelocities($0)
    }
  }
  
  @_specialize(where T == Double)
  @_specialize(where T == Float)
  public mutating func setVelocities<T: BinaryFloatingPoint>(
    _ buffer: UnsafeBufferPointer<SIMD3<T>>
  ) {
    ensureUniquelyReferenced()
    storage.eraseRarelyCachedProperties()
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
    storage.ensureAnchorVelocitiesValid()
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
  
  // WARNING: Call this every time the velocities change.
  func ensureAnchorVelocitiesValid() {
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

// MARK: - Properties

extension MM4RigidBody {
  /// The net angular velocity (in radians per picosecond) of the non-anchor
  /// atoms.
  public var angularVelocity: Quaternion<Float> {
    _read {
      storage.ensureAngularVelocityCached()
      yield storage.angularVelocity!
    }
    _modify {
      ensureUniquelyReferenced()
      storage.ensureCenterOfMassCached()
      storage.ensureAngularVelocityCached()
      let previous = storage.angularVelocity!
      yield &storage.angularVelocity!
      
      let next = storage.angularVelocity!
      guard next != previous else {
        return
      }
      guard let centerOfMass = storage.centerOfMass else {
        fatalError("This should never happen.")
      }
      
      let previousW = quaternionToVector(previous)
      let nextW = quaternionToVector(next)
      storage.withMasses(storage.nonAnchorMasses) { vMasses in
        for vID in 0..<storage.atoms.vectorCount {
          let rX = storage.vPositions[vID &* 3 &+ 0] - centerOfMass.x
          let rY = storage.vPositions[vID &* 3 &+ 1] - centerOfMass.y
          let rZ = storage.vPositions[vID &* 3 &+ 2] - centerOfMass.z
          var vX: MM4FloatVector = .zero
          var vY: MM4FloatVector = .zero
          var vZ: MM4FloatVector = .zero
          
          // Un-apply the previous bulk angular velocity.
          do {
            let w = previousW
            vX -= w.y * rZ - w.z * rY
            vY -= w.z * rX - w.x * rZ
            vZ -= w.x * rY - w.y * rX
          }
          
          // Apply the next bulk angular velocity.
          do {
            let w = nextW
            vX += w.y * rZ - w.z * rY
            vY += w.z * rX - w.x * rZ
            vZ += w.x * rY - w.y * rX
          }
          
          // Mask out the changes to velocity for anchors.
          let mass = vMasses[vID]
          vX.replace(with: MM4FloatVector.zero, where: mass .== 0)
          vY.replace(with: MM4FloatVector.zero, where: mass .== 0)
          vZ.replace(with: MM4FloatVector.zero, where: mass .== 0)
          
          // Write the new velocities as offsets relative to the existing ones.
          storage.vVelocities[vID &* 3 &+ 0] += vX
          storage.vVelocities[vID &* 3 &+ 1] += vY
          storage.vVelocities[vID &* 3 &+ 2] += vZ
        }
      }
      
      // Invalidate cached properties. Some could be restored, but err on the
      // side of simplicity for debugging.
      storage.eraseRarelyCachedProperties()
      storage.velocities = nil
      storage.ensureAnchorVelocitiesValid()
    }
  }
  
  /// The net linear velocity (in nanometers per picosecond) of the entire
  /// object.
  ///
  /// Every anchor's velocity is set to the rigid body's linear velocity.
  /// When importing velocities, all anchors must have the same velocity.
  public var linearVelocity: SIMD3<Float> {
    _read {
      storage.ensureLinearVelocityCached()
      yield storage.linearVelocity!
    }
    _modify {
      ensureUniquelyReferenced()
      storage.ensureLinearVelocityCached()
      let previous = storage.linearVelocity!
      yield &storage.linearVelocity!
      
      let velocityDifference = storage.linearVelocity! - previous
      guard any(velocityDifference .!= .zero) else {
        return
      }
      for vID in 0..<storage.atoms.vectorCount {
        storage.vVelocities[vID &* 3 &+ 0] += velocityDifference.x
        storage.vVelocities[vID &* 3 &+ 1] += velocityDifference.y
        storage.vVelocities[vID &* 3 &+ 2] += velocityDifference.z
      }
      if storage.atoms.count == 0 {
        storage.linearVelocity = .zero
      }
      
      // Invalidate cached properties. Some could be restored, but err on the
      // side of simplicity for debugging.
      storage.eraseRarelyCachedProperties()
      storage.velocities = nil
      storage.ensureAnchorVelocitiesValid()
    }
  }
}
