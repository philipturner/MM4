//
//  MM4RigidBody+Inertia.swift
//
//
//  Created by Philip Turner on 1/8/24.
//

extension MM4RigidBody {
  /// Symmetric matrix specifying the rigid body's moment of inertia.
  public var momentOfInertia: (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>) {
    fatalError("Not implemented.")
  }
  
  /// Rotate the object.
  /// - parameter angle: The angle to rotate, in radians.
  /// - parameter axis: The normalized vector to rotate around.
  ///
  /// If `axis` is not specified, the default value aligns with the current
  /// angular momentum. If the angular momentum is zero, you must enter zero
  /// for the angle.
  ///
  /// If the axis aligns with angular momentum, then angular momentum stays
  /// constant. Otherwise, the angular momentum is rotated along with the
  /// atom positions.
  public mutating func rotate(angle: Double, axis: SIMD3<Double>? = nil) {
    // TODO: Decompose the rotation into the diagonalized reference frame,
    // instead of adding a dependency on swift-numerics @ Quaternions. This
    // may simplify generation of the rotation matrix?
    //
    // Nvm - create a rotation matrix directly in the world space.
    
    ensureUniquelyReferenced()
    fatalError("Not implemented.")
  }
}

extension MM4RigidBodyStorage {
  
}

extension MM4RigidBodyStorage {
  func createMomentOfInertia() -> (SIMD3<Float>, SIMD3<Float>, SIMD3<Float>) {
    ensureCenterOfMassCached()
    guard let centerOfMass else {
      fatalError("This should never happen.")
    }
    guard atoms.count > 0 else {
      return (.zero, .zero, .zero)
    }
    
    var columns = (SIMD3<Double>.zero,
                   SIMD3<Double>.zero,
                   SIMD3<Double>.zero)
    withSegmentedLoop(chunk: 256) {
      var vXX: MM4FloatVector = .zero
      var vYY: MM4FloatVector = .zero
      var vZZ: MM4FloatVector = .zero
      var vXY: MM4FloatVector = .zero
      var vXZ: MM4FloatVector = .zero
      var vYZ: MM4FloatVector = .zero
      for vID in $0 {
        let x = vPositions[vID &* 3 &+ 0] - centerOfMass.x
        let y = vPositions[vID &* 3 &+ 1] - centerOfMass.y
        let z = vPositions[vID &* 3 &+ 2] - centerOfMass.z
        let mass = vMasses[vID]
        vXX.addProduct(mass, x * x)
        vYY.addProduct(mass, y * y)
        vZZ.addProduct(mass, z * z)
        vXY.addProduct(mass, x * y)
        vXZ.addProduct(mass, x * z)
        vYZ.addProduct(mass, y * z)
      }
      
      let XX = MM4DoubleVector(vXX).sum()
      let YY = MM4DoubleVector(vYY).sum()
      let ZZ = MM4DoubleVector(vZZ).sum()
      let XY = MM4DoubleVector(vXY).sum()
      let XZ = MM4DoubleVector(vXZ).sum()
      let YZ = MM4DoubleVector(vYZ).sum()
      columns.0 += SIMD3<Double>(YY + ZZ, -XY, -XZ)
      columns.1 += SIMD3<Double>(-XY, XX + ZZ, -YZ)
      columns.2 += SIMD3<Double>(-XZ, -YZ, XX + YY)
    }
    return (SIMD3<Float>(columns.0),
            SIMD3<Float>(columns.1),
            SIMD3<Float>(columns.2))
  }
}

extension MM4RigidBody {
  /// The net angular velocity (in radians per picosecond) of the entire
  /// object.
  public var angularVelocity: SIMD3<Double> {
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
      
      for vID in 0..<storage.atoms.vectorCount {
        let rX = storage.vPositions[vID &* 3 &+ 0] - centerOfMass.x
        let rY = storage.vPositions[vID &* 3 &+ 1] - centerOfMass.y
        let rZ = storage.vPositions[vID &* 3 &+ 2] - centerOfMass.z
        var vX: MM4FloatVector = .zero
        var vY: MM4FloatVector = .zero
        var vZ: MM4FloatVector = .zero
        
        // Un-apply the previous bulk angular velocity.
        do {
          let w = previous
          vX -= w.y * rZ - w.z * rY
          vY -= w.z * rX - w.x * rZ
          vZ -= w.x * rY - w.y * rX
        }
        
        // Apply the next bulk angular velocity.
        do {
          let w = next
          vX += w.y * rZ - w.z * rY
          vY += w.z * rX - w.x * rZ
          vZ += w.x * rY - w.y * rX
        }
        
        // Mask out the changes to velocity for anchors.
        let mass = storage.vMasses[vID]
        vX.replace(with: MM4FloatVector.zero, where: mass .== 0)
        vY.replace(with: MM4FloatVector.zero, where: mass .== 0)
        vZ.replace(with: MM4FloatVector.zero, where: mass .== 0)
        
        // Write the new velocities as offsets relative to the existing ones.
        storage.vVelocities[vID &* 3 &+ 0] += vX
        storage.vVelocities[vID &* 3 &+ 1] += vY
        storage.vVelocities[vID &* 3 &+ 2] += vZ
      }
      
      // Invalidate cached properties. Some could be restored, but err on the
      // side of simplicity for debugging.
      storage.eraseRarelyCachedProperties()
      storage.velocities = nil
    }
  }
}
