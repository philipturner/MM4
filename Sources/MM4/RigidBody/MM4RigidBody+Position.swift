//
//  MM4RigidBody+Position.swift
//  MM4
//
//  Created by Philip Turner on 11/20/23.
//

extension MM4RigidBody {
  /// The position (in nanometers) of each atom's nucleus.
  ///
  /// This is an ergonomic getter for positions. Behind the scenes, it
  /// automatically caches the results of swizzling all velocities into the
  /// non-vectorized representation.
  public var positions: [SIMD3<Float>] {
    _read {
      storage.ensurePositionsCached()
      yield storage.positions!
    }
  }
  
  public var vectorizedPositions: [MM4FloatVector] {
    _read {
      yield storage.vPositions
    }
    _modify {
      ensureUniquelyReferenced()
      storage.eraseRarelyCachedProperties()
      storage.centerOfMass = nil
      storage.positions = nil
      
      yield &storage.vPositions
      guard storage.vPositions.count == 3 * storage.atoms.vectorCount else {
        fatalError("Position buffer was not the correct size.")
      }
    }
  }
  
  @_specialize(where T == Double)
  @_specialize(where T == Float)
  public func getPositions<T: BinaryFloatingPoint>(
    _ buffer: UnsafeMutableBufferPointer<SIMD3<T>>
  ) {
    guard buffer.count == storage.atoms.count else {
      fatalError("Position buffer was not the correct size.")
    }
    let baseAddress = buffer.baseAddress.unsafelyUnwrapped
    
    for vID in 0..<storage.atoms.vectorCount {
      let x = storage.vPositions[vID &* 3 &+ 0]
      let y = storage.vPositions[vID &* 3 &+ 1]
      let z = storage.vPositions[vID &* 3 &+ 2]
      storage.swizzleFromVectorWidth((x, y, z), vID, baseAddress)
    }
  }
  
  public mutating func setPositions<T: BinaryFloatingPoint, U: Collection>(
    _ buffer: U
  ) where U.Element == SIMD3<T> {
    buffer.withContiguousStorageIfAvailable {
      setPositions($0)
    }
  }
  
  @_specialize(where T == Double)
  @_specialize(where T == Float)
  public mutating func setPositions<T: BinaryFloatingPoint>(
    _ buffer: UnsafeBufferPointer<SIMD3<T>>
  ) {
    ensureUniquelyReferenced()
    storage.eraseRarelyCachedProperties()
    storage.centerOfMass = nil
    storage.positions = nil
    
    guard buffer.count == storage.atoms.count else {
      fatalError("Position buffer was not the correct size.")
    }
    let baseAddress = buffer.baseAddress.unsafelyUnwrapped
    
    for vID in 0..<storage.atoms.vectorCount {
      let (x, y, z) = storage.swizzleToVectorWidth(vID, baseAddress)
      storage.vPositions[vID &* 3 &+ 0] = x
      storage.vPositions[vID &* 3 &+ 1] = y
      storage.vPositions[vID &* 3 &+ 2] = z
    }
  }
}

extension MM4RigidBodyStorage {
  func createPositions() -> [SIMD3<Float>] {
    let capacity = atoms.vectorCount * MM4VectorWidth
    let output = [SIMD3<Float>](unsafeUninitializedCapacity: capacity) {
      let baseAddress = $0.baseAddress.unsafelyUnwrapped
      for vID in 0..<atoms.vectorCount {
        let x = vPositions[vID &* 3 &+ 0]
        let y = vPositions[vID &* 3 &+ 1]
        let z = vPositions[vID &* 3 &+ 2]
        swizzleFromVectorWidth((x, y, z), vID, baseAddress)
      }
      $1 = atoms.count
    }
    return output
  }
  
  func createCenterOfMass() -> SIMD3<Float> {
    guard atoms.count > 0 else {
      return .zero
    }
    guard mass > 0 else {
      // I'm not sure any of the other properties can be computed in this
      // situation, so it is reasonable to throw a fatal error. It is better
      // than making up some behavior for the edge case, which require extensive
      // unit testing and complicate the code. The user can always use a custom
      // alternative to MM4RigidBody for parts that have every atom position
      // constrained.
      fatalError(
        "Could not create center of mass because all atoms had zero mass.")
    }
    var center: SIMD3<Double> = .zero
    withSegmentedLoop(chunk: 256) {
      var vCenterX: MM4FloatVector = .zero
      var vCenterY: MM4FloatVector = .zero
      var vCenterZ: MM4FloatVector = .zero
      for vID in $0 {
        let x = vPositions[vID &* 3 &+ 0]
        let y = vPositions[vID &* 3 &+ 1]
        let z = vPositions[vID &* 3 &+ 2]
        let mass = vMasses[vID]
        vCenterX.addProduct(mass, x)
        vCenterY.addProduct(mass, y)
        vCenterZ.addProduct(mass, z)
      }
      center.x += MM4DoubleVector(vCenterX).sum()
      center.y += MM4DoubleVector(vCenterY).sum()
      center.z += MM4DoubleVector(vCenterZ).sum()
    }
    return SIMD3<Float>(center) / mass
  }
  
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

// MARK: - Properties

extension MM4RigidBody {
  /// The total mass (in yoctograms) of all atoms.
  public var mass: Double {
    storage.mass
  }
  
  /// Center of mass.
  public var centerOfMass: SIMD3<Double> {
    _read {
      storage.ensureCenterOfMassCached()
      yield storage.centerOfMass!
    }
    _modify {
      ensureUniquelyReferenced()
      storage.ensureCenterOfMassCached()
      let previous = storage.centerOfMass!
      yield &storage.centerOfMass!
      
      let difference = storage.centerOfMass! - previous
      guard any(difference .!= .zero) else {
        return
      }
      for vID in 0..<storage.atoms.vectorCount {
        storage.vPositions[vID &* 3 &+ 0] += difference.x
        storage.vPositions[vID &* 3 &+ 1] += difference.y
        storage.vPositions[vID &* 3 &+ 2] += difference.z
      }
      if storage.atoms.count == 0 {
        storage.centerOfMass = .zero
      }
      
      // Invalidate cached properties. Some could be restored, but err on the
      // side of simplicity for debugging.
      storage.eraseRarelyCachedProperties()
      storage.positions = nil
    }
  }
  
  /// Symmetric matrix specifying the rigid body's moment of inertia.
  public var momentOfInertia: (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>) {
    storage.ensureMomentOfInertiaCached()
    return storage.momentOfInertia!
  }
  
  /// Rotate around the same axis as the current angular velocity.
  public mutating func rotate(angle: Double) {
    ensureUniquelyReferenced()
    fatalError("Not implemented.")
  }
}
