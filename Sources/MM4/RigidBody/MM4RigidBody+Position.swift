//
//  MM4RigidBody+Position.swift
//
//
//  Created by Philip Turner on 11/20/23.
//

extension MM4RigidBody {
  /// The position (in nanometers) of each atom's nucleus.
  public var positions: [SIMD3<Float>] {
    // _modify not supported b/c it requires very complex caching logic.
    // Workaround: use the exposed setPositions function.
    _read {
      storage.ensurePositionsCached()
      yield storage.positions!
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
}

extension MM4RigidBodyStorage {
  func createPositions() -> [SIMD3<Float>] {
    var output: [SIMD3<Float>] = []
    output.reserveCapacity(atoms.vectorCount * MM4VectorWidth)
    output.withUnsafeMutableBufferPointer { buffer in
      let baseAddress = buffer.baseAddress.unsafelyUnwrapped
      
      for vID in 0..<atoms.vectorCount {
        let x = vPositions[vID &* 3 &+ 0]
        let y = vPositions[vID &* 3 &+ 1]
        let z = vPositions[vID &* 3 &+ 2]
        swizzleFromVectorWidth((x, y, z), vID, baseAddress)
      }
    }
    return output
  }
  
  func createCenter(
    _ vMasses: UnsafePointer<MM4FloatVector>
  ) -> SIMD3<Double> {
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
    return center
  }
  
  func createCenterOfMass() -> SIMD3<Float> {
    guard atoms.count > 0 else {
      return .zero
    }
    if anchors.count == 0 {
      let center = withMasses(nonAnchorMasses, createCenter)
      return SIMD3<Float>(center / nonAnchorMass)
    } else {
      let center = withMasses(anchorMasses, createCenter)
      return SIMD3<Float>(center / anchorMass)
    }
  }
  
  func createMomentOfInertia() -> MM4MomentOfInertia {
    ensureCenterOfMassCached()
    guard let centerOfMass else {
      fatalError("This should never happen.")
    }
    guard atoms.count > 0 else {
      return MM4MomentOfInertia()
    }
    
    var momentOfInertia = MM4MomentOfInertia()
    withMasses(nonAnchorMasses) { vMasses in
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
        momentOfInertia.columns.0 += SIMD3<Double>(YY + ZZ, -XY, -XZ)
        momentOfInertia.columns.1 += SIMD3<Double>(-XY, XX + ZZ, -YZ)
        momentOfInertia.columns.2 += SIMD3<Double>(-XZ, -YZ, XX + YY)
      }
    }
    return momentOfInertia
  }
}

