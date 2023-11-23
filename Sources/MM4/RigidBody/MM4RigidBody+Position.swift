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
    _read { fatalError("Not implemented.") }
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
    
    guard buffer.count == atoms.count else {
      fatalError("Position buffer was not the correct size.")
    }
    let baseAddress = buffer.baseAddress.unsafelyUnwrapped
    
    for vID in 0..<atoms.vectorCount {
      var x: MM4FloatVector = .zero
      var y: MM4FloatVector = .zero
      var z: MM4FloatVector = .zero
      
      if vID == atoms.vectorCount &- 1 {
        let remaining = atoms.count &- vID &* MM4VectorWidth
        for lane in 0..<remaining {
          let position = baseAddress[vID &* MM4VectorWidth &+ lane]
          x[lane] = Float(position.x)
          y[lane] = Float(position.y)
          z[lane] = Float(position.z)
        }
      } else {
        for lane in 0..<MM4VectorWidth {
          let position = baseAddress[vID &* MM4VectorWidth &+ lane]
          x[lane] = Float(position.x)
          y[lane] = Float(position.y)
          z[lane] = Float(position.z)
        }
      }
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
    guard buffer.count == atoms.count else {
      fatalError("Position buffer was not the correct size.")
    }
    let baseAddress = buffer.baseAddress.unsafelyUnwrapped
    
    for vID in 0..<atoms.vectorCount {
      let x = storage.vPositions[vID &* 3 &+ 0]
      let y = storage.vPositions[vID &* 3 &+ 1]
      let z = storage.vPositions[vID &* 3 &+ 2]
      
      if vID == atoms.vectorCount &- 1 {
        let remaining = atoms.count &- vID &* MM4VectorWidth
        for lane in 0..<remaining {
          var position: SIMD3<T> = .zero
          position.x = T(x[lane])
          position.y = T(y[lane])
          position.z = T(z[lane])
          baseAddress[vID &* MM4VectorWidth &+ lane] = position
        }
      } else {
        for lane in 0..<MM4VectorWidth {
          var position: SIMD3<T> = .zero
          position.x = T(x[lane])
          position.y = T(y[lane])
          position.z = T(z[lane])
          baseAddress[vID &* MM4VectorWidth &+ lane] = position
        }
      }
    }
  }
}

extension MM4RigidBody {
  func withMasses<T>(_ closure: (UnsafePointer<MM4FloatVector>) -> T) -> T {
    self.masses.withUnsafeBufferPointer {
      let rawMasses = OpaquePointer($0.baseAddress)
      let vMasses = UnsafeRawPointer(rawMasses)!
        .assumingMemoryBound(to: MM4FloatVector.self)
      return closure(vMasses)
    }
  }
  
  func withSegmentedLoop(chunk: Int, _ closure: (Range<Int>) -> Void) {
    var loopEnd = 0
    while loopEnd < atoms.vectorCount {
      let loopStart = loopEnd
      loopEnd = min(loopEnd &+ chunk, atoms.vectorCount)
      closure(loopStart..<loopEnd)
    }
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
        let x = storage.vPositions[vID &* 3 &+ 0]
        let y = storage.vPositions[vID &* 3 &+ 1]
        let z = storage.vPositions[vID &* 3 &+ 2]
        vCenterX += x * vMasses[vID]
        vCenterY += y * vMasses[vID]
        vCenterZ += z * vMasses[vID]
      }
      center.x += MM4DoubleVector(vCenterX).sum()
      center.y += MM4DoubleVector(vCenterY).sum()
      center.z += MM4DoubleVector(vCenterZ).sum()
    }
    return center
  }
  
  func createCenterOfMass() -> SIMD3<Float> {
    if anchors.count == 0 {
      let center = withMasses(createCenter)
      return SIMD3<Float>(center / storage.mass)
    } else {
      let vMasses: UnsafeMutablePointer<MM4FloatVector> =
        .allocate(capacity: atoms.vectorCount)
      vMasses.initialize(repeating: .zero, count: atoms.vectorCount)
      defer { vMasses.deallocate() }
      
      let masses: UnsafeMutablePointer<Float> = .init(OpaquePointer(vMasses))
      var mass: Double = .zero
      for anchor in self.anchors {
        let value = self.masses[Int(anchor)]
        masses[Int(anchor)] = value
        mass += Double(value)
      }
      let center = createCenter(vMasses)
      return SIMD3<Float>(center / mass)
    }
  }
  
  func createAngularMass() -> MM4AngularMass {
    var angularMass = MM4AngularMass()
    guard anchors.count <= 1 else {
      return angularMass
    }
    
    ensureCenterOfMassCached()
    withMasses { vMasses in
      let centerOfMass = storage.centerOfMass!
      withSegmentedLoop(chunk: 256) {
        var vXX: MM4FloatVector = .zero
        var vYY: MM4FloatVector = .zero
        var vZZ: MM4FloatVector = .zero
        var vXY: MM4FloatVector = .zero
        var vXZ: MM4FloatVector = .zero
        var vYZ: MM4FloatVector = .zero
        for vID in $0 {
          let x = storage.vPositions[vID &* 3 &+ 0] - centerOfMass.x
          let y = storage.vPositions[vID &* 3 &+ 1] - centerOfMass.y
          let z = storage.vPositions[vID &* 3 &+ 2] - centerOfMass.z
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
        angularMass.columns.0 += SIMD3(YY + ZZ, -XY, -XZ)
        angularMass.columns.1 += SIMD3(-XY, XX + ZZ, -YZ)
        angularMass.columns.2 += SIMD3(-XZ, -YZ, XX + YY)
      }
    }
    return angularMass
  }
}
