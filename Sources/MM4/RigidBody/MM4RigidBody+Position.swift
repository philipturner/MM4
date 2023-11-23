//
//  MM4RigidBody+Position.swift
//
//
//  Created by Philip Turner on 11/20/23.
//

extension MM4RigidBody {
  /// The position (in nanometers) of each atom's nucleus.
  public var positions: [SIMD3<Float>] {
    // _modify not yet supported b/c it requires very complex caching logic.
    _read { fatalError("Not implemented.") }
  }
  
  @_specialize(where T == Double)
  @_specialize(where T == Float)
  mutating func setPositions<T: BinaryFloatingPoint>(
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
  func getPositions<T: BinaryFloatingPoint>(
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
  func createCenterOfMass() {
    var vAccumulatorX: MM4DoubleVector = .zero
    var vAccumulatorY: MM4DoubleVector = .zero
    var vAccumulatorZ: MM4DoubleVector = .zero
    masses.withUnsafeBufferPointer {
      let rawMasses = OpaquePointer($0.baseAddress)
      let vMasses = UnsafeRawPointer(rawMasses)!
        .assumingMemoryBound(to: MM4FloatVector.self)
      
      var vCenterX: MM4FloatVector = .zero
      var vCenterY: MM4FloatVector = .zero
      var vCenterZ: MM4FloatVector = .zero
      for vID in 0..<atoms.vectorCount {
        let x = storage.vPositions[vID &* 3 &+ 0]
        let y = storage.vPositions[vID &* 3 &+ 1]
        let z = storage.vPositions[vID &* 3 &+ 2]
        
        let vMass = vMasses[vID]
        vCenterX += x * vMass
        vCenterY += y * vMass
        vCenterZ += z * vMass
        
        if vID % 256 == 0 {
          vAccumulatorX += MM4DoubleVector(vCenterX)
          vAccumulatorY += MM4DoubleVector(vCenterY)
          vAccumulatorZ += MM4DoubleVector(vCenterZ)
          vCenterX = .zero
          vCenterY = .zero
          vCenterZ = .zero
        }
      }
      vAccumulatorX += MM4DoubleVector(vCenterX)
      vAccumulatorY += MM4DoubleVector(vCenterY)
      vAccumulatorZ += MM4DoubleVector(vCenterZ)
    }
    
    var center: SIMD3<Double> = .zero
    for lane in 0..<4 {
      center.x += Double(vAccumulatorX[lane])
      center.y += Double(vAccumulatorY[lane])
      center.z += Double(vAccumulatorZ[lane])
    }
    storage.centerOfMass = center / mass
  }
  
  func createAngularMass() {
    
  }
}
