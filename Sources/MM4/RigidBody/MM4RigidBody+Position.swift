//
//  MM4RigidBody+Position.swift
//
//
//  Created by Philip Turner on 11/20/23.
//

extension MM4RigidBody {
  public var positions: [SIMD3<Float>] {
    // We may want a similar API to that in MM4ForceField, which lazily caches
    // a version of the data that has been transformed. That way, the user
    // doesn't face O(n^2) computational complexity when treating it like an
    // array with per-element granularity.
    get { fatalError("Add wrapper class for lazily recomputing transformed positions.") }
    set { fatalError("Add wrapper class for lazily recomputing transformed positions.") }
  }
  
  @_specialize(where T == Double)
  @_specialize(where T == Float)
  mutating func setPositions<T: BinaryFloatingPoint>(
    _ buffer: UnsafeBufferPointer<SIMD3<T>>
  ) {
    guard buffer.count == atomCount else {
      fatalError("Position buffer was not the correct size.")
    }
    let baseAddress = buffer.baseAddress.unsafelyUnwrapped
    
    var vCenterX: MM4FloatVector = .zero
    var vCenterY: MM4FloatVector = .zero
    var vCenterZ: MM4FloatVector = .zero
    masses.withUnsafeBufferPointer {
      let rawMasses = OpaquePointer($0.baseAddress)
      let vMasses = UnsafeRawPointer(rawMasses)!
        .assumingMemoryBound(to: MM4FloatVector.self)
      
      for vID in 0..<atomVectorCount {
        var x: MM4FloatVector = .zero
        var y: MM4FloatVector = .zero
        var z: MM4FloatVector = .zero
        
        if vID == atomVectorCount - 1 {
          let remaining = atomCount - vID * MM4VectorWidth
          for lane in 0..<remaining {
            let position = baseAddress[vID * MM4VectorWidth + lane]
            x[lane] = Float(position.x)
            y[lane] = Float(position.y)
            z[lane] = Float(position.z)
          }
        } else {
          for lane in 0..<MM4VectorWidth {
            let position = baseAddress[vID * MM4VectorWidth + lane]
            x[lane] = Float(position.x)
            y[lane] = Float(position.y)
            z[lane] = Float(position.z)
          }
        }
        _positions[vID * 3 + 0] = x
        _positions[vID * 3 + 1] = y
        _positions[vID * 3 + 2] = z
        
        let vMass = vMasses[vID]
        vCenterX += x * vMass
        vCenterY += y * vMass
        vCenterZ += z * vMass
      }
    }
    
    var center: SIMD3<Double> = .zero
    for lane in 0..<4 {
      center.x += Double(vCenterX[lane])
      center.y += Double(vCenterY[lane])
      center.z += Double(vCenterZ[lane])
    }
    _centerOfMass.value = center / mass
  }
  
  @_specialize(where T == Double)
  @_specialize(where T == Float)
  func getPositions<T: BinaryFloatingPoint>(
    _ buffer: UnsafeMutableBufferPointer<SIMD3<T>>
  ) {
    guard buffer.count == atomCount else {
      fatalError("Position buffer was not the correct size.")
    }
    let baseAddress = buffer.baseAddress.unsafelyUnwrapped
    fatalError("Not implemented.")
  }
}
