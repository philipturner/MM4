//
//  MM4RigidBody+Position.swift
//
//
//  Created by Philip Turner on 11/20/23.
//

extension MM4RigidBody {
  /// The position (in nanometers) of each atom's nucleus.
  public var positions: [SIMD3<Float>] {
    _read { fatalError("Not implemented.") }
    _modify { fatalError("Not implemented.") }
  }
  
  @_specialize(where T == Double)
  @_specialize(where T == Float)
  mutating func setPositions<T: BinaryFloatingPoint>(
    _ buffer: UnsafeBufferPointer<SIMD3<T>>
  ) {
    ensureUniquelyReferenced()
    guard buffer.count == atoms.count else {
      fatalError("Position buffer was not the correct size.")
    }
    let baseAddress = buffer.baseAddress.unsafelyUnwrapped
    
    #if false
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
        vPositions[vID * 3 + 0] = x
        vPositions[vID * 3 + 1] = y
        vPositions[vID * 3 + 2] = z
        
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
    centerOfMass.value = center / mass
    #endif
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
    fatalError("Not implemented.")
  }
}
