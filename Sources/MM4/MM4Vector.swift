//
//  MM4Vector.swift
//  MM4  
//
//  Created by Philip Turner on 11/22/23.
//

// MARK: - Vector Types

let MM4VectorWidth: Int = 4

typealias MM4FloatVector = SIMD4<Float>
typealias MM4DoubleVector = SIMD4<Double>
typealias MM4Int8Vector = SIMD4<Int8>
typealias MM4Int16Vector = SIMD4<Int16>
typealias MM4Int32Vector = SIMD4<Int32>
typealias MM4Int64Vector = SIMD4<Int64>
typealias MM4UInt8Vector = SIMD4<UInt8>
typealias MM4UInt16Vector = SIMD4<UInt16>
typealias MM4UInt32Vector = SIMD4<UInt32>
typealias MM4UInt64Vector = SIMD4<UInt64>

// MARK: - Looping Over Vectors

extension MM4RigidBodyStorage {
  @_transparent
  func swizzleFromVectorWidth<T: BinaryFloatingPoint>(
    _ vectorized: (MM4FloatVector, MM4FloatVector, MM4FloatVector),
    _ vID: Int, _ baseAddress: UnsafeMutablePointer<SIMD3<T>>
  ) {
    let (x, y, z) = vectorized
    if vID == atoms.vectorCount &- 1 {
      let remaining = atoms.count &- vID &* MM4VectorWidth
      for lane in 0..<remaining {
        var vector3: SIMD3<T> = .zero
        vector3.x = T(x[lane])
        vector3.y = T(y[lane])
        vector3.z = T(z[lane])
        baseAddress[vID &* MM4VectorWidth &+ lane] = vector3
      }
    } else {
      for lane in 0..<MM4VectorWidth {
        var vector3: SIMD3<T> = .zero
        vector3.x = T(x[lane])
        vector3.y = T(y[lane])
        vector3.z = T(z[lane])
        baseAddress[vID &* MM4VectorWidth &+ lane] = vector3
      }
    }
  }
  
  @_transparent
  func swizzleToVectorWidth<T: BinaryFloatingPoint>(
    _ vID: Int, _ baseAddress: UnsafePointer<SIMD3<T>>
  ) -> (MM4FloatVector, MM4FloatVector, MM4FloatVector) {
    var x: MM4FloatVector = .zero
    var y: MM4FloatVector = .zero
    var z: MM4FloatVector = .zero
    
    if vID == atoms.vectorCount &- 1 {
      let remaining = atoms.count &- vID &* MM4VectorWidth
      for lane in 0..<remaining {
        let vector3 = baseAddress[vID &* MM4VectorWidth &+ lane]
        x[lane] = Float(vector3.x)
        y[lane] = Float(vector3.y)
        z[lane] = Float(vector3.z)
      }
    } else {
      for lane in 0..<MM4VectorWidth {
        let vector3 = baseAddress[vID &* MM4VectorWidth &+ lane]
        x[lane] = Float(vector3.x)
        y[lane] = Float(vector3.y)
        z[lane] = Float(vector3.z)
      }
    }
    return (x, y, z)
  }
  
  // Sum `chunk` number of scalars in FP32. Since the inner loop accumulates
  // every lane of the vector separately, we do not divide 'chunk' by the
  // vector width. We assume the user casts them to FP64 before reducing.
  @_transparent
  func withSegmentedLoop(chunk: Int, _ closure: (Range<Int>) -> Void) {
    var loopEnd = 0
    while loopEnd < atoms.vectorCount {
      let loopStart = loopEnd
      loopEnd = min(loopEnd &+ chunk, atoms.vectorCount)
      closure(loopStart..<loopEnd)
    }
  }
}
