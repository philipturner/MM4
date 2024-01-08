//
//  MM4RigidBodyStorage+Utilities.swift
//  MM4
//
//  Created by Philip Turner on 11/25/23.
//

// Source: https://stackoverflow.com/a/18504573
func invertMatrix3x3(
  _ col: (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>)
) -> (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>) {
  let determinant =
  col.0[0] * (col.1[1] * col.2[2] - col.2[1] * col.1[2]) -
  col.0[1] * (col.1[0] * col.2[2] - col.1[2] * col.2[0]) +
  col.0[2] * (col.1[0] * col.2[1] - col.1[1] * col.2[0])
  
  let result00 = col.1[1] * col.2[2] - col.2[1] * col.1[2]
  let result01 = col.0[2] * col.2[1] - col.0[1] * col.2[2]
  let result02 = col.0[1] * col.1[2] - col.0[2] * col.1[1]
  
  let result10 = col.1[2] * col.2[0] - col.1[0] * col.2[2]
  let result11 = col.0[0] * col.2[2] - col.0[2] * col.2[0]
  let result12 = col.1[0] * col.0[2] - col.0[0] * col.1[2]
  
  let result20 = col.1[0] * col.2[1] - col.2[0] * col.1[1]
  let result21 = col.2[0] * col.0[1] - col.0[0] * col.2[1]
  let result22 = col.0[0] * col.1[1] - col.1[0] * col.0[1]
  
  let column0 = SIMD3(result00, result01, result02) / determinant
  let column1 = SIMD3(result10, result11, result12) / determinant
  let column2 = SIMD3(result20, result21, result22) / determinant
  return (column0, column1, column2)
}

// diagonalizeMatrix3x3

// solveCubicPolynomial

extension MM4RigidBodyStorage {
//  @_transparent
//  func extractScalar(
//    _ scalarID: Int, _ array: [MM4FloatVector]
//  ) -> SIMD3<Float> {
//    let vID = scalarID / MM4VectorWidth
//    let lane = scalarID &- vID &* MM4VectorWidth
//    let x = array[vID &* 3 &+ 0][lane]
//    let y = array[vID &* 3 &+ 1][lane]
//    let z = array[vID &* 3 &+ 2][lane]
//    return SIMD3(x, y, z)
//  }
  
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
