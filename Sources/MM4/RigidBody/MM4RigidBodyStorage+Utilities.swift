//
//  MM4RigidBodyStorage+Utilities.swift
//
//
//  Created by Philip Turner on 11/25/23.
//

import Numerics

func quaternionToVector(_ quaternion: Quaternion<Float>) -> SIMD3<Float> {
  let angleAxis = quaternion.angleAxis
  let axis = angleAxis.axis
  
  if axis[0].isNaN || axis[1].isNaN || axis[2].isNaN {
    return .zero
  } else if angleAxis.length == 0 || angleAxis.angle.isNaN {
    return .zero
  } else {
    return angleAxis.angle * angleAxis.axis
  }
}

func vQuaternionPrepare(
  _ quaternion: Quaternion<Float>
) -> (real: Float, imaginary: SIMD3<Float>, p1: Float) {
  let real = quaternion.real
  let imaginary = quaternion.imaginary
  let p1 = real * real - (imaginary * imaginary).sum()
  return (real, imaginary, p1)
}

// Convert the quaternion to real-imaginary format beforehand to minimize
// overhead.
@inline(__always)
func vQuaternionAct(
  _ quaternion: (real: Float, imaginary: SIMD3<Float>, p1: Float),
  _ x: inout MM4FloatVector,
  _ y: inout MM4FloatVector,
  _ z: inout MM4FloatVector
) {
  let real = quaternion.real
  let imag = quaternion.imaginary
  
  let p1x = x * quaternion.p1
  let p1y = y * quaternion.p1
  let p1z = z * quaternion.p1
  
  let dot = imag.x * x + imag.y * y + imag.z * z
  let p2x = imag.x * dot
  let p2y = imag.y * dot
  let p2z = imag.z * dot
  
  // WARNING: `p3` is not multiplied by `real` yet.
  let p3x = imag.y * z - imag.z * y
  let p3y = imag.z * x - imag.x * z
  let p3z = imag.x * y - imag.y * x
  
  x = p1x.addingProduct(2, p2x.addingProduct(p3x, real))
  y = p1y.addingProduct(2, p2y.addingProduct(p3y, real))
  z = p1z.addingProduct(2, p2z.addingProduct(p3z, real))
}

extension MM4RigidBodyStorage {
  @inline(__always)
  func extractScalar(
    _ scalarID: Int, _ array: [MM4FloatVector]
  ) -> SIMD3<Float> {
    let vID = scalarID / MM4VectorWidth
    let lane = scalarID &- vID &* MM4VectorWidth
    let x = array[vID &* 3 &+ 0][lane]
    let y = array[vID &* 3 &+ 1][lane]
    let z = array[vID &* 3 &+ 2][lane]
    return SIMD3(x, y, z)
  }
  
  @inline(__always)
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
  
  @inline(__always)
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
  
  @inline(__always)
  func withMasses<T>(
    _ masses: [Float],
    _ closure: (UnsafePointer<MM4FloatVector>) -> T
  ) -> T {
    return masses.withUnsafeBufferPointer {
      let rawMasses = OpaquePointer($0.baseAddress)
      let vMasses = UnsafeRawPointer(rawMasses)!
        .assumingMemoryBound(to: MM4FloatVector.self)
      return closure(vMasses)
    }
  }
  
  @inline(__always)
  func withSegmentedLoop(chunk: Int, _ closure: (Range<Int>) -> Void) {
    var loopEnd = 0
    while loopEnd < atoms.vectorCount {
      let loopStart = loopEnd
      loopEnd = min(loopEnd &+ chunk, atoms.vectorCount)
      closure(loopStart..<loopEnd)
    }
  }
}

/// Moment of inertia.
///
/// This is hidden from the public API. One reason is that columns are stored in
/// double precision, while the public API should express most quantities in
/// single precision.
struct MM4MomentOfInertia {
  /// Symmetric matrix specifying the rigid body's moment of inertia.
  var columns: (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>)
  
  /// Initialize a moment of inertia with zero mass.
  init() {
    self.columns = (.zero, .zero, .zero)
  }
  
  /// The matrix is symmetric, but not exactly orthonormal. The inverse is not
  /// the same as the transpose.
  var inverse: (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>) {
    // Source: https://stackoverflow.com/a/18504573
    let col = columns
    let determinant =
    col.0[0] * (col.1[1] * col.2[2] - col.2[1] * col.1[2]) -
    col.0[1] * (col.1[0] * col.2[2] - col.1[2] * col.2[0]) +
    col.0[2] * (col.1[0] * col.2[1] - col.1[1] * col.2[0])
    let invdet = 1 / determinant
    
    let result00 = (col.1[1] * col.2[2] - col.2[1] * col.1[2]) * invdet
    let result01 = (col.0[2] * col.2[1] - col.0[1] * col.2[2]) * invdet
    let result02 = (col.0[1] * col.1[2] - col.0[2] * col.1[1]) * invdet
    
    let result10 = (col.1[2] * col.2[0] - col.1[0] * col.2[2]) * invdet
    let result11 = (col.0[0] * col.2[2] - col.0[2] * col.2[0]) * invdet
    let result12 = (col.1[0] * col.0[2] - col.0[0] * col.1[2]) * invdet
    
    let result20 = (col.1[0] * col.2[1] - col.2[0] * col.1[1]) * invdet
    let result21 = (col.2[0] * col.0[1] - col.0[0] * col.2[1]) * invdet
    let result22 = (col.0[0] * col.1[1] - col.1[0] * col.0[1]) * invdet
    
    let column0 = SIMD3(result00, result10, result20)
    let column1 = SIMD3(result01, result11, result21)
    let column2 = SIMD3(result02, result12, result22)
    return (column0, column1, column2)
  }
}
