//
//  MM4RigidBody+Affine.swift
//
//
//  Created by Philip Turner on 11/23/23.
//

import Numerics

/// Moment of inertia.
public struct MM4AngularMass {
  /// Symmetric matrix specifying the rigid body's moment of inertia.
  public var columns: (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>)
  
  /// Initialize a moment of inertia with zero mass.
  public init() {
    self.columns = (.zero, .zero, .zero)
  }
  
  /// The matrix is symmetric, but not exactly orthonormal. The inverse is not
  /// the same as the transpose.
  public var inverse: (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>) {
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

extension MM4RigidBody {
  /// If there is more than one anchor, the angular mass is zero.
  public var angularMass: MM4AngularMass {
    // no setter; instead use rotate()
    get { fatalError("Not implemented.") }
  }
  
  public var centerOfMass: SIMD3<Float> {
    _read { fatalError("Not implemented.") }
    _modify { fatalError("Not implemented.") }
  }
  
  public mutating func rotate(_ angle: Quaternion<Float>) {
    fatalError("Not implemented.")
  }
}
