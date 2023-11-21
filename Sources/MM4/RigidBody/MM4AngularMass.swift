//
//  MM4AngularMass.swift
//
//
//  Created by Philip Turner on 11/20/23.
//

/// Moment of inertia.
///
/// Source: [Stack Overflow](https://stackoverflow.com/a/18504573)
public struct MM4AngularMass {
  /// Symmetric matrix specifying the rigid body's moment of inertia.
  public var columns: (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>)
  
  /// Initialize a moment of inertia with zero mass.
  public init() {
    self.columns = (.zero, .zero, .zero)
  }
  
  /// Add an atom to the accumulator.
  mutating func append(mass: Double, relativePosition: SIMD3<Float>) {
    // From Wikipedia:
    // https://en.wikipedia.org/wiki/Rigid_body_dynamics#Mass_properties
    //
    // I_R = m * (I (S^T S) - S S^T)
    // where S is the column vector R - R_cm
    let positionSquared = relativePosition * relativePosition
    let STS = positionSquared[0] + positionSquared[1] + positionSquared[2]
    var column0 = SIMD3(STS, 0, 0)
    var column1 = SIMD3(0, STS, 0)
    var column2 = SIMD3(0, 0, STS)
    column0 -= relativePosition.x * relativePosition
    column1 -= relativePosition.y * relativePosition
    column2 -= relativePosition.z * relativePosition
    
    // Convert to FP64 before adding to the accumulator. The matrix is
    // symmetric, so it doesn't matter whether you mix up the rows and columns.
    columns.0 += mass * SIMD3<Double>(column0)
    columns.1 += mass * SIMD3<Double>(column1)
    columns.2 += mass * SIMD3<Double>(column2)
  }
  
  /// The matrix is symmetric, but not exactly orthonormal. Inversion is not a
  /// simple transpose operation.
  public var inverse: (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>) {
    // Source: https://stackoverflow.com/a/18504573
    //
    // double det = m(0, 0) * (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) -
    //              m(0, 1) * (m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0)) +
    //              m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));
    let determinant =
    columns.0[0]*(columns.1[1]*columns.2[2]-columns.2[1]*columns.1[2]) -
    columns.0[1]*(columns.1[0]*columns.2[2]-columns.1[2]*columns.2[0]) +
    columns.0[2]*(columns.1[0]*columns.2[1]-columns.1[1]*columns.2[0])
    let invdet = 1 / determinant
    
    // minv(0, 0) = (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) * invdet;
    // minv(0, 1) = (m(0, 2) * m(2, 1) - m(0, 1) * m(2, 2)) * invdet;
    // minv(0, 2) = (m(0, 1) * m(1, 2) - m(0, 2) * m(1, 1)) * invdet;
    let result00 = (columns.1[1]*columns.2[2]-columns.2[1]*columns.1[2])*invdet
    let result01 = (columns.0[2]*columns.2[1]-columns.0[1]*columns.2[2])*invdet
    let result02 = (columns.0[1]*columns.1[2]-columns.0[2]*columns.1[1])*invdet
    
    // minv(1, 0) = (m(1, 2) * m(2, 0) - m(1, 0) * m(2, 2)) * invdet;
    // minv(1, 1) = (m(0, 0) * m(2, 2) - m(0, 2) * m(2, 0)) * invdet;
    // minv(1, 2) = (m(1, 0) * m(0, 2) - m(0, 0) * m(1, 2)) * invdet;
    let result10 = (columns.1[2]*columns.2[0]-columns.1[0]*columns.2[2])*invdet
    let result11 = (columns.0[0]*columns.2[2]-columns.0[2]*columns.2[0])*invdet
    let result12 = (columns.1[0]*columns.0[2]-columns.0[0]*columns.1[2])*invdet
    
    // minv(2, 0) = (m(1, 0) * m(2, 1) - m(2, 0) * m(1, 1)) * invdet;
    // minv(2, 1) = (m(2, 0) * m(0, 1) - m(0, 0) * m(2, 1)) * invdet;
    // minv(2, 2) = (m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1)) * invdet;
    let result20 = (columns.1[0]*columns.2[1]-columns.2[0]*columns.1[1])*invdet
    let result21 = (columns.2[0]*columns.0[1]-columns.0[0]*columns.2[1])*invdet
    let result22 = (columns.0[0]*columns.1[1]-columns.1[0]*columns.0[1])*invdet
    
    let column0 = SIMD3(result00, result10, result20)
    let column1 = SIMD3(result01, result11, result21)
    let column2 = SIMD3(result02, result12, result22)
    return (column0, column1, column2)
  }
}
