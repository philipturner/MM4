//
//  MM4RigidBody+Affine.swift
//
//
//  Created by Philip Turner on 11/23/23.
//

import Numerics

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

// MARK: - Position

extension MM4RigidBody {
  /// Symmetric matrix specifying the rigid body's moment of inertia.
  ///
  /// If there is more than one anchor, this is the inertia of non-anchor atoms
  /// around the center of mass defined by anchors.
  public var momentOfInertia: (SIMD3<Float>, SIMD3<Float>, SIMD3<Float>) {
    // no setter; instead use rotate()
    get { fatalError("Not implemented.") }
  }
  
  /// Center of mass, treating anchors as astronomically larger than
  /// non-anchors.
  ///
  /// If there are any anchors, this is the mass-weighted average of the
  /// anchors.
  public var centerOfMass: SIMD3<Float> {
    _read { fatalError("Not implemented.") }
    _modify { fatalError("Not implemented.") }
  }
  
  /// Change the object's orientation by the specified 3D angle.
  public mutating func rotate(_ angle: Quaternion<Float>) {
    fatalError("Not implemented.")
  }
}

// MARK: - Velocity

// The setters for velocity have similar functionality to setters for position.

extension MM4RigidBody {
  // TODO: Add a one-line summary to each property's documentation.
  
  /// If the angular velocity is nonzero, the number of anchors cannot exceed 1.
  /// When importing velocities, if the number of anchors exceeds 1, the angular
  /// velocity is set to zero.
  public var angularVelocity: Quaternion<Float> {
    get { fatalError("Not implemented.") }
    set { fatalError("Not implemented.") }
  }
  
  /// The force is distributed evenly among all non-anchor atoms in the rigid
  /// body. If an anchor is explicitly selected as an external force target,
  /// there will be an error.
  public var externalForce: SIMD3<Float> {
    get { fatalError("Not implemented.") }
    set { fatalError("Not implemented.") }
  }
  
  /// Every anchor's velocity is set to the rigid body's linear velocity.
  /// When importing velocities, all anchors must have the same velocity.
  public var linearVelocity: SIMD3<Float> {
    get { fatalError("Not implemented.") }
    set { fatalError("Not implemented.") }
  }
}
