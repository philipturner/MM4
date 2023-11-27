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
    get {
      storage.ensureMomentOfInertiaCached()
      let momentOfInertia = storage.momentOfInertia!
      return (
        SIMD3<Float>(momentOfInertia.columns.0),
        SIMD3<Float>(momentOfInertia.columns.1),
        SIMD3<Float>(momentOfInertia.columns.2)
      )
    }
  }
  
  /// Center of mass, treating anchors as astronomically larger than
  /// non-anchors.
  ///
  /// If there are any anchors, this is the mass-weighted average of the
  /// anchors.
  public var centerOfMass: SIMD3<Float> {
    _read {
      storage.ensureCenterOfMassCached()
      yield storage.centerOfMass!
    }
    _modify {
      ensureUniquelyReferenced()
      storage.ensureCenterOfMassCached()
      storage.positions = nil
      let previous = storage.centerOfMass!
      yield &storage.centerOfMass!
      
      let difference = storage.centerOfMass! - previous
      guard any(difference .!= .zero) else {
        return
      }
      for vID in 0..<storage.atoms.vectorCount {
        storage.vPositions[vID &* 3 &+ 0] += difference.x
        storage.vPositions[vID &* 3 &+ 1] += difference.y
        storage.vPositions[vID &* 3 &+ 2] += difference.z
      }
    }
  }
  
  /// Change the object's orientation by the specified quaternion.
  public mutating func rotate(_ rotation: Quaternion<Float>) {
    ensureUniquelyReferenced()
    storage.ensureCenterOfMassCached()
    storage.positions = nil
    
    let centerOfMass = storage.centerOfMass!
    let w = rotation.angle * rotation.axis
    for vID in 0..<storage.atoms.vectorCount {
      var x = storage.vPositions[vID &* 3 &+ 0]
      var y = storage.vPositions[vID &* 3 &+ 1]
      var z = storage.vPositions[vID &* 3 &+ 2]
      let rX = x - centerOfMass.x
      let rY = y - centerOfMass.y
      let rZ = z - centerOfMass.z
      
      x += w.y * rZ - w.z * rY
      y += w.z * rX - w.x * rZ
      z += w.x * rY - w.y * rX
      storage.vPositions[vID &* 3 &+ 0] = x
      storage.vPositions[vID &* 3 &+ 1] = y
      storage.vPositions[vID &* 3 &+ 2] = z
    }
  }
}

// MARK: - Velocity

// The setters for velocity have similar functionality to setters for position.

extension MM4RigidBody {
  /// The net angular velocity (in radians per picosecond) of the non-anchor
  /// atoms.
  public var angularVelocity: Quaternion<Float> {
    _read {
      storage.ensureAngularVelocityCached()
      yield storage.angularVelocity!
    }
    _modify {
      ensureUniquelyReferenced()
      storage.ensureCenterOfMassCached()
      storage.ensureAngularVelocityCached()
      storage.velocities = nil
      let previous = storage.angularVelocity!
      yield &storage.angularVelocity!
      
      let next = storage.angularVelocity!
      guard next != previous else {
        return
      }
      guard let centerOfMass = storage.centerOfMass else {
        fatalError("This should never happen.")
      }
      
      let previousW = previous.angle * previous.axis
      let nextW = next.angle * next.axis
      storage.withMasses(storage.nonAnchorMasses) { vMasses in
        for vID in 0..<storage.atoms.vectorCount {
          let rX = storage.vPositions[vID &* 3 &+ 0] - centerOfMass.x
          let rY = storage.vPositions[vID &* 3 &+ 1] - centerOfMass.y
          let rZ = storage.vPositions[vID &* 3 &+ 2] - centerOfMass.z
          var vX: MM4FloatVector = .zero
          var vY: MM4FloatVector = .zero
          var vZ: MM4FloatVector = .zero
          
          // Un-apply the previous bulk angular velocity.
          do {
            let w = previousW
            vX -= w.y * rZ - w.z * rY
            vY -= w.z * rX - w.x * rZ
            vZ -= w.x * rY - w.y * rX
          }
          
          // Apply the next bulk angular velocity.
          do {
            let w = nextW
            vX += w.y * rZ - w.z * rY
            vY += w.z * rX - w.x * rZ
            vZ += w.x * rY - w.y * rX
          }
          
          // Mask out the changes to velocity for anchors.
          let mass = vMasses[vID]
          vX.replace(with: MM4FloatVector.zero, where: mass .== 0)
          vY.replace(with: MM4FloatVector.zero, where: mass .== 0)
          vZ.replace(with: MM4FloatVector.zero, where: mass .== 0)
          
          // Write the new velocities as offsets relative to the existing ones.
          storage.vVelocities[vID &* 3 &+ 0] += vX
          storage.vVelocities[vID &* 3 &+ 1] += vY
          storage.vVelocities[vID &* 3 &+ 2] += vZ
        }
      }
    }
  }
  
  /// The constant force (in piconewtons) exerted on the entire object.
  ///
  /// The force is distributed evenly among all non-anchor atoms in the rigid
  /// body. If an anchor is explicitly selected as a handle, there will be a
  /// fatal error.
  public var externalForce: SIMD3<Float> {
    _read {
      yield storage.externalForce
    }
    _modify {
      ensureUniquelyReferenced()
      yield &storage.externalForce
    }
  }
  
  /// The net linear velocity (in nanometers per picosecond) of the entire
  /// object.
  ///
  /// Every anchor's velocity is set to the rigid body's linear velocity.
  /// When importing velocities, all anchors must have the same velocity.
  public var linearVelocity: SIMD3<Float> {
    _read {
      storage.ensureLinearVelocityCached()
      yield storage.linearVelocity!
    }
    _modify {
      ensureUniquelyReferenced()
      storage.ensureLinearVelocityCached()
      storage.velocities = nil
      let previous = storage.linearVelocity!
      yield &storage.linearVelocity!
      
      let difference = storage.linearVelocity! - previous
      guard any(difference .!= .zero) else {
        return
      }
      for vID in 0..<storage.atoms.vectorCount {
        storage.vVelocities[vID &* 3 &+ 0] += difference.x
        storage.vVelocities[vID &* 3 &+ 1] += difference.y
        storage.vVelocities[vID &* 3 &+ 2] += difference.z
      }
    }
  }
}
