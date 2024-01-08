//
//  MM4RigidBodyStorage.swift
//  MM4
//
//  Created by Philip Turner on 11/22/23.
//

final class MM4RigidBodyStorage {
  // Reference frame; never mutated after initialization.
  var atoms: (count: Int, vectorCount: Int, nonAnchorCount: Int)
  var mass: Double = .zero
  var momentOfInertia: SIMD3<Double> = .zero
  var vMasses: [MM4FloatVector] = []
  var vPositions: [MM4FloatVector] = [] // relative positions in reference frame
  
  // Reference frame; updates are tracked over time, but never recomputed from
  // scratch after initialization.
  var centerOfMass: SIMD3<Double> = .zero
  var principalAxes: (
    SIMD3<Double>, SIMD3<Double>, SIMD3<Double>) = (.zero, .zero, .zero)
  
  // Bulk properties that are mutated instead of updating every atom's position
  // and velocity. When the atom velocities or forces are updated by the user,
  // some of these properties must be recomputed. The most common operation is
  // recomputing acceleration after force, center of mass, or rotation changes.
  //
  // You can track the caching of linear/angular acceleration separately, now
  // that forces are swizzled beforehand. Sort of like how center of mass and
  // moment of inertia are computed in two separate passes.
  var angularAcceleration: SIMD3<Double>?
  var angularMomentum: SIMD3<Double>?
  var linearAcceleration: SIMD3<Double>?
  var linearMomentum: SIMD3<Double>?
  var vForces: [MM4FloatVector]? // according to the outside observer
  var vVelocities: [MM4FloatVector]? // thermal velocities in reference frame
  
  // Allow fine-grained, O(n) access to positions, velocities, and forces.
  // Cache I/O to publicly presented positions, velocities, and forces.
  // Cache I/O to internal velocities and forces.
  // - When the object is being mutated, only one may be alive at any given
  //   moment: the public or internal copy.
  // - Otherwise, both may coexist. This is to minimize rounding error from
  //   constantly recomputing the other one.
  var cachedState: MM4State?
  
  init(parameters: MM4Parameters) {
    let atomCount = parameters.atoms.count
    let vectorCount = (atomCount + MM4VectorWidth - 1) / MM4VectorWidth
    self.atoms = (atomCount, vectorCount, .zero)
  }
  
  init(copying other: MM4RigidBodyStorage) {
    // Copy over the stored and cached properties.
  }
}
