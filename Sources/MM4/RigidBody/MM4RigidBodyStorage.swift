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
  var vVelocities: [MM4FloatVector] = []  // thermal velocities in reference frame
  
  // Reference frame; updates are tracked over time, but never recomputed from
  // scratch after initialization.
  var angularMomentum: SIMD3<Double> = .zero
  var centerOfMass: SIMD3<Double> = .zero
  var linearMomentum: SIMD3<Double> = .zero
  var principalAxes: (
    SIMD3<Double>, SIMD3<Double>, SIMD3<Double>) = (.zero, .zero, .zero)
  
  // Invalidate the accelerations whenever the center of mass, principal axes,
  // or forces change.
  var angularAcceleration: SIMD3<Double>?
  var forces: [SIMD3<Float>] = []
  var linearAcceleration: SIMD3<Double>?
  
  // Allow fine-grained, O(n) access to publicly presented positions and
  // velocities. You can invalidate these simply by deleting the array (setting
  // it to `nil`).
  var positions: [SIMD3<Float>]?
  var velocities: [SIMD3<Float>]?
  
  init(parameters: MM4Parameters) {
    let atomCount = parameters.atoms.count
    let vectorCount = (atomCount + MM4VectorWidth - 1) / MM4VectorWidth
    self.atoms = (atomCount, vectorCount, .zero)
  }
  
  init(copying other: MM4RigidBodyStorage) {
    // Copy over the stored and cached properties.
  }
}

extension MM4RigidBodyStorage {
  func ensurePositionsCached() {
    
  }
  
  func ensureVelocitiesCached() {
    
  }
}
