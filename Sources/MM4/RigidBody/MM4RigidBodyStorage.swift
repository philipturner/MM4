//
//  MM4RigidBody+Update.swift
//  MM4
//
//  Created by Philip Turner on 11/22/23.
//

final class MM4RigidBodyStorage {
  // Reference frame; never mutated after initialization.
  var atoms: (count: Int, vectorCount: Int)
  var mass: Double = .zero
  var momentOfInertia: SIMD3<Double> = .zero
  var vMasses: [MM4FloatVector] = []
  var vPositions: [MM4FloatVector] = [] // relative positions in reference frame
  var vVelocities: [MM4FloatVector] = []  // thermal velocities in reference frame
  
  // Reference frame; updates are tracked over time, but never recomputed from
  // scratch after initialization.
  var angularMomentum: SIMD3<Double> = .zero // in reference frame
  var centerOfMass: SIMD3<Double> = .zero
  var linearMomentum: SIMD3<Double> = .zero
  var principalAxes: (
    SIMD3<Double>, SIMD3<Double>, SIMD3<Double>
  ) = (.zero, .zero, .zero)
  
  // Invalidate the accelerations and forces whenever dependent properties
  // change. Lazily rematerialize them, if the dependent properties exist.
  var angularAcceleration: SIMD3<Double>? // in the global reference frame
  var forces: [SIMD3<Float>]?
  var linearAcceleration: SIMD3<Double>? // in the local reference frame
  
  // Allow fine-grained, O(n) access to publicly presented positions and
  // velocities. You can invalidate these simply by deleting the array (setting
  // it to `nil`).
  var positions: [SIMD3<Float>]?
  var velocities: [SIMD3<Float>]?
  
  init(descriptor: MM4RigidBodyDescriptor) {
    let parameters = descriptor.parameters!
    let atomCount = parameters.atoms.count
    let vectorCount = (atomCount + MM4VectorWidth - 1) / MM4VectorWidth
    self.atoms = (atomCount, vectorCount)
    
    guard let positions = descriptor.positions else {
      fatalError("Positions were not specified.")
    }
    
    // Swizzle the arrays into a format optimized for vector ALUs.
    createVectorizedMasses(parameters.atoms.masses)
    createVectorizedPositions(positions)
    createVectorizedVelocities(descriptor.velocities)
    
    // Linear Position
    mass = createMass()
    centerOfMass = createCenterOfMass(mass: mass)
    normalizeLinearPositions(centerOfMass: centerOfMass)
    
    // Angular Position
    let inertiaTensor = createInertiaTensor()
    (momentOfInertia, principalAxes) = diagonalize(matrix: inertiaTensor)
    normalizeOrientation(principalAxes: principalAxes)
    
    // Linear Momentum
    linearMomentum = createLinearMomentum()
    normalizeLinearVelocities(mass: mass, linearMomentum: linearMomentum)
    
    // Angular Momentum
    angularMomentum = createAngularMomentum()
    normalizeAngularVelocities(
      momentOfInertia: momentOfInertia,
      angularMomentum: angularMomentum)
    
    // Shift the momenta into the global reference frame.
    linearMomentum = gemv(matrix: principalAxes, vector: linearMomentum)
    angularMomentum = gemv(matrix: principalAxes, vector: angularMomentum)
  }
  
  init(copying other: MM4RigidBodyStorage) {
    // Copy over the stored and cached properties.
  }
}

extension MM4RigidBodyStorage {
  func ensurePositionsCached() {
    // compute the position by applying bulk CoM and MoI in the same function
  }
  
  func ensureVelocitiesCached() {
    // compute the velocity by applying bulk P and L in the same function
  }
  
  func ensureAccelerationsCached() {
    // compute both accelerations in the same function
  }
  
  func invalidatePositions() {
    
  }
  
  func invalidateVelocities() {
    
  }
  
  func invalidateForces() {
    
  }
}
