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
  var vPositions: [MM4FloatVector] = [] // relative to reference frame
  var vVelocities: [MM4FloatVector] = [] // relative to reference frame
  
  // Reference frame; updates are tracked, but never recomputed from scratch.
  var angularMomentum: SIMD3<Double> = .zero
  var centerOfMass: SIMD3<Double> = .zero
  var linearMomentum: SIMD3<Double> = .zero
  var principalAxes: (
    SIMD3<Double>, SIMD3<Double>, SIMD3<Double>
  ) = (.zero, .zero, .zero)
  
  // Invalidate the forces whenever dependent properties change.
  var forces: [SIMD3<Float>]?
  var netForce: SIMD3<Double>?
  var netTorque: SIMD3<Double>?
  
  // Allow fine-grained, O(n) access to publicly presented positions and
  // velocities.
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
    createVectorizedMasses(parameters.atoms.masses)
    createVectorizedPositions(positions)
    createVectorizedVelocities(descriptor.velocities)
    
    // Linear Position
    mass = createMass()
    centerOfMass = createCenterOfMass()
    normalizeLinearPositions(to: centerOfMass)
    
    // Linear Momentum
    linearMomentum = createLinearMomentum()
    normalizeLinearVelocities(to: linearMomentum / mass)
    
    // Angular Position
    let inertiaTensor = createInertiaTensor()
    (momentOfInertia, principalAxes) = diagonalize(matrix: inertiaTensor)
    normalizeOrientation(to: principalAxes)
    
    // Angular Momentum
    angularMomentum = createAngularMomentum()
    normalizeAngularVelocities(to: angularMomentum / momentOfInertia)
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
  
  func ensureForceAndTorqueCached() {
    // compute force and torque in the same function
  }
  
  func invalidatePositions() {
    
  }
  
  func invalidateVelocities() {
    
  }
  
  func invalidateForces() {
    
  }
}
