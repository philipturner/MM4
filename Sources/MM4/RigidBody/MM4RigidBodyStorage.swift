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
  
  init(descriptor: MM4RigidBodyDescriptor) throws {
    let masses = descriptor.masses!
    let atomCount = masses.count
    let vectorCount = (atomCount + MM4VectorWidth - 1) / MM4VectorWidth
    self.atoms = (atomCount, vectorCount)
    
    guard let positions = descriptor.positions else {
      fatalError("Positions were not specified.")
    }
    createVectorizedMasses(masses)
    createVectorizedPositions(positions)
    createVectorizedVelocities(descriptor.velocities)
    
    // Linear Position
    mass = createMass()
    guard mass > .leastNormalMagnitude else {
      throw MM4Error.defectiveInertiaTensor(
        (.zero, .zero, .zero), "Zero mass.")
    }
    centerOfMass = createCenterOfMass()
    normalizeLinearPositions(to: centerOfMass)
    
    // Linear Momentum
    linearMomentum = createLinearMomentum()
    normalizeLinearVelocities(to: linearMomentum / mass)
    
    // Angular Position
    let inertiaTensor = createInertiaTensor()
    let (Λ, Σ, failureReason) = diagonalize(matrix: inertiaTensor)
    
    guard let Λ, let Σ else {
      throw MM4Error.defectiveInertiaTensor(
        inertiaTensor, failureReason!)
    }
    (momentOfInertia, principalAxes) = (Λ, Σ)
    normalizeOrientation(to: principalAxes)
    
    // Angular Momentum
    angularMomentum = createAngularMomentum()
    normalizeAngularVelocities(to: angularMomentum / momentOfInertia)
  }
  
  init(copying other: MM4RigidBodyStorage) {
    // Copy only the sources of truth. Cached properties may be dirty. The
    // public API still shows the same semantics. If forces aren't 'nil', the
    // net force and torque will regenerate.
    atoms = other.atoms
    mass = other.mass
    momentOfInertia = other.momentOfInertia
    vMasses = other.vMasses
    vPositions = other.vPositions
    vVelocities = other.vVelocities
    
    angularMomentum = other.angularMomentum
    centerOfMass = other.centerOfMass
    linearMomentum = other.linearMomentum
    principalAxes = other.principalAxes
    
    forces = other.forces
  }
}

extension MM4RigidBodyStorage {
  func ensurePositionsCached() {
    guard positions == nil else {
      return
    }
    positions = Array(unsafeUninitializedCapacity: atoms.count) {
      $1 = atoms.count
      createPositions($0)
    }
  }
  
  func ensureVelocitiesCached() {
    guard velocities == nil else {
      return
    }
    velocities = Array(unsafeUninitializedCapacity: atoms.count) {
      $1 = atoms.count
      createVelocities($0)
    }
  }
  
  func ensureForceAndTorqueCached() {
    guard netForce == nil,
          netTorque == nil else {
      return
    }
    if let forces {
      forces.withUnsafeBufferPointer(setNetForceAndTorque(_:))
    }
  }
  
  func invalidatePositions() {
    positions = nil
    invalidateForces()
  }
  
  func invalidateVelocities() {
    velocities = nil
  }
  
  func invalidateForces() {
    forces = nil
    netForce = nil
    netTorque = nil
  }
}
