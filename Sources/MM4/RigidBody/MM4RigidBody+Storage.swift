//
//  MM4RigidBody+Storage.swift
//
//
//  Created by Philip Turner on 11/22/23.
//

import Numerics

// Sources of truth: anchors, mass, vPositions, vVelocities
// Positions may be decomposed repeatedly during operations like rotation
// Velocities may be decomposed repeatedly during operations like thermalization
final class MM4RigidBodyStorage {
  var anchors: Set<UInt32>
  
  var angularMass: MM4AngularMass?
  
  var angularVelocity: Quaternion<Float>?
  
  var centerOfMass: SIMD3<Double>?
  
  var thermalKineticEnergy: Double?
  
  var linearVelocity: SIMD3<Double>?
  
  var mass: Double
  
  var positions: [SIMD3<Float>]?
  
  var velocities: [SIMD3<Float>]?
  
  var vPositions: [MM4FloatVector]
  
  var vVelocities: [MM4FloatVector]
  
  init(atoms: MM4RigidBodyAtoms, parameters: MM4Parameters) {
    // Initialize stored properties.
    anchors = []
    mass = parameters.atoms.masses.reduce(Double(0)) {
      $0 + Double($1)
    }
    vPositions = Array(unsafeUninitializedCapacity: 3 * atoms.vectorCount) {
      $0.initialize(repeating: .zero)
      $1 = 3 * atoms.count
    }
    vVelocities = Array(unsafeUninitializedCapacity: 3 * atoms.vectorCount) {
      $0.initialize(repeating: .zero)
      $1 = 3 * atoms.count
    }
  }
  
  init(copying other: MM4RigidBodyStorage) {
    // Initialize stored properties, without copying cached ones.
    anchors = other.anchors
    mass = other.mass
    vPositions = other.vPositions
    vVelocities = other.vVelocities
  }
}

extension MM4RigidBody {
  /// Ensure all weak references point to the current storage object.
  mutating func ensureReferencesUpdated() {
    energy.kinetic.storage = storage
  }
  
  /// Ensure copy-on-write semantics.
  ///
  /// > WARNING: Call this before every mutating function.
  mutating func ensureUniquelyReferenced() {
    if !isKnownUniquelyReferenced(&storage) {
      storage = MM4RigidBodyStorage(copying: storage)
      ensureReferencesUpdated()
    }
  }
}
