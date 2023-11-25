//
//  MM4RigidBodyStorage.swift
//
//
//  Created by Philip Turner on 11/22/23.
//

import Numerics

final class MM4RigidBodyStorage {
  // Sources of truth.
  var anchors: Set<UInt32>
  var atoms: MM4RigidBodyAtoms
  var vPositions: [MM4FloatVector]
  var vVelocities: [MM4FloatVector]
  
  // Cached constants.
  var anchorMass: Double
  var anchorMasses: [Float]
  var nonAnchorMass: Double
  var nonAnchorMasses: [Float]
  
  // Frequently cached (special handling).
  var anchorVelocitiesValid: Bool?
  var centerOfMass: SIMD3<Float>?
  var positions: [SIMD3<Float>]?
  var velocities: [SIMD3<Float>]?
  
  // Rarely cached (frequently erased).
  var angularVelocity: Quaternion<Float>?
  var freeKineticEnergy: Double?
  var thermalKineticEnergy: Double?
  var linearVelocity: SIMD3<Float>?
  var momentOfInertia: MM4MomentOfInertia?
  
  init(
    anchors: Set<UInt32>,
    parameters: MM4Parameters
  ) {
    // Initialize stored properties.
    let atoms = MM4RigidBodyAtoms(count: parameters.atoms.count)
    self.anchors = anchors
    self.atoms = atoms
    vPositions = Array(unsafeUninitializedCapacity: 3 * atoms.vectorCount) {
      $0.initialize(repeating: .zero)
      $1 = 3 * atoms.count
    }
    vVelocities = Array(unsafeUninitializedCapacity: 3 * atoms.vectorCount) {
      $0.initialize(repeating: .zero)
      $1 = 3 * atoms.count
    }
    
    self.anchorMass = .zero
    self.anchorMasses = []
    self.nonAnchorMasses = parameters.atoms.masses
    if anchors.count > .zero {
      func generateArray() -> [Float] {
        let arraySize = atoms.vectorCount * MM4VectorWidth
        return Array(unsafeUninitializedCapacity: arraySize) {
          $0.initialize(repeating: .zero)
          $1 = atoms.count
        }
      }
      anchorMasses = generateArray()
      nonAnchorMasses = generateArray()
      
      for atomID in 0..<atoms.count {
        nonAnchorMasses[atomID] = parameters.atoms.masses[atomID]
      }
      for anchor in self.anchors {
        let mass = parameters.atoms.masses[Int(anchor)]
        anchorMass += Double(mass)
        anchorMasses[Int(anchor)] = mass
        nonAnchorMasses[Int(anchor)] = .zero
      }
      
    }
    self.nonAnchorMass = parameters.atoms.masses.reduce(Double(0)) {
      $0 + Double($1)
    }
    self.nonAnchorMass -= Double(anchorMass)
  }
  
  init(copying other: MM4RigidBodyStorage) {
    // Initialize stored properties, without copying cached ones.
    atoms = other.atoms
    anchors = other.anchors
    anchorMass = other.anchorMass
    anchorMasses = other.anchorMasses
    nonAnchorMass = other.nonAnchorMass
    nonAnchorMasses = other.nonAnchorMasses
    vPositions = other.vPositions
    vVelocities = other.vVelocities
  }
  
  func eraseFrequentlyCachedProperties() {
    // WARNING: Always ensure this reflects recently added properties.
    anchorVelocitiesValid = nil
    centerOfMass = nil
    positions = nil
    velocities = nil
  }
  
  func eraseRarelyCachedProperties() {
    // WARNING: Always ensure this reflects recently added properties.
    angularVelocity = nil
    freeKineticEnergy = nil
    thermalKineticEnergy = nil
    linearVelocity = nil
    momentOfInertia = nil
  }
  
  func ensureCenterOfMassCached() {
    if self.centerOfMass == nil {
      self.centerOfMass = createCenterOfMass()
    }
  }
  
  func ensureMomentOfInertiaCached() {
    if self.momentOfInertia == nil {
      self.momentOfInertia = createMomentOfInertia()
    }
  }
  
  func ensureLinearVelocityCached() {
    if self.linearVelocity == nil {
      self.linearVelocity = createLinearVelocity()
    }
  }
  
  func ensureAngularVelocityCached() {
    if self.angularVelocity == nil {
      self.angularVelocity = createAngularVelocity()
    }
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

