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
  var atoms: (count: Int, vectorCount: Int)
  var handles: Set<UInt32>
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
  var angularKineticEnergy: Double?
  var angularVelocity: Quaternion<Float>?
  var linearKineticEnergy: Double?
  var linearVelocity: SIMD3<Float>?
  var momentOfInertia: MM4MomentOfInertia?
  var thermalKineticEnergy: Double?
  
  init(
    anchors: Set<UInt32>,
    handles: Set<UInt32>,
    parameters: MM4Parameters
  ) {
    // Initialize stored properties.
    let atomCount = parameters.atoms.count
    let atomVectorCount = (atomCount + MM4VectorWidth - 1) / MM4VectorWidth
    self.anchors = anchors
    self.atoms = (atomCount, atomVectorCount)
    self.handles = handles
    self.vPositions = Array(unsafeUninitializedCapacity: 3 * atomVectorCount) {
      $0.initialize(repeating: .zero)
      $1 = 3 * atomCount
    }
    self.vVelocities = Array(unsafeUninitializedCapacity: 3 * atomVectorCount) {
      $0.initialize(repeating: .zero)
      $1 = 3 * atomCount
    }
    
    self.anchorMass = .zero
    self.anchorMasses = []
    self.nonAnchorMasses = parameters.atoms.masses
    if anchors.count > .zero {
      func generateArray() -> [Float] {
        let arraySize = atomVectorCount * MM4VectorWidth
        return Array(unsafeUninitializedCapacity: arraySize) {
          $0.initialize(repeating: .zero)
          $1 = atomCount
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
    handles = other.handles
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
    angularKineticEnergy = nil
    angularVelocity = nil
    linearKineticEnergy = nil
    linearVelocity = nil
    momentOfInertia = nil
    thermalKineticEnergy = nil
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
    if linearVelocity == nil {
      ensureAnchorVelocitiesValid()
      if anchors.count > 0 {
        let anchor = anchors.first!
        let velocity = extractScalar(Int(anchor), vVelocities)
        linearVelocity = velocity
      } else if atoms.count > 0 {
        linearVelocity = createLinearVelocity()
      } else {
        linearVelocity = .zero
      }
    }
  }
  
  func ensureAngularVelocityCached() {
    if angularVelocity == nil {
      ensureAnchorVelocitiesValid()
      if anchors.count <= 1, atoms.count > 0 {
        angularVelocity = createAngularVelocity()
      } else {
        angularVelocity = .zero
      }
    }
  }
  
  func ensureKineticEnergyCached() {
    if self.angularKineticEnergy == nil,
       self.linearKineticEnergy == nil,
       self.thermalKineticEnergy == nil {
      let angular = createAngularKineticEnergy()
      let linear = createLinearKineticEnergy()
      let total = createTotalKineticEnergy()
      let thermal = total - angular - linear
      self.angularKineticEnergy = angular
      self.linearKineticEnergy = linear
      self.thermalKineticEnergy = thermal
      
      let epsilon: Double = 1e-4
      if angular < epsilon * Double(atoms.count) {
        fatalError("Angular kinetic energy was too negative.")
      }
      if linear < epsilon * Double(atoms.count) {
        fatalError("Linear kinetic energy was too negative.")
      }
      if thermal < epsilon * Double(atoms.count) {
        fatalError("Thermal kinetic energy was too negative.")
      }
    } else if self.angularVelocity == nil ||
                self.linearKineticEnergy == nil ||
                self.thermalKineticEnergy == nil {
      fatalError(
        "Either all or none of the kinetic energies must be cached.")
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


