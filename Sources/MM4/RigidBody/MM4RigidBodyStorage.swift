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
  var externalForce: SIMD3<Float>
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
    self.externalForce = .zero
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
    
    func generateArray() -> [Float] {
      let arraySize = atomVectorCount * MM4VectorWidth
      return Array(unsafeUninitializedCapacity: arraySize) {
        $0.initialize(repeating: .zero)
        $1 = atomCount
      }
    }
    
    // We can't be 100% sure the backing array for masses will still
    // be padded to the vector width. In addition, it could be modified if any
    // anchors exist. The copy may not have the same capacity as the previous
    // array.
    nonAnchorMasses = generateArray()
    for atomID in 0..<atoms.count {
      nonAnchorMasses[atomID] = parameters.atoms.masses[atomID]
    }
    
    if anchors.count > .zero {
      anchorMasses = generateArray()
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
    externalForce = other.externalForce
    handles = other.handles
    vPositions = other.vPositions
    vVelocities = other.vVelocities
    
    anchorMass = other.anchorMass
    anchorMasses = other.anchorMasses
    nonAnchorMass = other.nonAnchorMass
    nonAnchorMasses = other.nonAnchorMasses
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
}

extension MM4RigidBodyStorage {
  func ensurePositionsCached() {
    if self.positions == nil {
      self.positions = createPositions()
    }
  }
  
  func ensureVelocitiesCached() {
    if self.velocities == nil {
      self.velocities = createVelocities()
    }
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
  
  func setThermalKineticEnergy(_ energy: Double) {
    ensureKineticEnergyCached()
    
    // Don't change the velocities unless thermal kinetic energy has changed.
    if energy != thermalKineticEnergy! {
      velocities = nil
      thermalKineticEnergy = energy
      createThermalVelocities()
    }
  }
}

extension MM4RigidBody {
  /// Ensure all weak references point to the current storage object.
  ///
  /// This function is underscored to prevent it from appearing in autocomplete.
  mutating func _ensureReferencesUpdated() {
    energy.kinetic.storage = storage
  }
  
  /// Ensure copy-on-write semantics.
  ///
  /// > WARNING: Call this before every mutating function.
  mutating func ensureUniquelyReferenced() {
    if !isKnownUniquelyReferenced(&storage) {
      storage = MM4RigidBodyStorage(copying: storage)
      _ensureReferencesUpdated()
    }
  }
}


