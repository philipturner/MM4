//
//  MM4RigidBodyStorage.swift
//
//
//  Created by Philip Turner on 11/22/23.
//

import Numerics

final class MM4RigidBodyStorage {
  // Sources of truth.
  var atoms: (count: Int, vectorCount: Int, nonAnchorCount: Int)
  var externalForces: [SIMD3<Float>]
  var mass: Double
  var vMasses: [MM4FloatVector]
  var vPositions: [MM4FloatVector]
  var vVelocities: [MM4FloatVector]
  
  // Frequently cached (special handling).
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
    parameters: MM4Parameters
  ) {
    let atomCount = parameters.atoms.count
    let vectorCount = (atomCount + MM4VectorWidth - 1) / MM4VectorWidth
    var nonAnchorCount = 0
    for mass in parameters.atoms.masses {
      if mass > 0 {
        nonAnchorCount &+= 1
      } else if mass < 0 {
        fatalError("Mass cannot be negative.")
      }
    }
    self.atoms = (atomCount, vectorCount, nonAnchorCount)
    self.externalForces = Array(repeating: .zero, count: atomCount)
    print("checkpoint 0")
    
    // Pad the arrays of masses, positions, and velocities, so the out-of-bounds
    // vector lanes aren't NAN.
    self.vMasses = Array(unsafeUninitializedCapacity: vectorCount) {
      $0.initialize(repeating: .zero)
      $1 = vectorCount
    }
    self.vPositions = Array(unsafeUninitializedCapacity: 3 * vectorCount) {
      $0.initialize(repeating: .zero)
      $1 = 3 * vectorCount
    }
    self.vVelocities = Array(unsafeUninitializedCapacity: 3 * vectorCount) {
      $0.initialize(repeating: .zero)
      $1 = 3 * vectorCount
    }
    
    // Copy the masses from 'MM4Parameters' into the special memory allocation
    // for vectorized masses.
    self.mass = .zero
    initializeMasses(parameters.atoms.masses)
  }
  
  private func initializeMasses(_ masses: [Float]) {
    masses.withContiguousStorageIfAvailable {
      let opaque = OpaquePointer($0.baseAddress!)
      let casted = UnsafePointer<MM4FloatVector>(opaque)
      
      var vMassAccumulator: MM4DoubleVector = .zero
      for vID in 0..<atoms.vectorCount {
        var vMass: MM4FloatVector = .zero
        
        if vID < atoms.vectorCount &- 1 {
          vMass = casted[vID]
        } else {
          let maxLanes = atoms.count - vID * MM4VectorWidth
          for laneID in 0..<maxLanes {
            let mass = masses[vID &* MM4VectorWidth &+ laneID]
            vMass[laneID] = mass
          }
        }
        
        vMassAccumulator += MM4DoubleVector(vMass)
        self.vMasses[vID] = vMass
      }
      self.mass = vMassAccumulator.sum()
    }
  }
  
  init(copying other: MM4RigidBodyStorage) {
    // Initialize stored properties, without copying cached ones.
    atoms = other.atoms
    externalForces = other.externalForces
    mass = other.mass
    vMasses = other.vMasses
    vPositions = other.vPositions
    vVelocities = other.vVelocities
  }
  
  func eraseFrequentlyCachedProperties() {
    // WARNING: Always ensure this reflects recently added properties.
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
    }  else if atoms.count == 0 {
      precondition(
        centerOfMass! == .zero,
        "Nonzero center of mass for empty rigid body.")
    }
  }
  
  func ensureMomentOfInertiaCached() {
    if self.momentOfInertia == nil {
      self.momentOfInertia = createMomentOfInertia()
    } else if atoms.count == 0 {
      precondition(
        momentOfInertia!.columns == (.zero, .zero, .zero),
        "Nonzero moment of inertia for empty rigid body.")
    }
  }
  
  func ensureLinearVelocityCached() {
    if linearVelocity == nil {
      if atoms.count > 0 {
        linearVelocity = createLinearVelocity()
      } else {
        linearVelocity = .zero
      }
    } else if atoms.count == 0 {
      precondition(
        linearVelocity! == .zero,
        "Nonzero linear velocity for empty rigid body.")
    }
  }
  
  func ensureAngularVelocityCached() {
    if angularVelocity == nil {
      if atoms.count > 0 {
        angularVelocity = createAngularVelocity()
      } else {
        angularVelocity = .zero
      }
    } else if atoms.count == 0 {
      precondition(
        angularVelocity! == .zero,
        "Nonzero angular velocity for empty rigid body.")
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
      if angular < -epsilon * Double(atoms.count) {
        fatalError("Angular kinetic energy was too negative.")
      }
      if linear < -epsilon * Double(atoms.count) {
        fatalError("Linear kinetic energy was too negative.")
      }
      if thermal < -epsilon * Double(atoms.count) {
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
  ///
  /// This function is underscored to prevent it from appearing in autocomplete.
  mutating func _ensureReferencesUpdated() {
    _energy.kinetic.storage = storage
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

// MARK: - Properties

extension MM4RigidBody {
  /// The constant force (in piconewtons) exerted on each atom.
  public var externalForces: [SIMD3<Float>] {
    _read {
      yield storage.externalForces
    }
    _modify {
      ensureUniquelyReferenced()
      yield &storage.externalForces
    }
  }
  
  /// The total mass (in yoctograms) of all atoms, excluding anchors.
  public var mass: Double {
    storage.mass
  }
}
