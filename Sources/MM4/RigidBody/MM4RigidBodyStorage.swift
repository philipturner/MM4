//
//  MM4RigidBodyStorage.swift
//  MM4
//
//  Created by Philip Turner on 11/22/23.
//

final class MM4RigidBodyStorage {
  // Sources of truth.
  var atoms: (count: Int, vectorCount: Int, nonAnchorCount: Int)
  var externalForces: [SIMD3<Float>]
  var mass: Float
  var vMasses: [MM4FloatVector]
  var vPositions: [MM4FloatVector]
  var vVelocities: [MM4FloatVector]
  
  // Frequently cached (special handling).
  var centerOfMass: SIMD3<Float>?
  var positions: [SIMD3<Float>]?
  var velocities: [SIMD3<Float>]?
  
  // Rarely cached (frequently erased).
  var angularVelocity: SIMD3<Float>?
  var linearVelocity: SIMD3<Float>?
  var momentOfInertia: (SIMD3<Float>, SIMD3<Float>, SIMD3<Float>)?
  
  init(parameters: MM4Parameters) {
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
      self.mass = Float(vMassAccumulator.sum())
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
    angularVelocity = nil
    linearVelocity = nil
    momentOfInertia = nil
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
        momentOfInertia! == (.zero, .zero, .zero),
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
}

extension MM4RigidBody {
  /// Ensure copy-on-write semantics.
  ///
  /// > WARNING: Call this before every mutating function.
  mutating func ensureUniquelyReferenced() {
    if !isKnownUniquelyReferenced(&storage) {
      storage = MM4RigidBodyStorage(copying: storage)
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
  public var mass: Float {
    storage.mass
  }
}
