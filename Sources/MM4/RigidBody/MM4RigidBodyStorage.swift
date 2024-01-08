//
//  MM4RigidBodyStorage.swift
//  MM4
//
//  Created by Philip Turner on 11/22/23.
//

final class MM4RigidBodyStorage {
  // Sources of truth, during the last time they were updated externally.
  // TODO: Define a different reference frame for positions, velocities, and
  // forces. In addition, define the transformation to the reference frame
  // defined by the diagonalized inertia tensor, and the current global
  // reference frame.
  var atoms: (count: Int, vectorCount: Int, nonAnchorCount: Int)
  var mass: Double
  var forces: [SIMD3<Float>] // according to the outside observer
  // Use case for forces (maybe?) write on single-core, compute/cache on
  // multi-core. Or can they be written by MM4ForceField on multi-core?
  // Or persisting forces across simulation cycles, even as positions change
  // (weird, but reasonable behavior for the API).
  // - SIMD3<Float>
  // - slightly more memory consumption than vectorized format, but swizzling
  //   only happens once in typical use cases. Plus, you can forward it directly
  //   to the user without extra caching logic.
  // - swizzling can also happen in parallel, while otherwise it would be a
  //   bottleneck for single-core CPU, potentially synchronized routines?
  //
  // - make the forcefield cache forces on single-core
  // - export to the rigid body on multi-core
  // - do swizzle in-place, just don't compute acceleration in-place yet
  //   - acceleration can change is positions change, but forces do not
  //   - invalidate acceleration when positions change, don't recompute
  //     on-the-spot as <s>positions are</s> bulk orientation is mutated
  // - multicore bandwidth bottleneck from copying; compact representation can
  //   provide a meaningful speedup now
  // - don't provide forces of the individual body; the user can just fetch
  //   from MM4ForceField at a specific range if they want
  var vMasses: [MM4FloatVector]
  var vPositions: [MM4FloatVector] // relative positions in reference frame
  var vVelocities: [MM4FloatVector] // thermal velocities in reference frame
  
  var lastUpdateCenterOfMass: SIMD3<Double> // ???
  var lastUpdatePrincipalAxes: SIMD3<Double> // ???
  
  // The diagonalized moment of inertia and principal axes. These are
  // updated immediately after positions, velocities, or forces are mutated.
  var angularVelocity: SIMD3<Double>
  var centerOfMass: SIMD3<Double>
  var linearVelocity: SIMD3<Double>
  var momentOfInertia: SIMD3<Double>
  var principalAxes: (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>)
  
  // Cached properties, according to the outside observer.
  var angularAcceleration: SIMD3<Double>?
  var linearAcceleration: SIMD3<Double>?
  var forces: [SIMD3<Float>]?
  var positions: [SIMD3<Float>]?
  var velocities: [SIMD3<Float>]?
  
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
