//
//  MM4ForceField+Properties.swift
//
//
//  Created by Philip Turner on 10/14/23.
//

import OpenMM

// Idea for adding torques:
// - Can't constrain "anchors" to have constant angular velocity, only constant
//   linear velocity.
// - A simulation can achieve near-constant angular velocity with a large
//   flywheel (which will be provided in the hardware catalog, and documented in
//   MM4).
// - Torques around a computed position:
//   - LocalCoordinatesSite - center of mass of a rigid body, exerts a
//     counterforce on the rigid body as a whole
//   - relative to a set of collinear anchor particles, which must all have the
//     same velocity
//     - check the collinearity constraint every time the user changes particle
//       velocities or modifies the anchors
// - Requires a new type (`MM4Torque`) that wraps a Swift `Quaternion` and an
//   enumeration with an associated value, which specifies the type of origin.

/// A data structure wrapping a system's energy.
public struct MM4ForceFieldEnergy {
  var forceField: MM4ForceField
  
  var _explosionThreshold: Double
  
  /// Whether to throw an error during an energy explosion.
  ///
  /// > Warning: Enabling this feature may significantly degrade performance.
  ///
  /// The default is `false`.
  ///
  /// Energy is tracked in low precision, as high precision is not needed
  /// to detect energy explosions.
  public var tracked: Bool = false
  
  init(forceField: MM4ForceField) {
    let atoms = forceField.system.parameters.atoms
    self.forceField = forceField
    self._explosionThreshold = 1e6 * (Double(atoms.count) / 1e4)
  }
  
  /// The threshold at which energy is considered to have exploded.
  ///
  /// The default is 1 million zJ per 10,000 atoms.
  public var explosionThreshold: Double {
    get { _explosionThreshold }
    set {
      guard newValue > 0 else {
        fatalError("Explosion threshold must be positive and nonzero.")
      }
      _explosionThreshold = newValue
    }
  }
  
  /// The system's total kinetic energy, in zeptojoules.
  public var kinetic: Double {
    forceField.ensureForcesAndEnergyCached()
    return forceField.cachedState.kineticEnergy!
  }
  
  /// The system's total potential energy, in zeptojoules.
  public var potential: Double {
    forceField.ensureForcesAndEnergyCached()
    return forceField.cachedState.potentialEnergy!
  }
}

extension MM4ForceField {
  /// The system's energy.
  ///
  /// To make the default behavior have high performance, energy is reported in
  /// low precision. To request a high-precision estimate, fetch it using an
  /// `MM4State`.
  public var energy: MM4ForceFieldEnergy {
    _energy
  }
  
  /// The net varying force (in piconewtons) exerted on each atom.
  public var forces: [SIMD3<Float>] {
    _read {
      ensureForcesAndEnergyCached()
      yield cachedState.forces!
    }
  }
  
  /// The position (in nanometers) of each atom's nucleus.
  public var positions: [SIMD3<Float>] {
    _read {
      ensurePositionsAndVelocitiesCached()
      yield cachedState.positions!
    }
    _modify {
      ensurePositionsAndVelocitiesCached()
      updateRecord.positions = true
      updateRecord.velocities = true
      yield &cachedState.positions!
    }
  }
  
  /// The bulk + thermal velocity (in nanometers per picosecond) of each atom.
  ///
  /// When thermalizing, the linear and angular momentum over every rigid body
  /// is conserved. Then, the thermal velocities are reinitialized. If you want
  /// more complex motion within the rigid body, fetch the thermalized
  /// velocities. Add the desired bulk velocity component to them, them set the
  /// new velocity values.
  public var velocities: [SIMD3<Float>] {
    _read {
      ensurePositionsAndVelocitiesCached()
      yield cachedState.velocities!
    }
    _modify {
      ensurePositionsAndVelocitiesCached()
      updateRecord.positions = true
      updateRecord.velocities = true
      yield &cachedState.velocities!
    }
  }
}

extension MM4ForceField {
  /// Indices of atoms that should be treated as having infinite mass in a
  /// simulation.
  ///
  /// An anchor's velocity does not vary due to thermal energy. Angular
  /// momentum is constrained according to the number of anchors present.
  /// - 0 anchors: conserve linear and angular momentum around center of mass.
  /// - 1 anchor: conserve linear and angular momentum around anchor.
  /// - collinear anchors: conserve linear momentum around average of
  ///   anchors, constrain angular momentum to the shared axis.
  /// - anchors form plane: conserve linear momentum around average of anchors,
  ///   force angular momentum to zero.
  public var anchors: Set<UInt32> {
    // _modify not supported b/c it requires very complex caching logic.
    // Workaround: import a new rigid body initialized with different anchors.
    _read {
      yield _anchors
    }
//    _modify {
//      updateRecord.anchors = true
//      yield &_anchors
//    }
  }
  
  /// The constant force (in piconewtons) exerted on each atom.
  ///
  /// The default value is zero for every atom.
  public var externalForces: [SIMD3<Float>] {
    _read {
      yield _externalForces
    }
    _modify {
      updateRecord.externalForces = true
      yield &_externalForces
    }
  }
  
  /// Atom indices for each rigid body.
  ///
  /// > Note: This is similar to molecules (`OpenMM_Context.getMolecules()`),
  /// with an additional restriction. The user must enter atoms for each
  /// molecule in one contiguous range of the atom list. Otherwise, the
  /// forcefield cannot initialize. See <doc:MM4ParametersDescriptor/bonds> for
  /// more details.
  ///
  /// Rigid bodies should have atoms laid out contiguously in memory, in Morton
  /// order. This format ensures spatial locality, which increases performance
  /// of nonbonded forces. Therefore, rigid bodies are contiguous ranges of the
  /// atom list.
  ///
  /// The set of rigid bodies must cover every atom in the system. No two ranges
  /// may overlap the same atom. If the array of rigid bodies is unspecified, it
  /// defaults to a range encompassing the entire system. This ensures the
  /// closed system's net momentum stays conserved.
  public var rigidBodyRanges: [Range<UInt32>] {
    _read {
      fatalError("TODO: New IR that stores rigid bodies in MM4ForceField")
    }
  }
}

