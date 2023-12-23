//
//  MM4ForceField+Properties.swift
//
//
//  Created by Philip Turner on 10/14/23.
//

import OpenMM

extension MM4ForceField {
  /// The constant force (in piconewtons) exerted on each atom.
  public var externalForces: [SIMD3<Float>] {
    _read {
      yield _externalForces
    }
    _modify {
      updateRecord.externalForces = true
      yield &_externalForces
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
  
  /// The linear velocity (in nanometers per picosecond) of each atom.
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
  /// The following rules apply to anchors:
  /// - They must all have zero external force.
  /// - They must all have the same velocity.
  ///
  /// An anchor's velocity does not vary due to thermal energy. Angular
  /// momentum is constrained according to the number of anchors present.
  /// - 0 anchors: conserve linear and angular momentum around center of mass.
  /// - 1 anchor: conserve linear and angular momentum around anchor.
  /// - multiple anchors: conserve momentum around average of anchors.
  ///   In the average, each anchor's weight is proportional to its atomic mass.
  public var anchors: Set<UInt32> {
    // _modify not supported b/c it requires very complex caching logic.
    // Workaround: import a new rigid body initialized with different anchors.
    _read {
      yield _anchors
    }
  }
  
  /// The net varying force (in piconewtons) exerted on each atom.
  public var forces: [SIMD3<Float>] {
    _read {
      ensureForcesAndEnergyCached()
      yield cachedState.forces!
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
  /// The set of rigid bodies must cover every atom in the system. No two ranges
  /// may overlap the same atom.
  ///
  /// Rigid bodies should have atoms laid out contiguously in memory, in Morton
  /// order. This format ensures spatial locality, which increases performance
  /// of nonbonded forces. Therefore, rigid bodies are contiguous ranges of the
  /// atom list.
  public var rigidBodyRanges: [Range<UInt32>] {
    _read {
      yield _rigidBodyRanges
    }
  }
}
