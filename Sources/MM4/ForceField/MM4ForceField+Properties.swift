//
//  MM4ForceField+Properties.swift
//
//
//  Created by Philip Turner on 10/14/23.
//

// Ergonomic APIs for accessing a force field's state. These delegate to
// a single-property instances of the batched functions, 'update' or 'state'.
// They are not intended to be used as a primary interface for interacting with
// the system (they have high overhead). Rather, the calling application's IR
// should handle most of the data processing. One should only transfer data
// in/out of OpenMM to set up the simulation.

// Use the following OpenMM functions to update the system.
// anchors
//   - setParticleParameters, updateParametersInContext
// externalForces
//   - setParticleParameters, updateParametersInContext
// positions
//   - setState
// velocities
//   - setState

// MARK: - Batched Functions

extension MM4ForceField {
  /// The net varying force (in piconewtons) exerted on each atom.
  ///
  /// > Note: This is a more ergonomic API, but less efficient than the batched
  /// functions <doc:MM4ForceField/update(descriptor:)> and
  /// <doc:MM4ForceField/state(descriptor:)>.
  public var forces: [SIMD3<Float>] {
    get { fatalError("Not implemented.") }
  }
  
  /// The system's total kinetic energy, in zeptojoules.
  ///
  /// > Note: This is a more ergonomic API, but less efficient than the batched
  /// functions <doc:MM4ForceField/update(descriptor:)> and
  /// <doc:MM4ForceField/state(descriptor:)>.
  public var kineticEnergy: Double {
    get { fatalError("Not implemented.") }
  }
  
  /// The position (in nanometers) of each atom's nucleus.
  ///
  /// > Note: This is a more ergonomic API, but less efficient than the batched
  /// functions <doc:MM4ForceField/update(descriptor:)> and
  /// <doc:MM4ForceField/state(descriptor:)>.
  public var positions: [SIMD3<Float>] {
    set { fatalError("Not implemented.") }
    get { fatalError("Not implemented.") }
  }
  
  /// The system's total potential energy, in zeptojoules.
  ///
  /// > Note: This is a more ergonomic API, but less efficient than the batched
  /// functions <doc:MM4ForceField/update(descriptor:)> and
  /// <doc:MM4ForceField/state(descriptor:)>.
  public var potentialEnergy: Double {
    get { fatalError("Not implemented.") }
  }
  
  /// The bulk + thermal velocity (in nanometers per picosecond) of each atom.
  ///
  /// > Note: This is a more ergonomic API, but less efficient than the batched
  /// functions <doc:MM4ForceField/update(descriptor:)> and
  /// <doc:MM4ForceField/state(descriptor:)>.
  ///
  /// When thermalizing, the linear and angular momentum over every rigid body
  /// is conserved. Then, the thermal velocities are reinitialized. If you want
  /// more complex motion within the rigid body, fetch the thermalized
  /// velocities. Add the desired bulk velocity component to them, them set the
  /// new velocity values.
  public var velocities: [SIMD3<Float>] {
    set { fatalError("Not implemented.") }
    get { fatalError("Not implemented.") }
  }
}

// MARK: - Simulation Setup

extension MM4ForceField {
  /// Indices of atoms that should ignore forces exerted on them.
  ///
  /// > Warning: Anchors may cause energy measurements to be
  /// nonsensical. This needs to be investigated further.
  ///
  /// This is implemented by treating particle like it has infinite mass (in
  /// OpenMM, by setting the mass to zero). Either set the particle's velocity
  /// to zero, or give all anchors in same rigid body the same linear velocity.
  /// Otherwise, the divergent paths will eventually rip the rigid body apart.
  /// The library currently doesn't check for adherence to this rule, but may
  /// enforce it in the future.
  public var anchors: [UInt32] {
    get { fatalError("Not implemented.") }
    set {
      // Reset every atom's mass to the one provided by the parameters. Then,
      // selectively set the anchors to zero. This may have more overhead, but
      // not that much more (~overhead of setting positions). It also reduces
      // the chance for bugs in a rarely tested edge case.
      fatalError("Not implemented.")
    }
  }
  
  // The force field could be optimized by deactivating the external force for
  // an incremental gain in performance. However, such fine-tuning optimizations
  // will come at a later date. The current "optimization" is that the force is
  // added each timestep (4 fs), but only computed every time the integrator
  // loop starts up.
  
  /// The constant force (in piconewtons) exerted on each atom.
  ///
  /// The default value is all zeroes for every particle, which may be used to
  /// deactivate the backing OpenMM force object.
  public var externalForces: [SIMD3<Float>] {
    get { fatalError("Not implemented.") }
    set { fatalError("Not implemented.") }
  }
  
  /// Atom indices for each rigid body.
  ///
  /// > Note: This is similar to molecules (`OpenMM_Context.getMolecules`),
  /// with an additional restriction. The user must enter the atoms from each
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
  public var rigidBodies: [Range<UInt32>] {
    get {
      // Create a new array with a different type. This isn't the fastest
      // approach, but the property should rarely be used in most cases. The
      // user should already have a data structure that separates the atoms into
      // rigid bodies during high-level operations.
      system.parameters.rigidBodies.map {
        let lowerBound = UInt32(truncatingIfNeeded: $0.lowerBound)
        let upperBound = UInt32(truncatingIfNeeded: $0.upperBound)
        return lowerBound..<upperBound
      }
    }
  }
}
