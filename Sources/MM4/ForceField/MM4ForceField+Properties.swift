//
//  MM4ForceField+Properties.swift
//
//
//  Created by Philip Turner on 10/14/23.
//

// Ergonomic APIs for accessing a force field's state. These delegate to
// a single-property instances of the batched functions, 'update' or 'state'.
//
// Use the following OpenMM functions to update the system.
// externalForces
//   - setParticleParameters, updateParametersInContext
//   - setIntegrationForceGroups, swapping integrators
// positions
//   - setState
// stationaryAtoms
//   - setParticleParameters, updateParametersInContext
//   - setIntegrationForceGroups, swapping integrators
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
  /// The constant force (in piconewtons) exerted on each atom.
  ///
  /// The default value is all zeroes for every particle, which disables the
  /// backing OpenMM force object.
  public var externalForces: [SIMD3<Float>] {
    set { fatalError("Not implemented.") }
    get { fatalError("Not implemented.") }
  }
  
  /// Atom indices for each rigid body.
  ///
  /// Rigid bodies should have atoms laid out contiguously in memory, in Morton
  /// order. This format ensures spatial locality, which increases performance
  /// of nonbonded forces. Therefore, rigid bodies must be entered as contiguous
  /// ranges of the atom list.
  ///
  /// The set of rigid bodies must cover every atom in the system. No two ranges
  /// may overlap the same atom. If the array of rigid bodies is unspecified, it
  /// defaults to a range encompassing the entire system. This ensures the
  /// closed system's net momentum stays conserved.
  public var rigidBodies: [Range<Int>] {
    set { fatalError("Not implemented.") }
    get { fatalError("Not implemented.") }
  }
  
  /// Whether each atom's absolute position should never change.
  ///
  /// This is implemented by setting the particle's mass and velocity to zero.
  ///
  /// > Warning: Stationary atoms may cause energy measurements to be
  /// nonsensical. This needs to be investigated further.
  public var stationaryAtoms: [Bool] {
    set { fatalError("Not implemented.") }
    get { fatalError("Not implemented.") }
  }
}
