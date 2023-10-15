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

extension MM4ForceField {
  /// The constant force (in piconewtons) exerted on each atom.
  ///
  /// > Note: This is a more ergonomic API, but less efficient than the batched
  /// functions <doc:MM4ForceField/update(descriptor:)> and
  /// <doc:MM4ForceField/state(descriptor:)>.
  public var externalForces: [SIMD3<Float>]? {
    set { fatalError("Not implemented.") }
    get { fatalError("Not implemented.") }
  }
  
  /// The net varying force (in piconewtons) exerted on each atom.
  ///
  /// > Note: This is a more ergonomic API, but less efficient than the batched
  /// functions <doc:MM4ForceField/update(descriptor:)> and
  /// <doc:MM4ForceField/state(descriptor:)>.
  public var forces: [SIMD3<Float>] {
    get { fatalError("Not implemented.") }
  }
  
  /// The system's total kinetic energy.
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
  
  /// The system's total potential energy.
  ///
  /// > Note: This is a more ergonomic API, but less efficient than the batched
  /// functions <doc:MM4ForceField/update(descriptor:)> and
  /// <doc:MM4ForceField/state(descriptor:)>.
  public var potentialEnergy: Double {
    get { fatalError("Not implemented.") }
  }
  
  /// Whether each atom's absolute position should never change.
  ///
  /// > Note: This is a more ergonomic API, but less efficient than the batched
  /// functions <doc:MM4ForceField/update(descriptor:)> and
  /// <doc:MM4ForceField/state(descriptor:)>.
  public var stationaryAtoms: [Bool]? {
    set { fatalError("Not implemented.") }
    get { fatalError("Not implemented.") }
  }
  
  /// The bulk + thermal velocity (in nanometers per picosecond) of each atom.
  ///
  /// > Note: This is a more ergonomic API, but less efficient than the batched
  /// functions <doc:MM4ForceField/update(descriptor:)> and
  /// <doc:MM4ForceField/state(descriptor:)>.
  public var velocities: [SIMD3<Float>] {
    set { fatalError("Not implemented.") }
    get { fatalError("Not implemented.") }
  }
}
