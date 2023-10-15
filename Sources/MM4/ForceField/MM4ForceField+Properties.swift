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
  /// The force (in piconewtons) exerted on each atom.
  ///
  /// This is a more ergonomic API, but less efficient than the batched
  /// functions <doc:MM4ForceField/update(descriptor:)> and
  /// <doc:MM4ForceField/state(descriptor:)>.
  public var externalForces: [SIMD3<Float>]? {
    set { fatalError("Not implemented.") }
    get { fatalError("Not implemented.") }
  }
  
  /// The position (in nanometers) of each atom's nucleus.
  ///
  /// This is a more ergonomic API, but less efficient than the batched
  /// functions <doc:MM4ForceField/update(descriptor:)> and
  /// <doc:MM4ForceField/state(descriptor:)>.
  public var positions: [SIMD3<Float>] {
    set { fatalError("Not implemented.") }
    get { fatalError("Not implemented.") }
  }
  
  /// Whether each atom's absolute position should never change.
  public var stationaryAtoms: [Bool]? {
    set { fatalError("Not implemented.") }
    get { fatalError("Not implemented.") }
  }
  
  /// The velocity (in nanometers per picosecond), of each atom.
  public var velocities: [SIMD3<Float>] {
    set { fatalError("Not implemented.") }
    get { fatalError("Not implemented.") }
  }
}
