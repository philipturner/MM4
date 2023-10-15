//
//  MM4ForceField+Update.swift
//
//
//  Created by Philip Turner on 10/14/23.
//

/// A configuration for an update to a simulated system.
public class MM4ForceFieldUpdateDescriptor {
  /// The constant force (in piconewtons) exerted on each atom.
  public var externalForces: [SIMD3<Float>]?
  
  /// The position (in nanometers) of each atom's nucleus.
  public var positions: [SIMD3<Float>]?
  
  /// Whether each atom's absolute position should never change.
  public var stationaryAtoms: [Bool]?
  
  /// The bulk velocity (in nanometers per picosecond), of each atom.
  ///
  /// If you set this, it is a good idea to
  /// [re-thermalize](<doc:MM4ForceField/thermalize(temperature:atoms:)>) the
  /// system.
  public var velocities: [SIMD3<Float>]?
  
  // In the future, "force field mode" may be updatable too.
}

extension MM4ForceField {
  /// Update the state of the system, while the simulation is being interrupted.
  public func update(descriptor: MM4ForceFieldUpdateDescriptor) {
    // Take the OpenMM object backing the `MM4System`, and call the respective
    // APIs to update it.
  }
}
