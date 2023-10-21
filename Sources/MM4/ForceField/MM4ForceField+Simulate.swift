//
//  MM4ForceField+Simulate.swift
//
//
//  Created by Philip Turner on 10/21/23.
//

extension MM4ForceField {
  /// Simulate the system's evolution for the specified time interval.
  ///
  /// - Parameter time: The time interval, in picoseconds.
  /// - Parameter maximumTimeStep: The largest time step that may be taken
  ///   during the simulation, in picoseconds.
  /// - throws: An error if the energy exploded.
  public func simulate(
    time: Double,
    maximumTimeStep: Double = 0.100 / 23 + 1e-8
  ) throws {
    // If the time doesn't divide evenly into 100 fs, use a temporary
    // integrator for the remainder, potentially with a slightly scaled-down
    // timestep.
    //
    // Dynamically switch between the integrators below. There should be
    // multiple OpenMM contexts surrounding the same system.
    // - 1-step custom MTS integrator (always used for minimization)
    // - 5-step custom MTS integrator
    // - 23-step custom MTS integrator
    //
    // Create a dictionary of OpenMM contexts, indexed by the number of steps
    // embedded in their integrators. 0 steps means the verlet integrator for
    // energy minimization. Dictionary elements should be lazily initialized
    // using an API for safely extracting them.
  }
}
