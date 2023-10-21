//
//  MM4ForceField+Properties.swift
//
//
//  Created by Philip Turner on 10/14/23.
//

// A set of actions that require O(n) compute-intensive procedures, either on
// the CPU or the GPU.

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
  
  /// Minimize the system's energy using the L-BFGS algorithm.
  ///
  /// This is one of few such algorithms with O(n) computational complexity. It
  /// is a limited-memory version of BFGS, an O(n^2) algorithm. BFGS, in turn,
  /// is an improvement on O(n^3) methods such as Newton's method.
  ///
  /// - Parameter tolerance: Accepted uncertainty in potential energy,
  ///   in zeptojoules.
  /// - Parameter maxIterations: Maximum number of force evaluations permitted
  ///   during the minimization. The default value, 0, puts no restrictions on
  ///   the number of evaluations.
  /// - throws: An error if the energy exploded.
  public func minimize(
    tolerance: Double = 10.0 * MM4ZJPerKJPerMol,
    maxIterations: Int = 0
  ) throws {
    // Use a different integrator that doesn't incorporate the fusion of
    // multiple timesteps. Also, zero out the system's bulk velocity during the
    // minimization. Restore the particles' velocities after the minimization is
    // finished.
  }
  
  /// Create random thermal velocities, while conserving the total (bulk)
  /// momentum of each rigid body.
  ///
  /// - Parameter temperature: The temperature to randomize thermal velocites
  ///   at, in kelvin.
  /// - Parameter rigidBodies: Indices of the rigid bodies to thermalize. If not
  ///   specified, it will thermalize the entire system.
  ///
  /// Thermalizing is recommended for any simulation that replicates macroscale
  /// conditions. The default is temperature 298.15 K, but other useful
  /// temperatures include liquid nitrogen (77.00 K) and liquid helium (4.15 K).
  ///
  /// Anchors have no velocity appended to them during thermalization. Before
  /// angular momentum is corrected around the true center of mass, it is first
  /// corrected around a center defined by any anchors. If a rigid body has only
  /// one anchor, that anchor's position is used. Otherwise, the center is the
  /// average position of all anchors in the rigid body.
  public func thermalize(
    temperature: Double = 298.15,
    rigidBodies: [Int]? = nil
  ) {
    
  }
}
