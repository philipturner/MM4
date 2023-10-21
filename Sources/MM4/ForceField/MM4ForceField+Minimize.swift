//
//  MM4ForceField+Minimize.swift
//
//
//  Created by Philip Turner on 10/21/23.
//

extension MM4ForceField {
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
}
