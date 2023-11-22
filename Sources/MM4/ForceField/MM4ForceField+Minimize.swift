//
//  MM4ForceField+Minimize.swift
//
//
//  Created by Philip Turner on 10/21/23.
//

import OpenMM

extension MM4ForceField {
  /// Minimize the system's potential energy, removing thermal potential energy.
  ///
  /// OpenMM uses the L-BFGS algorithm for gradient descent. It is an O(n)
  /// version of BFGS, an O(n^2) algorithm. BFGS improves upon O(n^3) methods
  /// such as Newton's method.
  ///
  /// - Parameter tolerance: Accepted uncertainty in potential energy,
  ///   in zeptojoules.
  /// - Parameter maxIterations: Maximum number of force evaluations permitted
  ///   during the minimization. The default value, 0, puts no restrictions on
  ///   the number of evaluations.
  /// - throws: <doc:MM4Error/energyDrift(_:)> if energy tracking is enabled.
  public func minimize(
    tolerance: Double = 10.0 * MM4ZJPerKJPerMol,
    maxIterations: Int = 0
  ) throws {
    flushUpdateRecord()
    invalidatePositionsAndVelocities()
    invalidateForcesAndEnergy()
    
    // Switch to an integrator that always reports the correct velocity.
    var integratorDescriptor = MM4IntegratorDescriptor()
    integratorDescriptor.start = true
    integratorDescriptor.end = true
    context.currentIntegrator = integratorDescriptor
    
    if trackingEnergy {
      // Record the current state.
      var stateDescriptor = MM4StateDescriptor()
      stateDescriptor.positions = true
      stateDescriptor.velocities = true
      let originalState = self.state(descriptor: stateDescriptor)
      
      // Check whether the system's energy will explode.
      func createEnergy() -> Double {
        var stateDescriptor = MM4StateDescriptor()
        stateDescriptor.energy = true
        
        let state = self.state(descriptor: stateDescriptor)
        return state.kineticEnergy! + state.potentialEnergy!
      }
      let startEnergy = createEnergy()
      context.step(1, timeStep: 1 * OpenMM_PsPerFs)
      let endEnergy = createEnergy()
      if abs(endEnergy - startEnergy) > thresholdEnergy {
        throw MM4Error.energyDrift(endEnergy - startEnergy)
      }
      
      // Restore the current state.
      self.positions = originalState.positions!
      self.velocities = originalState.velocities!
    }
    
    // Run the energy minimization.
    OpenMM_LocalEnergyMinimizer.minimize(
      context: context.context,
      tolerance: tolerance,
      maxIterations: maxIterations)
  }
}
