//
//  MM4ForceField+Simulate.swift
//  MM4
//
//  Created by Philip Turner on 10/21/23.
//

import OpenMM

extension MM4ForceField {
  /// The largest time step that may be taken during simulation, in picoseconds.
  /// Some steps may have a smaller duration.
  ///
  /// The default value is 4.35 fs.
  public var timeStep: Double {
    _read {
      yield _timeStep
    }
    _modify {
      yield &_timeStep
    }
  }
  
  /// Simulate the system's evolution for the specified time interval.
  ///
  /// - Parameter time: The time interval, in picoseconds.
  public func simulate(time: Double)  {
    if updateRecord.active() {
      flushUpdateRecord()
    }
    invalidatePositionsAndVelocities()
    invalidateForcesAndEnergy()
    
    // Internal state should be flushed the same way during this edge case as
    // when there is actually motion over time. It's less edge cases to check
    // and less complex behavior to reason about.
    if time == 0 {
      return
    }
    
    // Check whether the arguments are valid.
    guard time > 0, timeStep > 0 else {
      fatalError("Time or time step was invalid.")
    }
    
    // Create rough estimate of step count.
    var quotient = (time / timeStep).rounded(.down)
    var remainder = (time / timeStep) - quotient
    
    guard quotient >= 0, remainder >= 0 else {
      fatalError("This should never happen.")
    }
    
    // Correct for overshoot and undershoot from floating-point error.
    let epsilon: Double = 1e-4
    if remainder < epsilon {
      if quotient > 0 {
        quotient -= 1
        remainder += 1
      }
    } else if remainder > 1 + epsilon {
      fatalError("This should never happen.")
    }
    
    if quotient == 0 {
      print("block 1")
      var descriptor = MM4IntegratorDescriptor()
      descriptor.start = true
      descriptor.end = false
      context.currentIntegrator = descriptor
      context.step(1, timeStep: 0.002)
      
      descriptor.start = false
      descriptor.end = false
      context.currentIntegrator = descriptor
      context.step(9, timeStep: 0.002)
      
      descriptor.start = false
      descriptor.end = true
      context.currentIntegrator = descriptor
      context.step(1, timeStep: 0)
    } else {
      print("block 2")
      var descriptor = MM4IntegratorDescriptor()
      descriptor.start = true
      descriptor.end = false
      context.currentIntegrator = descriptor
      context.step(1, timeStep: timeStep)
      
      if quotient > 1 {
        descriptor.start = false
        descriptor.end = false
        context.currentIntegrator = descriptor
        context.step(Int(quotient - 1), timeStep: timeStep)
      }
      
      descriptor.start = false
      descriptor.end = true
      context.currentIntegrator = descriptor
      context.step(1, timeStep: remainder)
    }
  }
  
  /// Minimize the system's potential energy.
  ///
  /// OpenMM uses the L-BFGS algorithm for gradient descent. It is an O(n)
  /// version of BFGS, an O(n^2) algorithm. BFGS improves upon O(n^3) methods
  /// such as Newton's method.
  ///
  /// - Parameter tolerance: Accepted uncertainty in potential energy,
  ///   in zeptojoules. The default value is 10. This contrasts with the default
  ///   value for most OpenMM simulations, which is 16.6 zJ (10.0 kJ/mol).
  /// - Parameter maxIterations: Maximum number of force evaluations permitted
  ///   during the minimization. The default value, 0, puts no restrictions on
  ///   the number of evaluations.
  public func minimize(
    tolerance: Double = 10.0,
    maxIterations: Int = 0
  ) {
    if updateRecord.active() {
      flushUpdateRecord()
    }
    invalidatePositionsAndVelocities()
    invalidateForcesAndEnergy()
    
    // Switch to an integrator that always reports the correct velocity.
    var integratorDescriptor = MM4IntegratorDescriptor()
    integratorDescriptor.start = true
    integratorDescriptor.end = true
    context.currentIntegrator = integratorDescriptor
    
    // Run the energy minimization.
    //
    // The reporter doesn't do anything. You have to create a C++ class, which
    // is not possible through the OpenMM C API.
    OpenMM_LocalEnergyMinimizer.minimize(
      context: context.context,
      tolerance: tolerance,
      maxIterations: maxIterations,
      reporter: nil)
  }
}
