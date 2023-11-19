//
//  MM4ForceField+Simulate.swift
//
//
//  Created by Philip Turner on 10/21/23.
//

import OpenMM

extension MM4ForceField {
  /// Simulate the system's evolution for the specified time interval.
  ///
  /// - Parameter time: The time interval, in picoseconds.
  /// - Parameter timeStep: The largest time step that may be taken during the
  ///   simulation, in picoseconds. Some steps may have a smaller duration.
  /// - throws: <doc:MM4Error/energyDrift(_:)> if energy tracking is enabled.
  public func simulate(
    time: Double,
    timeStep: Double = 100 / 23 * OpenMM_PsPerFs
  ) throws {
    // Check whether the arguments are valid.
    guard time >= 0, timeStep > 0 else {
      fatalError("Time or time step was invalid.")
    }
    if time == 0 {
      return
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
    
    func createEnergy() -> Double {
      if trackingEnergy {
        var stateDescriptor = MM4StateDescriptor()
        stateDescriptor.energy = true
        
        let state = self.state(descriptor: stateDescriptor)
        return state.kineticEnergy! + state.potentialEnergy!
      } else {
        return 0
      }
    }
    let startEnergy = createEnergy()
    
    if quotient == 0 {
      var descriptor = MM4IntegratorDescriptor()
      descriptor.start = true
      descriptor.end = true
      context.currentIntegrator = descriptor
      context.step(1, timeStep: time)
    } else {
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
    
    let endEnergy = createEnergy()
    if abs(endEnergy - startEnergy) > thresholdEnergy {
      throw MM4Error.energyDrift(endEnergy - startEnergy)
    }
  }
}
