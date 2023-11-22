//
//  MM4ForceField+Simulate.swift
//
//
//  Created by Philip Turner on 10/21/23.
//

import OpenMM

/// A level of theory for the simulator.
public enum MM4LevelOfTheory: CaseIterable {
  /// The default time step is 4.348 femtoseconds.
  case molecularDynamics
  
  /// The default time step has yet to be determined.
  case rigidBodyDynamics
  
  /// The default time step in picoseconds.
  public var defaultTimeStep: Double {
    switch self {
    case .molecularDynamics:
      return 100 / 23 * OpenMM_PsPerFs
    case .rigidBodyDynamics:
      return 100 / 23 * OpenMM_PsPerFs
    }
  }
}

extension MM4ForceField {
  /// Simulate the system's evolution for the specified time interval.
  ///
  /// - Parameter time: The time interval, in picoseconds.
  /// - throws: <doc:MM4Error/energyDrift(_:)> if energy tracking is enabled.
  public func simulate(time: Double) throws {
    if updateRecord.active() {
      flushUpdateRecord()
      invalidatePositionsAndVelocities()
      invalidateForcesAndEnergy()
    }
    
    func createEnergy() -> Double {
      if trackingEnergy {
        return kineticEnergy + potentialEnergy
      } else {
        return 0
      }
    }
    let startEnergy = createEnergy()
    invalidatePositionsAndVelocities()
    invalidateForcesAndEnergy()
    
    guard _levelOfTheory.allSatisfy({ $0 == .molecularDynamics }) else {
      fatalError("Rigid body dynamics not implemented yet.")
    }
    if time == 0 {
      return
    }
    
    // Check whether the arguments are valid.
    let timeStep = self._timeStep[.molecularDynamics]!
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

extension MM4ForceField {
  /// The level of theory used to simulate each rigid body.
  ///
  /// The default value is `.molecularDynamics` for each rigid body. All rigid
  /// bodies must have the same level of theory, until hybrid simulation is
  /// supported.
  public func setLevelOfTheory(_ level: MM4LevelOfTheory, rigidBodyIndex: Int) {
    guard rigidBodyIndex >= 0, rigidBodyIndex < _levelOfTheory.count else {
      fatalError("Rigid body index out of bounds.")
    }
    _levelOfTheory[rigidBodyIndex] = level
  }
  
  /// The largest time step that may be taken during simulation, in picoseconds.
  /// Some steps may have a smaller duration.
  ///
  /// If not specified, the time step defaults to the value of
  /// [`defaultTimeStep`](<doc:MM4LevelOfTheory/defaultTimeStep>).
  public func setTimeStep(_ timeStep: Double, levelOfTheory: MM4LevelOfTheory) {
    _timeStep[levelOfTheory] = timeStep
  }
}
