//
//  MM4LevelOfTheory.swift
//
//
//  Created by Philip Turner on 11/25/23.
//

import OpenMM

/// A level of theory for the simulator.
///
/// > WARNING: Currently, only `.molecularMechanics` is supported.
///
/// Once an `MM4ForceField` is created, the level of theory cannot be changed.
/// The level in an imported or exported rigid body is simply ignored.
public enum MM4LevelOfTheory: CaseIterable, Hashable {
  /// The fourth molecular mechanics force field by Norman Allinger.
  ///
  /// The default time step is 4.348 femtoseconds.
  case molecularMechanics
  
  /// Simplified version of MM4 using only nonbonded forces.
  ///
  /// The default time step has yet to be determined.
  case rigidBodyMechanics
  
  /// The default time step in picoseconds.
  public var defaultTimeStep: Double {
    switch self {
    case .molecularMechanics:
      return 100 / 23 * OpenMM_PsPerFs
    case .rigidBodyMechanics:
      return 100 / 23 * OpenMM_PsPerFs
    }
  }
  
  public static var allCases: [MM4LevelOfTheory] {
    return [
      .molecularMechanics,
      .rigidBodyMechanics,
    ]
  }
}
