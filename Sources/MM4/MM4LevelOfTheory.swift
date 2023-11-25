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
/// Once an `MM4ForceField` is created, the level of theory can changed. A
/// change occurs when the imported or exported rigid body has a different level
/// of theory. The following changes are allowed:
/// - `.extendedTightBinding` to `.forceField`
/// - `.forceField` to `.extendedTightBinding`
/// - `.molecularMechanics` to `.rigidBodyMechanics`
/// - `.rigidBodyMechanics` to `.molecularMechanics`
public enum MM4LevelOfTheory: CaseIterable, Hashable {
  /// One of the algorithms from the xTB package.
  ///
  /// The default time step is 2.174 femtoseconds.
  case quantumMechanics(GFNLevelOfTheory)
  
  /// The fourth molecular mechanics force field by Norman Allinger.
  ///
  /// The default time step is 4.348 femtoseconds.
  case molecularMechanics
  
  /// Simplified version of MM4 using only nonbonded forces, with the same
  /// restrictions on bonding topology.
  ///
  /// The default time step has yet to be determined.
  case rigidBodyMechanics
  
  /// The default time step in picoseconds.
  public var defaultTimeStep: Double {
    switch self {
    case .quantumMechanics:
      return 50 / 23 * OpenMM_PsPerFs
    case .molecularMechanics:
      return 100 / 23 * OpenMM_PsPerFs
    case .rigidBodyMechanics:
      return 100 / 23 * OpenMM_PsPerFs
    }
  }
  
  public static var allCases: [MM4LevelOfTheory] {
    return [
      .quantumMechanics(.extendedTightBinding),
      .quantumMechanics(.forceField),
      .molecularMechanics,
      .rigidBodyMechanics,
    ]
  }
}

/// A level of theory for quantum mechanics.
///
/// If the MM4 level of theory is `.quantumMechanics`, you must specify the
/// GFN level of theory as well.
public enum GFNLevelOfTheory: CaseIterable, Hashable {
  /// A semi-empirical approximation to density functional theory.
  ///
  /// In the xTB package, this level of theory is named "GFN2-xTB".
  case extendedTightBinding
  
  /// An expensive force field supporting any bonding topology.
  ///
  /// In the xTB package, this level of theory is named "GFN-FF".
  case forceField
}
