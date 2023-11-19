//
//  MM4ForceField+Energy.swift
//
//
//  Created by Philip Turner on 11/18/23.
//

import OpenMM

extension MM4ForceField {
  /// The system's total kinetic energy, in zeptojoules.
  ///
  /// > Note: This is a more ergonomic API, but less efficient than the batched
  /// function <doc:MM4ForceField/state(descriptor:)>.
  public var kineticEnergy: Double {
    get {
      var descriptor = MM4StateDescriptor()
      descriptor.energy = true
      return state(descriptor: descriptor).kineticEnergy!
    }
  }
  
  /// The system's total potential energy, in zeptojoules.
  ///
  /// > Note: This is a more ergonomic API, but less efficient than the batched
  /// function <doc:MM4ForceField/state(descriptor:)>.
  public var potentialEnergy: Double {
    get {
      var descriptor = MM4StateDescriptor()
      descriptor.energy = true
      return state(descriptor: descriptor).potentialEnergy!
    }
  }
  
  /// Whether to track and debug energy explosions during simulation.
  ///
  /// > Warning: Enabling this feature may significantly degrade performance.
  ///
  /// For simulation, this feature is disabled by default. It is always enabled
  /// for minimization. In the future, this API may be deprecated. The
  /// replacement will allow control over whether energy tracking occurs during
  /// minimization.
  public var trackingEnergy: Bool {
    get { _trackingEnergy }
    set { _trackingEnergy = newValue }
  }
  
  /// The threshold for energy explosion is 1 million zJ @ 10,000 atoms. This
  /// implementation detail is not exposed to the public API yet. That fact may
  /// change if a significant need arises.
  var thresholdEnergy: Double {
    1e6 * Double(system.parameters.atoms.count) / 1e4
  }
}
