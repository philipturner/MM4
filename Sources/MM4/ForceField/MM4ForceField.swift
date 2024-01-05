//
//  MM4ForceField.swift
//
//
//  Created by Philip Turner on 9/10/23.
//

import OpenMM

/// A configuration for a force field.
///
/// > NOTE: This documentation page is still a draft. It may be inconsistent or
///   difficult to understand.
///
/// With this API, you have the set the positions, velocities, and external
/// forces manually. Their default values are all zero. This is the same
/// behavior as <doc:MM4RigidBody/init(parameters:)>.
///
/// We still deactivate certain forces by omitting parameter generation in
/// `MM4Parameters`. That was quite a powerful API design choice. But, there is
/// a new reason to include a force field descriptor. An esoteric use case where
/// you want to measure energy to beyond the (already reasonable) single
/// precision. Or, if you want to switch over to the 'Reference' implementation
/// to debug forces and/or energy. The user should have full control, at the
/// cost of being less ergonomic.
///
/// We design the ergonomic APIs so they should theoretically cover every single
/// use case. Less ergonomic APIs are provided when we foresee a critical issue
/// unique to million-atom systems, or other important use cases that would
/// otherwise be inaccessible. This includes energy conservation at the
/// precision of kT and modeling the shift in heat capacity from quantum
/// effects.
public struct MM4ForceFieldDescriptor {
  /// Required. The parameters to initialize internal forces with.
  public var parameters: MM4Parameters?
  
  /// Optional. The OpenMM platform to use for simulation.
  public var platform: OpenMM_Platform?
}

/// A force field simulator.
///
/// See the
/// [overview](https://philipturner.github.io/MM4/documentation/mm4#Overview)
/// for more information.
public class MM4ForceField {
  /// Delegate object wrapping the parameters and OpenMM system.
  var system: MM4System
  
  /// Delegate object wrapping the OpenMM context(s) and integrator.
  var context: MM4Context
  
  /// Caches I/O to OpenMM prevent `MM4RigidBody` exports from being O(n^2).
  var cachedState: MM4State
  
  /// Caches I/O to OpenMM prevent `MM4RigidBody` exports from being O(n^2).
  var updateRecord: MM4UpdateRecord
  
  /// Stores the system's energy.
  var _energy: MM4ForceFieldEnergy!
  
  /// Stores the external forces before reordering.
  var _externalForces: [SIMD3<Float>] = []
  
  /// Stores the time step, in picoseconds.
  var _timeStep: Double = 100 / 23 * OpenMM_PsPerFs
  
  /// Create a simulator using the specified configuration.
  public init(descriptor: MM4ForceFieldDescriptor) {
    guard let parameters = descriptor.parameters else {
      fatalError("Did not specify parameters in the force field descriptor.")
    }
    
    MM4Plugins.global.load()
    system = MM4System(parameters: parameters)
    context = MM4Context(system: system, platform: descriptor.platform)
    cachedState = MM4State()
    updateRecord = MM4UpdateRecord()
    
    _energy = MM4ForceFieldEnergy(forceField: self)
    _externalForces = Array(
      repeating: .zero, count: system.parameters.atoms.count)
  }
}
