//
//  MM4ForceField.swift
//  MM4
//
//  Created by Philip Turner on 9/10/23.
//

import OpenMM

/// A configuration for a force field.
public struct MM4ForceFieldDescriptor {
  /// Required. The cutoff to use for nonbonded interactions.
  ///
  /// The default value is 1.0 nm. This is 2.86σ for carbon and 2.45σ for
  /// silicon.
  ///
  /// > NOTE: This documentation page is still a draft. It may be inconsistent or
  ///   difficult to understand.
  ///
  /// Since germanium will rarely be used, use the 2.5σ cutoff for silicon. The
  /// slightly greater sigma for carbon allows greater accuracy in vdW forces
  /// for bulk diamond. 1.0 nm is also sufficient for charge-charge
  /// interactions.
  ///
  /// If the cutoff distance is `nil`, nonbonded forces are computed with the
  /// full O(n^2) cost. This may increase simulation speed for small-enough
  /// systems and large-enough GPUs.
  public var cutoffDistance: Float? = 1.0
  
  /// Required. The dielectric constant for the reaction field approximation to
  /// long-range electrostatic interactions.
  ///
  /// The default value is 5.7, the same as solid diamond. The actual
  /// value may fluctuate between 0 (vacuum) and greater than 10.0, for
  /// patches of moissanite or silicon. It is reasonable to assume the average
  /// is somewhere in the middle, at ~5.
  public var dielectricConstant: Float = 5.7
  
  /// Required. The integrator to use for simulation.
  ///
  /// The default value is `.verlet`.
  public var integrator: MM4IntegratorOptions = .verlet
  
  /// Required. The parameters that define internal forces.
  ///
  /// If there are multiple rigid bodies, you can combine their parameters with
  /// <doc:MM4Parameters/append(contentsOf:)>.
  public var parameters: MM4Parameters?
  
  /// Optional. The OpenMM platform to use for simulation.
  public var platform: OpenMM_Platform?
  
  /// Create a descriptor with the default properties.
  public init() {
    
  }
}

/// A force field simulator.
///
/// See the MM4
/// [overview](https://philipturner.github.io/MM4/documentation/mm4#Overview)
/// for an explanation of the force field.
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
  
  /// Stores the time step, in picoseconds.
  var _timeStep: Double
  
  /// Create a simulator using the specified configuration.
  public init(descriptor: MM4ForceFieldDescriptor) throws {
    guard let parameters = descriptor.parameters else {
      fatalError("No force field parameters were specified.")
    }
    
    system = MM4System(parameters: parameters, descriptor: descriptor)
    context = MM4Context(system: system, descriptor: descriptor)
    cachedState = MM4State()
    updateRecord = MM4UpdateRecord()
    
    switch descriptor.integrator {
    case .multipleTimeStep:
      _timeStep = 4.35 * OpenMM_PsPerFs
    case .verlet:
      _timeStep = 2.5 * OpenMM_PsPerFs
    }
    _energy = MM4ForceFieldEnergy(forceField: self)
  }
}
