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
  public var cutoffDistance: Float = 1.0
  
  /// Required. The dielectric constant for the reaction field approximation to
  /// long-range electrostatic interactions.
  ///
  /// The default value is 5.7, the same as solid diamond. The actual
  /// value may fluctuate between 0 (vacuum) and greater than 10.0, for
  /// patches of moissanite or silicon. It is reasonable to assume the average
  /// is somewhere in the middle, at ~5.
  public var dielectricConstant: Float = 5.7
  
  /// Optional. The parameters to initialize internal forces with.
  ///
  /// This parameter is mutually exclude with `rigidBodies`. You can either
  /// specify the parameters directly or enter them through the rigid bodies.
  /// You cannot do both. If `rigidBodies` is not specified, you must set
  /// `parameters`.
  public var parameters: MM4Parameters?
  
  /// Optional. The OpenMM platform to use for simulation.
  public var platform: OpenMM_Platform?
  
  /// Optional. The rigid bodies to initialize the system with.
  ///
  /// If you do not set the rigid bodies, you must set all positions, velocities,
  /// and external forces manually. Their default values are all zero. This is the
  /// same behavior as <doc:MM4RigidBody/init(parameters:)>.
  public var rigidBodies: [MM4RigidBody]?
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
    var parameters: MM4Parameters
    switch (descriptor.parameters, descriptor.rigidBodies) {
    case (.some(_), nil):
      parameters = descriptor.parameters!
    case (nil, .some(let rigidBodies)):
      parameters = MM4ForceField.createParameters(rigidBodies: rigidBodies)
    case (nil, nil):
      fatalError("Did not specify parameters or rigid bodies.")
    case (.some(_), .some(_)):
      fatalError("Specified both parameters and rigid bodies.")
    }
    
    MM4Plugins.global.load()
    system = MM4System(parameters: parameters, descriptor: descriptor)
    context = MM4Context(system: system, platform: descriptor.platform)
    cachedState = MM4State()
    updateRecord = MM4UpdateRecord()
    
    _energy = MM4ForceFieldEnergy(forceField: self)
    _externalForces = Array(
      repeating: .zero, count: system.parameters.atoms.count)
    
    if let rigidBodies = descriptor.rigidBodies {
      var atomCursor = 0
      for rigidBody in rigidBodies {
        let nextAtomCursor = atomCursor + rigidBody.parameters.atoms.count
        `import`(from: rigidBody, range: atomCursor..<nextAtomCursor)
        atomCursor = nextAtomCursor
      }
    }
  }
}
