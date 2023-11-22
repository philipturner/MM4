//
//  MM4ForceField.swift
//
//
//  Created by Philip Turner on 9/10/23.
//

/// A force field simulator.
///
/// See the [overview](``MM4``) for more information.
public class MM4ForceField {
  /// Delegate object wrapping the parameters and OpenMM system.
  var system: MM4System
  
  /// Delegate object wrapping the OpenMM context(s) and integrator.
  var context: MM4Context
  
  /// Caches I/O to OpenMM prevent `MM4RigidBody` exports from being O(n^2).
  var cachedState: MM4State
  
  /// Caches I/O to OpenMM prevent `MM4RigidBody` exports from being O(n^2).
  var updateRecord: MM4UpdateRecord
  
  /// Stores the anchor IDs before reordering.
  var _anchors: Set<UInt32> = []
  
  /// Stores the external forces before reordering.
  var _externalForces: [SIMD3<Float>] = []
  
  /// Stores the level of theory for each rigid body.
  var _levelOfTheory: [MM4LevelOfTheory] = []
  
  /// Stores the atom range for each rigid body.
  var _rigidBodyRanges: [Range<UInt32>] = []
  
  /// Stores the time step for each level of theory.
  var _timeStep: [MM4LevelOfTheory: Double] = [:]
  
  // .energy.tracked
  
  /// Whether to track and debug energy explosions.
  ///
  /// > Warning: Enabling this feature may significantly degrade performance.
  ///
  /// This feature is disabled by default.
  ///
  /// The energy is tracked in low precision, as high precision is not needed
  /// to detect energy explosions.
  public var trackingEnergy: Bool = false
  
  /// Create a simulator using the specified parameters and division into rigid
  /// bodies.
  public init(parameters: MM4Parameters, rigidBodyRanges: [Range<UInt32>]) {
    MM4Plugins.global.load()
    system = MM4System(parameters: parameters)
    context = MM4Context(system: system)
    cachedState = MM4State()
    updateRecord = MM4UpdateRecord()
    
    _externalForces = Array(
      repeating: .zero, count: system.parameters.atoms.count)
    
    _levelOfTheory = Array(
      repeating: .molecularDynamics, count: rigidBodyRanges.count)
    
    _rigidBodyRanges = rigidBodyRanges
    
    for level in MM4LevelOfTheory.allCases {
      _timeStep[level] = level.defaultTimeStep
    }
  }
}
