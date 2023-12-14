//
//  MM4ForceField.swift
//
//
//  Created by Philip Turner on 9/10/23.
//

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
  
  /// Stores the anchor IDs before reordering.
  var _anchors: Set<UInt32> = []
  
  /// Stores the system's energy.
  var _energy: MM4ForceFieldEnergy!
  
  /// Stores the external forces before reordering.
  var _externalForces: [SIMD3<Float>] = []
  
  /// Stores the time step for each level of theory.
  var _timeStep: [MM4LevelOfTheory: Double] = [:]
  
  // MARK: - Properties for Rigid Bodies
  
  /// Stores the level of theory for each rigid body.
  var _levelOfTheory: [MM4LevelOfTheory] = []
  
  /// Stores the anchor IDs separately for each rigid body.
  var _rigidBodyAnchors: [Set<UInt32>] = []
  
  /// Stores the external forces separately for each rigid body.
  var _rigidBodyExternalForces: [SIMD3<Float>] = []
  
  /// Stores the handles separately for each rigid body.
  var _rigidBodyHandles: [Set<UInt32>] = []
  
  /// Stores the atom range for each rigid body.
  var _rigidBodyRanges: [Range<UInt32>] = []
  
  /// Create a simulator using the specified parameters and division into rigid
  /// bodies.
  public init(parameters: MM4Parameters, rigidBodyRanges: [Range<UInt32>]) {
    MM4Plugins.global.load()
    system = MM4System(parameters: parameters)
    context = MM4Context(system: system)
    cachedState = MM4State()
    updateRecord = MM4UpdateRecord()
    
    // TODO: Find a way to fail when the user simulates gold atoms, without
    // setting them as anchors or using rigid body mechanics. The current
    // implementation treats them like noble gas particles.
    
    _energy = MM4ForceFieldEnergy(forceField: self)
    
    _externalForces = Array(
      repeating: .zero, count: system.parameters.atoms.count)
    
    for level in MM4LevelOfTheory.allCases {
      _timeStep[level] = level.defaultTimeStep
    }
    
    let rigidBodyCount = rigidBodyRanges.count
    _levelOfTheory = Array(
      repeating: parameters.levelOfTheory, count: rigidBodyCount)
    _rigidBodyAnchors = Array(repeating: [], count: rigidBodyCount)
    _rigidBodyExternalForces = Array(repeating: .zero, count: rigidBodyCount)
    _rigidBodyHandles = Array(repeating: [], count: rigidBodyCount)
    _rigidBodyRanges = rigidBodyRanges
  }
}
