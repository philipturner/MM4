//
//  MM4ForceField.swift
//
//
//  Created by Philip Turner on 9/10/23.
//

// There shouldn't need to be a force field descriptor. It should be possible
// to define its state entirely through the parameters. To have a simple
// implementation with just stretch/bend/nonbonded. Either the MM4Parameters
// was configured without torsions (motivating nonbonded exceptions to never be
// generated) or the user removed nonbondedExceptions14 manually.
//
// The user should specify the range where a rigid body's state is exported to.
// They can keep track of the ranges. In some cases, they might never need to
// preserve that information. Or they can decide that different atoms are
// classified as a single rigid body. The range they enter into the import/
// export functions must be self-consistent with the rigid body's atom count.

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
    _rigidBodyRanges = rigidBodyRanges
  }
}
