//
//  MM4ForceField+RigidBody.swift
//  MM4
//
//  Created by Philip Turner on 11/20/23.
//

extension MM4ForceField {
  /// Set the constant force (in piconewtons) exerted on each atom.
  ///
  /// The default value is zero for every atom.
  public func setExternalForces(_ forces: [SIMD3<Float>]) {
    guard forces.count == system.parameters.atoms.count else {
      fatalError("Number of external forces does not match atom count.")
    }
    
    let force = system.forces.external
    force.updateForces(forces, system: system)
    force.updateParametersInContext(context)
  }
  
  /// The net varying force (in piconewtons) exerted on each atom.
  public var forces: [SIMD3<Float>] {
    _read {
      ensureForcesAndEnergyCached()
      yield cachedState.forces!
    }
  }
  
  /// The position (in nanometers) of each atom's nucleus.
  public var positions: [SIMD3<Float>] {
    _read {
      ensurePositionsAndVelocitiesCached()
      yield cachedState.positions!
    }
    _modify {
      ensurePositionsAndVelocitiesCached()
      updateRecord.positions = true
      updateRecord.velocities = true
      yield &cachedState.positions!
    }
  }
  
  /// The velocity (in nanometers per picosecond) of each atom.
  public var velocities: [SIMD3<Float>] {
    _read {
      ensurePositionsAndVelocitiesCached()
      yield cachedState.velocities!
    }
    _modify {
      ensurePositionsAndVelocitiesCached()
      updateRecord.positions = true
      updateRecord.velocities = true
      yield &cachedState.velocities!
    }
  }
}

extension MM4ForceField {
  /// Write the force field's forces into the specified rigid body.
  /// - parameter rigidBody: The rigid body to update.
  /// - parameter range: The range of atoms in the force field that the rigid
  ///   body covers.
  public func exportForces(
    to rigidBody: inout MM4RigidBody,
    range: Range<Int>
  ) {
    // Cache the I/O accesses into OpenMM, otherwise this is O(n^2).
    if updateRecord.active() {
      flushUpdateRecord()
    }
    ensureForcesAndEnergyCached()
    guard range.startIndex >= 0,
          range.endIndex < system.parameters.atoms.count,
          range.count == rigidBody.parameters.atoms.count else {
      fatalError("Atom range was invalid.")
    }
    
    // Change the rigid body's forces.
    // TODO
  }
  
  /// Change the force field's positions to match the specified rigid body.
  /// - parameter rigidBody: The rigid body to update the force field's state
  ///   with.
  /// - parameter range: The range of atoms in the forcefield that the rigid
  ///   body covers.
  public func importPositions(
    from rigidBody: MM4RigidBody,
    range: Range<Int>
  ) {
    // Cache the I/O accesses into OpenMM, otherwise this is O(n^2).
    ensurePositionsAndVelocitiesCached()
    updateRecord.positions = true
    guard range.startIndex >= 0,
          range.endIndex < system.parameters.atoms.count,
          range.count == rigidBody.parameters.atoms.count else {
      fatalError("Atom range was invalid.")
    }
    
    // Change the force field's positions.
    cachedState.positions!.replaceSubrange(range, with: rigidBody.positions)
  }
  
  /// Change the force field's velocities to match the specified rigid body.
  /// - parameter rigidBody: The rigid body to update the force field's state
  ///   with.
  /// - parameter range: The range of atoms in the forcefield that the rigid
  ///   body covers.
  public func importVelocities(
    from rigidBody: MM4RigidBody,
    range: Range<Int>
  ) {
    // Cache the I/O accesses into OpenMM, otherwise this is O(n^2).
    ensurePositionsAndVelocitiesCached()
    updateRecord.velocities = true
    guard range.startIndex >= 0,
          range.endIndex < system.parameters.atoms.count,
          range.count == rigidBody.parameters.atoms.count else {
      fatalError("Atom range was invalid.")
    }
    
    // Change the force field's velocities.
    cachedState.velocities!.replaceSubrange(range, with: rigidBody.velocities)
  }
}
