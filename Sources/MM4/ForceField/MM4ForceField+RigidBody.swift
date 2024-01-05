//
//  MM4ForceField+RigidBody.swift
//
//
//  Created by Philip Turner on 11/20/23.
//

extension MM4ForceField {
  /// Create a simulator using the specified rigid bodies.
  public convenience init(rigidBodies: [MM4RigidBody]) {
    // Avoid a costly O(nlogn) series of reallocations while combining each
    // rigid body's parameters. If you're copying over the source code because
    // you need to use 'MM4ForceFieldDescriptor', you can exclude 90% of the
    // code here.
    var atomCapacity: Int = 0
    var bondCapacity: Int = 0
    var angleCapacity: Int = 0
    var torsionCapacity: Int = 0
    var ringCapacity: Int = 0
    var exception13Capacity: Int = 0
    var exception14Capacity: Int = 0
    
    var ranges: [Range<Int>] = []
    for rigidBody in rigidBodies {
      let oldAtomCapacity = atomCapacity
      
      let parameters = rigidBody.parameters
      atomCapacity += parameters.atoms.count
      bondCapacity += parameters.bonds.indices.count
      angleCapacity += parameters.angles.indices.count
      torsionCapacity += parameters.torsions.indices.count
      ringCapacity += parameters.rings.indices.count
      exception13Capacity += parameters.nonbondedExceptions13.count
      exception14Capacity += parameters.nonbondedExceptions14.count
      
      ranges.append(oldAtomCapacity..<atomCapacity)
    }
    
    // Elegantly handle the case where there are zero rigid bodies.
    var paramsDesc = MM4ParametersDescriptor()
    paramsDesc.atomicNumbers = []
    paramsDesc.bonds = []
    var parameters = try! MM4Parameters(descriptor: paramsDesc)
    parameters.atoms.reserveCapacity(atomCapacity)
    parameters.bonds.reserveCapacity(bondCapacity)
    parameters.angles.reserveCapacity(angleCapacity)
    parameters.torsions.reserveCapacity(torsionCapacity)
    parameters.rings.reserveCapacity(ringCapacity)
    parameters.nonbondedExceptions13.reserveCapacity(exception13Capacity)
    parameters.nonbondedExceptions14.reserveCapacity(exception14Capacity)
    parameters.atomsToBondsMap.reserveCapacity(bondCapacity)
    parameters.atomsToAtomsMap.reserveCapacity(atomCapacity)
    for rigidBody in rigidBodies {
      parameters.append(contentsOf: rigidBody.parameters)
    }
    
    var descriptor = MM4ForceFieldDescriptor()
    descriptor.parameters = parameters
    
    self.init(descriptor: descriptor)
    for (rigidBody, range) in zip(rigidBodies, ranges) {
      self.import(from: rigidBody, range: range)
    }
  }
  
  /// Write the force field's internal state to the specified rigid body.
  /// - parameter rigidBody: The rigid body to update.
  /// - parameter range: The range of atoms in the force field that the rigid
  ///   body covers.
  public func export(to rigidBody: inout MM4RigidBody, range: Range<Int>) {
    // Cache the I/O accesses into OpenMM, otherwise this is O(n^2).
    if updateRecord.active() {
      flushUpdateRecord()
    }
    ensurePositionsAndVelocitiesCached()
    guard range.startIndex >= 0,
          range.endIndex < system.parameters.atoms.count,
          range.count == rigidBody.parameters.atoms.count else {
      fatalError("Atom range was invalid.")
    }
    
    // Change the rigid body's external forces, positions, velocities.
    rigidBody.externalForces = Array(_externalForces[range])
    cachedState.positions![range].withContiguousStorageIfAvailable {
      rigidBody.setPositions($0)
    }
    cachedState.velocities![range].withContiguousStorageIfAvailable {
      rigidBody.setVelocities($0)
    }
  }
  
  /// Change the force field's internal state to match the specified rigid body.
  /// - parameter rigidBody: The rigid body to update the force field's state
  ///   with.
  /// - parameter range: The range of atoms in the forcefield that the rigid
  ///   body covers.
  public func `import`(from rigidBody: MM4RigidBody, range: Range<Int>) {
    // Cache the I/O accesses into OpenMM, otherwise this is O(n^2).
    ensurePositionsAndVelocitiesCached()
    updateRecord.externalForces = true
    updateRecord.positions = true
    updateRecord.velocities = true
    guard range.startIndex >= 0,
          range.endIndex < system.parameters.atoms.count,
          range.count == rigidBody.parameters.atoms.count else {
      fatalError("Atom range was invalid.")
    }
    
    // Change the force field's external forces, positions, and velocities.
    _externalForces.replaceSubrange(range, with: externalForces)
    cachedState.positions!.replaceSubrange(range, with: rigidBody.positions)
    cachedState.velocities!.replaceSubrange(range, with: rigidBody.velocities)
  }
}
