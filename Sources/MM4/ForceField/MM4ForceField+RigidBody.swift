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
    // rigid body's parameters.
    var atomCapacity: Int = 0
    var bondCapacity: Int = 0
    var angleCapacity: Int = 0
    var torsionCapacity: Int = 0
    var ringCapacity: Int = 0
    var exception13Capacity: Int = 0
    var exception14Capacity: Int = 0
    
    var ranges: [Range<UInt32>] = []
    for rigidBody in rigidBodies {
      let oldAtomCapacity = UInt32(atomCapacity)
      
      let parameters = rigidBody.parameters
      atomCapacity += parameters.atoms.count
      bondCapacity += parameters.bonds.indices.count
      angleCapacity += parameters.angles.indices.count
      torsionCapacity += parameters.torsions.indices.count
      ringCapacity += parameters.rings.indices.count
      exception13Capacity += parameters.nonbondedExceptions13.count
      exception14Capacity += parameters.nonbondedExceptions14.count
      
      ranges.append(oldAtomCapacity..<UInt32(atomCapacity))
    }
    
    // Elegantly handle the case where there are zero rigid bodies.
    var descriptor = MM4ParametersDescriptor()
    descriptor.atomicNumbers = []
    descriptor.bonds = []
    var parameters = try! MM4Parameters(descriptor: descriptor)
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
    self.init(parameters: parameters, rigidBodyRanges: ranges)
    for (index, rigidBody) in rigidBodies.enumerated() {
      self.import(from: rigidBody, index: index)
    }
  }
  
  /// Write the force field's internal state to the specified rigid body.
  ///
  /// This method does not overwrite the rigid body's anchors.
  ///
  /// If the exported rigid body has a different level of theory than when the
  /// force field was initialized, nothing happens. The discrepancy in level of
  /// theory is ignored, and the level currently used in the force field does
  /// not change.
  public func export(to rigidBody: inout MM4RigidBody, index: Int) {
    // Cache the I/O accesses into OpenMM, otherwise this is O(n^2).
    if updateRecord.active() {
      flushUpdateRecord()
    }
    ensurePositionsAndVelocitiesCached()
    guard index >= 0 && index < _rigidBodyRanges.count else {
      fatalError("Rigid body index out of bounds.")
    }
    
    // Change the rigid body's external forces, positions, velocities.
    let range = _rigidBodyRanges[index]
    _externalForces.withUnsafeBufferPointer { buffer in
      let baseAddress = buffer.baseAddress.unsafelyUnwrapped
      let pointer = UnsafeBufferPointer(
        start: baseAddress + Int(range.lowerBound),
        count: Int(range.upperBound - range.lowerBound))
      rigidBody.setExternalForces(pointer)
    }
    cachedState.positions!.withUnsafeBufferPointer { buffer in
      let baseAddress = buffer.baseAddress.unsafelyUnwrapped
      let pointer = UnsafeBufferPointer(
        start: baseAddress + Int(range.lowerBound),
        count: Int(range.upperBound - range.lowerBound))
      rigidBody.setPositions(pointer)
    }
    cachedState.velocities!.withUnsafeBufferPointer { buffer in
      let baseAddress = buffer.baseAddress.unsafelyUnwrapped
      let pointer = UnsafeBufferPointer(
        start: baseAddress + Int(range.lowerBound),
        count: Int(range.upperBound - range.lowerBound))
      rigidBody.setVelocities(pointer)
    }
  }
  
  /// Change the force field's internal state to match the specified rigid body.
  public func `import`(from rigidBody: MM4RigidBody, index: Int) {
    // Cache the I/O accesses into OpenMM, otherwise this is O(n^2).
    ensurePositionsAndVelocitiesCached()
    updateRecord.externalForces = true
    updateRecord.positions = true
    updateRecord.velocities = true
    guard index >= 0 && index < _rigidBodyRanges.count else {
      fatalError("Rigid body index out of bounds.")
    }
    
    // Assert that the imported anchor velocities are valid.
    rigidBody.storage.ensureAnchorVelocitiesValid()
    
    // Change the force field's external forces, positions, and velocities.
    let range = _rigidBodyRanges[index]
    _externalForces.withUnsafeMutableBufferPointer { buffer in
      let baseAddress = buffer.baseAddress.unsafelyUnwrapped
      let pointer = UnsafeMutableBufferPointer(
        start: baseAddress + Int(range.lowerBound),
        count: Int(range.upperBound - range.lowerBound))
      rigidBody.getExternalForces(pointer)
    }
    cachedState.positions!.withUnsafeMutableBufferPointer { buffer in
      let baseAddress = buffer.baseAddress.unsafelyUnwrapped
      let pointer = UnsafeMutableBufferPointer(
        start: baseAddress + Int(range.lowerBound),
        count: Int(range.upperBound - range.lowerBound))
      rigidBody.getPositions(pointer)
    }
    cachedState.velocities!.withUnsafeMutableBufferPointer { buffer in
      let baseAddress = buffer.baseAddress.unsafelyUnwrapped
      let pointer = UnsafeMutableBufferPointer(
        start: baseAddress + Int(range.lowerBound),
        count: Int(range.upperBound - range.lowerBound))
      rigidBody.getVelocities(pointer)
    }
    
    // Assert that the level of theory has not changed.
    guard _levelOfTheory[index] == rigidBody.parameters.levelOfTheory else {
      fatalError("Level of theory cannot change when importing a rigid body.")
    }
    
    // Check whether anchors have changed. If so, update the force field.
    let previousAnchors = _rigidBodyAnchors[index]
    if previousAnchors != rigidBody.anchors {
      updateRecord.anchors = true
      for anchor in previousAnchors {
        _anchors.remove(range.lowerBound + anchor)
      }
      for anchor in rigidBody.anchors {
        _anchors.insert(range.lowerBound + anchor)
      }
      _rigidBodyAnchors[index] = rigidBody.anchors
    }
  }
}
