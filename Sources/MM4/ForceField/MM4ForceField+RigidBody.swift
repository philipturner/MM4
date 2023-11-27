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
    var ranges: [Range<UInt32>] = []
    for rigidBody in rigidBodies {
      let oldAtomCapacity = UInt32(atomCapacity)
      
      let parameters = rigidBody.parameters
      atomCapacity += parameters.atoms.count
      bondCapacity += parameters.bonds.indices.count
      angleCapacity += parameters.angles.indices.count
      torsionCapacity += parameters.torsions.indices.count
      ringCapacity += parameters.rings.indices.count
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
    parameters.nonbondedExceptions13.reserveCapacity(atomCapacity)
    parameters.nonbondedExceptions14.reserveCapacity(atomCapacity)
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
  /// This method does not overwrite the rigid body's anchors or handles.
  ///
  /// If velocity is not the same across all anchors, there will be a fatal
  /// error. If the handles have unequal force or the anchors have nonzero force,
  /// the exported object's external force will be set to `.nan`.
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
    
    // Change the rigid body's positions and velocities.
    let range = _rigidBodyRanges[index]
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
    
    // Assert that the exported anchor velocities are valid.
    rigidBody.storage.ensureAnchorVelocitiesValid()
    
    // Check whether any anchors or handles changed. If so, update the object's
    // external force.
    let previousHandles = _rigidBodyHandles[index]
    if _rigidBodyAnchors[index] != rigidBody.anchors ||
        _rigidBodyExternalForces[index] != rigidBody.externalForce ||
        previousHandles != rigidBody.handles {
      
      if previousHandles.count > 0 {
        let atomID = range.lowerBound + previousHandles.first!
        let force = _externalForces[Int(atomID)]
        rigidBody.externalForce = force * Float(previousHandles.count)
      } else {
        rigidBody.externalForce = .zero
      }
      
      if any(_rigidBodyExternalForces[index] .!= .zero) {
        for anchor in rigidBody.anchors {
          if previousHandles.contains(anchor) {
            rigidBody.externalForce = .init(repeating: .nan)
          }
        }
        if rigidBody.handles != previousHandles {
          rigidBody.externalForce = .init(repeating: .nan)
        }
      }
    }
  }
  
  /// Change the force field's internal state to match the specified rigid body.
  public func `import`(from rigidBody: MM4RigidBody, index: Int) {
    // Cache the I/O accesses into OpenMM, otherwise this is O(n^2).
    ensurePositionsAndVelocitiesCached()
    updateRecord.positions = true
    updateRecord.velocities = true
    guard index >= 0 && index < _rigidBodyRanges.count else {
      fatalError("Rigid body index out of bounds.")
    }
    
    // Assert that the imported anchor velocities are valid.
    rigidBody.storage.ensureAnchorVelocitiesValid()
    
    // Change the force field's positions and velocities.
    let range = _rigidBodyRanges[index]
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
    
    // Check whether anchors or handles changed. If so, update the force field.
    let previousAnchors = _rigidBodyAnchors[index]
    if previousAnchors != rigidBody.anchors {
      for anchor in previousAnchors {
        _anchors.remove(range.lowerBound + anchor)
      }
      for anchor in rigidBody.anchors {
        _anchors.insert(range.lowerBound + anchor)
      }
      _rigidBodyAnchors[index] = rigidBody.anchors
    }
    do {
      let force = rigidBody.externalForce
      if force.x.isNaN || force.y.isNaN || force.z.isNaN {
        fatalError("Imported a rigid body with a NAN external forces. Perhaps the external force was invalid when you last exported data into the rigid body.")
      }
    }
    if _rigidBodyExternalForces[index] != rigidBody.externalForce ||
        _rigidBodyHandles[index] != rigidBody.handles {
      updateRecord.externalForces = true
      for atomID in range {
        _externalForces[Int(atomID)] = .zero
      }
      
      // If there are no handles, avoid a division by zero.
      var handleCount = Float(rigidBody.handles.count)
      if handleCount == 0 { handleCount = .greatestFiniteMagnitude }
      let perParticleForce = rigidBody.externalForce / handleCount
      for handle in rigidBody.handles {
        _externalForces[Int(handle)] = perParticleForce
      }
      _rigidBodyHandles[index] = rigidBody.handles
    }
  }
}
