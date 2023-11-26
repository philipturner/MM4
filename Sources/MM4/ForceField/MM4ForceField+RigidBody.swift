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
  
  /// This method does not overwrite the rigid body's anchors or handles.
  ///
  /// If velocity is not the same across all anchors, there will be a fatal
  /// error. The same is true for unequal forces across the handles, or nonzero forces
  /// on the anchors.
  public func export(to rigidBody: inout MM4RigidBody, index: Int) {
    // Cache the I/O accesses into OpenMM, otherwise this is O(n^2).
    // Assert that the object's anchor velocities are valid.
    fatalError("Not implemented.")
  }
  
  public func `import`(from rigidBody: MM4RigidBody, index: Int) {
    // Cache the I/O accesses into OpenMM, otherwise this is O(n^2).
    // Check whether anchors or handles changed.
    // Ensure the object's anchor velocities are valid.
    fatalError("Not implemented.")
  }
}
