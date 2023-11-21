//
//  MM4ForceField+RigidBody.swift
//
//
//  Created by Philip Turner on 11/20/23.
//

extension MM4ForceField {
  public convenience init(rigidBodies: [MM4RigidBody]) {
    fatalError("Not implemented.")
  }
  
  public func export(to rigidBody: inout MM4RigidBody, index: Int) {
    // Cache the I/O accesses into OpenMM, otherwise this is O(n^2).
    fatalError("Not implemented.")
  }
  
  public func `import`(from rigidBody: MM4RigidBody, index: Int) {
    // Cache the I/O accesses into OpenMM, otherwise this is O(n^2).
    fatalError("Not implemented.")
  }
}
