//
//  MM4RigidBody.swift
//  MM4
//
//  Created by Philip Turner on 11/19/23.
//

/// An enclosed group of covalently bonded atoms.
///
/// `MM4RigidBody` is an API to efficiently compute basic properties in rigid
/// body mechanics. It stores atoms in an internal representation optimized for
/// the vector units on modern CPUs.
public struct MM4RigidBody {
  /// The force field parameters cached for this rigid body.
  public let parameters: MM4Parameters
  
  /// The backing storage object.
  var storage: MM4RigidBodyStorage
  
  /// Create a rigid body using the specified configuration.
  public init(parameters: MM4Parameters) {
    self.parameters = parameters
    self.storage = MM4RigidBodyStorage(parameters: parameters)
  }
}
