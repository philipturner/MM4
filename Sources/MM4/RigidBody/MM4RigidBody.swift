//
//  MM4RigidBody.swift
//  MM4
//
//  Created by Philip Turner on 11/19/23.
//

public struct MM4RigidBodyDescriptor {
  /// Required.
  public var parameters: MM4Parameters?
  
  /// Required.
  ///
  /// After the rigid body is created, positions cannot be modified at the
  /// granularity of individual atoms. Doing so may deform the object, changing
  /// the inertial reference frame and violating some assumptions of rigid body
  /// mechanics.
  public var positions: [SIMD3<Float>]?
}

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
