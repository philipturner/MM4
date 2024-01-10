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
  
  /// Optional.
  ///
  /// The default value is zero for every atom.
  ///
  /// When the rigid body is created, atom velocities are decomposed into
  /// thermal velocities and bulk momenta. The thermal velocities cannot be
  /// modified, as that would require the momenta to be recomputed.
  public var velocities: [SIMD3<Float>]?
  
  public init() {
    
  }
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
  public init(descriptor: MM4RigidBodyDescriptor) throws {
    guard let parameters = descriptor.parameters else {
      fatalError("Parameters were not specified.")
    }
    self.parameters = parameters
    self.storage = try MM4RigidBodyStorage(descriptor: descriptor)
  }
  
  /// Ensures copy-on-write semantics. This is not exposed to the public API.
  ///
  /// > WARNING: Call this before every mutating function.
  mutating func ensureUniquelyReferenced() {
    if !isKnownUniquelyReferenced(&storage) {
      storage = MM4RigidBodyStorage(copying: storage)
    }
  }
}
