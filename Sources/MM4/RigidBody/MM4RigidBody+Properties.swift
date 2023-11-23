//
//  MM4RigidBody+Properties.swift
//
//
//  Created by Philip Turner on 11/22/23.
//

import Numerics

// MARK: - Energy

/// A data structure wrapping a rigid body's energy.
public struct MM4RigidBodyEnergy {
  /// The rigid body's kinetic energy.
  public var kinetic: MM4RigidBodyKineticEnergy = .init()
  
  init() { }
}

/// A data structure wrapping a rigid body's kinetic energy.
public struct MM4RigidBodyKineticEnergy {
  weak var storage: MM4RigidBodyStorage!
  
  init() { }
  
  /// Kinetic energy contribution from organized mechanical energy (linear
  /// velocity, angular velocity). Contributions from anchors are omitted.
  public var free: Double {
    fatalError("Not implemented.")
  }
  
  /// Kinetic energy contribution from disorganized thermal energy.
  /// Contributions from anchors are omitted.
  public var thermal: Double {
    get { fatalError("Not implemented.") }
    set { fatalError("Not implemented.") }
  }
}

// MARK: - Velocity

extension MM4RigidBody {
  /// If the angular velocity is nonzero, the number of anchors cannot exceed 1.
  /// When importing velocities, if the number of anchors exceeds 1, the angular
  /// velocity is set to zero.
  public var angularVelocity: Quaternion<Float> {
    get { fatalError("Not implemented.") }
    set { fatalError("Not implemented.") }
  }
  
  /// If the number of anchors exceeds 0, external force has no effect.
  public var externalForce: SIMD3<Float> {
    get { fatalError("Not implemented.") }
    set { fatalError("Not implemented.") }
  }
  
  /// If the number of anchors exceeds 1, external torque has no effect.
  ///
  /// Right now, external torque must be zero when simulating in the
  /// `.molecularDynamics` level of theory.
  public var externalTorque: Quaternion<Float> {
    get { fatalError("Not implemented.") }
    set { fatalError("Not implemented.") }
  }
  
  /// Every anchor's velocity is set to the rigid body's linear velocity.
  /// When importing velocities, all anchors must have the same velocity.
  public var velocity: SIMD3<Float> {
    get { fatalError("Not implemented.") }
    set { fatalError("Not implemented.") }
  }
}

// MARK: - Other Properties

extension MM4RigidBody {
  /// Indices of atoms that should be treated as having infinite mass.
  public var anchors: Set<UInt32> {
    _read { fatalError("Not implemented.") }
    _modify { fatalError("Not implemented.") }
  }
  
  /// The number of protons in each atom's nucleus.
  public var atomicNumbers: [UInt8] {
    parameters.atoms.atomicNumbers
  }
  
  /// Pairs of atom indices representing sigma bonds.
  public var bonds: [SIMD2<UInt32>] {
    parameters.bonds.indices
  }
  
  /// The object's total mass (in amu).
  public var mass: Double {
    fatalError("Not implemented.")
  }
  
  /// The mass (in amu) of each atom after hydrogen mass repartitioning.
  public var masses: [Float] {
    parameters.atoms.masses
  }
}
