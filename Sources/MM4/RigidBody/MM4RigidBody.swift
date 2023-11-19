//
//  MM4RigidBody.swift
//
//
//  Created by Philip Turner on 11/19/23.
//

import RealModule
import QuaternionModule

/*
 - MM4RigidBody structure that encapsulates all the rigid body calculations from thermalize(); MolecularRenderer will extend it with Lattice/Solid initializer that performs Morton reordering: rigid bodies are detached from specific simulators
 - MM4ForcesField initializer directly accepting rigid bodies, API to import/export to one rigid body at a time (but cache the I/O accesses into OpenMM)
   - export(to: inout RigidBody, index: Int)
   - import(from: RigidBody, index: Int)
 
 /// Means to extract atom positions and atomic numbers.
 public struct RigidBody {
   // atomicNumbers
   // bonds
   // relativePositions
   // relativeVelocities
   // - stores positions/velocities in a local coordinate space, internally
   //   projects to the global space using rigid body transforms
   
   // centerOfMass, rotationalInertia
   // position, rotation
   // velocity, angularVelocity
   
   // init(solid:) -> materializes Solid topology
   // init(lattice:) -> materializes Lattice topology
   // private init(...) -> reduces code duplication between Lattice and Solid
 }
 */

public struct MM4RigidBody {
  // MARK: - Properties summarizing topology
  
  public var anchors: [UInt32]
  public var atomicNumbers: [UInt8]
  public var bonds: [SIMD2<UInt32>]
  public var hydrogenMassRepartitioning: Float = 1.0
  public var positions: [SIMD3<Float>]
  
  // MARK: - Properties summarizing velocity
  
  /// If the angular velocity is nonzero, the number of anchors cannot exceed 1.
  /// When importing velocities, if the number of anchors exceeds 1, the angular
  /// velocity is set to zero.
  public var angularVelocity: Quaternion<Float> = .zero
  
  /// If the number of anchors exceeds 0, external force has no effect.
  public var externalForce: SIMD3<Float> = .zero
  
  /// If the number of anchors exceeds 1, external torque has no effect.
  ///
  /// Right now, external torque must be zero when simulating in the
  /// `.molecularDynamics` level of theory.
  public var externalTorque: Quaternion<Float> = .zero
  
  /// Every anchor's velocity is set to the rigid body's linear velocity.
  /// When importing velocities, all anchors must have the same velocity.
  public var linearVelocity: SIMD3<Float> = .zero
  
  /// Kinetic energy contributions from anchors are omitted when calculating
  /// temperature.
  public var temperature: Double = 0.0 // E = 3/2 NkT
  
  // MARK: - Computed properties
  
  public var masses: [Float] {
    get { fatalError() }
  }
  public var velocities: [Float] {
    get { fatalError() }
    set { fatalError() }
  }
}
