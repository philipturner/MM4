//
//  MM4RigidBody.swift
//
//
//  Created by Philip Turner on 11/19/23.
//

import Numerics

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

public struct MM4RigidBodyDescriptor {
  /// Required. The number of protons in each atom's nucleus.
  public var atomicNumbers: [UInt8]?
  
  /// Required. Pairs of atom indices representing sigma bonds.
  public var bonds: [SIMD2<UInt32>]?
  
  /// Optional. The amount of mass (in amu) to redistribute from a substituent
  /// atom to each covalently bonded hydrogen.
  ///
  /// If not specified, the default is 1 amu.
  public var hydrogenMassRepartitioning: Float?
  
  /// Required. The position (in nanometers) of each atom's nucleus.
  public var positions: [SIMD3<Float>]?
  
  public init() {
    
  }
}

public struct MM4RigidBody {
  // MARK: - Properties summarizing topology
  
  /// Indices of atoms that should be treated as having infinite mass in a
  /// simulation.
  ///
  /// An anchor's velocity does not vary due to thermal energy. Angular
  /// momentum is constrained according to the number of anchors present.
  /// - 0 anchors: conserve linear and angular momentum around center of mass.
  /// - 1 anchor: conserve linear and angular momentum around anchor.
  /// - collinear anchors: conserve linear momentum around average of
  ///   anchors, constrain angular momentum to the shared axis.
  /// - anchors form plane: conserve linear momentum around average of anchors,
  ///   force angular momentum to zero.
  public var anchors: [UInt32] = []
  
  /// The number of protons in each atom's nucleus.
  ///
  /// Due to some implementation details, this cannot be changed. To modify an
  /// atom's atomic number, create a different rigid body with a bonding
  /// topology reflecting the new valence.
  public internal(set) var atomicNumbers: [UInt8]
  
  /// Pairs of atom indices representing sigma bonds.
  ///
  /// This property is immutable because a future implementation may cache some
  /// force field parameters.
  public internal(set) var bonds: [SIMD2<UInt32>]
  
  /// The amount of mass (in amu) to redistribute from a substituent atom to
  /// each covalently bonded hydrogen.
  ///
  /// This property is immutable because a future implementation may cache some
  /// force field parameters.
  public internal(set) var hydrogenMassRepartitioning: Float
  
  /// The mass (in amu) of each atom after hydrogen mass repartitioning.
  ///
  /// This property is immutable because a future implementation may cache some
  /// force field parameters.
  public internal(set) var masses: [Float]
  
  /// A convenient method for accessing the atom count.
  var atomCount: Int { atomicNumbers.count }
  
  /// The high-performance storage format for atom positions.
  var _positions: [MM4FloatVector]
  
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
  
  /// Thermal energy in zeptojoules, calculated as total energy minus the
  /// contributions from bulk linear and angular velocity.
  ///
  /// Kinetic energy contributions from anchors are omitted when calculating
  /// temperature.
  public var thermalEnergy: Double = 0.0
  
  // MARK: - Initializer
  
  public init(descriptor: MM4RigidBodyDescriptor) throws {
    // Ensure the required descriptor properties were set.
    guard let descriptorAtomicNumbers = descriptor.atomicNumbers,
          let descriptorBonds = descriptor.bonds,
          let positions = descriptor.positions else {
      fatalError("Descriptor did not have the required properties.")
    }
    
    self.atomicNumbers = descriptorAtomicNumbers
    self.bonds = descriptorBonds
    self.hydrogenMassRepartitioning = descriptor.hydrogenMassRepartitioning ?? 1.0
    self.masses = descriptorAtomicNumbers.map { atomicNumber in
      MM4MassParameters.global.mass(atomicNumber: atomicNumber)
    }
    
    for bond in bonds {
      let atomicNumber1 = atomicNumbers[Int(bond[0])]
      let atomicNumber2 = atomicNumbers[Int(bond[1])]
      guard atomicNumber1 == 1 || atomicNumber2 == 1 else {
        continue
      }
      if atomicNumber1 == atomicNumber2 {
        fatalError("Hydrogen cannot be bonded to another hydrogen.")
      }
      
      let hydrogen = (atomicNumber1 == 1) ? bond[0] : bond[1]
      let nonHydrogen = (atomicNumber1 == 1) ? bond[1] : bond[0]
      masses[Int(hydrogen)] += hydrogenMassRepartitioning
      masses[Int(nonHydrogen)] -= hydrogenMassRepartitioning
    }
    
    fatalError("Call into position setter.")
  }
}
