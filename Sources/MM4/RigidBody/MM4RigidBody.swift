//
//  MM4RigidBody.swift
//
//
//  Created by Philip Turner on 11/19/23.
//

import Numerics

public struct MM4RigidBodyDescriptor {
  /// Required. The number of protons in each atom's nucleus.
  public var atomicNumbers: [UInt8]?
  
  /// Required. Pairs of atom indices representing sigma bonds.
  public var bonds: [SIMD2<UInt32>]?
  
  /// Required. The amount of mass (in amu) to redistribute from a substituent
  /// atom to each covalently bonded hydrogen.
  ///
  /// The default is 1 amu.
  public var hydrogenMassRepartitioning: Float = 1.0
  
  /// Required. The position (in nanometers) of each atom's nucleus.
  public var positions: [SIMD3<Float>]?
  
  public init() {
    
  }
}

public struct MM4RigidBody {
  // This type declaration is getting very unwieldy.
  //
  // TODO: Wrap all of these properties in delegate objects to de-clutter the
  // central type declaration. This approach makes it easier to debug, so it is
  // well worth the extra effort.
  
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
  public var anchors: Set<UInt32> = []
  
  /// The force field parameters cached for this rigid body.
  public internal(set) var parameters: MM4Parameters
  
  // MARK: - Properties hidden from the public API
  
  /// A convenient method for accessing the atom count.
  var atomCount: Int
  
  /// A convenient method for accessing the atom vector count.
  var atomVectorCount: Int
  
  /// A convenient method for addressing the end of the list.
  var atomVectorMask: SIMDMask<MM4UInt32Vector.MaskStorage>
  
  /// The high-performance storage format for center of mass.
  var centerOfMass: MM4CenterOfMass = .init()
  
  /// The high-performance storage format for atom positions.
  var vPositions: [MM4FloatVector]
  
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
  public var velocity: SIMD3<Float> = .zero
  
  /// Kinetic energy contribution from organized mechanical energy (linear
  /// velocity, angular velocity). Contributions from anchors are omitted.
  public var freeKineticEnergy: Double {
    fatalError("Not implemented.")
  }
  
  /// Kinetic energy contribution from disorganized thermal energy.
  /// Contributions from anchors are omitted.
  public var thermalKineticEnergy: Double = 0.0
  
  // MARK: - Initializer
  
  public init(descriptor: MM4RigidBodyDescriptor) throws {
    // Ensure the required descriptor properties were set.
    guard let descriptorAtomicNumbers = descriptor.atomicNumbers,
          let descriptorBonds = descriptor.bonds,
          let descriptorPositions = descriptor.positions else {
      fatalError("Descriptor did not have the required properties.")
    }
    
    var desc = MM4ParametersDescriptor()
    desc.atomicNumbers = descriptorAtomicNumbers
    desc.bonds = descriptorBonds
    desc.hydrogenMassRepartitioning = descriptor.hydrogenMassRepartitioning
    self.parameters = try MM4Parameters(descriptor: desc)
    self.centerOfMass.mass = parameters.atoms.masses.reduce(Double(0)) {
      $0 + Double($1)
    }
    
    self.atomCount = descriptorAtomicNumbers.count
    self.atomVectorCount = atomCount + MM4VectorWidth - 1
    self.atomVectorCount /= MM4VectorWidth
    
    let lastVectorStart = UInt32((atomVectorCount - 1) * MM4VectorWidth)
    var lastVector: MM4UInt32Vector = .init(repeating: lastVectorStart)
    for lane in 0..<MM4VectorWidth {
      lastVector[lane] &+= UInt32(lane)
    }
    self.atomVectorMask = lastVector .>= UInt32(atomCount)
    
    self.vPositions = Array(unsafeUninitializedCapacity: 3 * atomVectorCount) {
      $0.initialize(repeating: .zero)
      $1 = 3 * descriptorAtomicNumbers.count
    }
    descriptorPositions.withUnsafeBufferPointer {
      setPositions($0)
    }
  }
}
