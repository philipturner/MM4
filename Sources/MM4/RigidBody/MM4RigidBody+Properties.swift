//
//  MM4RigidBody+Properties.swift
//
//
//  Created by Philip Turner on 11/22/23.
//

extension MM4RigidBody {
  /// The number of protons in each atom's nucleus.
  public var atomicNumbers: [UInt8] {
    parameters.atoms.atomicNumbers
  }
  
  /// Pairs of atom indices representing sigma bonds.
  public var bonds: [SIMD2<UInt32>] {
    parameters.bonds.indices
  }
  
  /// The amount of mass (in amu) to redistribute from a substituent atom to
  /// each covalently bonded hydrogen.
  public var hydrogenMassRepartitioning: Float {
    parameters.hydrogenMassRepartitioning[0]
  }
  
  /// The object's total mass (in amu).
  public var mass: Double {
    centerOfMass.mass
  }
  
  /// The mass (in amu) of each atom after hydrogen mass repartitioning.
  public var masses: [Float] {
    parameters.atoms.masses
  }
}
