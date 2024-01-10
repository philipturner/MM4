//
//  MM4Parameters.swift
//  MM4
//
//  Created by Philip Turner on 9/10/23.
//

/// A configuration for a set of force field parameters.
public struct MM4ParametersDescriptor {
  /// Required. The number of protons in each atom's nucleus.
  public var atomicNumbers: [UInt8]?
  
  /// Required. Pairs of atom indices representing sigma bonds.
  public var bonds: [SIMD2<UInt32>]?
  
  /// Required. The forces to assign parameters for.
  ///
  /// The default value includes all available forces.
  ///
  /// Disabling certain forces may reduce the execution time required to
  /// generate parameters. For example, if torsion forces are excluded, the
  /// array of torsion parameters will be empty. `nonbondedExceptions14` will
  /// also be an empty array.
  public var forces: MM4ForceOptions = [
    .bend,
    .bendBend,
    .nonbonded,
    .stretch,
    .stretchBend,
    .stretchStretch,
    .torsion,
    .torsionBend,
    .torsionStretch,
  ]
  
  /// Required. The factor to multiply hydrogen mass by.
  ///
  /// If not specified, the default value gives hydrogens ~2 amu of mass.
  ///
  /// During hydrogen mass repartitioning, mass is added to hydrogens and
  /// removed from atoms covalently bonded to hydrogens. The resulting structure
  /// has the same total mass as before the transformation.
  public var hydrogenMassScale: Float = 2
  
  public init() {
    
  }
}

/// A set of force field parameters.
public struct MM4Parameters {
  /// Parameters for one atom.
  public var atoms: MM4Atoms = MM4Atoms()
  
  /// Parameters for a group of 2 atoms.
  public var bonds: MM4Bonds = MM4Bonds()
  
  /// Parameters for a group of 3 atoms.
  public var angles: MM4Angles = MM4Angles()
  
  /// Parameters for a group of 4 atoms.
  public var torsions: MM4Torsions = MM4Torsions()
  
  /// Parameters for a group of 5 atoms.
  public var rings: MM4Rings = MM4Rings()
  
  /// Atom pairs to be excluded from nonbonded and electrostatic interactions.
  public var nonbondedExceptions13: [SIMD2<UInt32>] = []
  
  /// Atom pairs that have reduced nonbonded and electrostatic interactions.
  public var nonbondedExceptions14: [SIMD2<UInt32>] = []
  
  /// Map from atoms to bonds that requires bounds checking.
  var atomsToBondsMap: [SIMD4<Int32>] = []
  
  /// Map from atoms to connected atoms that requires bounds checking.
  var atomsToAtomsMap: [SIMD4<Int32>] = []
  
  /// Create a set of parameters using the specified configuration.
  public init(descriptor: MM4ParametersDescriptor) throws {
    // Ensure the required descriptor properties were set.
    guard let descriptorAtomicNumbers = descriptor.atomicNumbers,
          let descriptorBonds = descriptor.bonds else {
      fatalError("Descriptor did not have the required properties.")
    }
    
    // Set the properties for conveniently iterating over the atoms.
    // Behavior should be well-defined when the atom count is zero.
    atoms.atomicNumbers = descriptorAtomicNumbers
    atoms.count = descriptorAtomicNumbers.count
    atoms.indices = 0..<descriptorAtomicNumbers.count
    bonds.indices = descriptorBonds.map { bond in
      return SIMD2(bond.min(), bond.max())
    }
    
    // Topology
    try createAtomsToBondsMap()
    try createAtomsToAtomsMap()
    try createTopology(forces: descriptor.forces)
    try createCenterTypes()
    
    // Atom Parameters
    try createAtomCodes()
    createMasses(hydrogenMassScale: descriptor.hydrogenMassScale)
    createNonbondedParameters(descriptor: descriptor)
    createNonbondedExceptions(forces: descriptor.forces)
    
    // Bond Parameters
    try createBondParameters(forces: descriptor.forces)
    try createAngleParameters(forces: descriptor.forces)
    try createTorsionParameters(forces: descriptor.forces)
    createElectronegativityEffectCorrections()
    createPartialCharges()
  }
  
  public mutating func append(contentsOf other: Self) {
    let atomOffset = UInt32(atoms.count)
    let bondOffset = UInt32(bonds.indices.count)
    atoms.append(contentsOf: other.atoms, atomOffset: atomOffset)
    bonds.append(contentsOf: other.bonds, atomOffset: atomOffset)
    angles.append(contentsOf: other.angles, atomOffset: atomOffset)
    torsions.append(contentsOf: other.torsions, atomOffset: atomOffset)
    rings.append(contentsOf: other.rings, atomOffset: atomOffset)
    
    nonbondedExceptions13 += other.nonbondedExceptions13.map {
      $0 &+ atomOffset
    }
    nonbondedExceptions14 += other.nonbondedExceptions14.map {
      $0 &+ atomOffset
    }
    atomsToBondsMap += other.atomsToBondsMap.map {
      let modified = $0 &+ Int32(truncatingIfNeeded: bondOffset)
      return $0.replacing(with: modified, where: $0 .>= 0)
    }
    atomsToAtomsMap += other.atomsToAtomsMap.map {
      let modified = $0 &+ Int32(truncatingIfNeeded: atomOffset)
      return $0.replacing(with: modified, where: $0 .>= 0)
    }
  }
}
