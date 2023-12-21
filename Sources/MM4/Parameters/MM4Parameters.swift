//
//  MM4Parameters.swift
//
//
//  Created by Philip Turner on 9/10/23.
//

/// A configuration for a set of force field parameters.
public struct MM4ParametersDescriptor {
  /// Required. The number of protons in each atom's nucleus.
  public var atomicNumbers: [UInt8]?
  
  /// Required. Pairs of atom indices representing sigma bonds.
  ///
  /// If the level of theory is `.rigidBodyMechanics`, then `bonds` must contain
  /// exclusively the bonds between hydrogen and non-hydrogen atoms. Such bonds
  /// are needed to properly implement hydrogen mass repartitioning and hydrogen
  /// reduction factors for vdW forces. MM4 does not parameterize any bonds that
  /// give a partial charge to hydrogen. Therefore, the generated partial
  /// charges are zero (even for bonds that are obviously polar). This rule
  /// mistreats the N-H bonds disallowed by `.molecularMechanics`, but accepted
  /// by rigid body mechanics. In addition, the MM4 default hydrogen reduction
  /// factor of 0.94 is assigned to every X-H bond except Si-H and Ge-H. You can
  /// override such parameters after the `MM4Parameters` object initializes.
  ///
  /// With rigid body mechanics, you may not specify any bonds between a pair of
  /// non-hydrogen atoms (e.g. C-C). This design choice removes an ambiguity
  /// about the source of truth for partial charges. To acquire accurate partial
  /// charges, you must generate them externally. Next, copy the charges
  /// into the `MM4Parameters` object after it initializes. For example, create
  /// a separate `MM4Parameters` with molecular mechanics and transfer over
  /// [`.atoms.parameters`](<doc:MM4Atoms/parameters>). Or, run the structure
  /// through the xTB program and record the partial charges.
  public var bonds: [SIMD2<UInt32>]?
  
  /// Required. The factor to multiply hydrogen mass by.
  ///
  /// If not specified, the default value gives hydrogens ~2 amu of mass.
  ///
  /// During hydrogen mass repartitioning, mass is added to hydrogens and
  /// removed from atoms covalently bonded to hydrogens. The resulting structure
  /// has the same total mass as before the transformation.
  public var hydrogenMassScale: Float = 2
  
  /// Required. The level of theory used for simulation.
  ///
  /// The default is `.molecularMechanics`.
  public var levelOfTheory: MM4LevelOfTheory = .molecularMechanics
  
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
  
  /// The level of theory used for simulation.
  public var levelOfTheory: MM4LevelOfTheory
  
  /// Atom pairs to be excluded from nonbonded and electrostatic interactions.
  var nonbondedExceptions13: [SIMD2<UInt32>] = []
  
  /// Atom pairs that have reduced nonbonded and electrostatic interactions.
  var nonbondedExceptions14: [SIMD2<UInt32>] = []
  
  /// Map from atoms to bonds that requires bounds checking.
  var atomsToBondsMap: [SIMD4<Int32>] = []
  
  /// Map from atoms to connected atoms that requires bounds checking.
  var atomsToAtomsMap: [SIMD4<Int32>] = []
  
  /// Create a set of parameters using the specified configuration.
  ///
  /// - throws: An error if there wasn't a parameter for a certain atom pair, or
  ///   the descriptor was invalid.
  ///
  /// This is a throwing initializer, allowing it to be used as a validation
  /// mechanism for structures that are potentially invalid. Enter the structure
  /// into the initializer, then try a different one if it fails. This removes
  /// the need to reimplement some of that logic in an automated search program.
  public init(descriptor: MM4ParametersDescriptor) throws {
    // Ensure the required descriptor properties were set.
    guard let descriptorAtomicNumbers = descriptor.atomicNumbers,
          let descriptorBonds = descriptor.bonds else {
      fatalError("Descriptor did not have the required properties.")
    }
    guard case .molecularMechanics = descriptor.levelOfTheory else {
      fatalError("Unsupported level of theory.")
    }
    self.levelOfTheory = descriptor.levelOfTheory
    
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
    try createTopology()
    try createCenterTypes()
    
    // Atom Parameters
    try createAtomCodes()
    createMasses(hydrogenMassScale: descriptor.hydrogenMassScale)
    createNonbondedParameters(hydrogenMassScale: descriptor.hydrogenMassScale)
    createNonbondedExceptions()
    
    // Bond Parameters
    try createBondParameters()
    try createAngleParameters()
    try createTorsionParameters()
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
      $0 &+ Int32(truncatingIfNeeded: bondOffset)
    }
    atomsToAtomsMap += other.atomsToAtomsMap.map {
      $0 &+ Int32(truncatingIfNeeded: atomOffset)
    }
  }
}
