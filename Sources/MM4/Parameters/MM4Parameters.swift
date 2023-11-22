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
  /// The bonding topology must arrange atoms from a single rigid body
  /// contiguously in memory. Otherwise, there will be an error when creating
  /// the parameters.
  public var bonds: [SIMD2<UInt32>]?
  
  /// Required. The amount of mass (in amu) to redistribute from a substituent
  /// atom to each covalently bonded hydrogen.
  ///
  /// If not specified, the default is 1 amu.
  ///
  /// Each array element corresponds to a different rigid body. See the note in
  /// <doc:MM4ParametersDescriptor/bonds> about ordering the atoms from each
  /// rigid body.
  public var hydrogenMassRepartitioning: Float = 1.0
  
  // TODO: Force the user to specify the different rigid bodies, providing
  // flexibility to covalently connect two different "rigid bodies". Or, a
  // better API for nano-parts that are partially elastic.
  //
  // Either way, the existing functionality for auto-detecting rigid bodies
  // needs to be erased.
  
  public init() {
    
  }
}

// TODO: Make an `MM4Parameters` mutating function that appends parameters
// from another instance.

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
  var nonbondedExceptions13: [SIMD2<UInt32>] = []
  
  /// Atom pairs that have reduced nonbonded and electrostatic interactions.
  var nonbondedExceptions14: [SIMD2<UInt32>] = []
  
  /// Map from atoms to bonds that can be efficiently traversed.
  var atomsToBondsMap: [SIMD4<Int32>] = []
  
  /// Map from atoms to connected atoms that can be efficiently traversed.
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
    let descriptorHMR = descriptor.hydrogenMassRepartitioning
    try createAtomCodes()
    createMasses(hydrogenMassRepartitioning: descriptorHMR)
    createVectorPadding()
    createNonbondedParameters(hydrogenMassRepartitioning: descriptorHMR)
    createNonbondedExceptions()
    
    // Bond Parameters
    try createBondParameters()
    try createAngleParameters()
    try createTorsionParameters()
    createPartialCharges()
  }
}
