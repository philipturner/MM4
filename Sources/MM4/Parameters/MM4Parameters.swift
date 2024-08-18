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
  /// The default value includes stretch, bend, stretch-bend, stretch-stretch,
  /// and nonbonded forces.
  ///
  /// Disabling certain forces may reduce the execution time required to
  /// generate parameters. For example, if torsion forces are excluded, the
  /// array of torsion parameters will be empty. `nonbondedExceptions14` will
  /// also be an empty array.
  public var forces: MM4ForceOptions = [
    .bend,
    .nonbonded,
    .stretch,
    .stretchBend,
    .stretchStretch,
  ]
  
  /// Required. The factor to multiply hydrogen mass by.
  ///
  /// If not specified, the default value gives hydrogens ~2 amu of mass.
  ///
  /// Hydrogen mass repartitioning (HMR) is a technique to maximize simulation
  /// speed. It makes hydrogens heaver and makes non-hydrogen atoms lighter,
  /// decreasing the C-H stretching frequency. The smaller frequency corresponds
  /// to a larger vibration period. Therefore, larger time steps can be used.
  public var hydrogenMassScale: Float = 2
  
  /// Create a descriptor with the default properties.
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
    Self.checkForces(descriptor.forces)
    
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
    createElectronegativityEffectCorrections()
    createPartialCharges()
    removeBondLengths(forces: descriptor.forces)
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
  
  // Validate that the specified combination of forces is okay. Otherwise,
  // cause a runtime crash.
  static func checkForces(_ forces: MM4ForceOptions) {
    let includeTorsions =
    forces.contains(.torsion) ||
    forces.contains(.torsionBend) ||
    forces.contains(.torsionStretch)
    
    // Use a conservative metric to determine whether angles are included. If
    // torsions are included but angles aren't, there's still some torsion
    // cross-terms that depend on equilibrium angle.
    let includeAngles =
    forces.contains(.bend) ||
    forces.contains(.bendBend) ||
    forces.contains(.stretchBend) ||
    forces.contains(.stretchStretch) ||
    includeTorsions
    
    // None of the other bonded forces will be stable without bond stretching.
    let includeBonds =
    forces.contains(.stretch) ||
    includeAngles ||
    includeTorsions
    
    if includeBonds {
      guard forces.contains(.stretch) else {
        fatalError(
          "The specified forces cannot be evaluated without '.stretch'.")
      }
    }
    if includeAngles {
      guard forces.contains(.bend) else {
        fatalError(
          "The specified forces cannot be evaluated without '.bend'.")
      }
    }
    if includeTorsions {
      guard forces.contains(.torsion) else {
        fatalError(
          "The specified forces cannot be evaluated without '.torsion'.")
      }
    }
  }
}
