//
//  MM4Parameters.swift
//  MM4
//
//  Created by Philip Turner on 9/10/23.
//

import Foundation

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
    let checkpoint0 = Date() // 0.0 ms 0.0 ms
    atoms.atomicNumbers = descriptorAtomicNumbers
    atoms.count = descriptorAtomicNumbers.count
    atoms.indices = 0..<descriptorAtomicNumbers.count
    bonds.indices = descriptorBonds.map { bond in
      return SIMD2(bond.min(), bond.max())
    }
    
    // Topology
    let checkpoint1 = Date() // 0.2 ms 0.2 ms
    try createAtomsToBondsMap()
    try createAtomsToAtomsMap()
    let checkpoint2 = Date() // 5.2 ms 8.2 ms
    try createTopology(forces: descriptor.forces)
    let checkpoint3 = Date() // 0.0 ms 0.0 ms
    try createCenterTypes()
    
    // Atom Parameters
    let checkpoint4 = Date() // 0.0 ms 0.0 ms
    try createAtomCodes()
    createMasses(hydrogenMassScale: descriptor.hydrogenMassScale)
    let checkpoint5 = Date() // 0.0 ms 0.0 ms
    createNonbondedParameters(descriptor: descriptor)
    let checkpoint6 = Date() // 1.7 ms 2.2 ms
    createNonbondedExceptions(forces: descriptor.forces)
    
    // Bond Parameters
    let checkpoint7 = Date() // 0.0 ms 0.7 ms
    try createBondParameters(forces: descriptor.forces)
    let checkpoint8 = Date() // 0.5 ms 3.2 ms
    try createAngleParameters(forces: descriptor.forces)
    let checkpoint9 = Date() // 1.0 ms 10.1 ms
    try createTorsionParameters(forces: descriptor.forces)
    let checkpoint10 = Date() // 0.7 ms 1.0 ms
    createElectronegativityEffectCorrections()
    let checkpoint11 = Date() //  0.0 ms 0.0 ms
    createPartialCharges()
    let checkpoint12 = Date()
    
    if atoms.count == 1514 {
      print("execution time:")
      let checkpoints = [
        checkpoint0, checkpoint1, checkpoint2, checkpoint3, checkpoint4,
        checkpoint5, checkpoint6, checkpoint7, checkpoint8, checkpoint9,
        checkpoint10, checkpoint11, checkpoint12
      ]
      let absoluteTimes = checkpoints.map { $0.timeIntervalSince1970 }
      let relativeTimes = (0..<checkpoints.count - 1).map { i in
        absoluteTimes[i + 1] - absoluteTimes[i]
      }
      for i in relativeTimes.indices {
        let ms = String(format: "%.1f", relativeTimes[i] * 1e3)
        print("- time interval \(i): \(ms) ms")
      }
    }
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
