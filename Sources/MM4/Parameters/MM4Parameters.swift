//
//  MM4Parameters.swift
//
//
//  Created by Philip Turner on 9/10/23.
//

/// A configuration for a set of force field parameters.
public class MM4ParametersDescriptor {
  /// Required. The number of protons in the atom's nucleus.
  public var atomicNumbers: [UInt8]?
  
  /// Required. Pairs of atom indices representing (potentially multiple)
  /// covalent bonds.
  ///
  /// The bonding topology must arrange atoms from a single rigid body
  /// contiguously in memory. Otherwise, there will be an error when creating
  /// the parameters.
  public var bonds: [SIMD2<UInt32>]?
  
  /// Optional. The bond order for each covalent bond, which may be fractional.
  ///
  /// If not specified, all covalent bonds are treated as sigma bonds.
  public var bondOrders: [Float]?
  
  /// Required. The amount of mass (in amu) to redistribute from a substituent
  /// atom to each covalently bonded hydrogen.
  ///
  /// The default is 1 amu.
  public var hydrogenMassRepartitioning: Double = 1.0
  
  public init() {
    
  }
}

// Owning classes like 'MM4Atoms' have immutable properties, but the structs
// themselves (e.g. 'MM4AtomParameters' are mutable. This creates
// consistent behavior with the C API, where a copy of the struct's data is
// returned. It is possible to modify the copy, which will not affect the source
// of truth.
//
// Another note, both the Swift and C APIs might have function call overhead to
// access a single element of any parameters array.

/// A set of force field parameters.
public class MM4Parameters {
  /// Parameters for one atom.
  public internal(set) var atoms: MM4Atoms = MM4Atoms()
  
  /// Parameters for a group of 2 atoms.
  public internal(set) var bonds: MM4Bonds = MM4Bonds()
  
  /// Parameters for a group of 3 atoms.
  public internal(set) var angles: MM4Angles = MM4Angles()
  
  /// Parameters for a group of 4 atoms.
  public internal(set) var torsions: MM4Torsions = MM4Torsions()
  
  /// Parameters for a group of 5 atoms.
  public internal(set) var rings: MM4Rings = MM4Rings()
  
  /// The amount of mass (in amu) redistributed from a substituent atom to each
  /// covalently bonded hydrogen.
  var hydrogenMassRepartitioning: Float = -1
  
  /// Atom pairs to be excluded from nonbonded and electrostatic interactions.
  var nonbondedExceptions13: [SIMD2<Int32>] = []
  
  /// Atom pairs that have reduced nonbonded and electrostatic interactions.
  var nonbondedExceptions14: [SIMD2<Int32>] = []
  
  /// Map from atoms to bonds that can be efficiently traversed.
  var atomsToBondsMap: UnsafeMutablePointer<SIMD4<Int32>>
  
  /// Map from bonds to atoms that can be efficiently traversed.
  var bondsToAtomsMap: UnsafeMutablePointer<SIMD2<Int32>>
  
  /// Map from atoms to connected atoms that can be efficiently traversed.
  var atomsToAtomsMap: UnsafeMutablePointer<SIMD4<Int32>>
  
  /// Grouping of atoms into rigid bodies.
  var rigidBodies: [Range<Int32>] = []
  
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
      fatalError("Descriptor did not have any atomic numbers or bonds.")
    }
    
    // Set the properties for conveniently iterating over the atoms.
    atoms.count = descriptorAtomicNumbers.count
    atoms.indices = 0..<descriptorAtomicNumbers.count
    
    // Check the bond orders.
    if let bondOrders = descriptor.bondOrders {
      precondition(
        bondOrders.allSatisfy { $0 == 1 }, "Pi bonds are not supported yet.")
    }
    
    // Compile the bonds into a map.
    bondsToAtomsMap = .allocate(capacity: descriptorBonds.count + 1)
    bondsToAtomsMap += 1
    bondsToAtomsMap[-1] = SIMD2(repeating: -1)
    for bondID in 0..<descriptorBonds.count {
      var bond = descriptorBonds[bondID]
      
      // Sort the indices in the bond, so the lower appears first.
      bond = SIMD2(bond.min(), bond.max())
      bondsToAtomsMap[bondID] = SIMD2(truncatingIfNeeded: bond)
      bonds.indices.append(SIMD2(truncatingIfNeeded: bond))
    }
    
    atoms.atomicNumbers = descriptorAtomicNumbers
    atomsToBondsMap = .allocate(capacity: atoms.atomicNumbers.count + 1)
    atomsToBondsMap += 1
    atomsToBondsMap[-1] = SIMD4(repeating: -1)
    for atomID in 0..<atoms.atomicNumbers.count {
      atomsToBondsMap[atomID] = SIMD4(repeating: -1)
    }
    
    for bondID in 0..<bonds.indices.count {
      let bond = bondsToAtomsMap[bondID]
      for j in 0..<2 {
        let atomID = Int(bond[j])
        var map = atomsToBondsMap[atomID]
        var succeeded = false
        for k in 0..<4 {
          if map[k] == -1 {
            map[k] = Int32(bondID)
            succeeded = true
            break
          }
        }
        if !succeeded {
          fatalError("An atom had more than 4 bonds.")
        }
        atomsToBondsMap[atomID] = map
      }
    }
    
    atomsToAtomsMap = .allocate(capacity: atoms.atomicNumbers.count + 1)
    atomsToAtomsMap += 1
    atomsToAtomsMap[-1] = SIMD4(repeating: -1)
    for atomID in 0..<atoms.atomicNumbers.count {
      let bondsMap = atomsToBondsMap[atomID]
      var atomsMap = SIMD4<Int32>(repeating: -1)
      for lane in 0..<4 {
        atomsMap[lane] = other(atomID: atomID, bondID: bondsMap[lane])
      }
      atomsToAtomsMap[atomID] = atomsMap
    }
    
    // Topology
    createRigidBodies()
    createTopology()
    createAtomCodes()
    createCenterTypes()
    
    // Per-Atom Parameters
    hydrogenMassRepartitioning = Float(descriptor.hydrogenMassRepartitioning)
    createMasses()
    createNonbondedParameters()
    createNonbondedExceptions()
    
    // Per-Bond Parameters
    createBondParameters()
    addElectrostaticCorrections()
    createPartialCharges()
    
    // Other Parameters
    createAngleParameters()
    createTorsionParameters()
  }
  
  deinit {
    (atomsToBondsMap - 1).deallocate()
    (bondsToAtomsMap - 1).deallocate()
    (atomsToAtomsMap - 1).deallocate()
  }
}

