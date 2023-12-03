//
//  MM4Parameters+Atoms.swift
//
//
//  Created by Philip Turner on 10/7/23.
//

/// Parameters for one atom.
public struct MM4Atoms {
  /// The number of protons in each atom's nucleus.
  public var atomicNumbers: [UInt8] = []
  
  /// The center type used to assign different parameters.
  public var centerTypes: [MM4CenterType?] = []
  
  /// The MM4 code for each atom in the system.
  public var codes: [MM4AtomCode] = []
  
  /// Convenient property for iterating over the atoms.
  public var count: Int = 0
  
  /// Convenient property for iterating over the atoms.
  public var indices: Range<Int> = 0..<0
  
  /// The mass (in yoctograms) of each atom after hydrogen mass repartitioning.
  public var masses: [Float] = []
  
  /// Each value corresponds to the atom at the same array index.
  public var parameters: [MM4AtomParameters] = []
  
  /// The smallest ring this is involved in.
  public var ringTypes: [UInt8] = []
  
  mutating func append(contentsOf other: Self, atomOffset: UInt32) {
    self.atomicNumbers += other.atomicNumbers
    self.centerTypes += other.centerTypes
    self.count += other.count
    self.indices = 0..<self.count
    self.masses += other.masses
    self.parameters += other.parameters
    self.ringTypes += other.ringTypes
  }
  
  mutating func reserveCapacity(_ minimumCapacity: Int) {
    atomicNumbers.reserveCapacity(minimumCapacity)
    centerTypes.reserveCapacity(minimumCapacity)
    codes.reserveCapacity(minimumCapacity)
    masses.reserveCapacity(minimumCapacity)
    parameters.reserveCapacity(minimumCapacity)
    ringTypes.reserveCapacity(minimumCapacity)
  }
}

/// MM4 codes for an element or an atom in a specific functional group.
public enum MM4AtomCode: UInt8, RawRepresentable {
  /// Carbon
  ///
  /// MM4 atom code: 1
  case alkaneCarbon = 1
  
  /// Hydrogen
  ///
  /// MM4 atom code: 5
  case hydrogen = 5
  
  /// Oxygen
  ///
  /// MM4 atom code: 6
  case oxygen = 6
  
  /// Nitrogen (trivalent)
  ///
  /// MM4 atom code: 8
  case nitrogen = 8
  
  /// Fluorine
  ///
  /// MM4 atom code: 11
  case fluorine = 11
  
  /// Sulfur
  ///
  /// MM4 atom code: 15
  case sulfur = 15
  
  /// Silicon
  ///
  /// MM4 atom code: 19
  case silicon = 19
  
  /// Phosphorus (trivalent)
  ///
  /// MM4 atom code: 25
  case phosphorus = 25
  
  /// Germanium
  ///
  /// MM4 atom code: 31
  case germanium = 31
  
  /// Carbon (5-ring)
  ///
  /// MM4 atom code: 123
  case cyclopentaneCarbon = 123
}

/// The number of hydrogens surrounding the carbon or silicon.
///
/// For non-carbon atoms, the center type reflects what it would be, if the atom
/// were replaced with a carbon, and its lone pairs replaced with hydrogens.
///
/// Methane carbons are disallowed in the simulator, as they will be simulated
/// very inaccurately. The fluorine parameters struggle with the edge case of
/// carbon tetrafluoride, and HMR creates the nonphysical carbon-8 isotope. In
/// general, primary sp3 carbons are also discouraged, for the same reasons
/// that methane is prohibited.
public enum MM4CenterType: UInt8 {
  case primary = 1
  case secondary = 2
  case tertiary = 3
  case quaternary = 4
}

/// Parameters for the vdW and electrostatic forces on an atom. Special values
/// are provided for hydrogen parameters, where applicable. Interactions
/// containing Si, Ge should be treated with MM3 heuristics:
/// - Dispersion factors for 1,4 nonbonded exceptions are eliminated.
/// - Hydrogen reduction factors for bonded hydrogens change from 0.94 to 0.923.
public struct MM4AtomParameters {
  /// Partial charge in units of proton charge.
  public var charge: Float
  
  /// Units: dimensionless
  public var dispersionFactor: Float
  
  /// Units:  kilocalorie / mole
  ///
  /// > WARNING: Convert aJ to kJ/mol.
  public var epsilon: (default: Float, hydrogen: Float)
  
  /// Units: dimensionless
  public var hydrogenReductionFactor: Float
  
  /// Units: angstrom
  ///
  /// > WARNING: Convert angstroms to nanometers.
  public var radius: (default: Float, hydrogen: Float)
}

extension MM4Parameters {
  /// - throws: `.missingParameter`, `.openValenceShell`
  mutating func createAtomCodes() throws {
    atoms.codes = try atoms.indices.map { atomID in
      let atomicNumber = atoms.atomicNumbers[atomID]
      let map = atomsToAtomsMap[atomID]
      var output: MM4AtomCode
      var valenceCount: Int
      var supportsHydrogen: Bool = false
      
      switch atomicNumber {
      case 1:
        output = .hydrogen
        valenceCount = 1
      case 6:
        let ringType = atoms.ringTypes[atomID]
        switch ringType {
        case 6:
          output = .alkaneCarbon
        case 5:
          output = .cyclopentaneCarbon
        default:
          fatalError("This should never happen.")
        }
        valenceCount = 4
        supportsHydrogen = true
      case 7:
        output = .nitrogen
        valenceCount = 3
      case 9:
        output = .fluorine
        valenceCount = 1
      case 14:
        output = .silicon
        valenceCount = 4
        supportsHydrogen = true
      case 15:
        output = .phosphorus
        valenceCount = 3
      case 16:
        output = .sulfur
        valenceCount = 2
      case 32:
        output = .germanium
        valenceCount = 4
      default:
        let address = createAddress(atomID)
        throw MM4Error.missingParameter([address])
      }
      
      if !supportsHydrogen {
        for lane in 0..<4 where map[lane] != -1 {
          let hydrogenID = Int(map[lane])
          if atoms.atomicNumbers[hydrogenID] == 1 {
            var addresses: [MM4Address] = []
            addresses.append(createAddress(atomID))
            addresses.append(createAddress(hydrogenID))
            throw MM4Error.missingParameter(addresses)
          }
        }
      }
      var valenceMask: SIMD4<Int32> = .zero
      valenceMask.replace(with: .one, where: map .!= -1)
      guard valenceMask.wrappedSum() == valenceCount else {
        let address = createAddress(atomID)
        let neighbors = createAddresses(map)
        throw MM4Error.openValenceShell(address, neighbors)
      }
      return output
    }
  }
  
  mutating func createMasses(hydrogenMassRepartitioning: Float) {
    atoms.masses = atoms.atomicNumbers.map { atomicNumber in
      /// Units: amu -> yg
      var mass = MM4MassParameters.global.mass(atomicNumber: atomicNumber)
      mass *= Float(MM4YgPerAmu)
      return mass
    }
    
    for atomID in atoms.indices
    where atoms.atomicNumbers[Int(atomID)] == 1 {
      atoms.masses[Int(atomID)] += hydrogenMassRepartitioning
      
      let map = atomsToBondsMap[Int(atomID)]
      let substituentID = Int(other(atomID: atomID, bondID: map[0]))
      atoms.masses[substituentID] -= hydrogenMassRepartitioning
    }
  }
  
  // Copying notes from the MM4_validator.swift script:
  //
  // The problems this script was built to solve, have been resolved through other
  // means. MM4 uses dipole-dipole interactions just like MM3; Jay was wrong about
  // that. Regarding the disfac=0.550, it won't change the shape of molecules
  // (internal structure) significantly. It definitely won't change the distance
  // between diamondoid surfaces, the primary area of concern regarding vdW.
  //
  // Assume any bona fide MM4 parameter had the torsions recalibrated against the
  // new disfac. That includes phosphorus, but doesn't stand for any phosphorus
  // torsions I stole from MM3. Accept the inaccuracy incurred here, just like the
  // inaccuracy from transferring Si angles/torsions to Ge.
  //
  // The hydrogen reduction factor is used **in addition** to mutations to X/H
  // radius. This combines in a wierd way to create R=0.92 for carbon, R=0.85
  // for nitrogen, oxygen, fluorine. Just use the following exceptions:
  //
  // - Si, Ge treat hydrogens like R=0.923
  //   - change the weight for the Si-H virtual sites
  // - Si, Ge treat other elements like the vdW parameters are from MM3
  //   - create a third conditional branch for the vdW energy function
  //   - particularly important: C-Si nonbonded interaction (0.08 A correction)
  // - Si, Ge treat other elements like disfac=1.000
  //   - change the parameter for 1,4 nonbonded exceptions
  //
  // Oxygen vdW parameters:
  // - Norman's actual parameters for solid O atoms likely didn't change by more
  //   than 0.01 A.
  // - Using the strange relationship that conflates the different vdW adjustments
  //   to create a single "R=0.85", I derived r=3.046 A, eps=0.0840.
  // - Si, P, S, Ge use the Hill function exactly instead of X/H vdW pairs.
  mutating func createNonbondedParameters(hydrogenMassRepartitioning: Float) {
    for atomID in atoms.indices {
      let atomicNumber = atoms.atomicNumbers[atomID]
      var epsilon: (default: Float, hydrogen: Float)
      var radius: (default: Float, hydrogen: Float)
      
      switch atomicNumber {
      case 1:
        // Set the hydrogen parameters to -1, so that a simple maximum
        // operation can be used to select the parameter for the heteroatom.
        // Furthermore, multiplication of the two atoms' hydrogen parameters can
        // be used as efficient logic for checking whether hydrogen exists. If
        // the result is negative, choose the hydrogen parameters (XOR gate).
        epsilon = (default: 0.017, hydrogen: -1)
        radius = (default: 1.640, hydrogen: -1)
      case 6:
        // If a carbon belongs to a rigid body with a certain amount of HMR,
        // assume all hydrogens it contacts have the same amount of HMR. This
        // may not be 100% accurate, but seems like the most logical course of
        // action. Most often, people will use the same amount of HMR for each
        // rigid body.
        let t = hydrogenMassRepartitioning
        let hydrogenRadius = t * (3.410 - 3.440) + 3.440
        epsilon = (default: 0.037, hydrogen: 0.024)
        radius = (default: 1.960, hydrogen: hydrogenRadius)
      case 7:
        epsilon = (default: 0.054, hydrogen: 0.110)
        radius = (default: 1.860, hydrogen: 3.110)
      case 8:
        epsilon = (default: 0.059, hydrogen: 0.084)
        radius = (default: 1.820, hydrogen: 3.046)
      case 9:
        epsilon = (default: 0.075, hydrogen: 0.092)
        radius = (default: 1.710, hydrogen: 2.870)
      case 14:
        epsilon = (default: 0.140, hydrogen: -1)
        radius = (default: 2.290, hydrogen: -1)
      case 15:
        epsilon = (default: 0.168, hydrogen: -1)
        radius = (default: 2.220, hydrogen: -1)
      case 16:
        epsilon = (default: 0.196, hydrogen: -1)
        radius = (default: 2.090, hydrogen: -1)
      case 32:
        epsilon = (default: 0.200, hydrogen: -1)
        radius = (default: 2.440, hydrogen: -1)
      default:
        fatalError("This should never happen.")
      }
      
      // Pre-compute the Hill function for hydrogen interactions.
      if atomicNumber != 1, epsilon.hydrogen < 0, radius.hydrogen < 0 {
        epsilon.hydrogen = (epsilon.default * 0.017).squareRoot()
        radius.hydrogen = radius.default + 1.640
      }
      
      // Use MM3 heuristics for silicon and germanium.
      var dispersionFactor: Float = 0.550
      var hydrogenReductionFactor: Float = 0.94
      if atomicNumber == 14 || atomicNumber == 32 {
        dispersionFactor = 1.000
        hydrogenReductionFactor = 0.923
      }
      
      atoms.parameters.append(
        MM4AtomParameters(
          charge: 0,
          dispersionFactor: dispersionFactor,
          epsilon: epsilon,
          hydrogenReductionFactor: hydrogenReductionFactor,
          radius: radius))
    }
  }
  
  mutating func createNonbondedExceptions() {
    // Create nonbonded exceptions.
    var nonbondedExceptions13Map: [SIMD2<UInt32>: Bool] = [:]
    var nonbondedExceptions14Map: [SIMD2<UInt32>: Bool] = [:]
    for torsion in torsions.indices {
      guard torsion[0] < torsion[3] else {
        fatalError("Torsion was not sorted.")
      }
      let pair = sortBond(SIMD2(torsion[0], torsion[3]))
      nonbondedExceptions14Map[pair] = true
    }
    for angle in angles.indices {
      guard angle[0] < angle[2] else {
        fatalError("Angle was not sorted.")
      }
      let pair = sortBond(SIMD2(angle[0], angle[2]))
      nonbondedExceptions13Map[pair] = true
      nonbondedExceptions14Map[pair] = nil
    }
    
    // Remove 1,2 interactions from erroneously being in either map. This can
    // often happen with 5-membered rings.
    for bond in bonds.indices {
      nonbondedExceptions13Map[bond] = nil
      nonbondedExceptions14Map[bond] = nil
    }
    nonbondedExceptions13 = nonbondedExceptions13Map.keys.map { $0 }
    nonbondedExceptions14 = nonbondedExceptions14Map.keys.map { $0 }
  }
}
