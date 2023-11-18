//
//  MM4Parameters+Atoms.swift
//
//
//  Created by Philip Turner on 10/7/23.
//

// MARK: - Functions for assigning per-atom parameters.

/// Parameters for one atom.
public struct MM4Atoms {
  /// The number of protons in the atom's nucleus.
  public internal(set) var atomicNumbers: [UInt8] = []
  
  /// The center type used to assign different parameters.
  public internal(set) var centerTypes: [MM4CenterType?] = []
  
  /// The MM4 code for each atom in the system.
  public internal(set) var codes: [MM4AtomCode] = []
  
  /// Convenient property for iterating over the atoms.
  public internal(set) var count: Int = 0
  
  /// Convenient property for iterating over the atoms.
  public internal(set) var indices: Range<Int> = 0..<0
  
  /// The mass of each atom after hydrogen mass repartitioning.
  public internal(set) var masses: [Float] = []
  
  /// Each value corresponds to the atom at the same array index.
  public internal(set) var parameters: [MM4AtomParameters] = []
  
  /// The smallest ring this is involved in.
  public internal(set) var ringTypes: [UInt8] = []
}

/// MM4 codes for an element or an atom in a specific functional group.
public enum MM4AtomCode: UInt8, RawRepresentable {
  /// Carbon (sp3)
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
  
  /// Carbon (sp3, 5-ring)
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

/// Parameters for the van der Waals force on a specific atom. Special values
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
  // Mass repartitioning is a quite easy task. It can be reimplemented from
  // scratch in 'RigidBody'.
  func createMasses() {
    atoms.masses = atoms.atomicNumbers.map { atomicNumber in
      MM4MassParameters.global.mass(atomicNumber: atomicNumber)
    }
    for atomID in 0..<atoms.atomicNumbers.count
    where atoms.atomicNumbers[atomID] == 1 {
      atoms.masses[atomID] += hydrogenMassRepartitioning
      
      // HMR affects any tetravalent atom.
      let map = atomsToBondsMap[atomID]
      guard map[0] != -1, map[1] == -1, map[2] == -1, map[3] == -1 else {
        fatalError("Hydrogen did not have exactly 1 bond.")
      }
      let substituentID = Int(other(atomID: atomID, bondID: map[0]))
      atoms.masses[substituentID] -= hydrogenMassRepartitioning
    }
  }
  
  func createAtomCodes() {
    atoms.codes = atoms.indices.map { atomID in
      let atomicNumber = atoms.atomicNumbers[atomID]
      let map = atomsToAtomsMap[atomID]
      var output: MM4AtomCode
      var valenceCount: Int
      var supportsHydrogen: Bool = false
      
      switch atomicNumber {
      case 1:
        output = .hydrogen
        valenceCount = 1
      case 5:
        let atomicNumbers = createAtomicNumbers(map: map)
        var nitrogenMask: SIMD4<UInt8> = .zero
        nitrogenMask.replace(with: .one, where: atomicNumbers .== 7)
        let nitrogenCount = nitrogenMask.wrappedSum()
        
        if nitrogenCount == 1 {
          valenceCount = 4
          fatalError("B-N dative bond not supported.")
        } else {
          fatalError("Boron had invalid number of nitrogen bonds.")
        }
      case 6:
        let ringType = atoms.ringTypes[atomID]
        switch ringType {
        case 6:
          output = .alkaneCarbon
        case 5:
          output = .cyclopentaneCarbon
        default:
          fatalError("Unsupported carbon ring type: \(ringType)")
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
        fatalError("Atomic number \(atomicNumber) not recognized.")
      }
      
      if !supportsHydrogen {
        for lane in 0..<4 where map[lane] != -1 {
          if atoms.atomicNumbers[Int(lane)] == 1 {
            fatalError(
              "Hydrogen not allowed for atomic number \(atomicNumber).")
          }
        }
      }
      var valenceMask: SIMD4<Int32> = .zero
      valenceMask.replace(with: .one, where: map .!= -1)
      guard valenceMask.wrappedSum() == valenceCount else {
        fatalError(
          "Valence count \(valenceMask.wrappedSum()) did not match expected \(valenceCount).")
      }
      return output
    }
  }
  
  func createCenterTypes() {
    let permittedAtomicNumbers: [UInt8] = [6, 7, 8, 14, 15, 16, 31]
    for atomID in atoms.indices {
      let atomicNumber = atoms.atomicNumbers[atomID]
      guard permittedAtomicNumbers.contains(atomicNumber) else {
        precondition(
          atomicNumber == 1 || atomicNumber == 9,
          "Atomic number \(atomicNumber) not recognized.")
        atoms.centerTypes.append(nil)
        continue
      }
      
      let map = atomsToAtomsMap[atomID]
      var otherElements: SIMD4<UInt8> = .zero
      for lane in 0..<4 {
        if map[lane] == -1 {
          otherElements[lane] = 1
        } else {
          let otherID = map[lane]
          otherElements[lane] = atoms.atomicNumbers[Int(otherID)]
        }
      }
      
      // In MM4, fluorine is treated like carbon when determining carbon types.
      // Allinger notes this may be a weakness of the forcefield. This idea has
      // been extended to encapsulate all non-hydrogen atoms.
      var matchMask: SIMD4<UInt8> = .zero
      matchMask.replace(with: .one, where: otherElements .!= 1)
      
      var carbonType: MM4CenterType
      switch matchMask.wrappedSum() {
      case 4:
        carbonType = .quaternary
      case 3:
        carbonType = .tertiary
      case 2:
        carbonType = .secondary
      case 1:
        carbonType = .primary
      default:
        fatalError("This should never happen.")
      }
      atoms.centerTypes.append(carbonType)
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
  func createNonbondedParameters() {
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
        let t = Float(hydrogenMassRepartitioning) - 0
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
        fatalError("Atomic number \(atomicNumber) not recognized.")
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
  
  func createNonbondedExceptions() {
    // Create nonbonded exceptions.
    var nonbondedExceptions13Map: [SIMD2<Int32>: Bool] = [:]
    var nonbondedExceptions14Map: [SIMD2<Int32>: Bool] = [:]
    for torsion in torsions.indices {
      guard torsion[0] < torsion[3] else {
        fatalError("Torsions were not sorted.")
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

