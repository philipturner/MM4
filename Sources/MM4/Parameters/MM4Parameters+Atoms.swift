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
  public var atomicNumbers: [UInt8] = []
  
  /// The center type used to assign different parameters.
  ///
  /// This is useful for debugging the assignment of parameters, to ensure the
  /// exact parameter as specified by the forcefield paper gets assigned.
  public var centerTypes: [MM4CenterType?] = []
  
  /// The MM4 code for each atom in the system.
  public var codes: [MM4AtomCode] = []
  
  /// Each value corresponds to the atom at the same array index.
  public var extendedParameters: [MM4AtomExtendedParameters?] = []
  
  /// The mass of each atom after hydrogen mass repartitioning.
  public var masses: [Float] = []
  
  /// Each value corresponds to the atom at the same array index.
  public var nonbondedParameters: [MM4NonbondedParameters] = []
  
  /// The smallest ring this is involved in.
  public var ringTypes: [UInt8] = []
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
  
  /// Nitrogen
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

/// Parameters for the van der Waals force on a specific atom, with an
/// alternative value for use in hydrogen interactions. This force does not
/// include electrostatic forces, which are handled separately in a bond-based
/// dipole-dipole interaction.
public struct MM4NonbondedParameters {
  /// Units:  kilocalorie / mole
  ///
  /// `heteroatom` includes carbon; the term was simply chosen as an antonym to
  /// hydrogen. Epsilons are computed using the geometric mean for heteroatoms,
  /// otherwise substitute directly with the hydrogen epsilon.
  ///
  /// > WARNING: Convert aJ to kJ/mol.
  public var epsilon: (heteroatom: Float, hydrogen: Float)
  
  /// Units: angstrom
  ///
  /// `heteroatom` includes carbon; the term was simply chosen as an antonym to
  /// hydrogen. Radii are computed using the arithmetic sum for heteroatoms,
  /// otherwise substitute directly with the hydrogen radius.
  ///
  /// > WARNING: Convert angstroms to nanometers.
  public var radius: (heteroatom: Float, hydrogen: Float)
}

public struct MM4AtomExtendedParameters {
  /// Partial charge in units of proton charge.
  public var charge: Float
}

extension MM4Parameters {
  func createMasses() {
    atoms.masses = atoms.atomicNumbers.map { atomicNumber in
      MM4MassParameters.global.mass(atomicNumber: atomicNumber)
    }
    for atomID in 0..<atoms.atomicNumbers.count
    where atoms.atomicNumbers[atomID] == 1 {
      atoms.masses[atomID] += hydrogenMassRepartitioning
      
      // HMR affects both carbon and silicon.
      let map = atomsToBondsMap[atomID]
      guard map[0] != -1, map[1] == -1, map[2] == -1, map[3] == -1 else {
        fatalError("Hydrogen did not have exactly 1 bond.")
      }
      let substituentID = Int(other(atomID: atomID, bondID: map[0]))
      atoms.masses[substituentID] -= hydrogenMassRepartitioning
    }
  }
  
  func createAtomCodes() {
    atoms.codes = atoms.atomicNumbers.indices.map { atomID in
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
    let permittedAtomicNumbers: [UInt8] = [6, 7, 14, 15, 16]
    for atomID in atoms.atomicNumbers.indices {
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
  
  // Two forces:
  // - Short-range vdW force between pairs of non-charged atoms smaller than
  //   the C-Si vdW radius.
  //   - Interaction group 1: small vdW radius and non-charged
  //   - Interaction group 2: any vdW radius and any charge status
  // - Electrostatics and either (a) larger vdW radii or (b) vdW interactions
  //   between pairs of charged particles.
  //   - Interaction group 1: large vdW radius and/or charged
  //   - Interaction group 2: large vdW radius and/or charged
  //
  // Near-term unoptimized alternative: two separate O(n^2) forces without any
  // cutoff, one for vdW and one for electrostatics
  func createNonbondedParameters() {
    for atomID in atoms.atomicNumbers.indices {
      let atomicNumber = atoms.atomicNumbers[atomID]
      var epsilon: (heteroatom: Float, hydrogen: Float)
      var radius: (heteroatom: Float, hydrogen: Float)
      
      switch atomicNumber {
      case 1:
        epsilon = (heteroatom: 0.017, hydrogen: 0.017)
        radius = (heteroatom: 1.640, hydrogen: 1.640)
      case 6:
        let t = Float(hydrogenMassRepartitioning) - 0
        let hydrogenRadius = t * (3.410 - 3.440) + 3.440
        epsilon = (heteroatom: 0.037, hydrogen: 0.024)
        radius = (heteroatom: 1.960, hydrogen: hydrogenRadius)
      case 7:
        // Parameters coming directly from the research papers:
        // hydrogen epsilon: 0.030 -> 0.110 = (1 / 80.7%)^6
        // hydrogen radius: 3.50 -> 3.110 = 88.9%
        epsilon = (heteroatom: 0.054, hydrogen: 0.110)
        radius = (heteroatom: 1.860, hydrogen: 3.110)
      case 8:
        // Extrapolating the missing MM4 parameters for oxygen:
        // - MM3 parameter for heteroatom epsilon/radius is likely the same. It
        //   is exactly the same for F, but slightly different for N. The
        //   hydrogen vdW paper specifically said "MM3/MM4 parameter",
        //   indicating that MM3(2000) used the same one as MM4.
        // - Tinker's parameters do not include the special C-H pairs outside of
        //   hydrogen bonding. But, I can extrapolate. The hydrogen vdW paper
        //   showed N having 2.8 A for MM4, but then changed that to 3.1 A for
        //   the final amines paper (~6 years later). Fluorine has 2.87 A, very
        //   close to the 2.9 A in hydrogen vdW (which used only 1 s.f.).
        // - Combine using the fluorine heuristic, decided upon at the end of
        //   the hydrogen vdW paper. Combine O + H with the usual method, but
        //   scale by 0.85 and 2.6.
        //
        // hydrogen epsilon: 0.032 -> 0.084 = (1 / 85.0%)^6
        // hydrogen radius: 3.46 -> 2.941 = 85.0%
        //   - 0.054 < 0.0590 < 0.075
        //   - 0.110 > 0.084 is NOT greater than 0.092
        //   - 1.860 > 1.8200 > 1.710
        //   - 3.110 > 2.941 > 2.870
        //
        // This attempt seems to have failed. Interpolate the scaling factors
        // between nitrogen and fluorine. The energy is rather large for N/O/F,
        // but potentially not from regular vdW interactions. This could be
        // emulating the absent induced polarization from electronegativity,
        // which is a good thing. Second-row elements and first-row carbon are
        // not very electronegative; it seems reasonable they should follow a
        // different trend for nonbonded hydrogen epsilons.
        //
        // hydrogen epsilon: 0.032 -> 0.098 = (1 / 82.9)^6
        // hydrogen radius: 3.46 -> 3.010 = 87.0%
        //   - 0.110 > 0.098 > 0.092
        //   - 3.110 > 3.010 > 2.870
        epsilon = (heteroatom: 0.059, hydrogen: 0.098)
        radius = (heteroatom: 1.820, hydrogen: 3.010)
      case 9:
        // Parameters coming directly from the research papers:
        // hydrogen epsilon: 0.036 -> 0.092 = (1 / 85.5%)^6
        // hydrogen radius: 3.35 -> 2.870 = 85.7%
        epsilon = (heteroatom: 0.075, hydrogen: 0.092)
        radius = (heteroatom: 1.710, hydrogen: 2.870)
      case 14:
        // For phosphorus and sulfur, we know the electronegative lone pair will
        // polarize hydrogens that get close. For silicon, it's more like
        // carbon. With carbon, the energy also decreased by 0.94x, when the
        // radius decreased that much. Instead of the formula used for H-N and
        // H-F, where epsilon scaled in the opposite direction to radius. I
        // think with silicon, the 0.94x correction might do more harm than
        // good.
        //
        // I will not change the energy of silicon. Carbon decreased it using
        // the factor 0.94x; more electronegative atoms increased it by an
        // amount roughly following a r^-6 law. Not changing at all is
        // treating silicon just like carbon, which seems reasonable.
        epsilon = (heteroatom: 0.140, hydrogen: 0.0488)
        radius = (heteroatom: 2.290, hydrogen: 3.930 * 0.94)
      case 15:
        // Use the MM3 parameters for phosphorus, as the MM4 paper doesn't
        // contain vdW parameters.
        var epsilonScale: Float = 1 / (0.94 * 0.94 * 0.94)
        epsilonScale *= epsilonScale // sixth power
        epsilon = (heteroatom: 0.168, hydrogen: 0.053 * epsilonScale)
        radius = (heteroatom: 2.220, hydrogen: 3.860 * 0.94)
      case 16:
        // Scale H-S vdW parameters by 0.94, as suggested for MM4.
        var epsilonScale: Float = 1 / (0.94 * 0.94 * 0.94)
        epsilonScale *= epsilonScale // sixth power
        epsilon = (heteroatom: 0.196, hydrogen: 0.0577 * epsilonScale)
        radius = (heteroatom: 2.090, hydrogen: 3.730 * 0.94)
      default:
        fatalError("Atomic number \(atomicNumber) not recognized.")
      }
      atoms.nonbondedParameters.append(
        MM4NonbondedParameters(epsilon: epsilon, radius: radius))
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
      nonbondedExceptions14Map[SIMD2(torsion[0], torsion[3])] = true
    }
    for angle in angles.indices {
      guard angle[0] < angle[2] else {
        fatalError("Angle was not sorted.")
      }
      nonbondedExceptions13Map[SIMD2(angle[0], angle[2])] = true
      nonbondedExceptions14Map[SIMD2(angle[0], angle[2])] = nil
    }
    nonbondedExceptions13 = nonbondedExceptions13Map.keys.map { $0 }
    nonbondedExceptions14 = nonbondedExceptions14Map.keys.map { $0 }
  }
}
