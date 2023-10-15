//
//  MM4Parameters+Torsions.swift
//
//
//  Created by Philip Turner on 10/7/23.
//

// MARK: - Functions for assigning per-torsion parameters.

/// Parameters for a group of 4 atoms.
public struct MM4Torsions {
  /// Each value corresponds to the torsion at the same array index.
  public var extendedParameters: [MM4TorsionExtendedParameters?] = []
  
  /// Groups of atom indices that form a torsion.
  public var indices: [SIMD4<Int32>] = []
  
  /// Map from a group of atoms to a torsion index.
  public var map: [SIMD4<Int32>: Int32] = [:]
  
  /// Each value corresponds to the torsion at the same array index.
  public var parameters: [MM4TorsionParameters] = []
  
  /// The smallest ring this is involved in.
  public var ringTypes: [UInt8] = []
}

/// Parameters for a torsion among hydrogen, carbon, and silicon atoms.
///
/// V1 term:
/// - zeroed out for X-C-C-H
/// - present for C-C-C-C
/// - present for 5-membered rings
/// - present for X-C-C-F
///
/// V3 term:
/// - present for X-C-C-H
/// - present for C-C-C-C
/// - present for 5-membered rings
/// - present for X-C-C-F
///
/// Vn term:
/// - 6 for some cases of X-C-C-H
/// - 2 for some cases of C-C-C-C
/// - zeroed out for 5-membered rings
/// - 2 for X-C-C-F
///
/// 1-term torsion-stretch:
/// - present for X-C-C-H
/// - present for C-C-C-C
/// - present for 5-membered rings
/// - zeroed out for X-C-C-F, to prioritize conciseness over performance
public struct MM4TorsionParameters {
  /// Units: kilocalorie / mole
  ///
  /// > WARNING: Convert kcal/mol to kJ/mol.
  public var V1: Float
  
  /// Units: kilocalorie / mole
  ///
  /// > WARNING: Convert kcal/mol to kJ/mol.
  public var Vn: Float
  
  /// Units: kilocalorie / mole
  ///
  /// > WARNING: Convert kcal/mol to kJ/mol.
  public var V3: Float
  
  /// The factor to multiply the angle with inside the cosine term for Vn.
  ///
  /// The value of `n` is most often 2. It must be an even integer.
  public var n: Float
  
  /// Units: kilocalorie / angstrom \* mole
  ///
  /// > WARNING: Convert kcal/mol to kJ/mol.
  public var Kts3: Float
}

/// Parameters for the various torsion forces unique to non-H/C/Si-containing
/// compounds (V4, V6, 3-term torsion-stretch, torsion-bend). This also includes
/// the bend-torsion-bend force, which is omitted from C-H torsions for
/// efficiency.
public struct MM4TorsionExtendedParameters {
  /// Units: kilocalorie / mole
  ///
  /// > WARNING: Convert kcal/mol to kJ/mol.
  public var V4: Float
  
  /// Units: kilocalorie / mole
  ///
  /// > WARNING: Convert kcal/mol to kJ/mol.
  public var V6: Float
  
  /// The V1-like term contributing to torsion-stretch stiffness.
  ///
  /// > WARNING: Convert kcal/mol to kJ/mol.
  public var Kts1: (left: Float, central: Float, right: Float)
  
  /// The V2-like term contributing to torsion-stretch stiffness.
  ///
  /// > WARNING: Convert kcal/mol to kJ/mol.
  public var Kts2: (left: Float, central: Float, right: Float)
  
  /// The V3-like term contributing to torsion-stretch stiffness.
  ///
  /// > WARNING: Convert kcal/mol to kJ/mol.
  public var Kts3: (left: Float, central: Float, right: Float)
  
  /// The V1-like term contributing to torsion-bend stiffness.
  ///
  /// > WARNING: Convert radians to degrees.
  ///
  /// > WARNING: Convert kcal/mol to kJ/mol.
  public var Ktb1: (left: Float, right: Float)
  
  /// The V2-like term contributing to torsion-bend stiffness.
  ///
  /// > WARNING: Convert radians to degrees.
  ///
  /// > WARNING: Convert kcal/mol to kJ/mol.
  public var Ktb2: (left: Float, right: Float)
  
  /// The V3-like term contributing to torsion-bend stiffness.
  ///
  /// > WARNING: Convert radians to degrees.
  ///
  /// > WARNING: Convert kcal/mol to kJ/mol.
  public var Ktb3: (left: Float, right: Float)
  
  /// Bend-torsion-bend constant.
  ///
  /// > WARNING: Convert radians to degrees.
  ///
  /// > WARNING: Convert kcal/mol to kJ/mol.
  public var Kbtb: Float
}

extension MM4Parameters {
  func createTorsionParameters() {
    for torsionID in torsions.indices.indices {
      let torsion = torsions.indices[torsionID]
      let ringType = torsions.ringTypes[torsionID]
      let codes = with5RingsRemoved {
        createAtomCodes(group: torsion, zero: SIMD4<UInt8>.zero)
      }
      if containsTwoNonCarbons(codes) {
        fatalError("Torsions may not contain two non-carbon atoms.")
      }
      let sortedCodes = sortTorsion(codes)
      
      let sortedBond = sortBond(SIMD2(torsion[1], torsion[2]))
      guard let bondID = bonds.map[sortedBond] else {
        fatalError("Could not fetch bond for torsion.")
      }
      let bondStiffness = bonds.parameters[Int(bondID)].stretchingStiffness
      
      // We need to double-check that we've been following this formula
      // consistently, for both angles and torsions: "Note that unless they are
      // explicitly in the table, parameters involving five-membered ring atom
      // types (122 and 123) are assigned parmeters which involve the regular
      // atom types (types 2 and 1). This is true for all parameters. I guess
      // this means parameters only contained in a 5-membered ring must
      // sometimes be dropped entirely, using the 6-membered ring parameters as
      // if zero atoms were 122 or 123.
      //
      // We should also clean up the parameters for torsion-stretch, assigning
      // them in a separate 'switch' statement. Don't try to list all
      // combinatorial permutations of 1 vs. 123.
      if Float.random(in: 0..<1) < 2 {
        fatalError("Fix the parameters for 5-membered rings.")
      }
      
      var V1: Float = 0.000
      var Vn: Float = 0.000
      var V3: Float
      var n: Float = 2
      var V4: Float?
      var V6: Float?
      
      var Ktb_l: SIMD3<Float>?
      var Ktb_r: SIMD3<Float>?
      var Kbtb: Float?
      
      // The formula from the MM4 alkene paper was ambiguous, specifying "-k":
      //   -k * Δl * Kts * (1 + cos(3ω))
      // The formula from the MM3 original paper was:
      //   11.995 * (Kts/2) * Δl * (1 + cos(3ω))
      // After running several parameters through Desmos, and comparing similar
      // ones (https://www.desmos.com/calculator/p5wqbw7tku), I think I have
      // identified a typo. "-k" was supposed to be "K_s^-1" or "1/K_s". This
      // would result in:
      // - C-C having ~33% less TS stiffness in MM4 than in MM3
      // - C-Csp2 (MM4) having ~17% less stiffness than Si-Csp2 (MM3)
      // - C-S (MM4) having nearly identical stiffness to C-Si (MM3)
      // - Central TS for H-C-C-F (MM4) having 42% less peak stiffness than C-S
      // - Central TS for C-C-C-F (MM4) having 2x larger peaks than C-S, and a
      //   new trough with ~2x the magnitude of the C-S peak
      // - Central TS for F-C-C-F (MM4) having 1% less peak stiffness than C-S,
      //   but a new trough with ~3.7x the magnitude of the C-S peak
      //
      // New formula:
      //   Δl * (Kts / Ks) * (1 + cos(3ω))
      
      var Kts: Float = 0.000
      var Kts_l: SIMD3<Float>?
      var Kts_c: SIMD3<Float>?
      var Kts_r: SIMD3<Float>?
      
      switch (sortedCodes[0], sortedCodes[1], sortedCodes[2], sortedCodes[3]) {
        // Carbon
      case (1, 1, 1, 1):
        (V1, Vn, V3) = (0.239, 0.024, 0.637)
        Kts = 0.660
      case (1, 1, 1, 5):
        V3 = 0.290
        Kts = 0.660
      case (5, 1, 1, 5), (5, 1, 123, 5):
        V3 = 0.260
        if sortedCodes[2] == 1 {
          Vn = 0.008
          n = 6
        }
        Kts = 0.660
      case (1, 123, 123, 1), (1, 123, 123, 123):
        (V1, V3) = (0.160, 0.550)
        Kts = 0.840
      case (5, 123, 123, 5):
        V3 = 0.300
        Kts = 0.840
      case (5, 123, 123, 123):
        V3 = 0.290
        Kts = 0.840
      case (123, 123, 123, 123):
        V1 = (ringType == 5) ? -0.150 : -0.120
        V3 = (ringType == 5) ? 0.160 : 0.550
        Kts = 0.840
      case (5, 1, 123, 123), (123, 1, 123, 5):
        V3 = 0.306
        Kts = 0.640
      case (1, 123, 123, 5):
        V3 = 0.306
        Kts = 0.840
        
        // Nitrogen
      case (1, 1, 1, 8):
        (V1, Vn, V3) = (1.139, 1.348, 1.582)
        (V4, V6) = (-0.140, 0.172)
        
        Ktb_l = SIMD3(0.007, 0.000, -0.005)
        Ktb_r = SIMD3(0.010, -0.004, -0.006)
        Kbtb = -0.05
        
        Kts_l = SIMD3(0.000, 0.000, 1.300)
        Kts_c = SIMD3(3.000, -3.100, 2.860)
        Kts_r = SIMD3(-0.600, 1.000, 1.500)
      case (5, 1, 1, 8):
        V3 = 0.455
        
        Ktb_r = SIMD3(-0.003, -0.003, 0.000)
        Kbtb = -0.05
      case (8, 1, 1, 8):
        (V1, Vn, V3) = (2.545, -2.520, 3.033)
        
        Ktb_l = SIMD3(-0.008, 0.008, -0.018)
        Ktb_r = SIMD3(-0.008, 0.008, -0.018)
        
        Kts_l = SIMD3(0.000, 5.000, 0.000)
        Kts_c = SIMD3(4.000, -5.560, 2.660)
        Kts_r = SIMD3(0.000, 5.000, 0.000)
      case (1, 1, 8, 1):
        (V1, Vn, V3) = (1.193, -0.337, 0.870)
        (V4, V6) = (0.228, -0.028)
        
        Ktb_l = SIMD3(-0.028, -0.020, 0.007)
        Ktb_r = SIMD3(0.003, -0.004, -0.009)
        
        Kts_l = SIMD3(5.000, 5.000, 0.000)
        Kts_c = SIMD3(1.000, 0.000, 4.580)
      case (5, 1, 8, 1):
        (V1, Vn, V3) = (0.072, -0.512, 0.562)
        
        Ktb_l = SIMD3(-0.010, -0.025, 0.000)
        Ktb_r = SIMD3(-0.003, -0.003, 0.000)
        Kbtb = 0.05
        
        Kts_l = SIMD3(3.100, 8.200, 0.000)
        Kts_c = SIMD3(0.000, -1.100, 0.580)
      case (5, 1, 8, 123):
        (Vn, V3) = (-0.450, 0.170)
        
        Ktb_l = SIMD3(-0.003, -0.003, 0.000)
        Kbtb = -0.12
      case (1, 8, 123, 5):
        (Vn, V3) = (0.550, 0.100)
        
        Ktb_r = SIMD3(-0.015, -0.028, 0.000)
      case (1, 8, 123, 123):
        (Vn, V3) = (-0.520, 0.180)
        
        Ktb_r = SIMD3(-0.003, -0.003, 0.000)
      case (123, 8, 123, 5):
        (V1, Vn, V3) = (0.072, -0.012, 0.597)
        
        Ktb_r = SIMD3(-0.015, -0.028, 0.000)
        
        Kts_c = SIMD3(0.000, 0.000, 0.580)
        Kts_r = SIMD3(0.000, 8.038, 0.000)
      case (5, 123, 123, 8):
        V3 = 0.374
        
        Ktb_r = SIMD3(-0.003, -0.003, 0.000)
        Kbtb = -0.09
      case (123, 8, 123, 123):
        if ringType == 5 {
          (V1, Vn, V3) = (1.150, -0.040, 0.860)
          Ktb_r = SIMD3(-0.008, -0.004, 0.000)
        } else {
          (Vn, V3) = (-0.520, 0.180)
          Ktb_r = SIMD3(-0.003, -0.003, 0.000)
        }
      case (8, 123, 123, 123):
        if ringType == 5 {
          V3 = 0.699
          Ktb_l = SIMD3(-0.008, -0.004, 0.000)
        } else {
          (Vn, V3) = (-0.850, 0.200)
          Ktb_l = SIMD3(-0.008, -0.004, 0.000)
        }
        
        Kts_l = SIMD3(1.250, 1.200, 0.000)
        Kts_c = SIMD3(4.000, -3.000, 2.660)
        Kts_r = SIMD3(1.500, 0.000, 2.500)
        
        // Fluorine
      case (1, 1, 1, 11):
        (V1, Vn, V3) = (-0.360, 0.380, 0.978)
        (V4, V6) = (0.240, 0.010)
        
        Ktb_l = SIMD3(0.000, -0.012, -0.009)
        Ktb_r = SIMD3(0.005, 0.004, 0.003)
        Kbtb = -0.06
        
        Kts_c = SIMD3(5.300, -4.800, 4.500)
      case (5, 1, 1, 11):
        (V1, Vn, V3) = (-0.460, 1.190, 0.420)
        (V4, V6) = (0.000, 0.000)
        
        Ktb_l = SIMD3(0.002, -0.022, 0.000)
        Ktb_r = SIMD3(0.000, 0.000, -0.001)
        Kbtb = 0.06
        
        Kts_c = SIMD3(0.000, -1.650, 2.400)
        Kts_r = SIMD3(0.000, 0.000, 0.550)
      case (11, 1, 1, 11):
        (V1, Vn, V3) = (-1.350, 0.305, 0.355)
        (V4, V6) = (0.000, 0.000)
        
        Ktb_l = SIMD3(0.000, -0.015, -0.003)
        Ktb_r = SIMD3(0.000, -0.015, -0.003)
        Kbtb = -0.06
        
        Kts_l = SIMD3(-5.550, 4.500, 0.000)
        Kts_c = SIMD3(4.800, -5.000, 1.559)
        Kts_r = SIMD3(-5.550, 4.500, 0.000)
        
        // Silicon
        //
        // For silicon, take the MM3 torsion-stretch constant, multiply by
        // 11.995 / 2, then multiply by bond stiffness.
      case (1, 1, 1, 19):
        if ringType == 5 {
          V3 = 0.850
          Kts = 0.840
        } else {
          Vn = 0.050
          V3 = 0.240
          Kts = 0.660
        }
      case (1, 1, 19, 1):
        if ringType == 5 {
          Vn = 0.800
          V3 = 0.000
        } else {
          V3 = 0.167
        }
        Kts = 0.036 * 11.995 / 2 * bondStiffness
      case (19, 1, 1, 19):
        V3 = 0.167
        Kts = 0.660
      case (19, 1, 19, 5):
        V3 = 0.167
        Kts = 0.036 * 11.995 / 2 * bondStiffness
      case (5, 1, 19, 1):
        V3 = 0.195
        Kts = 0.036 * 11.995 / 2 * bondStiffness
      case (5, 1, 19, 5):
        V3 = 0.177
        Kts = 0.036 * 11.995 / 2 * bondStiffness
      case (19, 1, 19, 1):
        V3 = 0.100
        Kts = 0.036 * 11.995 / 2 * bondStiffness
      case (1, 1, 19, 19):
        V3 = 0.300
        Kts = 0.036 * 11.995 / 2 * bondStiffness
      case (5, 1, 19, 19):
        V3 = 0.270
        Kts = 0.036 * 11.995 / 2 * bondStiffness
      case (1, 19, 19, 5):
        V3 = 0.127
        Kts = 0.012 * 11.995 / 2 * bondStiffness
      case (1, 19, 19, 1):
        V3 = 0.107
        Kts = 0.012 * 11.995 / 2 * bondStiffness
      case (1, 19, 19, 19):
        V3 = 0.350
        Kts = 0.012 * 11.995 / 2 * bondStiffness
      case (5, 19, 19, 5):
        V3 = 0.132
        Kts = 0.012 * 11.995 / 2 * bondStiffness
      case (5, 19, 19, 19):
        V3 = 0.070
        Kts = 0.012 * 11.995 / 2 * bondStiffness
      case (19, 19, 19, 19):
        V3 = (ringType == 5) ? 0.175 : 0.125
        Kts = 0.012 * 11.995 / 2 * bondStiffness
        
        // Phosphorus
        //
        // WARNING: Entering phosphorus in the reverse order that it appears in
        // the research paper. This may cause mistakes from confusing two
        // quantities. Place Ktb_r before Ktb_l in this section. Not all cases
        // have reversed order, but pay careful attention to the ones that do.
      case (5, 1, 25, 1):
        V3 = 0.300
        V6 = -0.050
        
        Ktb_r = SIMD3(0.000, 0.000, -0.001)
        Ktb_l = SIMD3(-0.001, -0.004, -0.001)
      case (5, 1, 1, 25):
        (Vn, V3) = (0.200, 0.360)
        
        Ktb_r = SIMD3(0.005, 0.000, 0.000)
      case (1, 1, 25, 1):
        (V1, Vn, V3) = (0.800, -0.400, 0.100)
        
        Ktb_r = SIMD3(-0.002, -0.003, -0.001)
        Ktb_l = SIMD3(-0.015, 0.007, -0.006)
        
        // Sulfur
        //
        // I am assigning torsion-stretch parameters from 6-ring sulfur torsions
        // to the equivalent 5-ring sulfur torsions, when I spot any missing
        // permutations of 5-ring vs. 6-ring carbons. There are likely some
        // cases that I missed, which will reveal themselves when encountered.
      case (5, 1, 15, 1), (5, 1, 15, 123), (1, 15, 123, 5):
        V3 = 0.540
        
        Ktb_l = SIMD3(0.012, 0.016, 0.000)
        Ktb_r = SIMD3(0.003, 0.000, 0.000)
        Kbtb = 0.135
        
        if sortedCodes[2] == 15 {
          Kts_l = SIMD3(0.000, 0.900, 0.000)
        } else {
          Kts_r = SIMD3(0.000, 0.900, 0.000)
        }
      case (1, 1, 15, 1), (1, 1, 15, 123), (123, 1, 15, 1), (123, 1, 15, 123):
        (V1, V3) = (0.410, 0.600)
        
        Ktb_l = SIMD3(0.006, 0.020, 0.000)
        Kbtb = 0.080
        
        Kts_l = SIMD3(0.000, 1.439, 0.000)
        Kts_c = SIMD3(0.000, 0.000, 1.559)
      case (15, 1, 1, 15), (15, 1, 123, 15), (15, 123, 123, 15):
        (V1, Vn, V3) = (0.461, 0.144, 1.511)
        
        // This torsion is symmetric, so setting only one of the angles to 0.003
        // seems suspicious. I will set both angles.
        Ktb_l = SIMD3(0.000, 0.003, 0.000)
        Ktb_r = SIMD3(0.000, 0.003, 0.000)
        Kbtb = 0.130
        
        Kts_l = SIMD3(1.919, 1.919, 0.000)
        Kts_r = SIMD3(1.919, 1.919, 0.000)
      case (5, 1, 1, 15):
        V3 = 0.460
        Kbtb = 0.050
      case
        (1, 1, 1, 15), (15, 1, 1, 123),
        (1, 1, 123, 15), (15, 1, 123, 1),
        (1, 123, 123, 15):
        (V1, Vn, V3) = (0.420, 0.100, 0.200)
        Kbtb = 0.090
        
        if sortedCodes[3] == 15 {
          Kts_l = SIMD3(4.798, 4.798, 0)
        } else {
          Kts_r = SIMD3(4.798, 4.798, 0)
        }
      case (5, 123, 123, 5):
        V3 = 0.200
        Kbtb = 0.020
      case (123, 15, 123, 1), (1, 15, 123, 123), (123, 15, 123, 123):
        if ringType == 5 && all(sortedCodes .== SIMD4(123, 15, 123, 123)) {
          (V1, Vn, V3) = (0.440, 0.300, 0.500)
        } else {
          (V1, V3) = (0.100, 0.200)
        }
      case (5, 1, 123, 15):
        (Vn, V3) = (0.330, 0.200)
        Kbtb = 0.050
      case (15, 123, 123, 123):
        if ringType == 5 {
          (V1, Vn, V3) = (0.040, 0.200, 0.300)
          Kbtb = 0.100
        } else {
          (V1, Vn, V3) = (0.520, 0.080, 0.250)
          Kbtb = 0.100
        }
        
        Kts_r = SIMD3(4.798, 4.798, 0)
      case (123, 15, 123, 5):
        V3 = 0.450
        Kbtb = 0.020
      case (15, 1, 15, 1), (15, 1, 15, 123),
        (1, 15, 123, 15), (123, 15, 123, 15):
        // Grabbing the S-C-S-C torsion parameters from MM3.
        (Vn, V3) = (-0.900, 0.300)
      default:
        fatalError("Unrecognized torsion: \(sortedCodes)")
      }
      
      if any(sortedCodes .!= codes) {
        swap(&Kts_l, &Kts_r)
        swap(&Ktb_l, &Ktb_r)
      }
      
      torsions.parameters.append(
        MM4TorsionParameters(
          V1: V1, Vn: Vn, V3: V3, n: n, Kts3: Kts))
      
      if V4 != nil || V6 != nil ||
          Ktb_l != nil || Ktb_r != nil || Kbtb != nil ||
          Kts_l != nil || Kts_c != nil || Kts_r != nil {
        torsions.extendedParameters.append(
          MM4TorsionExtendedParameters(
            V4: V4 ?? 0, V6: V6 ?? 0,
            Kts1: (Kts_l?[0] ?? 0, Kts_c?[0] ?? 0, Kts_r?[0] ?? 0),
            Kts2: (Kts_l?[1] ?? 0, Kts_c?[1] ?? 0, Kts_r?[1] ?? 0),
            Kts3: (Kts_l?[2] ?? 0, Kts_c?[2] ?? 0, Kts_r?[2] ?? 0),
            Ktb1: (Ktb_l?[0] ?? 0, Ktb_r?[0] ?? 0),
            Ktb2: (Ktb_l?[1] ?? 0, Ktb_r?[1] ?? 0),
            Ktb3: (Ktb_l?[2] ?? 0, Ktb_r?[2] ?? 0),
            Kbtb: Kbtb ?? 0))
      }
    }
  }
}
