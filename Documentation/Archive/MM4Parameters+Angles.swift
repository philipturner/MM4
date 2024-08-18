//
//  MM4Parameters+Angles.swift
//  MM4
//
//  Created by Philip Turner on 10/7/23.
//

/// Parameters for a group of 3 atoms.
public struct MM4Angles {
  /// Each value corresponds to the angle at the same array index.
  public var extendedParameters: [MM4AngleExtendedParameters?] = []
  
  /// Groups of atom indices that form an angle.
  public var indices: [SIMD3<UInt32>] = []
  
  /// Map from a group of atoms to an angle index.
  public var map: [SIMD3<UInt32>: UInt32] = [:]
  
  /// Each value corresponds to the angle at the same array index.
  public var parameters: [MM4AngleParameters] = []
  
  /// The smallest ring each angle is involved in.
  public var ringTypes: [UInt8] = []
  
  mutating func append(contentsOf other: Self, atomOffset: UInt32) {
    let angleOffset = UInt32(self.indices.count)
    self.extendedParameters += other.extendedParameters
    self.indices += other.indices.map {
      $0 &+ atomOffset
    }
    for key in other.map.keys {
      let value = other.map[key].unsafelyUnwrapped
      self.map[key &+ atomOffset] = value &+ angleOffset
    }
    self.parameters += other.parameters
    self.ringTypes += other.ringTypes
  }
}

/// Parameters for an angle between two bonds.
///
/// The parameters include bending stiffness and the multiplicative contribution
/// to bend-bend stiffness.
public struct MM4AngleParameters {
  /// Units: millidyne \* angstrom / radian^2
  ///
  /// > WARNING: Convert aJ to kJ/mol.
  public var bendBendStiffness: Float
  
  /// Units: millidyne \* angstrom / radian^2
  ///
  /// > WARNING: Convert aJ to kJ/mol.
  public var bendingStiffness: Float
  
  /// Units: degree
  ///
  /// > WARNING: Convert degrees to radians.
  public var equilibriumAngle: Float
  
  /// Units: millidyne / radian \* mole
  ///
  /// > WARNING: Convert aJ to kJ/mol.
  public var stretchBendStiffness: Float
}

/// Parameters for angle forces unique to highly electronegative elements.
public struct MM4AngleExtendedParameters {
  /// Stiffness for type 2 stretch-bend forces, affecting bonds not directly
  /// involved in this angle.
  ///
  /// > WARNING: Convert aJ to kJ/mol.
  public var stretchBendStiffness: Float
  
  /// > WARNING: Convert aJ to kJ/mol.
  public var stretchStretchStiffness: Float
}

extension MM4Parameters {
  /// - throws: `.missingParameter`
  mutating func createAngleParameters(forces: MM4ForceOptions) throws {
    for angleID in angles.indices.indices {
      let angle = angles.indices[angleID]
      let ringType = angles.ringTypes[angleID]
      let codes = with5RingsRemoved {
        createAtomCodes(group: angle, zero: SIMD3<UInt8>.zero)
      }
      
      func createAngleError() -> MM4Error {
        let map = SIMD3<Int32>(truncatingIfNeeded: angle)
        let addresses = createAddresses(SIMD4(map, -1))
        return MM4Error.missingParameter(addresses)
      }
      if containsTwoNonCarbons(codes) {
        throw createAngleError()
      }
      var sortedCodes = sortAngle(codes)
      
      // MARK: - Bend
      
      var bendingStiffnesses: SIMD3<Float>?
      var equilibriumAngles: SIMD3<Float>?
      var attemptCount = 2
      if !forces.contains(.bend) {
        bendingStiffnesses = .zero
        equilibriumAngles = .zero
        attemptCount = 0
      }
      
      var continueAttempt = false
      for attemptID in 0..<attemptCount {
        var codes: SIMD3<UInt8>
        if attemptID == 0 {
          codes = sortedCodes
        } else if continueAttempt {
          codes = sortedCodes.replacing(
            with: .one, where: sortedCodes .== 123)
          codes = sortAngle(codes)
        } else {
          continue
        }
        
        switch (codes[0], codes[1], codes[2]) {
          // Carbon
        case (1, 1, 1):
          bendingStiffnesses = SIMD3(repeating: 0.740)
          equilibriumAngles = SIMD3(109.500, 110.400, 111.800)
        case (1, 1, 5):
          bendingStiffnesses = SIMD3(0.590, 0.560, 0.600)
          equilibriumAngles = SIMD3(108.900, 109.470, 110.800)
        case (5, 1, 5):
          bendingStiffnesses = SIMD3(repeating: 0.540)
          equilibriumAngles = SIMD3(107.700, 107.800, 107.700)
        case (1, 1, 123):
          bendingStiffnesses = SIMD3(repeating: 0.740)
          equilibriumAngles = SIMD3(109.500, 110.500, 111.800)
        case (1, 123, 5):
          bendingStiffnesses = SIMD3(repeating: 0.560)
          equilibriumAngles = SIMD3(108.900, 109.470, 110.800)
        case (5, 1, 123):
          bendingStiffnesses = SIMD3(repeating: 0.560)
          equilibriumAngles = SIMD3(108.900, 109.470, 110.800)
        case (5, 123, 5):
          bendingStiffnesses = SIMD3(repeating: 0.620)
          equilibriumAngles = SIMD3(107.800, 107.800, 0.000)
        case (1, 123, 123):
          bendingStiffnesses = SIMD3(repeating: 0.740)
          equilibriumAngles = SIMD3(109.500, 110.500, 111.800)
        case (5, 123, 123):
          bendingStiffnesses = SIMD3(repeating: 0.580)
          equilibriumAngles = SIMD3(108.900, 109.470, 110.800)
        case (123, 123, 123):
          bendingStiffnesses = SIMD3(repeating: 0.740)
          if ringType == 5 {
            equilibriumAngles = SIMD3(108.300, 108.900, 109.000)
          } else {
            continueAttempt = true
          }
          
          // Nitrogen
        case (1, 1, 8):
          bendingStiffnesses = SIMD3(1.175, 1.165, 1.145)
          equilibriumAngles = SIMD3(106.4, 104.0, 104.6)
        case (5, 1, 8):
          bendingStiffnesses = SIMD3(0.850, 0.850, 1.110)
          equilibriumAngles = SIMD3(104.2, 105.0, 104.6)
        case (1, 8, 1):
          bendingStiffnesses = SIMD3(1.050, 0.970, .nan)
          equilibriumAngles = SIMD3(105.8, 106.6, .nan)
        case (1, 8, 123):
          bendingStiffnesses = SIMD3(0.880, 0.880, .nan)
          equilibriumAngles = SIMD3(109.4, 109.4, .nan)
        case (5, 123, 8):
          bendingStiffnesses = SIMD3(repeating: 0.500)
          equilibriumAngles = SIMD3(repeating: 109.4)
        case (8, 123, 123):
          bendingStiffnesses = SIMD3(repeating: 1.155)
          equilibriumAngles = SIMD3(repeating: 107.1)
          if ringType == 5 {
            equilibriumAngles![2] = 105.9
          }
        case (123, 8, 123):
          bendingStiffnesses = SIMD3(0.880, 0.880, .nan)
          if ringType == 5 {
            equilibriumAngles = SIMD3(105.2, 108.6, .nan)
          } else {
            continueAttempt = true
          }
        case (8, 1, 8):
          // This was taken from the MM3 parameters.
          bendingStiffnesses = SIMD3(repeating: 1.045)
          equilibriumAngles = SIMD3(repeating: 110.74)
          
          // Oxygen
        case (1, 1, 6):
          bendingStiffnesses = SIMD3(repeating: 1.275)
          equilibriumAngles = SIMD3(105.5, 106.2, 107.9)
        case (5, 1, 6):
          bendingStiffnesses = SIMD3(0.970, 0.870, 1.120)
          equilibriumAngles = SIMD3(106.9, 107.2, 106.6)
        case (6, 1, 6):
          bendingStiffnesses = SIMD3(repeating: 1.050)
          equilibriumAngles = SIMD3(108.0, 107.0, 107.1)
        case (1, 6, 1):
          bendingStiffnesses = SIMD3(repeating: 0.920)
          equilibriumAngles = SIMD3(repeating: 107.6)
        case (5, 123, 6):
          bendingStiffnesses = SIMD3(1.120, .nan, .nan)
          equilibriumAngles = SIMD3(106.5, .nan, .nan)
        case (6, 123, 6):
          bendingStiffnesses = SIMD3(repeating: 1.050)
          equilibriumAngles = SIMD3(110.0, 110.0, 107.1)
          if ringType == 5 {
            equilibriumAngles![2] = 107.7
          }
        case (6, 123, 123):
          bendingStiffnesses = SIMD3(repeating: 1.275)
          equilibriumAngles = SIMD3(105.5, 106.5, 107.9)
          if ringType == 5 {
            equilibriumAngles![1] = 105.5
            equilibriumAngles![2] = 105.9
          }
        case (123, 6, 123):
          if ringType == 5 {
            bendingStiffnesses = SIMD3(repeating: 0.920)
            equilibriumAngles = SIMD3(repeating: 110.0)
          } else {
            continueAttempt = true
          }
          
          // Fluorine
        case (1, 1, 11):
          bendingStiffnesses = SIMD3(repeating: 0.92)
          equilibriumAngles = SIMD3(106.90, 108.20, 109.30)
        case (5, 1, 11):
          bendingStiffnesses = SIMD3(0.82, 0.88, 0.98)
          equilibriumAngles = SIMD3(107.95, 107.90, 108.55)
        case (11, 1, 11):
          bendingStiffnesses = SIMD3(1.95, 2.05, 1.62)
          equilibriumAngles = SIMD3(104.30, 105.90, 108.08)
          
          // Silicon
        case (1, 1, 19):
          if ringType == 6 {
            bendingStiffnesses = SIMD3(repeating: 0.400)
            equilibriumAngles = SIMD3(109.00, 112.70, 111.50)
          } else {
            bendingStiffnesses = SIMD3(repeating: 0.550)
            equilibriumAngles = SIMD3(repeating: 107.20)
          }
        case (5, 1, 19):
          bendingStiffnesses = SIMD3(repeating: 0.540)
          equilibriumAngles = SIMD3(109.50, 110.00, 108.90)
        case (1, 19, 1):
          if ringType == 6 {
            bendingStiffnesses = SIMD3(repeating: 0.480)
            equilibriumAngles = SIMD3(109.50, 110.40, 109.20)
          } else {
            bendingStiffnesses = SIMD3(repeating: 0.650)
            equilibriumAngles = SIMD3(102.80, 103.80, 99.50)
          }
        case (19, 1, 19):
          // The wierd 120-degree like parameters are required to agree with
          // xTB results.
          bendingStiffnesses = SIMD3(repeating: 0.350)
          equilibriumAngles = SIMD3(109.50, 119.50, 117.00)
        case (1, 19, 5):
          bendingStiffnesses = SIMD3(repeating: 0.400)
          equilibriumAngles = SIMD3(109.30, 107.00, 110.00)
        case (5, 19, 5):
          bendingStiffnesses = SIMD3(repeating: 0.460)
          equilibriumAngles = SIMD3(106.50, 108.70, 109.50)
        case (1, 19, 19):
          bendingStiffnesses = SIMD3(repeating: 0.450)
          equilibriumAngles = SIMD3(repeating: 109.00)
        case (5, 19, 19):
          bendingStiffnesses = SIMD3(repeating: 0.350)
          equilibriumAngles = SIMD3(repeating: 109.40)
        case (19, 19, 19):
          if ringType == 6 {
            // Typo from the MM3 silicon paper and retained in the MM3(2000)
            // implementation donated to Tinker. Quaternary sp3 carbon has the
            // parameters 109.5-112.7-111.5, while sp3 silicon *should* have
            // something similar: 109.5-110.8-111.2. I think 118.00 was a typo
            // from the column four cells below: 19-22-22. Anything connected
            // to the other side of a cyclopropane carbon (60°) should have an
            // angle like 120°. This is not the first typo I have caught in
            // one of Allinger's research papers, see the note about the MM4
            // formula for the Torsion-Stretch cross-term.
            //
            // The stiffness does match up. Extrapolating the ratios of
            // 1-1-19 : 19-19-19 and 1-19-1 : 19-19-19 from 5-membered ring
            // variants, one gets 0.233 and 0.236 respectively for 19-19-19.
            // That is very close to 0.25, so I don't think that was messed
            // up.
            //
            // After doing deeper investigation with xTB, it turns out this
            // parameter is correct. 118.00 has lower RMS error than 109.50 and
            // 112.00.
            bendingStiffnesses = SIMD3(repeating: 0.250)
            equilibriumAngles = SIMD3(118, 110.80, 111.20)
          } else {
            bendingStiffnesses = SIMD3(repeating: 0.320)
            equilibriumAngles = SIMD3(repeating: 106.00)
          }
          
          // Phosphorus
        case (1, 25, 1):
          bendingStiffnesses = SIMD3(0.900, 0.725, .nan)
          equilibriumAngles = SIMD3(94.50, 97.90, .nan)
        case (5, 1, 25):
          bendingStiffnesses = SIMD3(0.475, 0.570, 0.717)
          equilibriumAngles = SIMD3(106.3, 108.3, 108.0)
        case (1, 1, 25):
          bendingStiffnesses = SIMD3(0.750, 0.825, 0.725)
          equilibriumAngles = SIMD3(107.05, 108.25, 109.55)
        case (25, 1, 25):
          // There is no angle parameter for P-C-P, so we'll have to examine
          // between nearby elements.
          //
          // | C  | N | O | F  |
          // | Si |   | S | Cl |
          //
          // Stiffness
          //
          // | 0.740 | 1.045 | 1.050 | 1.950 |
          // | 0.350 |       | 0.420 | 0.750 |
          //
          // Equilibrium Angle
          //
          // | 109.5 | 110.74 | 108.0 | 104.3 |
          // | 109.5 |        | 110.0 | 108.1 |
          //
          // Relationships between X-C-X and C-C-X parameters:
          // -  N: 1.175 -> 1.045, 106.4 -> 110.7
          // -  O: 1.275 -> 1.050, 105.5 -> 108.0
          // - Si: 0.400 -> 0.350, 109.0 -> 109.5
          // -  P: 0.750 ->      , 107.1 ->
          // -  S: 0.975 -> 0.420, 102.6 -> 110.00
          //
          // From C-C-X to X-C-X, angle always increases. The values for Si and
          // S are both within the narrow range 109.5-110.00. Using 109.5 for P
          // seems like the obvious choice.
          //
          // From the first row to the third row of the periodic table, there is
          // consistent pattern. It can be interpolated bilinearly to reach a
          // number of ~0.400. The nearby S parameter also jumps downward when
          // going from C-C-X to X-C-X. Silicon doesn't, but it's an exception
          // because it is electropositive (neutral C-C-C angle character).
          bendingStiffnesses = SIMD3(repeating: 0.4)
          equilibriumAngles = SIMD3(repeating: 109.5)
          
          // Sulfur
        case (5, 1, 15):
          bendingStiffnesses = SIMD3(repeating: 0.782)
          equilibriumAngles = SIMD3(108.9, 108.8, 105.8)
        case (1, 15, 1):
          bendingStiffnesses = SIMD3(0.920, .nan, .nan)
          equilibriumAngles = SIMD3(97.2, .nan, .nan)
        case (1, 1, 15):
          bendingStiffnesses = SIMD3(repeating: 0.975)
          equilibriumAngles = SIMD3(102.6, 105.7, 107.7)
        case (5, 123, 5):
          bendingStiffnesses = SIMD3(0.680, 0.680, .nan)
          equilibriumAngles = SIMD3(109.1, 107.5, .nan)
        case (1, 123, 15):
          bendingStiffnesses = SIMD3(repeating: 0.975)
          equilibriumAngles = SIMD3(102.6, 110.8, 107.7)
        case (123, 15, 123):
          bendingStiffnesses = SIMD3(0.920, .nan, .nan)
          equilibriumAngles = SIMD3(ringType == 5 ? 96.5 : 97.2, .nan, .nan)
        case (15, 123, 123):
          if ringType == 5 {
            bendingStiffnesses = SIMD3(repeating: 1.050)
            equilibriumAngles = SIMD3(108.0, 108.0, 108.5)
          } else {
            bendingStiffnesses = SIMD3(repeating: 0.975)
            equilibriumAngles = SIMD3(repeating: 106.2)
          }
        case (15, 1, 15), (15, 123, 15):
          // Grabbing the S-C-S angle parameters from MM3.
          bendingStiffnesses = SIMD3(repeating: 0.420)
          equilibriumAngles = SIMD3(repeating: 110.00)
          
          // Germanium
        case (1, 1, 31):
          if ringType == 5 {
            bendingStiffnesses = SIMD3(repeating: 0.680)
            equilibriumAngles = SIMD3(repeating: 105.5)
          } else {
            bendingStiffnesses = SIMD3(repeating: 0.450)
            equilibriumAngles = SIMD3(repeating: 109.3)
          }
        case (5, 1, 31):
          bendingStiffnesses = SIMD3(repeating: 0.420)
          equilibriumAngles = SIMD3(110.0, 111.9, 110.0)
        case (1, 31, 1):
          if ringType == 5 {
            bendingStiffnesses = SIMD3(repeating: 0.570)
            equilibriumAngles = SIMD3(repeating: 100.0)
          } else {
            bendingStiffnesses = SIMD3(repeating: 0.500)
            equilibriumAngles = SIMD3(109.5, 109.8, 110.5)
          }
        case (31, 1, 31):
          // There are no parameters for 31-1-31 in the MM3 forcefield. It
          // looks like silicon is a good first-order approximation for the
          // other angles. Therefore, I will use the silicon values here.
          //
          // xTB showed a germanium structure being very inaccurate when using
          // the large equilibrium angle parameter for sidewall silicon. I will
          // revert to just 109.5 for all types of Ge-C-Ge angles.
          //
          // The new values were interpolated halfway between the silicon
          // parameters and 109.5°. Then, the sidewall angle parameter was
          // refined to 112° to match GFN2-xTB results. 112° also made the
          // bridgehead parameter agree very well, in all metrics except nearby
          // H-Ge-H bond angle.
          //
          // Results of parameterizing germanium carbide, represented by
          // adamantane with all sidewall carbons replaced with Ge. This was
          // used to determine the bridgehead parameter (angle type 2). It is
          // not the original derivation (sidewall; angle type 3; 112°). It just
          // happened to converge at the same angle.
          //
          // | Ge-C-Ge equilibrium | 119.5 | 114.5 | 112   | 109.5 | GFN2-xTB |
          // | ------------------- | ----- | ----- | ----- | ----- | -------- |
          // | Ge-C-Ge angle       | 111.1 | 110.3 | 109.8 | 109.4 | 109.7    |
          // | C-Ge-C angle        | 106.0 | 107.9 | 108.8 | 109.6 | 109.0    |
          // | H-C-Ge angle        | 107.8 | 108.7 | 109.1 | 109.6 | 109.2    |
          // | H-Ge-C angle        | 110.9 | 110.6 | 110.4 | 110.2 | 110.2    |
          // | H-Ge-H angle        | 107.1 | 106.7 | 106.5 | 106.3 | 107.2    |
          bendingStiffnesses = SIMD3(repeating: 0.350)
          equilibriumAngles = SIMD3(109.50, 112, 112)
        case (1, 31, 5):
          bendingStiffnesses = SIMD3(repeating: 0.390)
          equilibriumAngles = SIMD3(110.2, 110.5, 111.5)
        case (5, 31, 5):
          bendingStiffnesses = SIMD3(repeating: 0.423)
          equilibriumAngles = SIMD3(107.5, 108.5, 109.5)
        case (1, 31, 31):
          bendingStiffnesses = SIMD3(repeating: 0.350)
          equilibriumAngles = SIMD3(repeating: 111.50)
        case (5, 31, 31):
          bendingStiffnesses = SIMD3(repeating: 0.350)
          equilibriumAngles = SIMD3(repeating: 114.5)
        case (31, 31, 31):
          // Silicon has some numbers similar to 112.50 when not in the
          // quaternary configuration: equilibriumAngles = SIMD3(109.50,
          // 110.80, 111.20). I need to do some benchmarks of solid germanium,
          // to see whether changing the type-1 interaction
          // constant will improve accuracy.
          //
          // The angle might need to be a variable of how many carbon atoms
          // are connected to the germanium. For now, I will set the type-1
          // constant to 109.5, which follows the pattern of carbon, and some
          // general trends in silicon parameters. I will primarily use this
          // for solid germanium or germanium carbide, where the equilibrium
          // is supposed to be at a tetrahedral conformation.
          bendingStiffnesses = SIMD3(repeating: 0.300)
          equilibriumAngles = SIMD3(repeating: 112.50)
          
          // Only set to 109.5 when simulating solid germanium.
          //          let map = atomsToAtomsMap[Int(angle[1])]
          //          var atomCodes: SIMD4<UInt8> = .zero
          //          for lane in 0..<4 {
          //            atomCodes[lane] = atoms.codes[Int(map[lane])].rawValue
          //          }
          //          if all(atomCodes .== 5 .| atomCodes .== 31) {
          //            equilibriumAngles![0] = 109.5
          //          }
          
          // This parameter change was reverted.
          
        default:
          continueAttempt = true
        }
      }
      guard let bendingStiffnesses,
            let equilibriumAngles else {
        throw createAngleError()
      }
      
      // MARK: - Bend-Bend, Stretch-Bend, Stretch-Stretch
      
      let originalCodes = sortedCodes
      sortedCodes.replace(with: .one, where: sortedCodes .== 123)
      sortedCodes = sortAngle(sortedCodes)
      
      var bendBendStiffness: Float
      var stretchBendStiffness: Float
      var stretchBendStiffness2: Float?
      var stretchStretchStiffness: Float?
      
      if forces.contains(.bendBend) ||
          forces.contains(.stretchBend) ||
          forces.contains(.stretchStretch) {
        if sortedCodes[0] == 5, sortedCodes[2] == 5 {
          bendBendStiffness = 0.000
          stretchBendStiffness = 0.000
        } else if any(sortedCodes .== 6) {
          // Oxygen
          if sortedCodes[1] == 1 {
            if sortedCodes[0] == 5 && sortedCodes[2] == 6 {
              bendBendStiffness = 0.20
            } else if any(sortedCodes .== 5) {
              bendBendStiffness = 0.24
            } else {
              bendBendStiffness = 0.30
            }
          } else {
            // If oxygen is in the center, it's impossible to have a bend-bend
            // interaction.
            bendBendStiffness = 0.00
          }
          
          switch (sortedCodes[0], sortedCodes[1], sortedCodes[2]) {
          case (1, 1, 6):
            if ringType == 5 && all(originalCodes .== SIMD3(6, 123, 123)) {
              stretchBendStiffness = 0.50
            } else {
              stretchBendStiffness = 0.02
            }
          case (5, 1, 6):
            stretchBendStiffness = 0.36
          case (1, 6, 1):
            if ringType == 5 && all(originalCodes .== SIMD3(123, 6, 123)) {
              stretchBendStiffness = 0.50
            } else {
              stretchBendStiffness = -0.12
            }
          default:
            throw createAngleError()
          }
        } else if any(sortedCodes .== 11) {
          // Fluorine
          guard sortedCodes[2] == 11 else {
            throw createAngleError()
          }
          switch sortedCodes[0] {
          case 1:
            bendBendStiffness = -0.10
            stretchBendStiffness = 0.160
            stretchBendStiffness2 = 0.000
            stretchStretchStiffness = 0.22
          case 5:
            bendBendStiffness = 0.00
            stretchBendStiffness = 0.160
            stretchBendStiffness2 = 0.000
            stretchStretchStiffness = -0.45
          case 11:
            bendBendStiffness = 0.09
            stretchBendStiffness = 0.140
            stretchBendStiffness2 = 0.275
            stretchStretchStiffness = 1.00
          default:
            throw createAngleError()
          }
        } else {
          switch sortedCodes[1] {
            // Carbon
          case 1:
            // Assume the MM4 paper's parameters for H-C-C/C-C-C also apply to
            // H-C-Si/C-C-Si/Si-C-Si.
            if all(sortedCodes .!= 5) {
              bendBendStiffness = 0.204
              stretchBendStiffness = (ringType == 5) ? 0.180 : 0.140
              
              // Exception for the nitrogen-containing bond C-C-N.
              if sortedCodes[0] == 1, sortedCodes[2] == 8 {
                stretchStretchStiffness = -0.10
              }
            } else {
              bendBendStiffness = 0.350
              stretchBendStiffness = 0.100
            }
            
            // Nitrogen
          case 8:
            bendBendStiffness = 0.204
            if codes[0] == 1, codes[2] == 1 {
              stretchBendStiffness = 0.04
            } else if codes[0] == 1, codes[0] == 123 {
              stretchBendStiffness = 0.30
            } else if codes[0] == 123, codes[2] == 123 {
              // The very large 0.30 parameter for 1-8-123 seems suspicious. I'm
              // going to set the default to 123-8-123 outside of a 5-membered
              // ring to that of 1-8-1. Often, the value inside the ring is
              // larger than in typical bonds. Not an order of magnitude
              // smaller.
              //
              // On second analysis, this parameter seems completely messed up.
              // It provides zero information, as the fallback would be 1-8-1,
              // not 1-8-123. Also, the table is malformatted there. Setting it
              // to 0.04, the same as the fallback, would be the safest choice.
              //
              // Perhaps the purpose was to clear up any confusion, as the
              // fallback rule was described in a different paper. In the
              // bend-bend section, a 1-8-123 parameter doesn't exist, so no
              // extra explanation is required to specify the 123-8-123 case.
              // Still, the presence of the 5-ring restriction is very
              // ambiguous.
              stretchBendStiffness = (ringType == 5) ? 0.04 : 0.04
            } else {
              throw createAngleError()
            }
            
            // Silicon
          case 19:
            if all(sortedCodes .!= 5) {
              bendBendStiffness = 0.24
              stretchBendStiffness = 0.10
            } else {
              bendBendStiffness = 0.30
              stretchBendStiffness = 0.06
            }
            
            // Phosphorus
          case 25:
            // There's no stretch-bend or bend-bend parameters in the phosphines
            // research paper. It seems some generic parameters were uniformly
            // applied to Si, P, and PO4 in MM3. They were not mentioned in the
            // MM4 paper, except a new bend-bend parameter for H-P-H. I don't
            // allow H-P-H angles in this forcefield.
            //
            // I will reuse the bend-bend and stretch-bend parameters from MM3.
            // They are the same as silicon.
            if all(sortedCodes .!= 5) {
              bendBendStiffness = 0.24
              stretchBendStiffness = 0.10
            } else {
              throw createAngleError()
            }
            
            // Sulfur
          case 15:
            // There cannot be a bend-bend interaction around a divalent sulfur.
            bendBendStiffness = 0.000
            if all(sortedCodes .== SIMD3(1, 15, 1)) {
              stretchBendStiffness = (ringType == 5) ? 0.280 : 0.150
            } else {
              throw createAngleError()
            }
            
            // Germanium
          case 31:
            // The parameters that Allinger created for MM3(2000) do not list
            // germanium under bend-bend parameters. However, silicon does have
            // a parameter. I will reuse that.
            bendBendStiffness = 0.240
            if all(sortedCodes .!= 5) {
              stretchBendStiffness = 0.450
            } else {
              stretchBendStiffness = 0.000
            }
            
          default:
            throw createAngleError()
          }
        }
      } else {
        bendBendStiffness = 0
        stretchBendStiffness = 0
      }
      
      if !forces.contains(.bendBend) {
        bendBendStiffness = 0
      }
      if !forces.contains(.stretchBend) {
        stretchBendStiffness = 0
        stretchBendStiffness2 = nil
      }
      if !forces.contains(.stretchStretch) {
        stretchStretchStiffness = nil
      }
      
      // MARK: - Identify Angle Type
      
      var angleType: Int
      if !forces.contains(.bend) {
        angleType = 1
      } else {
        var heavyAtomCount: Int = 0
        let map = atomsToAtomsMap[Int(angle[1])]
        for lane in 0..<4 where map[lane] != -1 {
          let otherID = UInt32(truncatingIfNeeded: map[lane])
          if all(angle .!= otherID),
             atoms.atomicNumbers[Int(otherID)] != 1 {
            heavyAtomCount &+= 1
          }
        }
        
        switch atoms.codes[Int(angle[1])] {
        case .alkaneCarbon, .cyclopentaneCarbon, .silicon, .germanium:
          switch heavyAtomCount {
          case 2: angleType = 1
          case 1: angleType = 2
          case 0: angleType = 3
          default: fatalError("Group IV atom had unexpected heavy atom count.")
          }
        case .nitrogen, .phosphorus:
          switch heavyAtomCount {
          case 1: angleType = 1
          case 0: fatalError("Group V atom was bonded to hydrogen.")
          default: fatalError("Group V atom had unexpected heavy atom count.")
          }
        case .oxygen, .sulfur:
          switch heavyAtomCount {
          case 0: angleType = 1
          default: fatalError("Group VI atom had unexpected heavy atom count.")
          }
        case .hydrogen, .fluorine:
          fatalError("Group VII atom cannot be the center of an angle.")
        }
      }
      
      guard !bendingStiffnesses[angleType - 1].isNaN,
            !equilibriumAngles[angleType - 1].isNaN else {
        print("Threw an angle error. angle type = \(angleType), bending stiffness = \(bendingStiffnesses), equilibrium angles = \(equilibriumAngles)")
        // Angle parameter was NaN.
        throw createAngleError()
      }
      angles.parameters.append(
        MM4AngleParameters(
          bendBendStiffness: bendBendStiffness,
          bendingStiffness: bendingStiffnesses[angleType - 1],
          equilibriumAngle: equilibriumAngles[angleType - 1],
          stretchBendStiffness: stretchBendStiffness))
      
      if stretchBendStiffness2 != nil || stretchStretchStiffness != nil {
        angles.extendedParameters.append(
          MM4AngleExtendedParameters(
            stretchBendStiffness: stretchBendStiffness2 ?? 0,
            stretchStretchStiffness: stretchStretchStiffness ?? 0))
      } else {
        angles.extendedParameters.append(nil)
      }
    }
  }
}
