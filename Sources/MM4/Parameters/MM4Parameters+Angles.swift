//
//  MM4Parameters+Angles.swift
//
//
//  Created by Philip Turner on 10/7/23.
//

// MARK: - Functions for assigning per-angle parameters.

/// Parameters for a group of 3 atoms.
public struct MM4Angles {
  /// Each value corresponds to the angle at the same array index.
  public var extendedParameters: [MM4AngleExtendedParameters?] = []
  
  /// Groups of atom indices that form an angle.
  public var indices: [SIMD3<Int32>] = []
  
  /// Map from a group of atoms to an angle index.
  public var map: [SIMD3<Int32>: Int32] = [:]
  
  /// Each value corresponds to the angle at the same array index.
  public var parameters: [MM4AngleParameters] = []
  
  /// The smallest ring this is involved in.
  public var ringTypes: [UInt8] = []
}

/// Parameters for an angle between two bonds, including bending stiffness
/// and multiplicative contribution to bend-bend stiffness.
public struct MM4AngleParameters {
  /// Units: millidyne^2 / radian^2
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

/// Parameters for the various angle forces unique to non-H/C/Si-containing
/// compounds.
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
  func createAngleParameters() {
    for angleID in angles.indices.indices {
      let angle = angles.indices[angleID]
      let ringType = angles.ringTypes[angleID]
      let codes = with5RingsRemoved {
        createAtomCodes(group: angle, zero: SIMD3<UInt8>.zero)
      }
      if containsTwoNonCarbons(codes) {
        fatalError("Angles may not contain two non-carbon atoms.")
      }
      var sortedCodes = sortAngle(codes)
      
      // MARK: - Bend
      
      var bendingStiffnesses: SIMD3<Float>?
      var equilibriumAngles: SIMD3<Float>?
      
      var continueAttempt = false
      for attemptID in 0..<2 {
        var codes: SIMD3<UInt8>
        if attemptID == 0 {
          codes = sortedCodes
        } else if continueAttempt {
          codes = sortedCodes.replacing(with: .one, where: sortedCodes .== 123)
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
          equilibriumAngles = SIMD3(108.300, 108.900, 109.000)
          
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
            // from the column four cells below: 19-22-22. Anything connected to
            // the other side of a cyclopropane carbon (60°) should have an angle
            // like 120°. This is not the first typo I have caught in one of
            // Allinger's research papers, see the note about the MM4 formula for
            // the Torsion-Stretch cross-term.
            //
            // The stiffness does match up. Extrapolating the ratios of 1-1-19 :
            // 19-19-19 and 1-19-1 : 19-19-19 from 5-membered ring variants, one
            // gets 0.233 and 0.236 respectively for 19-19-19. That is very close
            // to 0.25, so I don't think that was messed up.
            bendingStiffnesses = SIMD3(repeating: 0.250)
            equilibriumAngles = SIMD3(109.50, 110.80, 111.20)
          } else {
            bendingStiffnesses = SIMD3(repeating: 0.320)
            equilibriumAngles = SIMD3(repeating: 106.00)
          }
          
          // Phosphorus
        case (1, 25, 1):
          bendingStiffnesses = SIMD3(0.900, 0.725, .nan)
          equilibriumAngles = SIMD3(94.50, 97.90, .nan)
        case (1, 1, 25):
          bendingStiffnesses = SIMD3(0.750, 0.825, 0.725)
          equilibriumAngles = SIMD3(107.05, 108.25, 109.55)
          
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
        case (5, 31, 5):
          bendingStiffnesses = SIMD3(repeating: 0.423)
          equilibriumAngles = SIMD3(107.5, 108.5, 109.5)
        case (1, 31, 5):
          bendingStiffnesses = SIMD3(repeating: 0.390)
          equilibriumAngles = SIMD3(110.2, 110.5, 111.5)
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
        case (1, 1, 31):
          if ringType == 5 {
            bendingStiffnesses = SIMD3(repeating: 0.680)
            equilibriumAngles = SIMD3(repeating: 105.5)
          } else {
            bendingStiffnesses = SIMD3(repeating: 0.450)
            equilibriumAngles = SIMD3(repeating: 109.4)
          }
          // TODO: Add self-bonding parameters from Tinker. Extrapolate the
          // 31-1-31 parameter using silicon.
        default:
          break
        }
      }
      guard let bendingStiffnesses,
            let equilibriumAngles else {
        fatalError("Unrecognized angle: \(sortedCodes)")
      }
      
      // Factors in both the center type and the other atoms in the angle.
      var angleType: Int
      if sortedCodes[1] == 15 {
        angleType = 0
      } else {
        var matchMask: SIMD3<UInt8> = .zero
        matchMask.replace(with: .one, where: sortedCodes .== 5)
        let numHydrogens = Int(matchMask.wrappedSum())
        
        guard let centerType = atoms.centerTypes[Int(angle[1])] else {
          fatalError("Angle did not occur at tetravalent atom.")
        }
        switch centerType {
        case .quaternary:
          angleType = 1 - numHydrogens
        case .tertiary:
          angleType = 2 - numHydrogens
        case .secondary:
          angleType = 3 - numHydrogens
        case .primary:
          angleType = 4 - numHydrogens
        }
      }
      
      // MARK: - Bend-Bend, Stretch-Bend, Stretch-Stretch
      
      let originalCodes = sortedCodes
      sortedCodes.replace(with: .one, where: sortedCodes .== 123)
      sortedCodes = sortAngle(sortedCodes)
      
      var bendBendStiffness: Float
      var stretchBendStiffness: Float
      var stretchBendStiffness2: Float?
      var stretchStretchStiffness: Float?
      
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
          fatalError("Unrecognized oxygen angle: \(sortedCodes)")
        }
      } else if any(sortedCodes .== 11) {
        // Fluorine
        precondition(
          sortedCodes[2] == 11, "Unrecognized fluorine angle: \(sortedCodes)")
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
          fatalError("Unrecognized fluorine angle: \(sortedCodes)")
        }
      } else {
        switch sortedCodes[1] {
          // Carbon
        case 1:
          // Assume the MM4 paper's parameters for H-C-C/C-C-C also apply to
          // H-C-Si/C-C-Si/Si-C-Si.
          if any(sortedCodes .== 5) {
            bendBendStiffness = 0.350
            stretchBendStiffness = 0.100
          } else {
            bendBendStiffness = 0.204
            stretchBendStiffness = (ringType == 5) ? 0.180 : 0.140
            
            // Exception for the nitrogen-containing bond C-C-N.
            if sortedCodes[0] == 1, sortedCodes[2] == 8 {
              stretchStretchStiffness = -0.10
            }
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
            // ring to that of 1-8-1. Often, the value inside the ring is larger
            // than in typical bonds. Not an order of magnitude smaller.
            //
            // On second analysis, this parameter seems completely messed up.
            // It provides zero information, as the fallback would be 1-8-1, not
            // 1-8-123. Also, the table is malformatted there. Setting it to
            // 0.04, the same as the fallback, would be the safest choice.
            //
            // Perhaps the purpose was to clear up any confusion, as the
            // fallback rule was described in a different paper. In the
            // bend-bend section, a 1-8-123 parameter doesn't exist, so no
            // extra explanation is required to specify the 123-8-123 case.
            // Still, the presence of the 5-ring restriction is very ambiguous.
            stretchBendStiffness = (ringType == 5) ? 0.04 : 0.04
          } else {
            fatalError("Unrecognized nitrogen angle: \(sortedCodes)")
          }
          
          // Silicon
        case 19:
          if any(sortedCodes .== 5) {
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
          // applied to Si, P, and PO4 in MM3. They were removed from the MM4
          // paper, except a new bend-bend parameter for H-P-H. I don't allow
          // H-P-H angles in this forcefield.
          //
          // I assume this omission was intentional. The creators knew the
          // parameters existed, and they talked with Allinger about it. They
          // made a decision that the parameters weren't necessary, which is
          // generally good practice to avoid overfitting a forcefield.
          bendBendStiffness = 0
          stretchBendStiffness = 0
          
          // Sulfur
        case 15:
          // There cannot be a bend-bend interaction around a divalent sulfur.
          bendBendStiffness = 0.000
          if all(sortedCodes .== SIMD3(1, 15, 1)) {
            stretchBendStiffness = (ringType == 5) ? 0.280 : 0.150
          } else {
            fatalError("Unrecognized sulfur angle: \(sortedCodes)")
          }
          
          // Germanium
        case 31:
          // The parameters that Allinger created for MM3(2000) do not list
          // germanium under bend-bend parameters. Like sulfur, it is almost
          // certainly zero.
          bendBendStiffness = 0.000
          if any(sortedCodes .== 5) {
            stretchBendStiffness = 0.000
          } else {
            stretchBendStiffness = 0.450
          }
          
          
        default:
          fatalError("Unrecognized angle: \(sortedCodes)")
        }
      }
      
      guard !bendingStiffnesses[angleType - 1].isNaN,
            !equilibriumAngles[angleType - 1].isNaN else {
        fatalError("Angle parameter was NaN for angle: \(sortedCodes)")
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
