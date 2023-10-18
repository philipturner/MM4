//
//  MM4Force+Torsions.swift
//
//
//  Created by Philip Turner on 10/8/23.
//

import OpenMM

/// Torsion and torsion-stretch force.
class MM4TorsionForce: MM4Force {
  required init(system: MM4System) {
    // Eventually, we want to optimize this by fusing all the torsions
    // surrounding the same bond into a single invocation.
    //
    // Perhaps only do this for bonds that have no heteroatoms immediately next
    // to them. This should optimize bulk diamond and moissanite, without
    // introducing extra combination for highly mixed-element compounds. There
    // should be a separate force object for this.
    let force = OpenMM_CustomCompoundBondForce(numParticles: 4, energy: """
      torsion + torsionStretch;
      torsion = V1 * fourierExpansion1
              + Vn * fourierExpansionN
              + V3 * fourierExpansion3;
      torsionStretch = Kts3 * deltaLength * fourierExpansion3;
      
      deltaLength = distance(p2, p3) - equilibriumLength;
      fourierExpansion1 = 1 + cos(omega);
      fourierExpansionN = 1 - cos(n * omega);
      fourierExpansion3 = 1 + cos(3 * omega);
      omega = dihedral(p1, p2, p3, p4);
      """)
    force.addPerBondParameter(name: "V1")
    force.addPerBondParameter(name: "Vn")
    force.addPerBondParameter(name: "V3")
    force.addPerBondParameter(name: "n")
    force.addPerBondParameter(name: "Kts3")
    force.addPerBondParameter(name: "equilibriumLength")
    
    let particles = OpenMM_IntArray(size: 4)
    let array = OpenMM_DoubleArray(size: 6)
    let bonds = system.parameters.bonds
    let torsions = system.parameters.torsions
    for torsionID in torsions.indices.indices {
      let torsion = torsions.indices[torsionID]
      let parameters = torsions.parameters[torsionID]
      guard parameters.n.truncatingRemainder(dividingBy: 2) == 0 else {
        fatalError("'n' must be an even integer.")
      }
      if torsions.extendedParameters[torsionID] != nil {
        // Compute heteroatom V1/V2/V3 terms in the extended force.
        continue
      }
      
      // Units: kcal/mol -> kJ/mol
      //
      // WARNING: Divide all Vn torsion parameters by 2.
      array[0] = Double(parameters.V1) * OpenMM_KJPerKcal / 2
      array[1] = Double(parameters.Vn) * OpenMM_KJPerKcal / 2
      array[2] = Double(parameters.V3) * OpenMM_KJPerKcal / 2
      array[3] = Double(parameters.n)
      
      // Units: kcal/mol/angstrom -> kJ/mol/nm
      var Kts3 = Double(parameters.Kts3)
      Kts3 *= OpenMM_KJPerKcal
      Kts3 /= OpenMM_NmPerAngstrom
      array[4] = Kts3
      
      // Assume the bond is already sorted.
      let bond = SIMD2(torsion[1], torsion[2])
      guard let bondID = bonds.map[bond] else {
        fatalError("Invalid bond.")
      }
      let bondParameters = bonds.parameters[Int(bondID)]
      
      // No change in units for stretching stiffness.
      let stretchingStiffness = Double(bondParameters.stretchingStiffness)
      Kts3 /= stretchingStiffness
      
      // Units: angstrom -> nm
      var equilibriumLength = Double(bondParameters.equilibriumLength)
      equilibriumLength *= OpenMM_NmPerAngstrom
      array[5] = equilibriumLength
      
      let reorderedTorsion = system.reorder(torsion)
      for lane in 0..<4 {
        particles[lane] = reorderedTorsion[lane]
      }
      force.addBond(particles: particles, parameters: array)
    }
    super.init(forces: [force], forceGroup: 1)
  }
}

/// Torsion, torsion-stretch, torsion-bend, and bend-torsion-bend force.
class MM4TorsionExtendedForce: MM4Force {
  required init(system: MM4System) {
    // When fusing multiple carbon-like torsions into a single invocation: V4/V6
    // terms are very rare. There should be a separate force to handle only the
    // extended torsions with V4 or V6 terms.
    let force = OpenMM_CustomCompoundBondForce(numParticles: 4, energy: """
      torsion + torsionStretch + torsionBend + bendTorsionBend;
      torsion = V1 * fourierExpansion1
              + V2 * fourierExpansion2
              + V3 * fourierExpansion3
              + V4 * fourierExpansion4
              + V6 * fourierExpansion6;
      fourierExpansion4 = 1 - cos(4 * omega);
      fourierExpansion6 = 1 - cos(6 * omega);
      
      torsionStretch = torsionStretchLeft   * deltaLengthLeft
                     + torsionStretchCenter * deltaLengthCenter
                     + torsionStretchRight  * deltaLengthRight;
      torsionStretchLeft   = Kts1Left * fourierExpansion1
                           + Kts2Left * fourierExpansion2
                           + Kts3Left * fourierExpansion3;
      torsionStretchCenter = Kts1Center * fourierExpansion1
                           + Kts2Center * fourierExpansion2
                           + Kts3Center * fourierExpansion3;
      torsionStretchRight  = Kts1Right * fourierExpansion1
                           + Kts2Right * fourierExpansion2
                           + Kts3Right * fourierExpansion3;
      deltaLengthLeft   = distance(p1, p2) - equilibriumLengthLeft;
      deltaLengthCenter = distance(p2, p3) - equilibriumLengthCenter;
      deltaLengthRight  = distance(p3, p4) - equilibriumLengthRight;
      
      torsionBend = torsionBendLeft  * deltaThetaLeft
                  + torsionBendRight * deltaThetaRight;
      torsionBendLeft  = Ktb1Left * fourierExpansion1
                       + Ktb2Left * fourierExpansion2
                       + Ktb3Left * fourierExpansion3;
      torsionBendRight = Ktb1Right * fourierExpansion1
                       + Ktb2Right * fourierExpansion2
                       + Ktb3Right * fourierExpansion3;
      bendTorsionBend  = Kbtb * (fourierExpansion1 - 1)
                       * deltaThetaLeft * deltaThetaRight;
      deltaThetaLeft = angle(p1, p2, p3) - equilibriumAngleLeft;
      deltaThetaRight = angle(p2, p3, p4) - equilibriumAngleRight;
      
      fourierExpansion1 = 1 + cos(omega);
      fourierExpansion2 = 1 - cos(2 * omega);
      fourierExpansion3 = 1 + cos(3 * omega);
      omega = dihedral(p1, p2, p3, p4);
      """)
    force.addPerBondParameter(name: "V1")
    force.addPerBondParameter(name: "V2")
    force.addPerBondParameter(name: "V3")
    force.addPerBondParameter(name: "V4")
    force.addPerBondParameter(name: "V6")
    
    for position in ["Left", "Center", "Right"] {
      // Use the same technique to make code for appending the tuples for Kts/b
      // parameter more concise. Although that might require an inlined function
      // instead of a loop.
      force.addPerBondParameter(name: "Kts1\(position)")
      force.addPerBondParameter(name: "Kts2\(position)")
      force.addPerBondParameter(name: "Kts3\(position)")
    }
    for position in ["Left", "Center", "Right"] {
      force.addPerBondParameter(name: "equilibriumLength\(position)")
    }
    
    for position in ["Left", "Right"] {
      force.addPerBondParameter(name: "Ktb1\(position)")
      force.addPerBondParameter(name: "Ktb2\(position)")
      force.addPerBondParameter(name: "Ktb3\(position)")
    }
    force.addPerBondParameter(name: "Kbtb")
    for position in ["Left", "Right"] {
      force.addPerBondParameter(name: "equilibriumAngle\(position)")
    }
    
    let particles = OpenMM_IntArray(size: 4)
    let array = OpenMM_DoubleArray(
      size: 5 /*Vn*/ + 3 * 3 /*Kts*/ + 3 /*l_0*/
      + 2 * 3 /*Ktb*/ + 1 /*Kbtb*/ + 2 /*theta_0*/)
    let bonds = system.parameters.bonds
    let angles = system.parameters.angles
    let torsions = system.parameters.torsions
    for torsionID in torsions.indices.indices {
      guard let parameters = torsions.extendedParameters[torsionID] else {
        continue
      }
      let originalParameters = torsions.parameters[torsionID]
      let torsion = torsions.indices[torsionID]
      guard originalParameters.n == 2 else {
        fatalError("'n' must be two for extended torsions.")
      }
      
      // MARK: - Torsion
      
      // WARNING: Divide all Vn torsion parameters by 2.
      array[0] = Double(originalParameters.V1) * OpenMM_KJPerKcal / 2
      array[1] = Double(originalParameters.Vn) * OpenMM_KJPerKcal / 2
      array[2] = Double(originalParameters.V3) * OpenMM_KJPerKcal / 2
      array[3] = Double(parameters.V4) * OpenMM_KJPerKcal / 2
      array[4] = Double(parameters.V6) * OpenMM_KJPerKcal / 2
      
      // MARK: - Torsion-Stretch
      
      func createBondParams(_ index: Int) -> SIMD2<Float> {
        var bond = SIMD2(torsion[index], torsion[index + 1])
        bond = system.parameters.sortBond(bond)
        guard let bondID = bonds.map[bond] else {
          fatalError("Invalid bond.")
        }
        
        // No change in units for stretching stiffness.
        let bondParameters = bonds.parameters[Int(bondID)]
        return SIMD2(
          bondParameters.stretchingStiffness,
          bondParameters.equilibriumLength)
      }
      let bondParamsLeft = createBondParams(0)
      let bondParamsCenter = createBondParams(1)
      let bondParamsRight = createBondParams(2)
      
      // Units: kcal/mol/angstrom -> kJ/mol/nm
      var conversionFactors = SIMD3(
        repeating: OpenMM_KJPerKcal / OpenMM_NmPerAngstrom)
      let stiffnesses = SIMD3(
        bondParamsLeft[0],
        bondParamsRight[0],
        bondParamsCenter[0])
      conversionFactors /= SIMD3<Double>(stiffnesses)
      
      func addTorsionStretch(
        _ index: Int,
        _ tuple: (left: Float, central: Float, right: Float)
      ) {
        array[5 + index + 0] = Double(tuple.left) * conversionFactors[0]
        array[5 + index + 3] = Double(tuple.central) * conversionFactors[1]
        array[5 + index + 6] = Double(tuple.right) * conversionFactors[2]
      }
      addTorsionStretch(0, parameters.Kts1)
      addTorsionStretch(1, parameters.Kts2)
      addTorsionStretch(2, parameters.Kts3)
      
      // Units: angstrom -> nm
      array[14] = Double(bondParamsLeft[1]) * OpenMM_NmPerAngstrom
      array[15] = Double(bondParamsCenter[1]) * OpenMM_NmPerAngstrom
      array[16] = Double(bondParamsRight[1]) * OpenMM_NmPerAngstrom
      
      // MARK: - Torsion-Bend, Bend-Torsion-Bend
      
      // Bend-torsion-bend should remain unchanged (0.043828 converts directly
      // to attojoules, and should not be divided by 2. Same with torsion-bend,
      // which converts directly without dividing by 2.
      
      func createAngleParams(_ index: Int) -> Float {
        var angle = SIMD3(
          torsion[index], torsion[index + 1], torsion[index + 2])
        angle = system.parameters.sortAngle(angle)
        guard let angleID = angles.map[angle] else {
          fatalError("Invalid bond.")
        }
        
        let angleParameters = angles.parameters[Int(angleID)]
        return angleParameters.equilibriumAngle
      }
      let angleParamsLeft = createAngleParams(0)
      let angleParamsRight = createAngleParams(1)
      
      // Units: millidyne-angstrom/rad -> kJ/mol/rad
      func addTorsionBend(
        _ index: Int,
        _ tuple: (left: Float, right: Float)
      ) {
        array[17 + index + 0] = Double(tuple.left) * MM4KJPerMolPerAJ
        array[17 + index + 3] = Double(tuple.right) * MM4KJPerMolPerAJ
      }
      addTorsionBend(0, parameters.Ktb1)
      addTorsionBend(1, parameters.Ktb2)
      addTorsionBend(2, parameters.Ktb3)
      
      // Units: millidyne-angstrom/rad^2 -> kJ/mol/rad^2
      array[23] = Double(parameters.Kbtb) * MM4KJPerMolPerAJ
      
      // Units: degree -> rad
      array[24] = Double(angleParamsLeft) * OpenMM_RadiansPerDegree
      array[25] = Double(angleParamsRight) * OpenMM_RadiansPerDegree
      
      let reorderedTorsion = system.reorder(torsion)
      for lane in 0..<4 {
        particles[lane] = reorderedTorsion[lane]
      }
      force.addBond(particles: particles, parameters: array)
    }
    super.init(forces: [force], forceGroup: 2)
  }
}
