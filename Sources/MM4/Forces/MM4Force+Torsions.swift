//
//  MM4Force+Torsions.swift
//  MM4
//
//  Created by Philip Turner on 10/8/23.
//

import OpenMM

/// Torsion and torsion-stretch force.
class MM4TorsionForce: MM4Force {
  required init(system: MM4System, descriptor: MM4ForceFieldDescriptor) {
    // Eventually, we want to optimize this by fusing all the torsions
    // surrounding the same bond into a single invocation.
    func createForce() -> OpenMM_CustomCompoundBondForce {
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
      return force
    }
    var force: OpenMM_CustomCompoundBondForce!
    
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
      if parameters.V1 == 0,
         parameters.Vn == 0,
         parameters.V3 == 0,
         parameters.Kts3 == 0 {
        continue
      }
      
      // Units: kcal/mol -> kJ/mol
      //          kJ/mol -> zJ
      //
      // WARNING: Divide all Vn torsion parameters by 2.
      let unitConversionFactor = OpenMM_KJPerKcal * MM4ZJPerKJPerMol
      array[0] = Double(parameters.V1) * unitConversionFactor / 2
      array[1] = Double(parameters.Vn) * unitConversionFactor / 2
      array[2] = Double(parameters.V3) * unitConversionFactor / 2
      array[3] = Double(parameters.n)
      
      // Units: kcal/mol/angstrom -> kJ/mol/angstrom
      //          kJ/mol/angstrom -> kJ/mol/nm
      //                kJ/mol/nm -> zJ/nm
      var Kts3 = Double(parameters.Kts3)
      Kts3 *= OpenMM_KJPerKcal
      Kts3 /= OpenMM_NmPerAngstrom
      Kts3 *= MM4ZJPerKJPerMol
      array[4] = Kts3
      
      // Assume the bond is already sorted.
      let bond = SIMD2(torsion[1], torsion[2])
      guard let bondID = bonds.map[bond] else {
        fatalError("Invalid bond.")
      }
      
      // Units: angstrom -> nm
      let bondParameters = bonds.parameters[Int(bondID)]
      var equilibriumLength = Double(bondParameters.equilibriumLength)
      equilibriumLength *= OpenMM_NmPerAngstrom
      array[5] = equilibriumLength
      
      let reorderedTorsion = system.reorder(torsion)
      for lane in 0..<4 {
        particles[lane] = reorderedTorsion[lane]
      }
      
      if force == nil {
        force = createForce()
      }
      force.addBond(particles: particles, parameters: array)
    }
    super.init(forces: [force], forceGroup: 1)
  }
}

/// Torsion, torsion-stretch, torsion-bend, and bend-torsion-bend force.
class MM4TorsionExtendedForce: MM4Force {
  required init(system: MM4System, descriptor: MM4ForceFieldDescriptor) {
    // When fusing multiple carbon-like torsions into a single invocation: V4/V6
    // terms are very rare. There should be a separate force to handle only the
    // extended torsions with V4 or V6 terms.
    func createForce() -> OpenMM_CustomCompoundBondForce {
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
      return force
    }
    var force: OpenMM_CustomCompoundBondForce!
    
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
      
      do {
        // Units: kcal/mol -> kJ/mol
        //          kJ/mol -> zJ
        //
        // WARNING: Divide all Vn torsion parameters by 2.
        let unitConversionFactor = OpenMM_KJPerKcal * MM4ZJPerKJPerMol
        array[0] = Double(originalParameters.V1) * unitConversionFactor / 2
        array[1] = Double(originalParameters.Vn) * unitConversionFactor / 2
        array[2] = Double(originalParameters.V3) * unitConversionFactor / 2
        array[3] = Double(parameters.V4) * unitConversionFactor / 2
        array[4] = Double(parameters.V6) * unitConversionFactor / 2
      }
      
      // MARK: - Torsion-Stretch
      
      func addTorsionStretch(
        _ index: Int,
        _ tuple: (`left`: Float, central: Float, `right`: Float)
      ) {
        // Units: kcal/mol/angstrom -> kJ/mol/angstrom
        //          kJ/mol/angstrom -> kJ/mol/nm
        //                kJ/mol/nm -> zJ/nm
        var params = SIMD3<Double>(
          SIMD3(tuple.left, tuple.central, tuple.right))
        params *= OpenMM_KJPerKcal
        params /= OpenMM_NmPerAngstrom
        params *= MM4ZJPerKJPerMol
        array[5 + index + 0] = params[0]
        array[5 + index + 3] = params[1]
        array[5 + index + 6] = params[2]
      }
      addTorsionStretch(0, parameters.Kts1)
      addTorsionStretch(1, parameters.Kts2)
      addTorsionStretch(2, parameters.Kts3)
      
      func createEquilibriumLength(_ index: Int) -> Double {
        var bond = SIMD2(torsion[index], torsion[index + 1])
        bond = system.parameters.sortBond(bond)
        guard let bondID = bonds.map[bond] else {
          fatalError("Invalid bond.")
        }
        
        // Units: angstrom -> nm
        let parameters = bonds.parameters[Int(bondID)]
        var equilibriumLength = Double(parameters.equilibriumLength)
        equilibriumLength *= OpenMM_NmPerAngstrom
        return equilibriumLength
      }
      array[14] = createEquilibriumLength(0)
      array[15] = createEquilibriumLength(1)
      array[16] = createEquilibriumLength(2)
      
      // MARK: - Torsion-Bend, Bend-Torsion-Bend
      
      // Bend-torsion-bend should remain unchanged (0.043828 converts directly
      // to attojoules, and should not be divided by 2. Same with torsion-bend,
      // which converts directly without dividing by 2.
      
      func addTorsionBend(
        _ index: Int,
        _ tuple: (left: Float, right: Float)
      ) {
        // Units: millidyne-angstrom/rad -> kJ/mol/rad
        //                    kJ/mol/rad -> zJ/rad
        var params = SIMD2<Double>(
          SIMD2(tuple.left, tuple.right))
        params *= MM4KJPerMolPerAJ
        params *= MM4ZJPerKJPerMol
        array[17 + index + 0] = params[0]
        array[17 + index + 3] = params[1]
      }
      addTorsionBend(0, parameters.Ktb1)
      addTorsionBend(1, parameters.Ktb2)
      addTorsionBend(2, parameters.Ktb3)
      
      do {
        // Units: millidyne-angstrom/rad^2 -> kJ/mol/rad^2
        //                    kJ/mol/rad^2 -> zJ/rad^2
        var Kbtb = Double(parameters.Kbtb)
        Kbtb *= MM4KJPerMolPerAJ
        Kbtb *= MM4ZJPerKJPerMol
        array[23] = Kbtb
      }
      
      func createEquilibriumAngle(_ index: Int) -> Double {
        var angle = SIMD3(
          torsion[index], torsion[index + 1], torsion[index + 2])
        angle = system.parameters.sortAngle(angle)
        guard let angleID = angles.map[angle] else {
          fatalError("Invalid bond.")
        }
        
        // Units: degree -> rad
        let parameters = angles.parameters[Int(angleID)]
        var equilibriumAngle = Double(parameters.equilibriumAngle)
        equilibriumAngle *= OpenMM_RadiansPerDegree
        return equilibriumAngle
      }
      array[24] = createEquilibriumAngle(0)
      array[25] = createEquilibriumAngle(1)
      
      let reorderedTorsion = system.reorder(torsion)
      for lane in 0..<4 {
        particles[lane] = reorderedTorsion[lane]
      }
      
      if force == nil {
        force = createForce()
      }
      force.addBond(particles: particles, parameters: array)
    }
    super.init(forces: [force], forceGroup: 1)
  }
}
