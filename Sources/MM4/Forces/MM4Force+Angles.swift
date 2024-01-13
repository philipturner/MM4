//
//  MM4Force+Angles.swift
//  MM4
//
//  Created by Philip Turner on 10/8/23.
//

import Foundation
import OpenMM

/// Angle bend and stretch-bend force.
class MM4BendForce: MM4Force {
  required init(system: MM4System, descriptor: MM4ForceFieldDescriptor) {
    // https://www.desmos.com/calculator/shl9ovintw
    //
    // Leaving the formula as-is doesn't create a monotonically increasing
    // function until 180° (red). Green is what I previously had with the old
    // MM4. Scaling the X axis into degrees (blue), it does seem like the
    // correct function. There is a correction in the MM3(2000) parameters to
    // use 0.000000022 for the sextic term, so I will enter it there. I don't
    // see how the higher-order terms contribute at all to angle stiffness; the
    // curve looks exactly like that for the quadratic version.
    //
    // Except:
    // - correcting the 9e-10 parameter to 2.2e-8 completely changed the red
    //   curve's behavior. Look at the graph here:
    //   https://www.desmos.com/calculator/vaoiwzbxpi
    // - now, the red curve does monotonically increase toward infinity, and
    //   the shape is very obviously not a perfect parabola
    // - the plot still has much too low overall magnitude. Originally, I
    //   created the (180/pi)^2 pre-multiplier hack because the angle force
    //   was so low, it didn't even exist. Scaling the angle from radians to
    //   degrees beforehand will have the same effect as making the curve
    //   steeper; the same energy maximum (rise) must be reached with much
    //   shorter run (180 -> pi)
    // - converting 0.021914 to 71.94, the same value as the bond stretching
    //   coefficient, gets us the correct magnitude
    // - scale the remaining, relative multipliers by increasing powers of
    //   the (180/pi) correction factor
    //
    // New energy formula
    // - https://www.desmos.com/calculator/p4fyruchbi
    // - x-axis is radians
    // - y-axis is attojoules
    // - the y-axis will be multiplied by the verbatim parameter for bending
    //   stiffness, which is 0.3-0.4 for carbon
    // - 0.15 attojoules with one radian of bending seems very reasonable
    //
    // Other notes:
    // - ignore the units in the MM4 formula; they actually return something in
    //   Kcal/mol
    // - 143.88 is the conversion factor from aJ -> kcal/mol
    // - 71.94 = 143.88 / 2
    //
    func createForce() -> OpenMM_CustomCompoundBondForce {
      let correction = 180 / Float.pi
      // let bendingStiffness = /*71.94*/ 1.00 * bendingStiffness
      let cubicTerm = 0.014 * correction
      let quarticTerm = 5.6e-5 * pow(correction, 2)
      let quinticTerm = 7.0e-7 * pow(correction, 3)
      let sexticTerm = 2.2e-8 * pow(correction, 4)
      
      let force = OpenMM_CustomCompoundBondForce(numParticles: 3, energy: """
      bend + stretchBend;
      bend = bendingStiffness * deltaTheta^2 * (
        1
        - \(cubicTerm) * deltaTheta
        + \(quarticTerm) * deltaTheta^2
        - \(quinticTerm) * deltaTheta^3
        + \(sexticTerm) * deltaTheta^4
      );
      stretchBend = stretchBendStiffness * deltaTheta * (
        deltaLengthLeft + deltaLengthRight
      );
      
      deltaTheta = angle(p1, p2, p3) - equilibriumAngle;
      deltaLengthLeft = distance(p1, p2) - equilibriumLengthLeft;
      deltaLengthRight = distance(p3, p2) - equilibriumLengthRight;
      """)
      force.addPerBondParameter(name: "bendingStiffness")
      force.addPerBondParameter(name: "equilibriumAngle")
      force.addPerBondParameter(name: "stretchBendStiffness")
      force.addPerBondParameter(name: "equilibriumLengthLeft")
      force.addPerBondParameter(name: "equilibriumLengthRight")
      return force
    }
    var force: OpenMM_CustomCompoundBondForce!
    
    let particles = OpenMM_IntArray(size: 3)
    let array = OpenMM_DoubleArray(size: 5)
    let bonds = system.parameters.bonds
    let angles = system.parameters.angles
    for angleID in angles.indices.indices {
      let angle = angles.indices[angleID]
      let parameters = angles.parameters[angleID]
      guard parameters.bendingStiffness != 0 else {
        precondition(
          parameters.bendBendStiffness == 0,
          "Cannot have bend-bend force without bend force.")
        precondition(
          parameters.stretchBendStiffness == 0,
          "Cannot have stretch-bend force without bend force.")
        continue
      }
      
      // Units: millidyne-angstrom/rad^2 -> kJ/mol/rad^2
      //                    kJ/mol/rad^2 -> zJ/rad^2
      //
      // WARNING: 143 needs to be divided by 2 before it becomes 71.94.
      var bendingStiffness = Double(parameters.bendingStiffness)
      bendingStiffness *= MM4KJPerMolPerAJ
      bendingStiffness /= 2
      bendingStiffness *= MM4ZJPerKJPerMol
      
      // Units: degree -> rad
      var equilibriumAngle = Double(parameters.equilibriumAngle)
      equilibriumAngle *= OpenMM_RadiansPerDegree
      
      // Units: millidyne-angstrom/rad^2 -> kJ/mol/rad^2
      //                    kJ/mol/rad^2 -> zJ/rad^2
      //
      // This part does not need to be divided by 2; it was never divided by
      // 2 in the first place (2.5118 was used instead of 1.2559).
      var stretchBendStiffness = Double(parameters.stretchBendStiffness)
      stretchBendStiffness *= MM4KJPerMolPerAJ
      stretchBendStiffness *= MM4ZJPerKJPerMol
      
      // Units: angstrom -> nm
      let bondLeft = system.parameters.sortBond(SIMD2(angle[0], angle[1]))
      let bondRight = system.parameters.sortBond(SIMD2(angle[1], angle[2]))
      
      @_transparent
      func createLength(_ bond: SIMD2<UInt32>) -> Double {
        guard let bondID = bonds.map[bond] else {
          fatalError("Invalid bond.")
        }
        let parameters = bonds.parameters[Int(bondID)]
        var equilibriumLength = Double(parameters.equilibriumLength)
        equilibriumLength *= OpenMM_NmPerAngstrom
        return equilibriumLength
      }
      
      let reorderedAngle = system.reorder(angle)
      for lane in 0..<3 {
        particles[lane] = reorderedAngle[lane]
      }
      array[0] = bendingStiffness
      array[1] = equilibriumAngle
      array[2] = stretchBendStiffness
      array[3] = createLength(bondLeft)
      array[4] = createLength(bondRight)
      
      if force == nil {
        force = createForce()
      }
      force.addBond(particles: particles, parameters: array)
    }
    super.init(forces: [force], forceGroup: 2)
  }
}

/// Angle bend-bend force.
class MM4BendBendForce: MM4Force {
  required init(system: MM4System, descriptor: MM4ForceFieldDescriptor) {
    var includeBendBend = false
    for params in system.parameters.angles.parameters {
      if params.bendBendStiffness != 0 {
        includeBendBend = true
      }
    }
    guard includeBendBend else {
      super.init(forces: [], forceGroup: 2)
      return
    }
    
    // Sequence of indices for both generating energy expressions and fetching
    // angles during per-atom iteration.
    let indexSequence: [SIMD2<Int>] = [
      SIMD2(0, 1), SIMD2(1, 2), SIMD2(2, 0),
      SIMD2(0, 3), SIMD2(1, 3), SIMD2(2, 3)
    ]
    
    // Make a different optimized bend-bend force for trivalent and tetravalent
    // atoms.
    func createForce(
      valenceCount: Int, angleCount: Int
    ) -> OpenMM_CustomCompoundBondForce {
      var energy: String = ""
      var interactions: [SIMD2<Int>] = []
      for i in 0..<angleCount {
        for j in (i + 1)..<angleCount {
          interactions.append(SIMD2(i, j))
        }
      }
      energy += "-("
      for (index, interaction) in interactions.enumerated() {
        energy += "bendBend\(interaction[0])\(interaction[1])"
        if index < interactions.count - 1 {
          energy += " + "
        } else {
          energy += ");\n"
        }
      }
      for interaction in interactions {
        energy += "bendBend\(interaction[0])\(interaction[1])"
        energy += " = term\(interaction[0]) * term\(interaction[1]);\n"
      }
      
      // The center particle is always first, so shift the index sequence by +1
      // before entering into the energy expression.
      for i in 0..<angleCount {
        energy += "term\(i) = bendBendStiffness\(i) * deltaTheta\(i);\n"
        energy += "deltaTheta\(i) = "
        
        let sequence = indexSequence[i]
        energy += "angle(p\(1 + sequence[0]), p0, p\(1 + sequence[1]))"
        energy += " - equilibriumAngle\(i);\n";
      }
      return OpenMM_CustomCompoundBondForce(
        numParticles: 1 + valenceCount, energy: energy)
    }
    
    // These is no bend-bend force for divalent atoms.
    var trivalentForce: OpenMM_CustomCompoundBondForce!
    var tetravalentForce: OpenMM_CustomCompoundBondForce!
    
    let particleArrays = [
      OpenMM_IntArray(size: 1 + 3),
      OpenMM_IntArray(size: 1 + 4)
    ]
    let parameterArrays = [
      OpenMM_DoubleArray(size: 3 * 2),
      OpenMM_DoubleArray(size: 6 * 2)
    ]
    let atoms = system.parameters.atoms
    let angles = system.parameters.angles
    
    for atomID in atoms.indices {
      let atomicNumber = atoms.atomicNumbers[atomID]
      var valenceCount: Int
      switch atomicNumber {
      case 1: valenceCount = 1 // H
      case 6: valenceCount = 4 // C
      case 7: valenceCount = 3 // N
      case 8: valenceCount = 2 // O
      case 9: valenceCount = 1 // F
      case 14: valenceCount = 4 // Si
      case 15: valenceCount = 3 // P
      case 16: valenceCount = 2 // S
      case 32: valenceCount = 4 // Ge
      default: fatalError("Atomic number not recognized: \(atomicNumber)")
      }
      if valenceCount < 3 {
        continue
      }
      
      let angleCount = (valenceCount == 3) ? 3 : 6
      let arrayIndex = (valenceCount == 3) ? 0 : 1
      let particles = particleArrays[arrayIndex]
      let map = system.parameters.atomsToAtomsMap[atomID]
      particles[0] = Int(system.reorderedIndices[atomID])
      
      for i in 0..<valenceCount {
        let reorderedIndex = system.reorderedIndices[Int(map[i])]
        particles[1 + i] = Int(reorderedIndex)
      }
      
      let array = parameterArrays[arrayIndex]
      var includeParticles = false
      for i in 0..<angleCount {
        let sequence = indexSequence[i]
        var angle = SIMD3(
          UInt32(truncatingIfNeeded: map[sequence[0]]),
          UInt32(truncatingIfNeeded: atomID),
          UInt32(truncatingIfNeeded: map[sequence[1]]))
        angle = system.parameters.sortAngle(angle)
        guard let angleID = angles.map[angle] else {
          fatalError("Angle did not exist.")
        }
        
        // Units: millidyne-angstrom/rad^2 -> kJ/mol/rad^2
        //                    kJ/mol/rad^2 -> zJ/rad^2
        let parameters = angles.parameters[Int(angleID)]
        var bendBendStiffness = Double(parameters.bendBendStiffness)
        bendBendStiffness *= MM4KJPerMolPerAJ
        bendBendStiffness *= MM4ZJPerKJPerMol
        array[2 * i + 0] = bendBendStiffness
        
        // Units: degree -> rad
        var equilibriumAngle = Double(parameters.equilibriumAngle)
        equilibriumAngle *= OpenMM_RadiansPerDegree
        array[2 * i + 1] = equilibriumAngle
        
        if bendBendStiffness != 0 {
          includeParticles = true
        }
      }
      guard includeParticles else {
        continue
      }
      
      if arrayIndex == 0 {
        if trivalentForce == nil {
          trivalentForce = createForce(valenceCount: 3, angleCount: 3)
        }
        trivalentForce.addBond(particles: particles, parameters: array)
      } else {
        if tetravalentForce == nil {
          tetravalentForce = createForce(valenceCount: 4, angleCount: 6)
        }
        tetravalentForce.addBond(particles: particles, parameters: array)
      }
    }
    super.init(forces: [trivalentForce, tetravalentForce], forceGroup: 2)
  }
}

/// Type 2 stretch-bend and stretch-stretch force.
class MM4BendExtendedForce: MM4Force {
  required init(system: MM4System, descriptor: MM4ForceFieldDescriptor) {
    func createForce() -> OpenMM_CustomCompoundBondForce {
      let force = OpenMM_CustomCompoundBondForce(numParticles: 5, energy: """
      stretchBend + stretchStretch;
      stretchBend = stretchBendStiffness * deltaTheta * (
        deltaLength4 + deltaLength5
      );
      stretchStretch = stretchStretchStiffness * (
        deltaLength2 * deltaLength3
      );
      
      deltaTheta = angle(p2, p1, p3) - equilibriumAngle;
      deltaLength2 = distance(p1, p2) - equilibriumLength2;
      deltaLength3 = distance(p1, p3) - equilibriumLength3;
      deltaLength4 = distance(p1, p4) - equilibriumLength4;
      deltaLength5 = distance(p1, p5) - equilibriumLength5;
      """)
      force.addPerBondParameter(name: "stretchBendStiffness")
      force.addPerBondParameter(name: "stretchStretchStiffness")
      force.addPerBondParameter(name: "equilibriumAngle")
      force.addPerBondParameter(name: "equilibriumLength2")
      force.addPerBondParameter(name: "equilibriumLength3")
      force.addPerBondParameter(name: "equilibriumLength4")
      force.addPerBondParameter(name: "equilibriumLength5")
      return force
    }
    var force: OpenMM_CustomCompoundBondForce!
    
    // Iterate by angles instead of by atoms this time.
    let particles = OpenMM_IntArray(size: 5)
    let array = OpenMM_DoubleArray(size: 7)
    let bonds = system.parameters.bonds
    let angles = system.parameters.angles
    for angleID in angles.indices.indices {
      let angle = angles.indices[angleID]
      guard let parameters = angles.extendedParameters[angleID] else {
        continue
      }
      let originalParameters = angles.parameters[angleID]
      
      // Units: aJ/... -> kJ/mol/...
      //    kJ/mol/... -> zJ/...
      var stretchBendStiffness = Double(parameters.stretchBendStiffness)
      stretchBendStiffness *= MM4KJPerMolPerAJ
      stretchBendStiffness *= MM4ZJPerKJPerMol
      array[0] = stretchBendStiffness
      
      // Units: aJ/... -> kJ/mol/...
      //    kJ/mol/... -> zJ/...
      var stretchStretchStiffness = Double(parameters.stretchStretchStiffness)
      stretchStretchStiffness *= MM4KJPerMolPerAJ
      stretchStretchStiffness *= MM4ZJPerKJPerMol
      array[1] = stretchStretchStiffness
      
      // Units: degree -> rad
      var equilibriumAngle = Double(originalParameters.equilibriumAngle)
      equilibriumAngle *= OpenMM_RadiansPerDegree
      array[2] = equilibriumAngle
      
      let map = system.parameters.atomsToAtomsMap[Int(angle[1])]
      var sortedMap: SIMD4<Int32> = .init(repeating: -1)
      var sortedStart: Int = 0
      var sortedEnd: Int = 3
      for lane in 0..<4 {
        guard map[lane] != -1 else {
          fatalError("Unexpected behavior creating bend extended force.")
        }
        if angle[0] == map[lane] || angle[2] == map[lane] {
          sortedMap[sortedStart] = map[lane]
          sortedStart += 1
        } else {
          sortedMap[sortedEnd] = map[lane]
          sortedEnd -= 1
        }
      }
      guard sortedStart == 2, sortedEnd == 2 else {
        // There should be two atoms from the angle, two not from the angle.
        fatalError("Unexpected number of atoms matched angle.")
      }
      
      let particleID0 = Int(angle[1])
      particles[0] = Int(system.reorderedIndices[particleID0])
      for i in 0..<4 {
        let reorderedIndex = system.reorderedIndices[Int(map[i])]
        particles[1 + i] = Int(reorderedIndex)
        
        
        let bondPart1: UInt32 = angle[1]
        let bondPart2: UInt32 = UInt32(truncatingIfNeeded: sortedMap[i])
        var bond = SIMD2<UInt32>(bondPart1, bondPart2)
        bond = system.parameters.sortBond(bond)
        guard let bondID = bonds.map[bond] else {
          fatalError("Bond was not in the map.")
        }
        
        // Units: angstrom -> nm
        let parameters = bonds.parameters[Int(bondID)]
        var equilibriumLength = Double(parameters.equilibriumLength)
        equilibriumLength *= OpenMM_NmPerAngstrom
        array[2 + i] = equilibriumLength
      }
      
      if force == nil {
        force = createForce()
      }
      force.addBond(particles: particles, parameters: array)
    }
    super.init(forces: [force], forceGroup: 2)
  }
}
