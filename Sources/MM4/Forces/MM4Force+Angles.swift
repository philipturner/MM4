//
//  MM4Force+Angles.swift
//
//
//  Created by Philip Turner on 10/8/23.
//

import OpenMM

/// Angle bend and stretch-bend force.
class MM4BendForce: MM4Force {
  required init(system: MM4System) {
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
    let correction = 180 / Float.pi
    // let bendingStiffness = /*71.94*/ 1.00 * bendingStiffness
    let cubicTerm = 0.014 * correction
    let quarticTerm = 5.6e-5 * correction * correction
    let quinticTerm = 7.0e-7 * correction * correction * correction
    let sexticTerm = 2.2e-8  * correction * correction * correction * correction
    
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
      deltaLengthRight = distance(p1, p2) - equilibriumLengthRight;
      """)
    force.addPerBondParameter(name: "bendingStiffness")
    force.addPerBondParameter(name: "equilibriumAngle")
    force.addPerBondParameter(name: "stretchBendStiffness")
    force.addPerBondParameter(name: "equilibriumLengthLeft")
    force.addPerBondParameter(name: "equilibriumLengthRight")
    
    let particles = OpenMM_IntArray(size: 3)
    let array = OpenMM_DoubleArray(size: 5)
    let bonds = system.parameters.bonds
    let angles = system.parameters.angles
    for angleID in angles.indices.indices {
      let angle = angles.indices[angleID]
      let parameters = angles.parameters[angleID]
      
      // Units: millidyne-angstrom/rad^2 -> kJ/mol/rad^2
      //
      // WARNING: 143 needs to be divided by 2 before it becomes 71.94.
      var bendingStiffness = Double(parameters.bendingStiffness)
      bendingStiffness *= MM4_KJPerMolPerAJ
      bendingStiffness /= 2
      
      // Units: degree -> rad
      var equilibriumAngle = Double(parameters.equilibriumAngle)
      equilibriumAngle *= OpenMM_RadiansPerDegree
      
      // Units: millidyne-angstrom/rad^2 -> kJ/mol/rad^2
      //
      // This part does not need to be divided by 2; it was never divided by
      // 2 in the first place (2.5118 was used instead of 1.2559).
      var stretchBendStiffness = Double(parameters.stretchBendStiffness)
      stretchBendStiffness *= MM4_KJPerMolPerAJ
      
      // Units: angstrom -> nm
      let bondLeft = system.parameters.sortBond(SIMD2(angle[0], angle[1]))
      let bondRight = system.parameters.sortBond(SIMD2(angle[1], angle[2]))
      
      @inline(__always)
      func createLength(_ bond: SIMD2<Int32>) -> Double {
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
      force.addBond(particles: particles, parameters: array)
    }
    super.init(force: force, forceGroup: 2)
  }
}

/// Angle bend-bend force.
///
/// This object may be deleted during a future optimization.
class MM4BendBendForce: MM4Force {
  
}

class MM4BendExtendedForc: MM4Force {
  
}
