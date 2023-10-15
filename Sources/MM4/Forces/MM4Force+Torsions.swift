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
    let force = OpenMM_CustomCompoundBondForce(numParticles: 4, energy: """
      torsion + torsionStretch;
      torsion = V1 * fourierExpansion1
              + Vn * fourierExpansionN
              + V3 * fourierExpansion3;
      torsionStretch = Kts3 * deltaLength * fourierExpansion3;
      fourierExpansion1 = 1 + cos(omega);
      fourierExpansionN = 1 - cos(n * omega);
      fourierExpansion3 = 1 + cos(3 * omega);
      omega = dihedral(p1, p2, p3, p4);
      deltaLength = distance(p2, p3) - equilibriumLength;
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
      
      // Units: kcal/mol -> kJ/mol
      //
      // WARNING: Divide all Vn torsion parameters by 2.
      var V1 = Double(parameters.V1)
      var Vn = Double(parameters.Vn)
      var V3 = Double(parameters.V3)
      V1 *= OpenMM_KJPerKcal / 2
      Vn *= OpenMM_KJPerKcal / 2
      V3 *= OpenMM_KJPerKcal / 2
      
      // Units: kcal/mol/angstrom -> kJ/mol/nm
      var Kts3 = Double(parameters.Kts3)
      Kts3 *= OpenMM_KJPerKcal
      Kts3 /= OpenMM_NmPerAngstrom
      
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
      
      let reorderedTorsion = system.reorder(torsion)
      for lane in 0..<4 {
        particles[lane] = reorderedTorsion[lane]
      }
      array[0] = V1
      array[1] = Vn
      array[2] = V3
      array[3] = Double(parameters.n)
      array[4] = Kts3
      array[5] = equilibriumLength
      force.addBond(particles: particles, parameters: array)
    }
    super.init(force: force, forceGroup: 1)
  }
}

class MM4TorsionExtendedForce: MM4Force {
  // When fusing multiple carbon-like torsions into a single invocation:
  //
  // Non-carbon torsions only have V4/V6 when there's a carbon in the middle.
  // Therefore, we could do an optimization that computes Vn/V3 as a substitute
  // for V4/V6. Then we don't need to compute "ghost torsions" to fill in the
  // nonexistent 3rd and 4th valences for certain elements.
  //
  // At least, most often. The following exception requires adding 1 ghost atom
  // to a regular carbon torsion:
  // - 1, 1, 8, 1
  
  // WARNING: Divide all Vn torsion parameters by 2.
}
