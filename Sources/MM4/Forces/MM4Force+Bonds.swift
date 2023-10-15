//
//  MM4Force+Bonds.swift
//
//
//  Created by Philip Turner on 10/8/23.
//

import OpenMM

/// Morse bond stretch force.
class MM4StretchForce: MM4Force {
  required init(system: MM4System) {
    // Using "beta" instead of "alpha", as it's the character used in
    // Nanosystems 3.3.3(a).
    //
    // beta = sqrt(ks / 2De)
    let force = OpenMM_CustomBondForce(energy: """
      stretch;
      stretch = potentialWellDepth * (
        1 - exp(-beta * deltaLength)
      )^2 - 1);
      deltaLength = r - equilibriumLength;
      """)
    force.addPerBondParameter(name: "potentialWellDepth")
    force.addPerBondParameter(name: "beta")
    force.addPerBondParameter(name: "equilibriumLength")
    
    let array = OpenMM_DoubleArray(size: 3)
    let bonds = system.parameters.bonds
    for bondID in bonds.indices.indices {
      // Pre-multiply constants in formulas as much as possible. For example,
      // the "beta" constant for bond stretch is multiplied by
      // 'OpenMM_AngstromsPerNm'. This reduces the amount of computation during
      // force execution.
      let bond = bonds.indices[bondID]
      let parameters = bonds.parameters[bondID]
      
      // Units: millidyne-angstrom -> kJ/mol
      var potentialWellDepth = Double(parameters.potentialWellDepth)
      potentialWellDepth *= MM4_KJPerMolPerAJ
      
      // Units: angstrom^-1 -> nm^-1
      var beta = Double(
        parameters.stretchingStiffness / (2 * parameters.stretchingStiffness)
      ).squareRoot()
      beta /= OpenMM_NmPerAngstrom
      
      // Units: angstrom -> nm
      var equilibriumLength = Double(parameters.equilibriumLength)
      equilibriumLength *= OpenMM_NmPerAngstrom
      
      let particles = system.reorder(bond)
      array[0] = potentialWellDepth
      array[1] = beta
      array[2] = equilibriumLength
      force.addBond(particles: particles, parameters: array)
    }
    super.init(force: force, forceGroup: 2)
  }
}
