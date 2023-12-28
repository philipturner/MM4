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
      potentialWellDepth * ((
        1 - exp(-beta * (r - equilibriumLength))
      )^2 - 1);
      """)
    force.addPerBondParameter(name: "potentialWellDepth")
    force.addPerBondParameter(name: "beta")
    force.addPerBondParameter(name: "equilibriumLength")
    var forceActive = false
    
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
      potentialWellDepth *= MM4KJPerMolPerAJ
      
      // Units: angstrom^-1 -> nm^-1
      var beta = Double(
        parameters.stretchingStiffness / (2 * parameters.potentialWellDepth)
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
      forceActive = true
    }
    super.init(forces: [force], forcesActive: [forceActive], forceGroup: 2)
  }
}
