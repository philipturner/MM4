//
//  MM4Force+Nonbonded.swift
//
//
//  Created by Philip Turner on 10/8/23.
//

import Foundation
import OpenMM

// The MVP forcefield is minimally optimized. It has a cutoff using 2.5
// sigma for the largest atom supported (germanium). Keeping the cutoff constant
// will minimize variation between different structures during debugging. It
// also masks out non-charged particles from the electrostatic force.
// Electrostatic should be cheaper to compute than vdW, per force evaluation.
//
// For now, not segregating vdW into two different atom groups.
class MM4NonbondedForce: MM4Force {
  required init(system: MM4System) {
    let force = OpenMM_CustomNonbondedForce(energy: """
      epsilon * (
        -2.25 * (length / r)^6 +
        1.84e5 * exp(-12.00 * (r / length))
      );
      epsilon = select(isHydrogenBond, heteroatomEpsilon, hydrogenEpsilon);
      radius = select(isHydrogenBond, heteroatomRadius, hydrogenRadius);
      
      isHydrogenBond = step(hydrogenEpsilon1 * hydrogenEpsilon2);
      heteroatomEpsilon = sqrt(epsilon1 * epsilon2);
      hydrogenEpsilon = max(hydrogenEpsilon1, hydrogenEpsilon2);
      heteroatomRadius = radius1 + radius2;
      hydrogenRadius = max(hydrogenRadius1, hydrogenRadius2);
      """)
    force.addPerParticleParameter(name: "epsilon")
    force.addPerParticleParameter(name: "hydrogenEpsilon")
    force.addPerParticleParameter(name: "radius")
    force.addPerParticleParameter(name: "hydrogenRadius")
    
    let germaniumRadius = 2.440 * OpenMM_NmPerAngstrom
    let cutoff = germaniumRadius * 2.5 * OpenMM_SigmaPerVdwRadius
    force.nonbondedMethod = .noCutoff
    force.useSwitchingFunction = true
    force.cutoffDistance = cutoff
    force.switchingDistance = cutoff * pow(1.0 / 3, 1.0 / 6)
    
    let array = OpenMM_DoubleArray(size: 3)
    let atoms = system.parameters.atoms
    for atomID in atoms.atomicNumbers.indices {
      let parameters = atoms.nonbondedParameters[atomID]
      
      // Units: kcal/mol -> kJ/mol
      let (epsilon, hydrogenEpsilon) = parameters.epsilon
      array[0] = Double(epsilon) * OpenMM_KJPerKcal
      array[1] = Double(hydrogenEpsilon) * OpenMM_KJPerKcal
      
      // Units: angstrom -> nm
      let (radius, hydrogenRadius) = parameters.radius
      array[2] = Double(radius) * OpenMM_NmPerAngstrom
      array[3] = Double(hydrogenRadius) * OpenMM_NmPerAngstrom
      
      force.addParticle(parameters: array)
    }
    
    let bonds = system.parameters.bonds
    let bondPairs = OpenMM_BondArray(size: bonds.indices.count)
    for bondID in bonds.indices.indices {
      let bond = bonds.indices[bondID]
      bondPairs[bondID] = system.reorder(bond)
    }
    force.createExclusionsFromBonds(bondPairs, bondCutoff: 3)
    super.init(forces: [force], forceGroup: 2)
  }
}

class MM4NonbondedExceptionForce: MM4Force {
  required init(system: MM4System) {
    // It seems like "disfac" was the dispersion factor, similar to the DISP-14
    // keyword in Tinker. Keep the Pauli repulsion force the same though.
    let force = OpenMM_CustomBondForce(energy: """
      epsilon * (
        \(-2.25 * 0.550) * (equilibriumLength / r)^6 +
        1.84e5 * exp(-12.00 * (r / equilibriumLength))
      );
      """)
    force.addPerBondParameter(name: "epsilon")
    force.addPerBondParameter(name: "equilibriumLength")
    
    let array = OpenMM_DoubleArray(size: 2)
    let exceptions = system.parameters.nonbondedExceptions14
    let atoms = system.parameters.atoms
    for exception in exceptions {
      let parameters1 = atoms.nonbondedParameters[Int(exception[0])]
      let parameters2 = atoms.nonbondedParameters[Int(exception[1])]
      
      var epsilon: Float
      var equilibriumLength: Float
      if parameters1.epsilon.hydrogen * parameters2.epsilon.hydrogen < 0 {
        epsilon = max(parameters1.epsilon.hydrogen,
                      parameters2.epsilon.hydrogen)
        equilibriumLength = max(parameters1.radius.hydrogen,
                                parameters2.radius.hydrogen)
      } else {
        epsilon = sqrt(parameters1.epsilon.heteroatom *
                       parameters2.epsilon.heteroatom)
        equilibriumLength = parameters1.radius.heteroatom +
        /**/                parameters2.radius.heteroatom
      }
      
      // Units: kcal/mol -> kJ/mol, angstrom -> nm
      array[0] = Double(epsilon) * OpenMM_KJPerKcal
      array[1] = Double(equilibriumLength) * OpenMM_NmPerAngstrom
      
      let particles = system.reorder(exception)
      force.addBond(particles: particles, parameters: array)
    }
    super.init(forces: [force], forceGroup: 1)
  }
}
