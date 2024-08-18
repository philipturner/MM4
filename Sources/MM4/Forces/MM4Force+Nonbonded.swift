//
//  MM4Force+Nonbonded.swift
//  MM4
//
//  Created by Philip Turner on 10/8/23.
//

import Foundation
import OpenMM

class MM4NonbondedForce: MM4Force {
  required init(system: MM4System, descriptor: MM4ForceFieldDescriptor) {
    var includeNonbonded = false
    for params in system.parameters.atoms.parameters {
      if params.epsilon.default != 0 || params.epsilon.hydrogen != 0 {
        includeNonbonded = true
      }
    }
    guard includeNonbonded else {
      super.init(forces: [], forceGroup: 1)
      return
    }
    
    // WARNING
    //
    // The hydrogens needs to be shifted toward C/Si/Ge by a factor of 0.94.
    // Run the forcefield through initially without that modification. Only
    // after the code is thoroughly debugged, add the virtual particles and
    // reorder them. The new vdW force will be projected onto both the hydrogen
    // and the non-hydrogen. How does this affect partial charges? Is the bond
    // dipole modified to accomodate the different length?
    //
    // Luckily, there are no polarized atom-hydrogen bonds in this MM4
    // implementation. No dipole-dipole interactions or projected charge-charge
    // interactions can involve a virtual sites. Hydrogens have 0 partial
    // charge. Computing the coulomb interaction on their virtual sites creates
    // zero energy, removing the need to account for the position being
    // different.
    
    // An issue arises with this functional form: the extremely
    // malformed case where hydrogens approach closer than the vdW radius. This
    // may happen when a (100) crystal surface has hydrogen collisions. The
    // hydrogens often fall into an attractive region of the potential,
    // whose minimum is r = 0.
    //
    // To overcome the issue, the functional form becomes piecewise in this
    // extreme region. The region never appears in any realistic minimized
    // structure. This graph (https://www.desmos.com/calculator/5eypewn7fl)
    // shows the following potentials (in order):
    // - vdW H-H
    // - close-range section of clamped vdW H-H
    //   - freezes the attractive term at the value with r = 0.5 * r_vdW
    // - Morse C-H
    //
    // The clamping doesn't occur until the hydrogens are in an extremely steep
    // region of the vdW potential. The graph is in zeptojoules. While the
    // minimum of the vdW well is ~0.1 zJ, the crossover with the clamped region
    // occurs on the scale of attojoules. This can be seen by changing the 1000x
    // multiplier for every equation to 1x.
    
    let force = OpenMM_CustomNonbondedForce(energy: """
      epsilon * (
        -2.25 * (min(2, radius / safe_r))^6 +
        1.84e5 * exp(-12.00 * (safe_r / radius))
      );
      safe_r = max(0.001, r);
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
    
    if let cutoffDistance = descriptor.cutoffDistance {
      force.nonbondedMethod = .cutoffNonPeriodic
      force.cutoffDistance = Double(cutoffDistance)
      
      force.useSwitchingFunction = true
      force.switchingDistance = Double(cutoffDistance * pow(1.0 / 3, 1.0 / 6))
    } else {
      force.nonbondedMethod = .noCutoff
    }
    
    
    let array = OpenMM_DoubleArray(size: 4)
    let atoms = system.parameters.atoms
    for atomID in atoms.indices {
      let parameters = atoms.parameters[Int(atomID)]
      
      // Units: kcal/mol -> kJ/mol
      //          kJ/mol -> zJ
      let (epsilon, hydrogenEpsilon) = parameters.epsilon
      array[0] = Double(epsilon) * OpenMM_KJPerKcal * MM4ZJPerKJPerMol
      array[1] = Double(hydrogenEpsilon) * OpenMM_KJPerKcal * MM4ZJPerKJPerMol
      
      // Units: angstrom -> nm
      let (radius, hydrogenRadius) = parameters.radius
      array[2] = Double(radius) * OpenMM_NmPerAngstrom
      array[3] = Double(hydrogenRadius) * OpenMM_NmPerAngstrom
      force.addParticle(parameters: array)
      
      // Give the original hydrogens zero vdW energy.
      if atoms.atomicNumbers[atomID] == 1 {
        array[0] = 0
        array[1] = 0
        force.addParticle(parameters: array)
      }
    }
    
    system.createExceptions(force: force)
    super.init(forces: [force], forceGroup: 1)
  }
}
