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

/// This force only computes the correction to vdW, while the electrostatic
/// exception force computes the correction to partial charges. Either method,
/// using 1,3 or 1,4 exceptions, wouldn't change whether these interactions fall
/// on the diagonal. It also wouldn't change the compute cost due to divergence.
/// The version here may actually decrease compute cost a little, as the
/// exp(-12) term is omitted.
class MM4NonbondedExceptionForce: MM4Force {
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
    
    // It seems like "disfac" was the dispersion factor, similar to the DISP-14
    // keyword in Tinker. Keep the Pauli repulsion force the same though.
    func createForce() -> OpenMM_CustomBondForce {
      let dispersionFactor: Double = 0.550
      let correction = dispersionFactor - 1
      let force = OpenMM_CustomBondForce(energy: """
      epsilon * (
        \(-2.25 * correction) * (min(2, radius / r))^6
      );
      """)
      force.addPerBondParameter(name: "epsilon")
      force.addPerBondParameter(name: "radius")
      return force
    }
    
    // Separate force for (Si, Ge) to (H, C, Si, Ge) interactions.
    func createLegacyForce() -> OpenMM_CustomBondForce {
      let legacyForce = OpenMM_CustomBondForce(energy: """
      legacyEpsilon * (
        -2.25 * (min(2, legacyRadius / r))^6 +
        1.84e5 * exp(-12.00 * (r / legacyRadius))
      ) - epsilon * (
        -2.25 * (min(2, radius / r))^6 +
        1.84e5 * exp(-12.00 * (r / radius))
      );
      """)
      legacyForce.addPerBondParameter(name: "epsilon")
      legacyForce.addPerBondParameter(name: "radius")
      legacyForce.addPerBondParameter(name: "legacyEpsilon")
      legacyForce.addPerBondParameter(name: "legacyRadius")
      return legacyForce
    }
    
    var force: OpenMM_CustomBondForce!
    var legacyForce: OpenMM_CustomBondForce!
    
    // Ideally, one would test how different choices affect the accuracy of
    // molecular structures, similar to the testing of electrostatic exceptions:
    // - dispersionFactor = 0.550
    // - dispersionFactor = 1.000
    // - applying the dispersion factor to the Pauli repulsion force
    let array = OpenMM_DoubleArray(size: 4)
    let exceptions = system.parameters.nonbondedExceptions14
    let atoms = system.parameters.atoms
    for exception in exceptions {
      let parameters1 = atoms.parameters[Int(exception[0])]
      let parameters2 = atoms.parameters[Int(exception[1])]
      
      let epsilon: Float
      let radius: Float
      if parameters1.epsilon.hydrogen * parameters2.epsilon.hydrogen < 0 {
        epsilon = max(parameters1.epsilon.hydrogen,
                      parameters2.epsilon.hydrogen)
        radius = max(parameters1.radius.hydrogen,
                     parameters2.radius.hydrogen)
      } else {
        epsilon = sqrt(parameters1.epsilon.default *
                       parameters2.epsilon.default)
        radius = parameters1.radius.default +
        /**/     parameters2.radius.default
      }
      
      // Units: kcal/mol -> zJ, angstrom -> nm
      array[0] = Double(epsilon) * OpenMM_KJPerKcal * MM4ZJPerKJPerMol
      array[1] = Double(radius) * OpenMM_NmPerAngstrom
      
      let atomicNumber1 = atoms.atomicNumbers[Int(exception[0])]
      let atomicNumber2 = atoms.atomicNumbers[Int(exception[0])]
      let atomicNumbers = SIMD2(atomicNumber1, atomicNumber2)
      let hydrogenMask = atomicNumbers .== 1
      let carbonMask = atomicNumbers .== 6
      let siliconMask = atomicNumbers .== 14
      let germaniumMask = atomicNumbers .== 32
      
      // Use MM3 parameters for torsions involving Si, Ge. The parameters for
      // an exception with an H bonded to C may be slightly off. They use the
      // R=0.94 distance that MM3 was not parameterized with. The expected
      // impact is small, dwarfed by the difference in C-Si and C-Ge exceptions.
      // In addition, correcting the hydrogens would complicate the code.
      var selectedForce: OpenMM_CustomBondForce
      if any(siliconMask .| germaniumMask),
         all(hydrogenMask .| carbonMask .| siliconMask .| germaniumMask) {
        // Change the 8-bit masks to 32-bit for selecting FP32 numbers.
        let hydrogenMask = SIMD2<Int32>(truncatingIfNeeded: atomicNumbers) .== 1
        let carbonMask = SIMD2<Int32>(truncatingIfNeeded: atomicNumbers) .== 6
        
        // Change the H/C epsilons to match MM3.
        var epsilons = SIMD2(parameters1.epsilon.default,
                             parameters2.epsilon.default)
        epsilons.replace(with: 0.020, where: hydrogenMask)
        epsilons.replace(with: 0.027, where: carbonMask)
        
        // Change the H/C radii to match MM3.
        var radii = SIMD2(parameters1.radius.default,
                          parameters2.radius.default)
        radii.replace(with: 1.62, where: hydrogenMask)
        radii.replace(with: 2.04, where: carbonMask)
        
        // Create alternative vdW parameters using the Hill function.
        let legacyEpsilon = sqrt(epsilons[0] * epsilons[1])
        let legacyRadius = radii[0] + radii[1]
        
        // Units: kcal/mol -> zJ, angstrom -> nm
        array[2] = Double(legacyEpsilon) * OpenMM_KJPerKcal * MM4ZJPerKJPerMol
        array[3] = Double(legacyRadius) * OpenMM_NmPerAngstrom
        
        if legacyForce == nil {
          legacyForce = createLegacyForce()
        }
        selectedForce = legacyForce
      } else {
        if force == nil {
          force = createForce()
        }
        selectedForce = force
      }
      
      let particles = system.virtualSiteReorder(exception)
      selectedForce.addBond(particles: particles, parameters: array)
    }
    super.init(forces: [force, legacyForce], forceGroup: 1)
  }
}
