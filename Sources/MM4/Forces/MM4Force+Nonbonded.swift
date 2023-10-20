//
//  MM4Force+Nonbonded.swift
//
//
//  Created by Philip Turner on 10/8/23.
//

import Foundation
import OpenMM

// The MVP forcefield isn't maximally optimized. It has a cutoff using 2.5
// sigma for the largest atom supported (germanium), which also happens to be a
// reasonable cutoff for electrostatics with RF. In the future, the vdW force
// will be segregated into two different sets of atom groups. The more expensive
// set contains contiguous tiles with large or charged atom pairs.
class MM4NonbondedForce: MM4Force {
  required init(system: MM4System) {
    let vacuumPermittivity: Double = 8.8541878128e-12
    let coulombConstant = 1 / (4 * Double.pi * vacuumPermittivity)
    var scaleFactor = coulombConstant
    
    let coulombsPerElementaryCharge: Double = 801088317 / 5e27
    scaleFactor *= coulombsPerElementaryCharge * coulombsPerElementaryCharge
    
    let metersPerNm: Double = 1e-9
    scaleFactor /= metersPerNm
    
    // Units: J -> aJ -> kJ/mol
    let AJPerJ: Double = 1e18
    scaleFactor *= AJPerJ * MM4KJPerMolPerAJ
    
    // Minized code from OpenMM reaction field kernel. Note that the r^-3
    // scaling force is corrected to r^-2 in the calling kernel.
    /*
     {
#ifdef USE_CUTOFF
   unsigned int includeInteraction = (!isExcluded && r2 < CUTOFF_SQUARED);
   real tempForce = 0.0f;
#if HAS_COULOMB
 #ifdef USE_CUTOFF
   const real prefactor = ONE_4PI_EPS0*CHARGE1*CHARGE2;
   tempForce += prefactor*(invR - 2.0f*REACTION_FIELD_K*r2);
   tempEnergy += includeInteraction ? prefactor*(invR + REACTION_FIELD_K*r2 - REACTION_FIELD_C) : 0;
 #else
   const real prefactor = ONE_4PI_EPS0*CHARGE1*CHARGE2*invR;
   tempForce += prefactor;
   tempEnergy += includeInteraction ? prefactor : 0;
 #endif
#endif
   dEdR += includeInteraction ? tempForce*invR*invR : 0;
}
     */
    let force = OpenMM_CustomNonbondedForce(energy: """
      vanDerWaals + electrostatic;
      vanDerWaals = epsilon * (
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
    
    let array = OpenMM_DoubleArray(size: 4)
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

// ========================================================================== //
//
// Create a counterforce at short range that corrects for 1/2 of a
// dipole falling on the 1-3 border. In the paper about fluorine, Allinger
// mentioned some 1-4 fluorines on fluoroethane having a repulsive effect
// from their dipoles. This suggests 1-4 interactions ("sclfac = 1.000"?) are
// included, but perhaps not the 1-3 interactions.
//
// Create a different force to handle the dipole-dipole interactions along
// every possible path between two bonds, which wasn't fully handled when
// recording the 1-3/1-4 nonbonded exceptions. In other words, enumerate every
// torsion out there and create a nonbonded exception force for it.
//
// Find all the bonds on atom 1 that are present in a torsion. Repeat the same
// process for atom 4. Then, compute an O(n^2) algorithm on all of the dipole
// permutations ('n' in number of dipoles, often n = 2-3). A final force
// computes a charge-charge interaction with the partial charge generated by
// the involved bonds omitted.
//
// ========================================================================== //
//
// More efficient alternative: instead of projecting onto dipoles, factor in
// only the charges on 1-3 atoms that result from the dipoles present in 1-4
// atoms. This has an effect similar to scaling those atoms' charge-charge
// interactions by 0.5 in conventional force fields. It reproduces the same
// dynamics to a first-order approximation, but avoids something potentially as
// costly as torsion forces.
//
// Also, from the paper on induced dipoles:
//
// > It was later found that distances among interacting
// > centers, d, need to be empirically adjusted to give
// > better results in the induced dipole moments cal-
// > culations. They are adjusted as follows: (a) for all
// > nongeminal cases, 1,4 interactions and higher, all
// > distances are increased by a factor of 1.3,
// >     d = 1.3 * d
// > The need of using this scaling factor can be ratio-
// > nalized by the fact that we use midpoints of bonds,
// > instead of atoms, as interacting centers. Two ap-
// > proaches are basically the same when two dipoles
// > are far away from each other, however, the scaling
// > factor 1.3 is needed when two dipoles are very close
// > to each other. This scaling factor indicates that the
// > real source of polarization in the molecules might
// > better be considered to be at the atoms instead of
// > at the bonds. (b) For geminal cases, such as the 1,3
// > interaction in the moiety of x-Z-y, if one and ONLY
// > one of the two atoms attached to Z (x or y) is hy-
// > hdrogen, then the distances are increased by a factor
// > of 3.0,
// >     d = 3.0 ∗ d
// > otherwise the distances are increased by a factor
// > of 2.0,
// >     d = 2.0 ∗ d
// > The reason for these rather large scaling factors for
// > geminal cases is that the idealized dipole interaction
// > tensor is no longer valid when two dipoles are lo-
// > cated so close to one another.
//
// Summary: when two dipoles are part of the same angle, the dipole-dipole
// interaction model no longer applies **with respect to induced polarization**.
// At the distance of 1,4, there is still distortion from modeling atoms as two
// dipoles. This suggests that a partial charge model, where the partial charge
// centers are close to the nuclei, is qualitatively better than the
// dipole-dipole interaction formula **with respect to induced polarization**.
//
// The electrostatic exception force is not involving any costly induced
// dipoles. But, the insight may justify using the partial charges not only for
// performance, but also for correctness. This would be a bit different than
// existing forcefields that omit 1,3 and scale 1,4 by a half. It is like
// scaling 1,3 by a half and 1,4 by the full amount. However, both MM3 and MM4
// seem to be parameterized with the dipole-dipole approach with zero
// attenuation for 1,4 interactions.
//
// Update: dipole-dipole forces will be used initially, rather than
//         charge-charge forces
//
// ========================================================================== //

// There is a separate force for vdW/charge and another force accepting every
// permutation of bonds present in the 1,4 atoms' torsions.
class MM4NonbondedExceptionForce: MM4Force {
  required init(system: MM4System) {
    // Equation for dipole-dipole interaction:
    // https://janheyda.files.wordpress.com/2015/08/electrostatics-multipoles.pdf
    //
    // debye = 0.02081943 electron-nm
    // charge-charge = q1 q2 / 4 pi epsilon_0
    // dipole-dipole = mu1 mu2 / 4 pi * (
    //   normalize(p1 - p2) * normalize(p4 - p3) + ...
    //
    
    // It seems like "disfac" was the dispersion factor, similar to the DISP-14
    // keyword in Tinker. Keep the Pauli repulsion force the same though.
    let force = OpenMM_CustomBondForce(energy: """
      vanDerWaals + electrostatic;
      vanDerWaals = epsilon * (
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
