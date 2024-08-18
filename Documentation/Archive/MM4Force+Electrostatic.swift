//
//  MM4+Electrostatic.swift
//  MM4
//
//  Created by Philip Turner on 10/19/23.
//

import Foundation
import OpenMM

/// To isolate the electrostatic force from the switching function, we create
/// two forces that end up fused into the same group.
class MM4ElectrostaticForce: MM4Force {
  static var prefactor: Double {
    let vacuumPermittivity: Double = 8.8541878128e-12
    let coulombConstant = 1 / (4 * Double.pi * vacuumPermittivity)
    var prefactor = coulombConstant
    
    let coulombsPerElementaryCharge: Double = 801088317 / 5e27
    prefactor *= coulombsPerElementaryCharge * coulombsPerElementaryCharge
    
    let metersPerNm: Double = 1e-9
    prefactor /= metersPerNm
    
    // Units: J -> aJ
    //       aJ -> kJ/mol
    //   kJ/mol -> zJ
    let AJPerJ: Double = 1e18
    prefactor *= AJPerJ
    prefactor *= MM4KJPerMolPerAJ
    prefactor *= MM4ZJPerKJPerMol
    return prefactor
  }
  
  static func reactionFieldConstants(
    cutoffDistance: Float, dielectricConstant: Float
  ) -> (K: Float, C: Float) {
    // Testing whether there's a sharp cutoff as stated on Wikipedia:
    // https://www.desmos.com/calculator/zgmywq7xui
    //
    // No; the function has an exact minimum at the cutoff distance. Force acts
    // just like if energy were conserved. OpenMM adds a correction term (C) to
    // make the potential energy be zero at the cutoff.
    
    // Source (from OpenMM C++ code):
    //
    // double reactionFieldK = pow(force.getCutoffDistance(), -3.0)
    // * (force.getReactionFieldDielectric()-1.0)
    // / (2.0*force.getReactionFieldDielectric()+1.0);
    //
    // double reactionFieldC = (1.0 / force.getCutoffDistance())
    // * (3.0*force.getReactionFieldDielectric())
    // / (2.0*force.getReactionFieldDielectric()+1.0);
    let reactionFieldK = pow(cutoffDistance, -3) * (dielectricConstant - 1)
    /**/               / (2 * dielectricConstant + 1)
    let reactionFieldC = (1 / cutoffDistance) * (3 * dielectricConstant)
    /**/               / (2 * dielectricConstant + 1)
    return (reactionFieldK, reactionFieldC)
  }
  
  required init(system: MM4System, descriptor: MM4ForceFieldDescriptor) {
    var includeElectrostatic = false
    for params in system.parameters.atoms.parameters {
      if params.charge != 0 {
        includeElectrostatic = true
      }
    }
    guard includeElectrostatic else {
      super.init(forces: [], forceGroup: 1)
      return
    }
    
    // Minimized code from OpenMM reaction field kernel. Note that the r^-3
    // scaling force term becomes r^-2 in the calling kernel.
    /*
    {
    #ifdef USE_CUTOFF
      unsigned int includeInteraction = (!isExcluded && r2 < CUTOFF_SQUARED);
      real tempForce = 0.0f;
      const real prefactor = ONE_4PI_EPS0*CHARGE1*CHARGE2;
      tempForce += prefactor*(invR - 2.0f*REACTION_FIELD_K*r2);
      tempEnergy += includeInteraction ? prefactor*(invR + REACTION_FIELD_K*r2 - REACTION_FIELD_C) : 0;
      dEdR += includeInteraction ? tempForce*invR*invR : 0;
    }
     */
    
    let prefactor = MM4ElectrostaticForce.prefactor
    var K, C: Float
    if let cutoffDistance = descriptor.cutoffDistance {
      (K, C) = MM4ElectrostaticForce.reactionFieldConstants(
        cutoffDistance: cutoffDistance,
        dielectricConstant: descriptor.dielectricConstant)
    } else {
      (K, C) = (.zero, .zero)
    }
    
    let force = OpenMM_CustomNonbondedForce(energy: """
      \(prefactor) * charge1 * charge2 * (
        1 / r + \(K) * r^2 - \(C)
      );
      """)
    force.addPerParticleParameter(name: "charge")
    
    if let cutoffDistance = descriptor.cutoffDistance {
      // Force is discontinuous at the cutoff
      // (https://www.desmos.com/calculator/a9my66qauu). However, adding a
      // switching function like with the vdW force makes it worse. It
      // introduces a jump in force as the energy drops faster than before.
      // There's not much wrong with discontinuous forces, but there's
      // definitely something wrong with discontinuous energies.
      force.nonbondedMethod = .cutoffNonPeriodic
      force.cutoffDistance = Double(cutoffDistance)
    } else {
      force.nonbondedMethod = .noCutoff
    }
    
    let array = OpenMM_DoubleArray(size: 1)
    let atoms = system.parameters.atoms
    for atomID in atoms.indices {
      let parameters = atoms.parameters[Int(atomID)]
      
      // Units: elementary charge
      array[0] = Double(parameters.charge)
      force.addParticle(parameters: array)
      
      // Give the original hydrogens zero charge.
      if atoms.atomicNumbers[atomID] == 1 {
        array[0] = 0
        force.addParticle(parameters: array)
      }
    }
    
    system.createExceptions(force: force)
    super.init(forces: [force], forceGroup: 1)
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
// > drogen, then the distances are increased by a factor
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

/// In the future, an optimization could fuse the O(bonds^2) interactions into a
/// single kernel invocation. This would also pre-compute the O(bonds^2) partial
/// charge interactions to "undo" the regular electrostatic force.
class MM4ElectrostaticExceptionForce: MM4Force {
  required init(system: MM4System, descriptor: MM4ForceFieldDescriptor) {
    var includeElectrostaticException = false
    for params in system.parameters.atoms.parameters {
      if params.charge != 0 {
        includeElectrostaticException = true
      }
    }
    if system.parameters.nonbondedExceptions14.count == 0 {
      includeElectrostaticException = false
    }
    guard includeElectrostaticException else {
      super.init(forces: [], forceGroup: 1)
      return
    }
    
    // Equation for dipole-dipole interaction:
    // https://janheyda.files.wordpress.com/2015/08/electrostatics-multipoles.pdf
    //
    // debye = 0.02081943 electron-nm
    // charge-charge = q1 q2 / 4 pi epsilon_0
    // dipole-dipole = mu1 mu2 / 4 pi epsilon_0 * (
    //   normalize(p1 - p2) * normalize(p4 - p3) + ...
    
    // For each bond-bond interaction, compute the projected charge onto each
    // 1,4 atom and undo it.
    let prefactor = MM4ElectrostaticForce.prefactor
    var K, C: Float
    if let cutoffDistance = descriptor.cutoffDistance {
      (K, C) = MM4ElectrostaticForce.reactionFieldConstants(
        cutoffDistance: cutoffDistance,
        dielectricConstant: descriptor.dielectricConstant)
    } else {
      (K, C) = (.zero, .zero)
    }
    
    // It is currently unknown whether MM4 includes the 1-2,3-4 dipole-dipole
    // interaction. This must be resolved through testing.
    #if true
    let force = OpenMM_CustomCompoundBondForce(numParticles: 4, energy: """
      \(prefactor) * -chargeCharge;
      chargeCharge = chargeChargeProduct * (
        1 / r14 + \(K) * r14^2 - \(C)
      );
      r14 = distance(p1, p4);
      """)
    #else
    let force = OpenMM_CustomCompoundBondForce(numParticles: 4, energy: """
      \(prefactor) * (dipoleDipole - chargeCharge);
      
      dipoleDipole = dipoleDipoleProduct * invLength_muij^3 * lengthScale * (
        (xi * xj + yi * yj + zi * zj) - 3 * invLength_muij^2 * (
          xi * x_torsion + yi * y_torsion + zi * z_torsion
        ) * (
          x_torsion * xj + y_torsion * yj + z_torsion * zj
        )
      );
      x_torsion = xi - xj;
      y_torsion = yi - yj;
      z_torsion = zi - zj;
      lengthScale = invLength_mui * invLength_muj;
      invLength_mui = 1 / sqrt(xi^2 + yi^2 + zi^2);
      invLength_muj = 1 / sqrt(xj^2 + yj^2 + zj^2);
      invLength_muij = 1 / sqrt(xij^2 + yij^2 + zij^2);
      
      xij = x12 - x34;
      yij = y12 - y34;
      zij = z12 - z34;
      xi = x1 - x12;
      yi = y1 - y12;
      zi = z1 - z12;
      xj = x4 - x34;
      yj = y4 - y34;
      zj = z4 - z34;
      x12 = (x1 + x2) / 2;
      y12 = (y1 + y2) / 2;
      z12 = (z1 + z2) / 2;
      x34 = (x3 + x4) / 2;
      y34 = (y3 + y4) / 2;
      z34 = (z3 + z4) / 2;
      
      chargeCharge = chargeChargeProduct * (
        1 / r14 + \(K) * r14^2 - \(C)
      );
      r14 = distance(p1, p4);
      """)
    #endif
    force.addPerBondParameter(name: "dipoleDipoleProduct")
    force.addPerBondParameter(name: "chargeChargeProduct")
    
    // Collect all the torsions that exist between a 1,4 pair.
    var pairsToTorsionsMap: [SIMD2<UInt32>: SIMD8<Int32>] = [:]
    let bonds = system.parameters.bonds
    let torsions = system.parameters.torsions
    for torsionID in torsions.indices.indices {
      let torsion = torsions.indices[torsionID]
      var pair = SIMD2(torsion[0], torsion[3])
      pair = system.parameters.sortBond(pair)
      
      guard var map = pairsToTorsionsMap[pair] else {
        pairsToTorsionsMap[pair] = SIMD8(
          Int32(torsionID), -1, -1, -1, -1, -1, -1, 1)
        continue
      }
      guard map[8] < 7 else {
        fatalError("Pair to torsion map was full.")
      }
      map[Int(map[8])] = Int32(torsionID)
      map[8] += 1
      pairsToTorsionsMap[pair] = map
    }
    
    // Create force instances for each exception.
    let particles = OpenMM_IntArray(size: 4)
    let array = OpenMM_DoubleArray(size: 2)
    var arrayLeft: [SIMD2<UInt32>] = []
    var arrayRight: [SIMD2<UInt32>] = []
    for exception in system.parameters.nonbondedExceptions14 {
      guard let map = pairsToTorsionsMap[exception] else {
        fatalError("No torsions found for 1,4 exception.")
      }
      arrayLeft.removeAll(keepingCapacity: true)
      arrayRight.removeAll(keepingCapacity: true)
      
      for lane in 0..<Int(map[8]) {
        // 'bondLeft' and 'bondRight' are not necessarily sorted in numerical
        // order, just in a convenient form where the inner atom goes first.
        //
        // WARNING: Before entering particles into the OpenMM kernel, swap the
        // two indices in 'bondLeft'.
        var torsion = torsions.indices[Int(map[lane])]
        if torsion[1] == exception[1] || torsion[2] == exception[0] {
          torsion = SIMD4(torsion[3], torsion[2], torsion[1], torsion[0])
        }
        let bondLeft = SIMD2(torsion[1], torsion[0])
        let bondRight = SIMD2(torsion[2], torsion[3])
        
        func append(_ bond: SIMD2<UInt32>, array: inout [SIMD2<UInt32>]) {
          for i in 0..<array.count {
            if all(array[i] .== bond) {
              break
            }
            if i + 1 == array.count {
              array.append(bond)
              break
            }
          }
        }
        append(bondLeft, array: &arrayLeft)
        append(bondRight, array: &arrayRight)
      }
      
      /// - Returns: Dipole in e-nm, partial charges in e.
      func project(_ bond: SIMD2<UInt32>) -> (
        dipoleMoment: Float, charge: SIMD2<Float>
      )? {
        let sorted = system.parameters.sortBond(bond)
        guard let bondID = bonds.map[sorted],
              let parameters = bonds.extendedParameters[Int(bondID)] else {
          return nil
        }
        var partialCharges = system.parameters.projectDipole(
          parameters.dipoleMoment, bondID: bondID)
        if any(bond .!= sorted) {
          partialCharges = SIMD2(partialCharges[1], partialCharges[0])
        }
        
        // Units: elementary charge * angstrom -> elementary charge * nm
        var dipoleMoment = parameters.dipoleMoment * Float(MM4EAngstromPerDebye)
        dipoleMoment *= Float(OpenMM_NmPerAngstrom)
        return (dipoleMoment, partialCharges)
      }
      
      for bondLeft in arrayLeft {
        guard let (dipoleLeft, chargesLeft) = project(bondLeft) else {
          continue
        }
        for bondRight in arrayRight {
          guard let (dipoleRight, chargesRight) = project(bondRight) else {
            continue
          }
          
          let dipoleDipoleProduct = dipoleLeft * dipoleRight
          let chargeChargeProduct = chargesLeft[1] * chargesRight[1]
          guard dipoleDipoleProduct > 0 || chargeChargeProduct > 0 else {
            continue
          }
          array[0] = Double(dipoleDipoleProduct)
          array[1] = Double(chargeChargeProduct)
          
          // WARNING: Before entering particles into the OpenMM kernel, swap the
          // two indices in 'bondLeft'.
          var original: SIMD4<UInt32> = .zero
          original[0] = bondLeft[1]
          original[1] = bondLeft[0]
          original[2] = bondRight[0]
          original[3] = bondRight[1]
          
          let reordered = system.virtualSiteReorder(original)
          for lane in 0..<4 {
            particles[lane] = reordered[lane]
          }
          force.addBond(particles: particles, parameters: array)
        }
      }
    }
    super.init(forces: [force], forceGroup: 1)
  }
}
