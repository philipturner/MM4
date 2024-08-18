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
