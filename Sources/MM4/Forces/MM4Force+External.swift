//
//  MM4Force+External.swift
//
//
//  Created by Philip Turner on 10/14/23.
//

import OpenMM

class MM4ExternalForce: MM4Force {
  required init(system: MM4System) {
    let force = OpenMM_CustomExternalForce(energy: """
      x * fx + y * fy + z * fz;
      """)
    force.addPerParticleParameter(name: "fx")
    force.addPerParticleParameter(name: "fy")
    force.addPerParticleParameter(name: "fz")
    
    let array = OpenMM_DoubleArray(size: 3)
    array[0] = 0
    array[1] = 0
    array[2] = 0
    
    let atoms = system.parameters.atoms
    for atomID in atoms.indices {
      force.addParticle(atomID, parameters: array)
    }
    
    super.init(forces: [force], forceGroup: 1)
  }
  
  /// > WARNING: Reorder the forces before entering into this object.
  func updateForces(_ forces: [SIMD3<Float>], context: MM4Context) {
    guard let forceObject = self.forces[0] as? OpenMM_CustomExternalForce else {
      fatalError("External force was not an external force.")
    }
    let array = OpenMM_DoubleArray(size: 3)
    
    for atomID in forces.indices {
      // Units: zJ/nm -> kJ/mol/nm
      var force = SIMD3<Double>(forces[atomID])
      force *= MM4KJPerMolPerZJ
      
      // Force is the negative gradient of potential energy.
      force = -force
      
      for lane in 0..<3 {
        array[lane] = force[lane]
      }
      forceObject.setParticleParameters(
        index: atomID, particle: atomID, parameters: array)
    }
    forceObject.updateParametersInContext(context.context)
    
    // TODO: Create a single context/integrator accepting a global parameter for
    // number of loop iterations. Add while blocks to the OpenMM bindings.
    fatalError("Rewrite the Integrator implementation from scratch.")
  }
}
