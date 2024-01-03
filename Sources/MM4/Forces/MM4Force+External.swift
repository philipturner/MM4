//
//  MM4Force+External.swift
//
//
//  Created by Philip Turner on 10/14/23.
//

import OpenMM

class MM4ExternalForce: MM4ForceGroup {
  required init(system: MM4System) {
    let force = OpenMM_CustomExternalForce(energy: """
      x * fx + y * fy + z * fz;
      """)
    force.addPerParticleParameter(name: "fx")
    force.addPerParticleParameter(name: "fy")
    force.addPerParticleParameter(name: "fz")
    var forceActive = false
    
    let array = OpenMM_DoubleArray(size: 3)
    array[0] = 0
    array[1] = 0
    array[2] = 0
    
    for reorderedID in system.reorderedIndices {
      force.addParticle(Int(reorderedID), parameters: array)
      forceActive = true
    }
    super.init(forces: [force], forcesActive: [forceActive], forceGroup: 1)
  }
  
  /// Do not reorder the forces before entering into this function.
  func updateForces(_ forces: [SIMD3<Float>], system: MM4System) {
    let forceObject = self.forces[0] as! OpenMM_CustomExternalForce
    let array = OpenMM_DoubleArray(size: 3)
    
    for (originalID, reorderedID) in system.reorderedIndices.enumerated() {
      // Units: zJ/nm -> kJ/mol/nm
      var force = SIMD3<Double>(forces[originalID])
      force *= MM4KJPerMolPerZJ
      
      // Force is the negative gradient of potential energy.
      force = -force
      
      for lane in 0..<3 {
        array[lane] = force[lane]
      }
      forceObject.setParticleParameters(
        index: originalID, particle: Int(reorderedID), parameters: array)
    }
  }
  
  /// This must be called every time the forces change.
  func updateParametersInContext(_ context: MM4Context) {
    let forceObject = self.forces[0] as! OpenMM_CustomExternalForce
    forceObject.updateParametersInContext(context.context)
  }
}
