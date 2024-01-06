//
//  MM4Force+External.swift
//  MM4
//
//  Created by Philip Turner on 10/14/23.
//

import OpenMM

class MM4ExternalForce: MM4Force {
  required init(system: MM4System, descriptor: MM4ForceFieldDescriptor) {
    // There is no need to convert from kJ/mol to zJ here.
    let force = OpenMM_CustomExternalForce(energy: """
      x * slope_x + y * slope_y + z * slope_z;
      """)
    force.addPerParticleParameter(name: "slope_x")
    force.addPerParticleParameter(name: "slope_y")
    force.addPerParticleParameter(name: "slope_z")
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
      // Force is the negative gradient of potential energy.
      let slope = SIMD3<Double>(-forces[originalID])
      for lane in 0..<3 {
        array[lane] = slope[lane]
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
