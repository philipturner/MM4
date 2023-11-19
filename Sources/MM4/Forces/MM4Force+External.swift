//
//  MM4Force+External.swift
//
//
//  Created by Philip Turner on 10/14/23.
//

import OpenMM

class MM4ExternalForce: MM4Force {
  /// Whether each particle participates in the extended force. This masks out
  /// virtual sites.
  var inclusionMask: [Bool] = []
  
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
    
    inclusionMask.reserveCapacity(system.reorderedIndices.count)
    for reorderedID in system.reorderedIndices.indices {
      if reorderedID >= system.virtualSiteCount {
        let originalID = Int(system.originalIndices[reorderedID])
        let atomicNumber = system.parameters.atoms.atomicNumbers[originalID]
        if atomicNumber == 1 {
          inclusionMask.append(false)
          continue
        }
      }
      
      inclusionMask.append(true)
      force.addParticle(reorderedID, parameters: array)
    }
    super.init(forces: [force], forceGroup: 1)
  }
  
  /// Reorder the forces before entering into this function. Forces
  /// corresponding to virtual sites are ignored.
  func updateForces(_ forces: [SIMD3<Float>], context: MM4Context) {
    let forceObject = self.forces[0]
    guard let forceObject = forceObject as? OpenMM_CustomExternalForce else {
      fatalError("External force was not an external force.")
    }
    let array = OpenMM_DoubleArray(size: 3)
    
    for reorderedID in forces.indices {
      guard inclusionMask[reorderedID] else {
        continue
      }
      
      // Units: zJ/nm -> kJ/mol/nm
      var force = SIMD3<Double>(forces[reorderedID])
      force *= MM4KJPerMolPerZJ
      
      // Force is the negative gradient of potential energy.
      force = -force
      
      for lane in 0..<3 {
        array[lane] = force[lane]
      }
      forceObject.setParticleParameters(
        index: reorderedID, particle: reorderedID, parameters: array)
    }
    forceObject.updateParametersInContext(context.context)
  }
}
