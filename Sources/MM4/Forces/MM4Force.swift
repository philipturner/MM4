//
//  MMForce.swift
//
//
//  Created by Philip Turner on 10/14/23.
//

import OpenMM

class MM4Force {
  /// The OpenMM object containing the force.
  var forces: [OpenMM_Force]
  
  /// How many times to evaluate this force per timestep.
  var forceGroup: Int
  
  init(forces: [OpenMM_Force], forceGroup: Int) {
    for force in forces {
      force.forceGroup = forceGroup
    }
    self.forces = forces
    self.forceGroup = forceGroup
  }
  
  required init(system: MM4System) {
    fatalError("Not implemented.")
  }
}
