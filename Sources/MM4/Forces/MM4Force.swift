//
//  MMForce.swift
//
//
//  Created by Philip Turner on 10/14/23.
//

import OpenMM

class MM4Force {
  /// The OpenMM object containing the force.
  var force: OpenMM_Force
  
  /// How many times to evaluate this force per timestep.
  var forceGroup: Int
  
  init(force: OpenMM_Force, forceGroup: Int) {
    self.force = force
    self.forceGroup = forceGroup
  }
  
  required init(system: MM4System) {
    fatalError("Not implemented.")
  }
}
