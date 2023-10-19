//
//  MM4Integrator.swift
// 
//
//  Created by Philip Turner on 10/3/23.
//

import OpenMM

/// Stores an OpenMM context generated from an integrator variant.
class MM4Context {
  var context: OpenMM_Context
  
  init(system: MM4System, integrator: MM4Integrator) {
    fatalError("Not implemented.")
  }
}

extension MM4ForceField {
  func context(descriptor: MM4IntegratorDescriptor) -> MM4Context {
    if let context = contexts[descriptor] {
      return context
    }
    
    let integrator = MM4Integrator(descriptor: descriptor)
    let context = MM4Context(system: system, integrator: integrator)
    contexts[descriptor] = context
    return context
  }
}
