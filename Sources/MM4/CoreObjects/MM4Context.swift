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
    self.context = OpenMM_Context(
      system: system.system, integrator: integrator.integrator)
  }
}

extension MM4ForceField {
//  func context(descriptor: MM4IntegratorDescriptor) -> MM4Context {
//    if let context = contexts[descriptor] {
//      return context
//    }
//    
//    let integrator = MM4Integrator(descriptor: descriptor)
//    let context = MM4Context(system: system, integrator: integrator)
//    contexts[descriptor] = context
//    return context
//  }
  
//  func switchContext(_ context: MM4Context) {
//    if latestContext === context {
//      // There is no need to switch contexts.
//      return
//    }
//    
//    let state = latestContext.context.state(types: [.positions, .velocities])
//    context.context.state = state
//    latestContext = context
//    
//    // Update any forces whose parameters change based on user input.
//    system.forces.external.updateForces(context: context)
//  }
}

