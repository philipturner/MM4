//
//  MM4Integrator.swift
//  MM4
//
//  Created by Philip Turner on 10/3/23.
//

import OpenMM

/// Stores an OpenMM context generated from an integrator variant.
class MM4Context {
  var compoundIntegrator: OpenMM_CompoundIntegrator
  var context: OpenMM_Context
  var integrators: [MM4IntegratorDescriptor: Int] = [:]
  
  init(system: MM4System, platform: OpenMM_Platform?) {
//    self.compoundIntegrator = OpenMM_VerletIntegrator(stepSize: 0)
    self.compoundIntegrator = OpenMM_CompoundIntegrator()
    
    for start in [false, true] {
      for end in [false, true] {
        var descriptor = MM4IntegratorDescriptor()
        descriptor.start = start
        descriptor.end = end
        
        let integrator = MM4Integrator(descriptor: descriptor)
        integrator.integrator.transfer()
        let index = compoundIntegrator.addIntegrator(integrator.integrator)
        integrators[descriptor] = index
      }
    }
    
    if let platform {
      self.context = OpenMM_Context(
        system: system.system,
        integrator: compoundIntegrator,
        platform: platform)
    } else {
      self.context = OpenMM_Context(
        system: system.system,
        integrator: compoundIntegrator)
    }
  }
  
  var currentIntegrator: MM4IntegratorDescriptor {
    get { fatalError("Not implemented.") }
    set {
      guard let index = integrators[newValue] else {
        fatalError("This should never happen.")
      }
      compoundIntegrator.currentIntegrator = index
    }
  }
  
  /// Modeled after how the OpenMM `integrator.step` API is typically used -
  /// without an argument label for steps.
  func step(_ steps: Int, timeStep: Double) {
    compoundIntegrator.stepSize = timeStep
    compoundIntegrator.step(steps)
  }
}
