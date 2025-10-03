//
//  MM4Integrator.swift
//  MM4
//
//  Created by Philip Turner on 10/3/23.
//

import OpenMM

struct MM4ContextDescriptor {
  var integratorOptions: MM4IntegratorOptions?
  var platform: OpenMM_Platform?
  var system: MM4System?
}

/// Encapsulates an OpenMM context and the various integrators.
class MM4Context {
  var compoundIntegrator: OpenMM_CompoundIntegrator?
  var context: OpenMM_Context
  var customIntegrators: [MM4CustomIntegratorDescriptor: Int] = [:]
  var integrator: OpenMM_Integrator
  var verletIntegrator: OpenMM_VerletIntegrator?
  
  init(descriptor: MM4ContextDescriptor) {
    guard let integratorOptions = descriptor.integratorOptions,
          let platform = descriptor.platform,
          let system = descriptor.system else {
      fatalError("Descriptor was incomplete.")
    }
    
    switch integratorOptions {
    case .multipleTimeStep:
      self.compoundIntegrator = OpenMM_CompoundIntegrator()
      self.integrator = compoundIntegrator!
      
      for start in [false, true] {
        for end in [false, true] {
          var descriptor = MM4CustomIntegratorDescriptor()
          descriptor.start = start
          descriptor.end = end
          
          let integrator = MM4CustomIntegrator(descriptor: descriptor)
          integrator.integrator.transfer()
          let index = compoundIntegrator!.addIntegrator(integrator.integrator)
          customIntegrators[descriptor] = index
        }
      }
    case .verlet:
      self.verletIntegrator = OpenMM_VerletIntegrator(stepSize: 0)
      self.integrator = verletIntegrator!
    }
    
    self.context = OpenMM_Context(
      system: system.system,
      integrator: integrator,
      platform: platform)
  }
  
  var currentIntegrator: MM4CustomIntegratorDescriptor {
    get { fatalError("Not implemented.") }
    set {
      guard let index = customIntegrators[newValue],
            let compoundIntegrator else {
        fatalError("This should never happen.")
      }
      compoundIntegrator.currentIntegrator = index
    }
  }
  
  /// Modeled after how the OpenMM `integrator.step` API is typically used -
  /// without an argument label for steps.
  func step(_ steps: Int, timeStep: Double) {
    integrator.stepSize = timeStep
    integrator.step(steps)
  }
}
