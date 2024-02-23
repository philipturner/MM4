//
//  MM4Integrator.swift
//  MM4
//
//  Created by Philip Turner on 10/3/23.
//

import OpenMM

/// Encapsulates an OpenMM context and the various integrators.
class MM4Context {
  var compoundIntegrator: OpenMM_CompoundIntegrator?
  var context: OpenMM_Context
  var customIntegrators: [MM4CustomIntegratorDescriptor: Int] = [:]
  var integrator: OpenMM_Integrator
  var verletIntegrator: OpenMM_VerletIntegrator?
  
  init(system: MM4System, descriptor: MM4ForceFieldDescriptor) {
    switch descriptor.integrator {
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
    
    if let platform = descriptor.platform {
      self.context = OpenMM_Context(
        system: system.system,
        integrator: integrator,
        platform: platform)
    } else {
      self.context = OpenMM_Context(
        system: system.system,
        integrator: integrator)
    }
    
    print()
    print("serialized system")
    print("=================")
    print(OpenMM_XmlSerializer.serializeSystem(system.system))
    print("=================")
    
    print()
    print("serialized integrator")
    print("=====================")
    print(OpenMM_XmlSerializer.serializeIntegrator(integrator))
    print("=====================")
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
