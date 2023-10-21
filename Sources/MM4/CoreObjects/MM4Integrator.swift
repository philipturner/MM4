//
//  MM4Integrator.swift
//  
//
//  Created by Philip Turner on 10/14/23.
//

import OpenMM

/// A configuration for an integrator.
class MM4IntegratorDescriptor: Hashable {
  /// The number of force evaluations fused together to reduce redundant
  /// computation of forces. Setting this to one creates a Verlet integrator.
  var fusedTimeSteps: Int = 1
  
  init() {
    
  }
  
  static func == (
    lhs: MM4IntegratorDescriptor,
    rhs: MM4IntegratorDescriptor
  ) -> Bool {
    guard lhs.fusedTimeSteps == rhs.fusedTimeSteps else {
      return false
    }
    return true
  }
  
  func hash(into hasher: inout Hasher) {
    hasher.combine(fusedTimeSteps)
  }
}

class MM4Integrator {
  var integrator: OpenMM_CustomIntegrator
  
  /// Create an integrator using the specified configuration.
  init(descriptor: MM4IntegratorDescriptor) {
    self.integrator = OpenMM_CustomIntegrator(stepSize: 1 * OpenMM_PsPerFs)
    
    integrator.addPerDofVariable(name: "force0", initialValue: 0)
    integrator.addComputePerDof(variable: "force0", expression: "f0")
    for i in 0..<descriptor.fusedTimeSteps {
      if i == 0 {
        integrator.addComputePerDof(variable: "v", expression: """
          v + 0.5 * (dt / \(descriptor.fusedTimeSteps)) * (force0 + f1) / m
          """)
        integrator.addComputePerDof(variable: "v", expression: """
          v + 0.25 * (dt / \(descriptor.fusedTimeSteps)) * f2 / m
          """)
      } else {
        integrator.addComputePerDof(variable: "v", expression: """
          v + 0.5 * (dt * 2 / \(descriptor.fusedTimeSteps)) * (force0 + f1) / m
          """)
        integrator.addComputePerDof(variable: "v", expression: """
          v + 0.25 * (dt * 2 / \(descriptor.fusedTimeSteps)) * f2 / m
          """)
      }
      
      integrator.addComputePerDof(variable: "x", expression: """
        x + 0.5 * (dt / \(descriptor.fusedTimeSteps)) * v
        """)
      integrator.addComputePerDof(variable: "v", expression: """
        v + 0.5 * (dt / \(descriptor.fusedTimeSteps)) * f2 / m
        """)
      integrator.addComputePerDof(variable: "x", expression: """
        x + 0.5 * (dt / \(descriptor.fusedTimeSteps)) * v
        """)
      
      if i + 1 == descriptor.fusedTimeSteps {
        integrator.addComputePerDof(variable: "v", expression: """
          v + 0.25 * (dt / \(descriptor.fusedTimeSteps)) * f2 / m
          """)
        integrator.addComputePerDof(variable: "v", expression: """
          v + 0.5 * (dt / \(descriptor.fusedTimeSteps)) * (force0 + f1) / m
          """)
      }
    }
  }
  
  /// Modeled after how the OpenMM `integrator.step` API is typically used -
  /// without an argument label for steps.
  func step(_ steps: Int, timeStep: Double) {
    integrator.stepSize = timeStep
    integrator.step(steps)
  }
}
