//
//  MM4Integrator.swift
//  MM4  
//
//  Created by Philip Turner on 10/14/23.
//

import OpenMM

/// A configuration for an integrator.
struct MM4IntegratorDescriptor: Hashable {
  /// Whether to correct velocities for the start of leapfrog integration
  /// intervals.
  var start: Bool = false
  
  /// Whether to correct velocities for the end of leapfrog integration
  /// intervals.
  var end: Bool = false
  
  init() {
    
  }
  
  static func == (
    lhs: MM4IntegratorDescriptor,
    rhs: MM4IntegratorDescriptor
  ) -> Bool {
    guard lhs.start == rhs.start,
          lhs.end == rhs.end else {
      return false
    }
    return true
  }
  
  func hash(into hasher: inout Hasher) {
    hasher.combine(start)
    hasher.combine(end)
  }
}

class MM4Integrator {
  var integrator: OpenMM_CustomIntegrator
  
  /// Create an integrator using the specified configuration.
  init(descriptor: MM4IntegratorDescriptor) {
    self.integrator = OpenMM_CustomIntegrator(stepSize: 0)
    
    if descriptor.start || true {
      integrator.addComputePerDof(variable: "v", expression: """
        v + 0.5 * dt * f1 / m
        """)
      integrator.addComputePerDof(variable: "v", expression: """
        v + 0.25 * dt * f2 / m
        """)
    } else {
      integrator.addComputePerDof(variable: "v", expression: """
        v + 1.0 * dt * f1 / m
        """)
      integrator.addComputePerDof(variable: "v", expression: """
        v + 0.5 * dt * f2 / m
        """)
    }
    
    integrator.addComputePerDof(variable: "x", expression: """
      x + 0.5 * dt * v
      """)
    integrator.addConstrainPositions()
    
    integrator.addComputePerDof(variable: "v", expression: """
      v + 0.5 * dt * f2 / m
      """)
    integrator.addComputePerDof(variable: "x", expression: """
      x + 0.5 * dt * v
      """)
    integrator.addConstrainPositions()
    
    if descriptor.end || true {
      integrator.addComputePerDof(variable: "v", expression: """
        v + 0.25 * dt * f2 / m
        """)
      integrator.addComputePerDof(variable: "v", expression: """
        v + 0.5 * dt * f1 / m
        """)
    }
  }
}
