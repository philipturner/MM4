//
//  MM4Integrator.swift
//  MM4  
//
//  Created by Philip Turner on 10/14/23.
//

import OpenMM

/// Options for customizing the integrator used for simulation.
public enum MM4IntegratorOptions {
  /// Executes bonded forces at twice the rate of nonbonded forces.
  ///
  /// The default time step is 4.35 fs.
  ///
  /// This integrator uses a multiple time-stepping (MTS) scheme. Cheaper bonded
  /// forces, such as bond-stretch and bond-bend, are only stable at ~2 fs
  /// without constraints. Expensive forces like torsions, nonbonded, and
  /// electrostatic can execute at ~4 fs.
  ///
  /// MTS has a smaller coefficient to O(n) scaling than the
  /// Verlet integrator. However, it has a larger O(1) prefactor. Only use
  /// it for large systems that are not latency-bound.
  case multipleTimeStep
  
  /// Executes all forces at the same rate.
  ///
  /// The default time step is 2.5 fs.
  ///
  /// > WARNING: The velocity Verlet integrator does not report correct
  ///   energies. For precise energy measurements in energy-conserving
  ///   simulations, consider using the `.multipleTimeStep` integrator.
  case verlet
}

/// A configuration for a custom integrator.
struct MM4CustomIntegratorDescriptor: Hashable {
  /// Whether to correct velocities for the start of leapfrog integration
  /// intervals.
  var start: Bool = false
  
  /// Whether to correct velocities for the end of leapfrog integration
  /// intervals.
  var end: Bool = false
  
  init() {
    
  }
  
  static func == (
    lhs: MM4CustomIntegratorDescriptor,
    rhs: MM4CustomIntegratorDescriptor
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

class MM4CustomIntegrator {
  var integrator: OpenMM_CustomIntegrator
  
  /// Create an integrator using the specified configuration.
  init(descriptor: MM4CustomIntegratorDescriptor) {
    self.integrator = OpenMM_CustomIntegrator(stepSize: 0)
    
    integrator.addConstrainPositions()
    if descriptor.start {
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
    
    if descriptor.end {
      integrator.addConstrainPositions()
      integrator.addComputePerDof(variable: "v", expression: """
        v + 0.5 * dt * f1 / m
        """)
      integrator.addComputePerDof(variable: "v", expression: """
        v + 0.25 * dt * f2 / m
        """)
    }
  }
}
