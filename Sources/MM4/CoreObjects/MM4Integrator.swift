//
//  MM4Integrator.swift
//  
//
//  Created by Philip Turner on 10/14/23.
//

/// A configuration for an integrator.
class MM4IntegratorDescriptor: Hashable {
  /// The number of force evaluations fused together to reduce redundant
  /// computation of forces. This must be one for a Verlet integrator.
  var fusedTimeSteps: Int = 1
  
  /// This must be false for a Verlet integrator.
  var multipleTimeSteps: Bool = false
  
  init() {
    
  }
  
  static func == (
    lhs: MM4IntegratorDescriptor,
    rhs: MM4IntegratorDescriptor
  ) -> Bool {
    guard lhs.fusedTimeSteps == rhs.fusedTimeSteps,
          lhs.multipleTimeSteps == rhs.multipleTimeSteps else {
      return false
    }
    return true
  }
  
  func hash(into hasher: inout Hasher) {
    hasher.combine(fusedTimeSteps)
    hasher.combine(multipleTimeSteps)
  }
}

class MM4Integrator {
  /// Create an integrator using the specified configuration.
  init(descriptor: MM4IntegratorDescriptor) {
    
  }
}
