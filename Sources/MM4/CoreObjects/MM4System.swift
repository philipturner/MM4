//
//  MM4System.swift
//  MM4
//
//  Created by Philip Turner on 10/14/23.
//

import OpenMM

/// Encapsulates an OpenMM system and the associated force objects.
///
/// This object takes ownership of the `parameters` passed in.
class MM4System {
  /// The forces used by the system.
  var forces: MM4Forces!
  
  /// The location where the parameters are owned.
  var parameters: MM4Parameters
  
  /// Indices may eventually be rearranged for performance.
  var reorderedIndices: [UInt32] = []
  
  /// The backing OpenMM system object.
  var system: OpenMM_System
  
  init(parameters: MM4Parameters, descriptor: MM4ForceFieldDescriptor) {
    // Initialize base properties.
    self.system = OpenMM_System()
    self.parameters = parameters
    
    // Create virtual sites.
    self.createReorderedIndices()
    self.createMasses()
    self.createVirtualSites()
    
    // Create force objects.
    self.forces = MM4Forces(system: self, descriptor: descriptor)
    forces.addForces(to: system)
  }
}

extension MM4System {
  @_transparent
  func reorder(_ indices: SIMD2<UInt32>) -> SIMD2<Int> {
    var output: SIMD2<UInt32> = .zero
    for i in 0..<indices.scalarCount {
      output[i] = reorderedIndices[Int(indices[i])]
    }
    return SIMD2(truncatingIfNeeded: output)
  }
  
  @_transparent
  func reorder(_ indices: SIMD3<UInt32>) -> SIMD3<Int> {
    var output: SIMD3<UInt32> = .zero
    for i in 0..<indices.scalarCount {
      output[i] = reorderedIndices[Int(indices[i])]
    }
    return SIMD3(truncatingIfNeeded: output)
  }
  
  @_transparent
  func reorder(_ indices: SIMD4<UInt32>) -> SIMD4<Int> {
    var output: SIMD4<UInt32> = .zero
    for i in 0..<indices.scalarCount {
      output[i] = reorderedIndices[Int(indices[i])]
    }
    return SIMD4(truncatingIfNeeded: output)
  }
}
