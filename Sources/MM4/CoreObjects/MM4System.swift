//
//  MM4System.swift
//  MM4
//
//  Created by Philip Turner on 10/14/23.
//

import OpenMM
import QuartzCore

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
  
  /// The number of virtual sites in the system.
  var virtualSiteCount: Int = 0
  
  init(parameters: MM4Parameters, descriptor: MM4ForceFieldDescriptor) {
    let checkpoint0 = CACurrentMediaTime()
    
    // Initialize base properties.
    self.system = OpenMM_System()
    self.parameters = parameters
    
    let checkpoint1 = CACurrentMediaTime()
    
    // Create virtual sites.
    self.createReorderedIndices()
    self.createMasses()
    self.createVirtualSites()
    
    let checkpoint2 = CACurrentMediaTime()
    
    // Create force objects.
    self.forces = MM4Forces(system: self, descriptor: descriptor)
    let checkpoint3 = CACurrentMediaTime()
    forces.addForces(to: system)
    let checkpoint4 = CACurrentMediaTime()
    
    print()
    print("MM4System.init")
    print("0 to 1: \(checkpoint1 - checkpoint0)")
    print("1 to 2: \(checkpoint2 - checkpoint1)")
    print("2 to 3: \(checkpoint3 - checkpoint2)")
    print("3 to 4: \(checkpoint4 - checkpoint3)")
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
