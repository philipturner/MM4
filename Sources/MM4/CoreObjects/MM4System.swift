//
//  MM4System.swift
//
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
  
  /// Map from reordered indices to original indices.
  var originalIndices: [Int32] = []
  
  /// The location where the parameters are owned.
  var parameters: MM4Parameters
  
  /// Indices may eventually be rearranged for performance.
  var reorderedIndices: [Int32] = []
  
  /// The backing OpenMM system object.
  var system: OpenMM_System
  
  /// The number of virtual sites.
  var virtualSiteCount: Int = 0
  
  /// Whether each atom is a virtual site.
  var virtualSiteMask: [Bool] = []
  
  init(parameters: MM4Parameters) {
    // Initialize base properties.
    self.system = OpenMM_System()
    self.parameters = parameters
    
    // Create virtual sites.
    self.createReorderedIndices()
    self.createMasses()
    self.createVirtualSites()
    self.createVirtualSiteMask()
    
    // Create force objects.
    self.forces = MM4Forces(system: self)
    forces.addForces(to: system)
  }
}

extension MM4System {
  @inline(__always)
  func reorder(_ indices: SIMD2<Int32>) -> SIMD2<Int> {
    var output: SIMD2<Int32> = .zero
    for i in 0..<indices.scalarCount {
      output[i] = reorderedIndices[Int(output[i])]
    }
    return SIMD2(truncatingIfNeeded: output)
  }
  
  @inline(__always)
  func reorder(_ indices: SIMD3<Int32>) -> SIMD3<Int> {
    var output: SIMD3<Int32> = .zero
    for i in 0..<indices.scalarCount {
      output[i] = reorderedIndices[Int(output[i])]
    }
    return SIMD3(truncatingIfNeeded: output)
  }
  
  @inline(__always)
  func reorder(_ indices: SIMD4<Int32>) -> SIMD4<Int> {
    var output: SIMD4<Int32> = .zero
    for i in 0..<indices.scalarCount {
      output[i] = reorderedIndices[Int(output[i])]
    }
    return SIMD4(truncatingIfNeeded: output)
  }
}
