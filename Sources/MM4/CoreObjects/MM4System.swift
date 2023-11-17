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
  /// Concise method for fetching the atom count. This should not be exposed to
  /// the public API.
  var atomCount: Int { parameters.atoms.atomicNumbers.count }
  
  /// Bond pairs using the reordered indices.
  var bondPairs: OpenMM_BondArray
  
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
  
  /// Map from reordered indices to potential virtual site indices.
  var virtualSiteIndices: [Int32] = []
  
  init(parameters: MM4Parameters) {
    // Initialize base properties.
    self.bondPairs = OpenMM_BondArray(size: parameters.bonds.indices.count)
    self.system = OpenMM_System()
    self.parameters = parameters
    
    // Create virtual sites.
    self.createReorderedIndices(parameters: parameters)
    self.createMasses(parameters: parameters)
    self.createVirtualSites(parameters: parameters)
    
    // Create force objects.
    self.createBondPairs(parameters: parameters)
    self.forces = MM4Forces(system: self)
    forces.addForces(to: system)
  }
}

extension MM4System {
  private func createBondPairs(parameters: MM4Parameters) {
    let bonds = parameters.bonds
    for bondID in bonds.indices.indices {
      let bond = bonds.indices[bondID]
      bondPairs[bondID] = reorder(bond)
    }
  }
  
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
