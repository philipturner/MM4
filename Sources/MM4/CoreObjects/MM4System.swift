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
  /// Bond pairs using the reordered indices.
  var bondPairs: OpenMM_BondArray
  
  /// The forces used by the system.
  var forces: MM4Forces!
  
  /// The location where the parameters are owned.
  var parameters: MM4Parameters
  
  /// Indices may eventually be rearranged for performance.
  var reorderedIndices: [Int32]
  
  /// The backing OpenMM system object.
  var system: OpenMM_System
  
  init(parameters: MM4Parameters) {
    self.parameters = parameters
    
    self.system = OpenMM_System()
    for mass in parameters.atoms.masses {
      system.addParticle(mass: Double(mass))
    }
    
    // Store a mapping from current indices to reversed indices in the force
    // objects. Eventually, separate the atoms into two groups of "small" vs
    // "large" atoms, creating different zones of internally contiguous tiles.
    //
    // The index reversing is a litmus test, to ensure the code is aware of the
    // index reordering at every step.
    self.reorderedIndices = parameters.atoms.atomicNumbers.indices.map {
      Int32(parameters.atoms.atomicNumbers.count - 1 - $0)
    }
    
    let bonds = parameters.bonds
    self.bondPairs = OpenMM_BondArray(size: bonds.indices.count)
    for bondID in bonds.indices.indices {
      let bond = bonds.indices[bondID]
      bondPairs[bondID] = reorder(bond)
    }
    
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
