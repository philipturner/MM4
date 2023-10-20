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
  var parameters: MM4Parameters
  
  var reorderedIndices: [Int32]
  
  var bondPairs: OpenMM_BondArray
  
  init(parameters: MM4Parameters) {
    self.parameters = parameters
    
    // Store a mapping from current indices to reversed indices in the force
    // objects. Eventually, separate the atoms into two groups of "small" vs
    // "large" atoms, creating different zones of internally contiguous tiles.
    reorderedIndices = parameters.atoms.atomicNumbers.indices.map {
      Int32(parameters.atoms.atomicNumbers.count - 1 - $0)
    }
    
    let bonds = parameters.bonds
    bondPairs = OpenMM_BondArray(size: bonds.indices.count)
    for bondID in bonds.indices.indices {
      let bond = bonds.indices[bondID]
      bondPairs[bondID] = reorder(bond)
    }
    
    // Forces:
    // - [x] Angles
    // - [x] Bonds
    // - [x] Electrostatic
    // - [x] External
    // - [x] Nonbonded
    // - [x] Torsions
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
