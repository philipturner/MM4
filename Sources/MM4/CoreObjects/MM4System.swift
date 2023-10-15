//
//  MM4System.swift
//
//
//  Created by Philip Turner on 10/14/23.
//

/// Encapsulates an OpenMM system and the associated force objects.
///
/// This object takes ownership of the `parameters` passed in.
class MM4System {
  var parameters: MM4Parameters
  
  var reorderedIndices: [Int32]
  
  init(parameters: MM4Parameters) {
    self.parameters = parameters
    
    // Store a mapping from current indices to reversed indices in the force
    // objects. Eventually, separate the atoms into two groups of "small" vs
    // "large" atoms, creating different zones of internally contiguous tiles.
    reorderedIndices = parameters.atoms.atomicNumbers.indices.map {
      Int32(parameters.atoms.atomicNumbers.count - 1 - $0)
    }
  }
}

extension MM4System {
  @inline(__always)
  func reorder<T: SIMD>(_ indices: T) -> T
  where T.Scalar == Int32 {
    var output = indices
    for i in 0..<T.scalarCount {
      output[i] = reorderedIndices[Int(output[i])]
    }
    return output
  }
}
