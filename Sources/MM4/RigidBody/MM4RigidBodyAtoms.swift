//
//  MM4RigidBodyAtoms.swift
//
//
//  Created by Philip Turner on 11/25/23.
//

/// Wrapper for ergonomically vectorizing computations.
struct MM4RigidBodyAtoms {
  /// A convenient method for accessing the atom count.
  var count: Int
  
  /// A convenient method for accessing the atom vector count.
  var vectorCount: Int
  
  /// A convenient method for addressing the end of the list.
  var vectorMask: SIMDMask<MM4UInt32Vector.MaskStorage>
  
  init(count: Int) {
    self.count = count
    self.vectorCount = count + MM4VectorWidth - 1
    self.vectorCount /= MM4VectorWidth
    
    let lastVectorStart = UInt32((vectorCount - 1) * MM4VectorWidth)
    var lastVector: MM4UInt32Vector = .init(repeating: lastVectorStart)
    for lane in 0..<MM4VectorWidth {
      lastVector[lane] &+= UInt32(lane)
    }
    self.vectorMask = lastVector .>= UInt32(count)
  }
}
