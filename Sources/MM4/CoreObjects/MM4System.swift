//
//  MM4System.swift
//
//
//  Created by Philip Turner on 10/14/23.
//

/// Encapsulates an OpenMM system and the associated force objects.
class MM4System {
  init(parameters: MM4Parameters) {
    // Store a mapping from current indices to reversed indices in the force
    // objects. Eventually, separate the atoms into two groups of "small" vs
    // "large" atoms, creating different zones of internally contiguous tiles.
  }
}
