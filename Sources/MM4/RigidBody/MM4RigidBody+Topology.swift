//
//  MM4RigidBody+Topology.swift
//  
//
//  Created by Philip Turner on 11/20/23.
//

// positions, center of mass, moment of inertia

extension MM4RigidBody {
  public var positions: [SIMD3<Float>] {
    get { fatalError() }
    set { fatalError() }
  }
  
  // functions for writing/reading the computed properties from
  // pre-allocated memory or slices of larger buffers
}
