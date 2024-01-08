//
//  MM4RigidBody+Force.swift
//
//
//  Created by Philip Turner on 1/7/24.
//

extension MM4RigidBody {
  public var forces: [SIMD3<Float>] {
    _read {
      fatalError("Not implemented.")
    }
    _modify {
      fatalError("Not implemented.")
    }
  }
  
  public var netForce: SIMD3<Double> {
    fatalError("Not implemented.")
  }
  
  public var netTorque: SIMD3<Double> {
    fatalError("Not implemented.")
  }
}
