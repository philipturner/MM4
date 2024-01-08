//
//  MM4RigidBody+Acceleration.swift
//
//
//  Created by Philip Turner on 1/7/24.
//

// Instead of net force and net torque, expose an API for linear and angular
// accelerations. In the reference frame, the moment of inertia is
// diagonalized and does not need to be inverted. This property may reduce the
// amount of rounding error when computing accelerations.

extension MM4RigidBody {
  public var forces: [SIMD3<Float>] {
    fatalError("Not implemented.")
  }
  
  public func setForces() {
    fatalError("Not implemented.")
  }
  
  public var angularAcceleration: SIMD3<Double> {
    fatalError("Not implemented.")
  }
  
  public var linearAcceleration: SIMD3<Double> {
    fatalError("Not implemented.")
  }
}
