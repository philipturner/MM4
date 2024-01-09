//
//  MM4RigidBody+Force.swift
//
//
//  Created by Philip Turner on 1/7/24.
//

extension MM4RigidBody {
  /// The forces (in piconewtons) on each atom.
  ///
  /// Force is a function of position and the surrounding environment. When the
  /// rigid body changes position, the forces change as well. Therefore, the
  /// forces are invalidated (turned into `nil`). To restore the forces,
  /// you must assign them to the `forces` property again.
  public var forces: [SIMD3<Float>]? {
    _read {
      fatalError("Not implemented.")
    }
    _modify {
      fatalError("Not implemented.")
    }
  }
  
  public var netForce: SIMD3<Double>? {
    fatalError("Not implemented.")
  }
  
  // Make a similar warning to that in 'angularMomentum'.
  public var netTorque: SIMD3<Double>? {
    fatalError("Not implemented.")
  }
}
