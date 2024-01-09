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

// TODO: Function for creating net force and torque separately. Until we start
// aggressively optimizing for performance, fusion would be a premature
// optimization that makes debugging harder.
//
// Or make them together. Start with a function that computes torque. It should
// be easy to make it accumulate force as well.

// create

extension MM4RigidBodyStorage {
  func createNetTorque(
    _ buffer: UnsafeBufferPointer<SIMD3<Float>>
  ) -> SIMD3<Double> {
    guard buffer.count == atoms.count else {
      fatalError("Force buffer was not the correct size.")
    }
    let baseAddress = buffer.baseAddress.unsafelyUnwrapped
    
    var netTorque: SIMD3<Double> = .zero
    withSegmentedLoop(chunk: 256) {
      var vTorqueX: MM4FloatVector = .zero
      var vTorqueY: MM4FloatVector = .zero
      var vTorqueZ: MM4FloatVector = .zero
      for vID in $0 {
        let (x, y, z) = swizzleToVectorWidth(vID, baseAddress)
        
      }
    }
  }
}
