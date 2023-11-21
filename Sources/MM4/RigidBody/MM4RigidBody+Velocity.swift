//
//  MM4RigidBody+Velocity.swift
//
//
//  Created by Philip Turner on 11/20/23.
//

extension MM4RigidBody {
  // Handle decomposition of accepted velocities
  // Delegate the rotational velocity part to MM4RigidBody+Rotation
  // Delegate the thermal energy part to MM4RigidBody+Temperature
  
  public var velocities: [SIMD3<Float>] {
    get { fatalError() }
    set { fatalError() }
  }
  
  @_specialize(where T == Double)
  @_specialize(where T == Float)
  func setVelocities<T: BinaryFloatingPoint>(
    _ buffer: UnsafeBufferPointer<SIMD3<T>>
  ) {
    guard buffer.count == atomCount else {
      fatalError("Velocity buffer was not the correct size.")
    }
  }
  
  @_specialize(where T == Double)
  @_specialize(where T == Float)
  func getVelocities<T: BinaryFloatingPoint>(
    _ buffer: UnsafeMutableBufferPointer<SIMD3<T>>
  ) {
    guard buffer.count == atomCount else {
      fatalError("Position buffer was not the correct size.")
    }
  }
}
