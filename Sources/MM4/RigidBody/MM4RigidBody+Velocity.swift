//
//  MM4RigidBody+Velocity.swift
//
//
//  Created by Philip Turner on 11/20/23.
//

extension MM4RigidBody {
  /// The bulk + thermal velocity (in nanometers per picosecond) of each atom.
  public var velocities: [SIMD3<Float>] {
    // _modify not yet supported b/c it requires very complex caching logic.
    _read { fatalError("Not implemented.") }
  }
  
  @_specialize(where T == Double)
  @_specialize(where T == Float)
  func setVelocities<T: BinaryFloatingPoint>(
    _ buffer: UnsafeBufferPointer<SIMD3<T>>
  ) {
    guard buffer.count == atoms.count else {
      fatalError("Velocity buffer was not the correct size.")
    }
  }
  
  @_specialize(where T == Double)
  @_specialize(where T == Float)
  func getVelocities<T: BinaryFloatingPoint>(
    _ buffer: UnsafeMutableBufferPointer<SIMD3<T>>
  ) {
    guard buffer.count == atoms.count else {
      fatalError("Position buffer was not the correct size.")
    }
  }
}
