//
//  MM4RigidBody+Position.swift
//
//
//  Created by Philip Turner on 11/20/23.
//

// positions, center of mass

// cache the center of mass when positions are changed

extension MM4RigidBody {
  public var positions: [SIMD3<Float>] {
    get { fatalError() }
    set { fatalError() }
  }
  
  @_specialize(where T == Double)
  @_specialize(where T == Float)
  mutating func setPositions<T: BinaryFloatingPoint>(
    _ buffer: UnsafeBufferPointer<SIMD3<T>>
  ) {
    guard buffer.count == atomCount else {
      fatalError("Position buffer was not the correct size.")
    }
  }
  
  @_specialize(where T == Double)
  @_specialize(where T == Float)
  func getPositions<T: BinaryFloatingPoint>(
    _ buffer: UnsafeMutableBufferPointer<SIMD3<T>>
  ) {
    guard buffer.count == atomCount else {
      fatalError("Position buffer was not the correct size.")
    }
  }
}
