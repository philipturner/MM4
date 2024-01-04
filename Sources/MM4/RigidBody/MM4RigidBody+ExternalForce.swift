//
//  MM4RigidBody+ExternalForce.swift
//
//
//  Created by Philip Turner on 12/22/23.
//

import Numerics

extension MM4RigidBody {
  /// The constant force (in piconewtons) exerted on each atom.
  public var externalForces: [SIMD3<Float>] {
    // _modify not supported b/c it requires checking whether an anchor was
    // changed to something invalid. Implementing this in an ergonomic way would
    // result in an O(n^2) implementation, as the atom that was changed during
    // _modify is unknown.
    _read {
      yield storage.externalForces
    }
  }
  
  @_specialize(where T == Double)
  @_specialize(where T == Float)
  public mutating func setExternalForces<T: BinaryFloatingPoint>(
    _ array: Array<SIMD3<T>>
  ) {
    array.withUnsafeBufferPointer {
      setExternalForces($0)
    }
  }
  
  @_specialize(where T == Double)
  @_specialize(where T == Float)
  public mutating func setExternalForces<T: BinaryFloatingPoint>(
    _ buffer: UnsafeBufferPointer<SIMD3<T>>
  ) {
    ensureUniquelyReferenced()
    storage.eraseRarelyCachedProperties()
    
    guard buffer.count == storage.atoms.count else {
      fatalError("Velocity buffer was not the correct size.")
    }
    let baseAddress = buffer.baseAddress.unsafelyUnwrapped
    
    for atomID in 0..<storage.atoms.count {
      storage.externalForces[atomID] = SIMD3(baseAddress[atomID])
    }
    storage.ensureAnchorExternalForcesValid()
  }
  
  @_specialize(where T == Double)
  @_specialize(where T == Float)
  public func getExternalForces<T: BinaryFloatingPoint>(
    _ buffer: UnsafeMutableBufferPointer<SIMD3<T>>
  ) {
    guard buffer.count == storage.atoms.count else {
      fatalError("Velocity buffer was not the correct size.")
    }
    let baseAddress = buffer.baseAddress.unsafelyUnwrapped
    
    for atomID in 0..<storage.atoms.count {
      baseAddress[atomID] = SIMD3(storage.externalForces[atomID])
    }
  }
}

extension MM4RigidBodyStorage {
  // WARNING: Call this every time the external forces change.
  func ensureAnchorExternalForcesValid() {
    var valid = true
    for anchor in anchors {
      let externalForce = externalForces[Int(anchor)]
      if externalForce != .zero {
        valid = false
      }
    }
    guard valid else {
      fatalError("Anchor external forces are invalid.")
    }
  }
}
