//
//  MM4RigidBody+Velocity.swift
//  MM4
//
//  Created by Philip Turner on 11/20/23.
//

extension MM4RigidBody {
  /// The velocity (in nanometers per picosecond) of each atom.
  public var velocities: [SIMD3<Float>] {
    _read {
      storage.ensureVelocitiesCached()
      yield storage.velocities!
    }
  }
  
  /// The net linear momentum, in yoctogram-nanometers per picosecond.
  public var linearMomentum: SIMD3<Double> {
    _read {
      fatalError("Not implemented.")
    }
    _modify {
      ensureUniquelyReferenced()
      fatalError("Not implemented.")
    }
  }
  
  /// The net angular momentum, in yoctogram-radians per picosecond.
  public var angularMomentum: SIMD3<Double> {
    _read {
      fatalError("Not implemented.")
    }
    _modify {
      ensureUniquelyReferenced()
      fatalError("Not implemented.")
    }
  }
}

extension MM4RigidBodyStorage {
  // Use this function to import velocities into the force field.
  func createVelocities(
    _ buffer: UnsafeMutableBufferPointer<SIMD3<Float>>
  ) {
    guard buffer.count == atoms.count else {
      fatalError("Velocity buffer was not the correct size.")
    }
    let baseAddress = buffer.baseAddress.unsafelyUnwrapped
    
    // TODO: Transform into the global reference frame.
    for vID in 0..<atoms.vectorCount {
      let x = vVelocities[vID &* 3 &+ 0]
      let y = vVelocities[vID &* 3 &+ 1]
      let z = vVelocities[vID &* 3 &+ 2]
      swizzleFromVectorWidth((x, y, z), vID, baseAddress)
    }
  }
  
  func createVectorizedVelocities(_ velocities: [SIMD3<Float>]?) {
    vVelocities = Array(repeating: .zero, count: 3 * atoms.vectorCount)
    guard let velocities else {
      return
    }
    precondition(
      velocities.count == atoms.count,
      "Initial velocities array has incorrect size.")
    
    velocities.withContiguousStorageIfAvailable {
      let baseAddress = $0.baseAddress.unsafelyUnwrapped
      for vID in 0..<atoms.vectorCount {
        let (x, y, z) = swizzleToVectorWidth(vID, baseAddress)
        vVelocities[vID &* 3 &+ 0] = x
        vVelocities[vID &* 3 &+ 1] = y
        vVelocities[vID &* 3 &+ 2] = z
      }
    }
  }
  
  func createLinearMomentum() -> SIMD3<Double> {
    var momentum: SIMD3<Double> = .zero
    withSegmentedLoop(chunk: 256) {
      var vMomentumX: MM4FloatVector = .zero
      var vMomentumY: MM4FloatVector = .zero
      var vMomentumZ: MM4FloatVector = .zero
      for vID in $0 {
        let x = vVelocities[vID &* 3 &+ 0]
        let y = vVelocities[vID &* 3 &+ 1]
        let z = vVelocities[vID &* 3 &+ 2]
        let mass = vMasses[vID]
        vMomentumX.addProduct(mass, x)
        vMomentumY.addProduct(mass, y)
        vMomentumZ.addProduct(mass, z)
      }
      momentum.x += MM4DoubleVector(vMomentumX).sum()
      momentum.y += MM4DoubleVector(vMomentumY).sum()
      momentum.z += MM4DoubleVector(vMomentumZ).sum()
    }
    return momentum
  }
  
  func createAngularVelocity() -> SIMD3<Float> {
    var momentum: SIMD3<Double> = .zero
    ensureCenterOfMassCached()
    guard let centerOfMass else {
      fatalError("This should never happen.")
    }
    withSegmentedLoop(chunk: 256) {
      var vMomentumX: MM4FloatVector = .zero
      var vMomentumY: MM4FloatVector = .zero
      var vMomentumZ: MM4FloatVector = .zero
      for vID in $0 {
        let rX = vPositions[vID &* 3 &+ 0] - centerOfMass.x
        let rY = vPositions[vID &* 3 &+ 1] - centerOfMass.y
        let rZ = vPositions[vID &* 3 &+ 2] - centerOfMass.z
        let vX = vVelocities[vID &* 3 &+ 0]
        let vY = vVelocities[vID &* 3 &+ 1]
        let vZ = vVelocities[vID &* 3 &+ 2]
        let mass = vMasses[vID]
        vMomentumX.addProduct(mass, rY * vZ - rZ * vY)
        vMomentumY.addProduct(mass, rZ * vX - rX * vZ)
        vMomentumZ.addProduct(mass, rX * vY - rY * vX)
      }
      momentum.x += MM4DoubleVector(vMomentumX).sum()
      momentum.y += MM4DoubleVector(vMomentumY).sum()
      momentum.z += MM4DoubleVector(vMomentumZ).sum()
    }
    
    precondition(
      momentOfInertia .!= .zero, "Moment of inertia was not yet initialized.")
    
    
    
    
    fatalError(
      "Angular velocity cannot be created because moment of inertia and reference frame are unknown.")
    
    ensureMomentOfInertiaCached()
    guard let momentOfInertia else {
      fatalError("This should never happen.")
    }
    let inverse = invertMatrix3x3(momentOfInertia)
    let velocityX = inverse.0 * Float(momentum.x)
    let velocityY = inverse.1 * Float(momentum.y)
    let velocityZ = inverse.2 * Float(momentum.z)
    return velocityX + velocityY + velocityZ
  }
}


