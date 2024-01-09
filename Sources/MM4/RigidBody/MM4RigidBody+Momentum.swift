//
//  MM4RigidBody+Velocity.swift
//  MM4
//
//  Created by Philip Turner on 11/20/23.
//

// Velocity and Momentum

// MARK: - Public API

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

// MARK: - Linear Properties

extension MM4RigidBodyStorage {
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
    var linearMomentum: SIMD3<Double> = .zero
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
      linearMomentum.x += MM4DoubleVector(vMomentumX).sum()
      linearMomentum.y += MM4DoubleVector(vMomentumY).sum()
      linearMomentum.z += MM4DoubleVector(vMomentumZ).sum()
    }
    return linearMomentum
  }
  
  func normalizeLinearVelocities(mass: Double, linearMomentum: SIMD3<Double>) {
    let linearVelocity = linearMomentum / mass
    for vID in 0..<atoms.vectorCount {
      vVelocities[vID &* 3 &+ 0] -= Float(linearVelocity.x)
      vVelocities[vID &* 3 &+ 1] -= Float(linearVelocity.y)
      vVelocities[vID &* 3 &+ 2] -= Float(linearVelocity.z)
    }
  }
}

// MARK: - Angular Properties

// These functions assume the center of mass is (0, 0, 0), and both position
// and velocity have the same orientation.

extension MM4RigidBodyStorage {
  func createAngularMomentum() -> SIMD3<Double> {
    var angularMomentum: SIMD3<Double> = .zero
    withSegmentedLoop(chunk: 256) {
      var vMomentumX: MM4FloatVector = .zero
      var vMomentumY: MM4FloatVector = .zero
      var vMomentumZ: MM4FloatVector = .zero
      for vID in $0 {
        let rX = vPositions[vID &* 3 &+ 0]
        let rY = vPositions[vID &* 3 &+ 1]
        let rZ = vPositions[vID &* 3 &+ 2]
        let vX = vVelocities[vID &* 3 &+ 0]
        let vY = vVelocities[vID &* 3 &+ 1]
        let vZ = vVelocities[vID &* 3 &+ 2]
        let mass = vMasses[vID]
        vMomentumX.addProduct(mass, rY * vZ - rZ * vY)
        vMomentumY.addProduct(mass, rZ * vX - rX * vZ)
        vMomentumZ.addProduct(mass, rX * vY - rY * vX)
      }
      angularMomentum.x += MM4DoubleVector(vMomentumX).sum()
      angularMomentum.y += MM4DoubleVector(vMomentumY).sum()
      angularMomentum.z += MM4DoubleVector(vMomentumZ).sum()
    }
    return angularMomentum
  }
  
  func normalizeAngularVelocities(
    angularMomentum: SIMD3<Double>,
    inertiaTensor: (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>)
  ) {
    
    
    /*
    let inverse = invertMatrix3x3(momentOfInertia)
    let velocityX = inverse.0 * Float(momentum.x)
    let velocityY = inverse.1 * Float(momentum.y)
    let velocityZ = inverse.2 * Float(momentum.z)
    return velocityX + velocityY + velocityZ
    
    for vID in 0..<storage.atoms.vectorCount {
      let rX = storage.vPositions[vID &* 3 &+ 0] - centerOfMass.x
      let rY = storage.vPositions[vID &* 3 &+ 1] - centerOfMass.y
      let rZ = storage.vPositions[vID &* 3 &+ 2] - centerOfMass.z
      var vX: MM4FloatVector = .zero
      var vY: MM4FloatVector = .zero
      var vZ: MM4FloatVector = .zero
      
      // Un-apply the previous bulk angular velocity.
      do {
        let w = previous
        vX -= w.y * rZ - w.z * rY
        vY -= w.z * rX - w.x * rZ
        vZ -= w.x * rY - w.y * rX
      }
      
      // Apply the next bulk angular velocity.
      do {
        let w = next
        vX += w.y * rZ - w.z * rY
        vY += w.z * rX - w.x * rZ
        vZ += w.x * rY - w.y * rX
      }
      
      // Mask out the changes to velocity for anchors.
      let mass = storage.vMasses[vID]
      vX.replace(with: MM4FloatVector.zero, where: mass .== 0)
      vY.replace(with: MM4FloatVector.zero, where: mass .== 0)
      vZ.replace(with: MM4FloatVector.zero, where: mass .== 0)
      
      // Write the new velocities as offsets relative to the existing ones.
      storage.vVelocities[vID &* 3 &+ 0] += vX
      storage.vVelocities[vID &* 3 &+ 1] += vY
      storage.vVelocities[vID &* 3 &+ 2] += vZ
    }
     */
  }
  
  // Transforms velocities into the global reference frame.
  // TODO: Transform into the global reference frame.
  func createVelocities(
    _ buffer: UnsafeMutableBufferPointer<SIMD3<Float>>
  ) {
    guard buffer.count == atoms.count else {
      fatalError("Velocity buffer was not the correct size.")
    }
    let baseAddress = buffer.baseAddress.unsafelyUnwrapped
    
    for vID in 0..<atoms.vectorCount {
      let x = vVelocities[vID &* 3 &+ 0]
      let y = vVelocities[vID &* 3 &+ 1]
      let z = vVelocities[vID &* 3 &+ 2]
      swizzleFromVectorWidth((x, y, z), vID, baseAddress)
    }
  }
}


