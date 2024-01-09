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
  public var momentum: SIMD3<Double> {
    get {
      fatalError("Not implemented.")
    }
    set {
      ensureUniquelyReferenced()
      fatalError("Not implemented.")
    }
  }
  
  /// The net angular momentum, in yoctogram-radians per picosecond.
  ///
  /// In contrast to linear momentum, angular momentum cannot be projected onto
  /// the global reference frame. It is not parallel to the axis of rotation.
  /// Rather, it is the magnitude of angular momentum with respect to each
  /// eigenpair of the inertia tensor.
  public var angularMomentum: SIMD3<Double> {
    get {
      fatalError("Not implemented.")
    }
    set {
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
  
  func normalizeLinearVelocities(to linearVelocity: SIMD3<Double>) {
    for vID in 0..<atoms.vectorCount {
      vVelocities[vID &* 3 &+ 0] -= Float(linearVelocity.x)
      vVelocities[vID &* 3 &+ 1] -= Float(linearVelocity.y)
      vVelocities[vID &* 3 &+ 2] -= Float(linearVelocity.z)
    }
  }
}

// MARK: - Angular Properties

// These functions assume the linear position, linear velocity, and orientation
// are already normalized.

extension MM4RigidBodyStorage {
  func createAngularMomentum() -> SIMD3<Double> {
    var angularMomentum: SIMD3<Double> = .zero
    withSegmentedLoop(chunk: 256) {
      var vMomentumX: MM4FloatVector = .zero
      var vMomentumY: MM4FloatVector = .zero
      var vMomentumZ: MM4FloatVector = .zero
      for vID in $0 {
        // w = r x v
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
  
  // L = I w
  // w = I^{-1} L
  // w = (Σ Λ ΣT)^{-1} L
  // w = (ΣT)^{-1} Λ^{-1} Σ^{-1} L
  // w = Σ (1/Λ) ΣT L
  //
  // If already transformed into the eigenbasis, Σ is the identity matrix.
  // L / Λ is the formula to compute angular velocity.
  func normalizeAngularVelocities(to angularVelocity: SIMD3<Double>) {
    let w = SIMD3<Float>(angularVelocity)
    
    for vID in 0..<atoms.vectorCount {
      // v = w x r
      let rX = vPositions[vID &* 3 &+ 0]
      let rY = vPositions[vID &* 3 &+ 1]
      let rZ = vPositions[vID &* 3 &+ 2]
      vVelocities[vID &* 3 &+ 0] -= w.y * rZ - w.z * rY
      vVelocities[vID &* 3 &+ 1] -= w.z * rX - w.x * rZ
      vVelocities[vID &* 3 &+ 2] -= w.x * rY - w.y * rX
    }
  }
  
  func createVelocities(
    _ buffer: UnsafeMutableBufferPointer<SIMD3<Float>>
  ) {
    guard buffer.count == atoms.count else {
      fatalError("Velocity buffer was not the correct size.")
    }
    let baseAddress = buffer.baseAddress.unsafelyUnwrapped
    
    let Σ = (
      SIMD3<Float>(principalAxes.0),
      SIMD3<Float>(principalAxes.1),
      SIMD3<Float>(principalAxes.2))
    let w̄ = SIMD3<Float>(angularMomentum / momentOfInertia)
    let v̄ = SIMD3<Float>(linearMomentum / mass)
    
    for vID in 0..<atoms.vectorCount {
      let rX = vPositions[vID &* 3 &+ 0]
      let rY = vPositions[vID &* 3 &+ 1]
      let rZ = vPositions[vID &* 3 &+ 2]
      var vX = vVelocities[vID &* 3 &+ 0]
      var vY = vVelocities[vID &* 3 &+ 1]
      var vZ = vVelocities[vID &* 3 &+ 2]
      
      // Σ * (v + w̄ x r) + v̄
      vX += w̄.y * rZ - w̄.z * rY
      vY += w̄.z * rX - w̄.x * rZ
      vZ += w̄.x * rY - w̄.y * rX
      var x = Σ.0.x * vX + Σ.1.x * vY + Σ.2.x * vZ
      var y = Σ.0.y * vX + Σ.1.y * vY + Σ.2.y * vZ
      var z = Σ.0.z * vX + Σ.1.z * vY + Σ.2.z * vZ
      x += v̄.x
      y += v̄.y
      z += v̄.z
      
      swizzleFromVectorWidth((x, y, z), vID, baseAddress)
    }
  }
}
