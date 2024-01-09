//
//  MM4RigidBody+Position.swift
//  MM4
//
//  Created by Philip Turner on 11/20/23.
//

import QuaternionModule

// MARK: - Public API

extension MM4RigidBody {
  /// The position (in nanometers) of each atom's nucleus.
  public var positions: [SIMD3<Float>] {
    storage.ensurePositionsCached()
    return storage.positions!
  }
  
  /// The total mass (in yoctograms).
  public var mass: Double {
    storage.mass
  }
  
  /// The center of mass, in nanometers.
  ///
  /// When the center of mass is modified, forces are invalidated.
  public var centerOfMass: SIMD3<Double> {
    _read {
      yield storage.centerOfMass
    }
    _modify {
      ensureUniquelyReferenced()
      storage.invalidatePositions()
      yield &storage.centerOfMass
    }
  }
  
  /// The eigenvalues of the inertia tensor, in descending order.
  public var momentOfInertia: SIMD3<Double> {
    storage.momentOfInertia
  }
  
  /// The eigenvectors of the inertia tensor.
  public var principalAxes: (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>) {
    storage.principalAxes
  }
  
  /// Rotate the object.
  /// - parameter angle: The angle to rotate, in radians.
  /// - parameter axis: The normalized vector to rotate around.
  ///
  /// If `axis` is not specified, the default value aligns with the current
  /// angular momentum. If the angular momentum is zero, you must enter zero
  /// for the angle.
  public mutating func rotate(angle: Double) {
    // TODO: Delete the comment. The angular momentum is never mutated
    // because it's always relative to the principal axes.
    //
    // If you want to add an 'axis' argument in the future (with a default value
    // of NAN, that will not be a breaking API change). However, I see no
    // reason to include an axis argument at the moment. We don't do rotations,
    // momentum, or torque with respect to axes in the global coordinate space.
    // That would complicate the function quite a bit.
    //
    // The intended use case is not constructing a scene. People can do that
    // with quaternion rotations. The use case is integration during a rigid
    // body dynamics simulation. This restricts what kind of data
    // transformations you can perform. One restriction is that rotations must
    // occur in the direction the object is actually rotating.
    ensureUniquelyReferenced()
    storage.invalidatePositions()
    storage.invalidateVelocities()
    
    let magnitude = (angularMomentum * angularMomentum).sum()
    if magnitude < .leastNormalMagnitude else {
      precondition(angle.magnitude < .leastNormalMagnitude, "Angular momentum was zero, but rotation angle was nonzero: \(angle)")
    }
    
    if (angularMomentum * angularMomentum).sum() < .leastNormalMagnitude {
      guard angle.magnitude
    }
    
    if let axis {
      fatalError("'axis' argument not supported yet.")
    } else {
      
    }
    
    // Create a rotation from a quaternion. Use it to generate a rotation
    // matrix.
  }
}

// MARK: - Linear Properties

extension MM4RigidBodyStorage {
  func createVectorizedMasses(_ masses: [Float]) {
    vMasses = Array(repeating: .zero, count: atoms.vectorCount)
    
    masses.withContiguousStorageIfAvailable {
      let baseAddress = $0.baseAddress.unsafelyUnwrapped
      let vAddress = UnsafePointer<MM4FloatVector>(OpaquePointer(baseAddress))
      
      withSegmentedLoop(chunk: 256) {
        for vID in $0 {
          var vMass: MM4FloatVector = .zero
          if vID == atoms.vectorCount &- 1 {
            let remaining = atoms.count &- vID &* MM4VectorWidth
            for lane in 0..<remaining {
              vMass[lane] = baseAddress[vID &* MM4VectorWidth &+ lane]
            }
          } else {
            vMass = vAddress[vID]
          }
          vMasses[vID] = vMass
        }
      }
    }
  }
  
  func createVectorizedPositions(_ positions: [SIMD3<Float>]?) {
    vPositions = Array(repeating: .zero, count: 3 * atoms.vectorCount)
    guard let positions else {
      fatalError("Initial positions must be specified.")
    }
    precondition(
      positions.count == atoms.count,
      "Initial velocities array has incorrect size.")
    
    positions.withContiguousStorageIfAvailable {
      let baseAddress = $0.baseAddress.unsafelyUnwrapped
      for vID in 0..<atoms.vectorCount {
        let (x, y, z) = swizzleToVectorWidth(vID, baseAddress)
        vPositions[vID &* 3 &+ 0] = x
        vPositions[vID &* 3 &+ 1] = y
        vPositions[vID &* 3 &+ 2] = z
      }
    }
  }
  
  func createMass() -> Double {
    var mass: Double = .zero
    withSegmentedLoop(chunk: 256) {
      var vMassAccumulator: MM4FloatVector = .zero
      for vID in $0 {
        vMassAccumulator += vMasses[vID]
      }
      mass += MM4DoubleVector(vMassAccumulator).sum()
    }
    return mass
  }
  
  func createCenterOfMass() -> SIMD3<Double> {
    var centerOfMass: SIMD3<Double> = .zero
    withSegmentedLoop(chunk: 256) {
      var vAccumulatorX: MM4FloatVector = .zero
      var vAccumulatorY: MM4FloatVector = .zero
      var vAccumulatorZ: MM4FloatVector = .zero
      for vID in $0 {
        let x = vPositions[vID &* 3 &+ 0]
        let y = vPositions[vID &* 3 &+ 1]
        let z = vPositions[vID &* 3 &+ 2]
        let mass = vMasses[vID]
        vAccumulatorX.addProduct(mass, x)
        vAccumulatorY.addProduct(mass, y)
        vAccumulatorZ.addProduct(mass, z)
      }
      centerOfMass.x += MM4DoubleVector(vAccumulatorX).sum()
      centerOfMass.y += MM4DoubleVector(vAccumulatorY).sum()
      centerOfMass.z += MM4DoubleVector(vAccumulatorZ).sum()
    }
    
    precondition(mass > 0, "Mass must be positive.")
    return centerOfMass / mass
  }
  
  func normalizeLinearPositions(to centerOfMass: SIMD3<Double>) {
    for vID in 0..<atoms.vectorCount {
      vPositions[vID &* 3 &+ 0] -= Float(centerOfMass.x)
      vPositions[vID &* 3 &+ 1] -= Float(centerOfMass.y)
      vPositions[vID &* 3 &+ 2] -= Float(centerOfMass.z)
    }
  }
}

// MARK: - Angular Properties

// These functions assume the linear position is already normalized.

extension MM4RigidBodyStorage {
  func createInertiaTensor() -> (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>) {
    var inertiaTensor: (
      SIMD3<Double>, SIMD3<Double>, SIMD3<Double>
    ) = (.zero, .zero, .zero)
    
    withSegmentedLoop(chunk: 256) {
      var vXX: MM4FloatVector = .zero
      var vYY: MM4FloatVector = .zero
      var vZZ: MM4FloatVector = .zero
      var vXY: MM4FloatVector = .zero
      var vXZ: MM4FloatVector = .zero
      var vYZ: MM4FloatVector = .zero
      for vID in $0 {
        let x = vPositions[vID &* 3 &+ 0]
        let y = vPositions[vID &* 3 &+ 1]
        let z = vPositions[vID &* 3 &+ 2]
        let mass = vMasses[vID]
        vXX.addProduct(mass, x * x)
        vYY.addProduct(mass, y * y)
        vZZ.addProduct(mass, z * z)
        vXY.addProduct(mass, x * y)
        vXZ.addProduct(mass, x * z)
        vYZ.addProduct(mass, y * z)
      }
      
      let XX = MM4DoubleVector(vXX).sum()
      let YY = MM4DoubleVector(vYY).sum()
      let ZZ = MM4DoubleVector(vZZ).sum()
      let XY = MM4DoubleVector(vXY).sum()
      let XZ = MM4DoubleVector(vXZ).sum()
      let YZ = MM4DoubleVector(vYZ).sum()
      inertiaTensor.0 += SIMD3<Double>(YY + ZZ, -XY, -XZ)
      inertiaTensor.1 += SIMD3<Double>(-XY, XX + ZZ, -YZ)
      inertiaTensor.2 += SIMD3<Double>(-XZ, -YZ, XX + YY)
    }
    return inertiaTensor
  }
  
  func normalizeOrientation(
    to principalAxes: (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>)
  ) {
    let Σ = (
      SIMD3<Float>(principalAxes.0),
      SIMD3<Float>(principalAxes.1),
      SIMD3<Float>(principalAxes.2))
    
    for vID in 0..<atoms.vectorCount {
      let rX = vPositions[vID &* 3 &+ 0]
      let rY = vPositions[vID &* 3 &+ 1]
      let rZ = vPositions[vID &* 3 &+ 2]
      let r = (rX, rY, rZ)
      vPositions[vID &* 3 &+ 0] = dot(vector: Σ.0, scalars: r)
      vPositions[vID &* 3 &+ 1] = dot(vector: Σ.1, scalars: r)
      vPositions[vID &* 3 &+ 2] = dot(vector: Σ.2, scalars: r)
      
      let vX = vVelocities[vID &* 3 &+ 0]
      let vY = vVelocities[vID &* 3 &+ 1]
      let vZ = vVelocities[vID &* 3 &+ 2]
      let v = (vX, vY, vZ)
      vVelocities[vID &* 3 &+ 0] = dot(vector: Σ.0, scalars: v)
      vVelocities[vID &* 3 &+ 1] = dot(vector: Σ.1, scalars: v)
      vVelocities[vID &* 3 &+ 2] = dot(vector: Σ.2, scalars: v)
    }
  }
  
  func createPositions(
    _ buffer: UnsafeMutableBufferPointer<SIMD3<Float>>
  ) {
    guard buffer.count == atoms.count else {
      fatalError("Position buffer was not the correct size.")
    }
    let baseAddress = buffer.baseAddress.unsafelyUnwrapped
    
    let Σ = (
      SIMD3<Float>(principalAxes.0),
      SIMD3<Float>(principalAxes.1),
      SIMD3<Float>(principalAxes.2))
    let r̄ = SIMD3<Float>(centerOfMass)
    
    for vID in 0..<atoms.vectorCount {
      let rX = vPositions[vID &* 3 &+ 0]
      let rY = vPositions[vID &* 3 &+ 1]
      let rZ = vPositions[vID &* 3 &+ 2]
      
      // Σ * r + r̄
      var x = Σ.0.x * rX + Σ.1.x * rY + Σ.2.x * rZ
      var y = Σ.0.y * rX + Σ.1.y * rY + Σ.2.y * rZ
      var z = Σ.0.z * rX + Σ.1.z * rY + Σ.2.z * rZ
      x += r̄.x
      y += r̄.y
      z += r̄.z
      
      swizzleFromVectorWidth((x, y, z), vID, baseAddress)
    }
  }
}
