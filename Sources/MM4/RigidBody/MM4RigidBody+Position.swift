//
//  MM4RigidBody+Position.swift
//  MM4
//
//  Created by Philip Turner on 11/20/23.
//

// Position and Inertia

// MARK: - Public API

extension MM4RigidBody {
  /// The position (in nanometers) of each atom's nucleus.
  public var positions: [SIMD3<Float>] {
    _read {
      storage.ensurePositionsCached()
      yield storage.positions!
    }
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
      fatalError("Not implemented.")
    }
    _modify {
      ensureUniquelyReferenced()
      fatalError("Not implemented.")
    }
  }
  
  /// The eigenvalues of the inertia tensor, in descending order.
  public var momentOfInertia: SIMD3<Double> {
    fatalError("Not implemented.")
  }
  
  /// The eigenvectors of the inertia tensor.
  public var principalAxes: (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>) {
    fatalError("Not implemented.")
  }
  
  /// Rotate the object.
  /// - parameter angle: The angle to rotate, in radians.
  /// - parameter axis: The normalized vector to rotate around.
  ///
  /// If `axis` is not specified, the default value aligns with the current
  /// angular momentum. If the angular momentum is zero, you must enter zero
  /// for the angle.
  ///
  /// If the axis aligns with angular momentum, angular momentum stays
  /// constant. If it does not align, angular momentum is rotated along with the
  /// atom positions. In either case, forces are invalidated because the atom
  /// positions have changed.
  public mutating func rotate(angle: Double, axis: SIMD3<Double>? = nil) {
    // Create a rotation from a quaternion. Use it to generate a rotation
    // matrix.
    
    ensureUniquelyReferenced()
    fatalError("Not implemented.")
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
  
  func createCenterOfMass(mass: Double) -> SIMD3<Double> {
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
  
  func normalizeLinearPositions(centerOfMass: SIMD3<Double>) {
    for vID in 0..<atoms.vectorCount {
      vPositions[vID &* 3 &+ 0] -= Float(centerOfMass.x)
      vPositions[vID &* 3 &+ 1] -= Float(centerOfMass.y)
      vPositions[vID &* 3 &+ 2] -= Float(centerOfMass.z)
    }
  }
}

// MARK: - Angular Properties

// These functions assume the center of mass is (0, 0, 0), and both position
// and velocity have the same orientation.

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
    principalAxes: (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>)
  ) {
    // Measure the similarity between the eigenvector and the scalars.
    @_transparent
    func dot(
      _ vector: SIMD3<Double>,
      _ scalars: (MM4FloatVector, MM4FloatVector, MM4FloatVector)
    ) -> MM4FloatVector {
      Float(vector.x) * scalars.0 +
      Float(vector.y) * scalars.1 +
      Float(vector.z) * scalars.2
    }
    
    for vID in 0..<atoms.vectorCount {
      let rX = vPositions[vID &* 3 &+ 0]
      let rY = vPositions[vID &* 3 &+ 1]
      let rZ = vPositions[vID &* 3 &+ 2]
      vPositions[vID &* 3 &+ 0] = dot(principalAxes.0, (rX, rY, rZ))
      vPositions[vID &* 3 &+ 1] = dot(principalAxes.1, (rX, rY, rZ))
      vPositions[vID &* 3 &+ 2] = dot(principalAxes.2, (rX, rY, rZ))
      
      let vX = vVelocities[vID &* 3 &+ 0]
      let vY = vVelocities[vID &* 3 &+ 1]
      let vZ = vVelocities[vID &* 3 &+ 2]
      vVelocities[vID &* 3 &+ 0] = dot(principalAxes.0, (vX, vY, vZ))
      vVelocities[vID &* 3 &+ 1] = dot(principalAxes.1, (vX, vY, vZ))
      vVelocities[vID &* 3 &+ 2] = dot(principalAxes.2, (vX, vY, vZ))
    }
  }
  
  func createPositions(
    _ buffer: UnsafeMutableBufferPointer<SIMD3<Float>>
  ) {
    guard buffer.count == atoms.count else {
      fatalError("Position buffer was not the correct size.")
    }
    let baseAddress = buffer.baseAddress.unsafelyUnwrapped
    
    // TODO: Transform into the global reference frame.
    for vID in 0..<atoms.vectorCount {
      let x = vPositions[vID &* 3 &+ 0]
      let y = vPositions[vID &* 3 &+ 1]
      let z = vPositions[vID &* 3 &+ 2]
      swizzleFromVectorWidth((x, y, z), vID, baseAddress)
    }
  }
}
