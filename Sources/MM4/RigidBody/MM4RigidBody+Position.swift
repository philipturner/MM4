//
//  MM4RigidBody+Position.swift
//  MM4
//
//  Created by Philip Turner on 11/20/23.
//

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
  public var centerOfMass: SIMD3<Double> {
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
  // Use this function to import positions into the force field.
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
  
  func createVectorizedMasses(_ masses: [Float]) {
    vMasses = Array(repeating: .zero, count: atoms.vectorCount)
    guard mass == 0, atoms.nonAnchorCount == 0 else {
      fatalError("Vectorized masses were already initialized.")
    }
    
    masses.withContiguousStorageIfAvailable {
      let baseAddress = $0.baseAddress.unsafelyUnwrapped
      let vAddress = UnsafePointer<MM4FloatVector>(OpaquePointer(baseAddress))
      
      withSegmentedLoop(chunk: 256) {
        var vMassAccumulator: MM4FloatVector = .zero
        var vNonAnchorCount: MM4FloatVector = .zero
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
          
          vMassAccumulator += vMass
          vMass.replace(with: .one, where: vMass .!= .zero)
          vNonAnchorCount += vMass
        }
        mass += MM4DoubleVector(vMassAccumulator).sum()
        atoms.nonAnchorCount += Int(vNonAnchorCount.sum())
      }
    }
  }
  
  func createCenterOfMass() {
    precondition(
      centerOfMass == .zero, "Center of mass was already initialized.")
    guard atoms.nonAnchorCount > 0 else {
      return
    }
    
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
    
    precondition(mass > 0, "Mass was not yet initialized.")
    centerOfMass /= mass
  }
}
