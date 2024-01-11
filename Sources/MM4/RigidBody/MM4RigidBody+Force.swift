//
//  MM4RigidBody+Force.swift
//
//
//  Created by Philip Turner on 1/7/24.
//

extension MM4RigidBody {
  /// The force (in piconewtons) on each atom.
  ///
  /// Force is a function of position and the surrounding environment. When the
  /// rigid body changes position, the forces change as well. Therefore, the
  /// forces are invalidated (turned into `nil`). To restore the forces,
  /// you must assign them to the `forces` property again.
  public var forces: [SIMD3<Float>]? {
    _read {
      yield storage.forces
    }
    _modify {
      ensureUniquelyReferenced()
      storage.netForce = nil
      storage.netTorque = nil
      yield &storage.forces
    }
  }
  
  /// The derivative of linear momentum with respect to time, in piconewtons.
  public var netForce: SIMD3<Double>? {
    storage.ensureForceAndTorqueCached()
    return storage.netForce
  }
  
  /// The derivative of angular momentum with respect to time, in piconewtons
  /// per nanometer.
  public var netTorque: SIMD3<Double>? {
    storage.ensureForceAndTorqueCached()
    return storage.netTorque
  }
}

extension MM4RigidBodyStorage {
  func setNetForceAndTorque(
    _ buffer: UnsafeBufferPointer<SIMD3<Float>>
  ) {
    guard buffer.count == atoms.count else {
      fatalError("Force buffer was not the correct size.")
    }
    let baseAddress = buffer.baseAddress.unsafelyUnwrapped
    
    let Σ = (
      SIMD3<Float>(principalAxes.0),
      SIMD3<Float>(principalAxes.1),
      SIMD3<Float>(principalAxes.2))
    
    var netForce: SIMD3<Double> = .zero
    var netTorque: SIMD3<Double> = .zero
    withSegmentedLoop(chunk: 256) {
      var vForceX: MM4FloatVector = .zero
      var vForceY: MM4FloatVector = .zero
      var vForceZ: MM4FloatVector = .zero
      var vTorqueX: MM4FloatVector = .zero
      var vTorqueY: MM4FloatVector = .zero
      var vTorqueZ: MM4FloatVector = .zero
      for vID in $0 {
        var (fX, fY, fZ) = swizzleToVectorWidth(vID, baseAddress)
        vForceX += fX
        vForceY += fY
        vForceZ += fZ
        let rX = vPositions[vID &* 3 &+ 0]
        let rY = vPositions[vID &* 3 &+ 1]
        let rZ = vPositions[vID &* 3 &+ 2]
        
        // r x (ΣT * f)
        let f = (fX, fY, fZ)
        fX = dot(vector: Σ.0, scalars: f)
        fY = dot(vector: Σ.1, scalars: f)
        fZ = dot(vector: Σ.2, scalars: f)
        vTorqueX += rY * fZ - rZ * fY
        vTorqueY += rZ * fX - rX * fZ
        vTorqueZ += rX * fY - rY * fX
      }
      netForce.x += MM4DoubleVector(vForceX).sum()
      netForce.y += MM4DoubleVector(vForceY).sum()
      netForce.z += MM4DoubleVector(vForceZ).sum()
      netTorque.x += MM4DoubleVector(vTorqueX).sum()
      netTorque.y += MM4DoubleVector(vTorqueY).sum()
      netTorque.z += MM4DoubleVector(vTorqueZ).sum()
    }
    self.netForce = netForce
    self.netTorque = netTorque
  }
}
