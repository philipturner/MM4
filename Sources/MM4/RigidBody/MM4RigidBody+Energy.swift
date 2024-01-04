//
//  MM4RigidBody+Energy.swift
//
//
//  Created by Philip Turner on 12/22/23.
//

/// The energy of a rigid body.
public struct MM4RigidBodyEnergy {
  /// The kinetic energy of the rigid body.
  public var kinetic: MM4RigidBodyKineticEnergy = .init()
  
  init() { }
}

/// The kinetic energy of a rigid body.
public struct MM4RigidBodyKineticEnergy {
  weak var storage: MM4RigidBodyStorage!
  
  init() { }
  
  /// Kinetic energy contribution from organized mechanical energy, present in
  /// the bulk angular velocity.
  public var angular: Double {
    storage.ensureKineticEnergyCached()
    return storage.angularKineticEnergy!
  }
  
  /// Kinetic energy contribution from organized mechanical energy, present in
  /// the bulk linear velocity.
  public var linear: Double {
    storage.ensureKineticEnergyCached()
    return storage.linearKineticEnergy!
  }
  
  /// Kinetic energy contribution from disorganized thermal energy.
  public var thermal: Double {
    storage.ensureKineticEnergyCached()
    return storage.thermalKineticEnergy!
  }
}

extension MM4RigidBody {
  /// The rigid body's energy.
  public var energy: MM4RigidBodyEnergy {
    return _energy
  }
}

extension MM4RigidBodyStorage {
  // Linear kinetic energy from non-anchor atoms.
  func createLinearKineticEnergy() -> Double {
    ensureLinearVelocityCached()
    guard let linearVelocity else {
      fatalError("This should never happen.")
    }
    
    let v = SIMD3<Double>(linearVelocity)
    return mass * (v * v).sum() / 2
  }
  
  // Angular kinetic energy about an angular mass defined by non-anchor atoms.
  func createAngularKineticEnergy() -> Double {
    ensureMomentOfInertiaCached()
    ensureAngularVelocityCached()
    guard let momentOfInertia,
          let angularVelocity else {
      fatalError("This should never happen.")
    }
    
    let I = momentOfInertia
    let w = SIMD3<Double>(quaternionToVector(angularVelocity))
    let velocityX = I.columns.0 * w.x
    let velocityY = I.columns.1 * w.y
    let velocityZ = I.columns.2 * w.z
    let Iw = velocityX + velocityY + velocityZ
    return (w * Iw).sum() / 2
  }
  
  // Total translational kinetic energy from non-anchor atoms.
  func createTotalKineticEnergy() -> Double {
    // This function should return the same value as separating the bulk and
    // thermal velocities, computing the energies separately, then summing.
    
    // Anchors should never be included in kinetic energy, either free or
    // thermal. For thermal, the kinetic energy from anchor velocities will
    // contribute zero to the total energy. This should let velocity rescaling
    // work properly.
    var kinetic: Double = .zero
    withMasses { vMasses in
      withSegmentedLoop(chunk: 256) {
        var vKineticX: MM4FloatVector = .zero
        var vKineticY: MM4FloatVector = .zero
        var vKineticZ: MM4FloatVector = .zero
        for vID in $0 {
          let x = vVelocities[vID &* 3 &+ 0]
          let y = vVelocities[vID &* 3 &+ 1]
          let z = vVelocities[vID &* 3 &+ 2]
          let mass = vMasses[vID]
          vKineticX.addProduct(mass, x * x)
          vKineticY.addProduct(mass, y * y)
          vKineticZ.addProduct(mass, z * z)
        }
        kinetic += MM4DoubleVector(vKineticX).sum()
        kinetic += MM4DoubleVector(vKineticY).sum()
        kinetic += MM4DoubleVector(vKineticZ).sum()
      }
    }
    return kinetic / 2
  }
}
