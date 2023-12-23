//
//  MM4RigidBody+Energy.swift
//
//
//  Created by Philip Turner on 12/22/23.
//

/// The energy of a rigid body.
public struct MM4RigidBodyEnergy {
  /// The kinetic energy of the rigid body's non-anchor atoms.
  public var kinetic: MM4RigidBodyKineticEnergy = .init()
  
  init() { }
}

/// The kinetic energy of a rigid body's non-anchor atoms.
public struct MM4RigidBodyKineticEnergy {
  weak var storage: MM4RigidBodyStorage!
  
  init() { }
  
  /// Kinetic energy contribution from organized mechanical energy, present in
  /// the bulk angular velocity.
  ///
  /// If there are more than one anchors, the free kinetic energy from
  /// angular momentum should ideally be zero. However, the system can
  /// experience angular displacements from the orientation strictly enforced by
  /// anchors.
  public var angular: Double {
    storage.ensureKineticEnergyCached()
    return storage.angularKineticEnergy!
  }
  
  /// Kinetic energy contribution from organized mechanical energy, present in
  /// the bulk linear velocity.
  ///
  /// If there are any anchors, the free kinetic energy from the linear
  /// velocity is treated specially. It equals the total mass of non-anchor
  /// atoms, combined with the anchors' linear velocity.
  public var linear: Double {
    storage.ensureKineticEnergyCached()
    return storage.linearKineticEnergy!
  }
  
  /// Kinetic energy contribution from disorganized thermal energy.
  ///
  /// Contributions from anchors are omitted.
  public var thermal: Double {
    _read {
      storage.ensureKineticEnergyCached()
      yield storage.thermalKineticEnergy!
    }
    _modify {
      storage.ensureKineticEnergyCached()
      var energy = storage.thermalKineticEnergy!
      yield &energy
      storage.setThermalKineticEnergy(energy)
    }
  }
}

extension MM4RigidBody {
  /// The rigid body's energy.
  public var energy: MM4RigidBodyEnergy {
    _read {
      yield _energy
    }
    _modify {
      ensureUniquelyReferenced()
      yield &_energy
    }
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
    return nonAnchorMass * (v * v).sum() / 2
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
    withMasses(nonAnchorMasses) { vMasses in
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
