//
//  MM4RigidBody+Properties.swift
//
//
//  Created by Philip Turner on 11/22/23.
//

import Numerics

// MARK: - Energy

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

// MARK: - Other Properties

extension MM4RigidBody {
  /// Indices of atoms that should be treated as having infinite mass.
  public var anchors: Set<UInt32> {
    // _modify not supported b/c it requires very complex caching logic.
    // Workaround: import a new rigid body initialized with different anchors.
    _read { yield storage.anchors }
  }
  
  /// The number of protons in each atom's nucleus.
  public var atomicNumbers: [UInt8] {
    parameters.atoms.atomicNumbers
  }
  
  /// Indices of atoms where external force is applied.
  public var handles: Set<UInt32> {
    // _modify not supported b/c it requires very complex caching logic.
    // Workaround: import a new rigid body initialized with different targets.
    _read { yield storage.handles }
  }
  
  /// Pairs of atom indices representing sigma bonds.
  public var bonds: [SIMD2<UInt32>] {
    parameters.bonds.indices
  }
  
  /// The total mass (in yoctograms) of all atoms, excluding anchors.
  public var mass: Double {
    storage.nonAnchorMass
  }
  
  /// The mass (in yoctograms) of each atom.
  ///
  /// This is different than the masses in `MM4Parameters`. `MM4RigidBody`
  /// zeroes out the mass for each anchor, while `MM4Parameters` does not.
  public var masses: [Float] {
    storage.nonAnchorMasses
  }
}
