//
//  MM4ForceField+Properties.swift
//  MM4
//
//  Created by Philip Turner on 10/14/23.
//

import OpenMM

extension MM4ForceField {
  /// The constant force (in piconewtons) exerted on each atom.
  public var externalForces: [SIMD3<Float>] {
    _read {
      yield _externalForces
    }
    _modify {
      updateRecord.externalForces = true
      yield &_externalForces
    }
  }
  
  /// The net varying force (in piconewtons) exerted on each atom.
  public var forces: [SIMD3<Float>] {
    _read {
      ensureForcesAndEnergyCached()
      yield cachedState.forces!
    }
  }
  
  /// The position (in nanometers) of each atom's nucleus.
  public var positions: [SIMD3<Float>] {
    _read {
      ensurePositionsAndVelocitiesCached()
      yield cachedState.positions!
    }
    _modify {
      ensurePositionsAndVelocitiesCached()
      updateRecord.positions = true
      updateRecord.velocities = true
      yield &cachedState.positions!
    }
  }
  
  /// The linear velocity (in nanometers per picosecond) of each atom.
  ///
  /// When thermalizing, the linear and angular momentum over every rigid body
  /// is conserved. Then, the thermal velocities are reinitialized. If you want
  /// more complex motion within the rigid body, fetch the thermalized
  /// velocities. Add the desired bulk velocity component to them, them set the
  /// new velocity values.
  public var velocities: [SIMD3<Float>] {
    _read {
      ensurePositionsAndVelocitiesCached()
      yield cachedState.velocities!
    }
    _modify {
      ensurePositionsAndVelocitiesCached()
      updateRecord.positions = true
      updateRecord.velocities = true
      yield &cachedState.velocities!
    }
  }
}
