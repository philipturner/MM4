//
//  MM4ForceField+RigidBody.swift
//  MM4
//
//  Created by Philip Turner on 11/20/23.
//

/// A data structure wrapping a system's energy.
///
/// We always report energy, time, and unit conversion constants in double
/// precision. Potential energy can sometimes have a massive absolute value
/// that dwarfs the magnitude of relative values. Rigid body bulk properties
/// (only accessed by the CPU) are stored in FP64 to minimize drift from
/// accumulated rounding error. Everything else is in FP32. Force
/// field parameters are in single precision because we don't know their values
/// with a lot of certainty. In addition, cached parameters may consume a large
/// amount of memory.
///
/// Some internal summation algorithms temporarily employ mixed
/// precision. Both the CPU code for `MM4RigidBody` and GPU code for the OpenMM
/// Metal plugin use block summation. The majority of the compute cost happens
/// on FP32, either to harness more GFLOPS and reduce register pressure
/// (CPU) or provide compatibility with a specific hardware vendor (GPU). The
/// rounding error is capped to whatever error occurs within a block
/// (~1000 elements). Outside of the block, the number is casted to
/// FP64 and summed with no additional rounding error. This prevents unexpected
/// issues popping up with large systems that exceed 1000 atoms. The mantissa
/// of the quantity is still ~24 bits; we just prevent it from degrading to
/// ~16 bits.
public struct MM4ForceFieldEnergy {
  var forceField: MM4ForceField
  
  init(forceField: MM4ForceField) {
    self.forceField = forceField
  }
  
  /// The system's total kinetic energy, in zeptojoules.
  public var kinetic: Double {
    forceField.ensureForcesAndEnergyCached()
    return forceField.cachedState.kineticEnergy!
  }
  
  /// The system's total potential energy, in zeptojoules.
  public var potential: Double {
    forceField.ensureForcesAndEnergyCached()
    return forceField.cachedState.potentialEnergy!
  }
}

extension MM4ForceField {
  /// The system's energy.
  public var energy: MM4ForceFieldEnergy {
    _energy
  }
  
  /// The constant force (in piconewtons) exerted on each atom.
  ///
  /// The default value is zero for every atom.
  ///
  /// > WARNING: There is no getter for this property. You must assign an entire
  ///   array to it at once. Otherwise, there will be a runtime crash.
  public var externalForces: [SIMD3<Float>] {
    get {
      fatalError("You cannot retrieve the values of external forces.")
    }
    set {
      guard newValue.count == system.parameters.atoms.count else {
        fatalError("Number of external forces does not match atom count.")
      }
      
      let force = system.forces.external
      force.updateForces(newValue, system: system)
      force.updateParametersInContext(context)
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
  
  /// The velocity (in nanometers per picosecond) of each atom.
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
