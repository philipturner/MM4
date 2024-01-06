//
//  MM4ForceField+Energy.swift
//  MM4
//
//  Created by Philip Turner on 12/22/23.
//

import OpenMM

/// A data structure wrapping a system's energy.
///
/// > NOTE: This documentation page is still a draft. It may be inconsistent or
///   difficult to understand.
///
/// We always report energy, time, and unit conversion constants in double
/// precision. Potential energy can sometimes have a massive absolute value
/// that dwarfs the magnitude of relative values.
/// Everything else, including rigid body bulk properties susceptible to the
/// same degree of rounding error, are reported in single precision. Force
/// field parameters and similar quantities (heat capacity) are in single
/// precision because we don't know their value with a lot of certainty.
///
/// However, some internal summation algorithms temporarily employ mixed
/// precision. Both the CPU code for `MM4RigidBody` and GPU code for the OpenMM
/// Metal plugin use block summation. The majority of the compute cost happens
/// on FP32, either to improve execution speed and reduce register pressure
/// (CPU) or provide compatibility with a specific hardware vendor (GPU). The
/// rounding error is capped to whatever error occurs within a block
/// (~1000 elements). Outside of the block, the number is casted to
/// FP64 and summed with no additional rounding error. This prevents unexpected
/// issues popping up with large systems that exceed 1000 elements. The mantissa
/// of the quantity is still ~32 bits; we just prevent it from degrading to
/// something like ~16 bits.
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
  
  /// Minimize the system's potential energy, removing thermal potential energy.
  ///
  /// OpenMM uses the L-BFGS algorithm for gradient descent. It is an O(n)
  /// version of BFGS, an O(n^2) algorithm. BFGS improves upon O(n^3) methods
  /// such as Newton's method.
  ///
  /// - Parameter tolerance: Accepted uncertainty in potential energy,
  ///   in zeptojoules. The default value is 10. This contrasts with the default
  ///   value for most OpenMM simulations, which is 16.6 zJ (10.0 kJ/mol).
  /// - Parameter maxIterations: Maximum number of force evaluations permitted
  ///   during the minimization. The default value, 0, puts no restrictions on
  ///   the number of evaluations.
  public func minimize(
    tolerance: Double = 10.0,
    maxIterations: Int = 0
  ) {
    if updateRecord.active() {
      flushUpdateRecord()
    }
    invalidatePositionsAndVelocities()
    invalidateForcesAndEnergy()
    
    // Switch to an integrator that always reports the correct velocity.
    var integratorDescriptor = MM4IntegratorDescriptor()
    integratorDescriptor.start = true
    integratorDescriptor.end = true
    context.currentIntegrator = integratorDescriptor
    
    // Run the energy minimization.
    //
    // The reporter doesn't do anything. You have to create a C++ class, which
    // is not possible through the OpenMM C API.
    OpenMM_LocalEnergyMinimizer.minimize(
      context: context.context,
      tolerance: tolerance,
      maxIterations: maxIterations,
      reporter: nil)
  }
}
