//
//  MM4State.swift
//
//
//  Created by Philip Turner on 9/10/23.
//

/// A configuration for a frame of a simulation.
public class MM4StateDescriptor {
  /// Required. Whether to report the system's total kinetic and potential
  /// energy.
  ///
  /// The default is `false`.
  public var energy: Bool = false
  
  /// Required. Whether to report the force exerted on each atom.
  ///
  /// The default is `false`.
  public var forces: Bool = false
  
  /// Required. Whether to report each atom's position.
  ///
  /// The default is `false`.
  public var positions: Bool = false
  
  /// Required. Whether to report each atom's velocity.
  ///
  /// The default is `false`.
  public var velocities: Bool = false
  
  public init() {
    
  }
}

/// A frame of a simulation.
public class MM4State {
  /// The net varying force (in piconewtons) exerted on each atom.
  ///
  /// This is converted from kJ/mol/nm to piconewtons.
  public internal(set) var forces: [SIMD3<Float>]?
  
  /// The system's total kinetic energy, in zeptojoules.
  public internal(set) var kineticEnergy: Double?
  
  /// The position (in nanometers) of each atom's nucleus.
  public internal(set) var positions: [SIMD3<Float>]?
  
  /// The system's total potential energy, in zeptojoules.
  public internal(set) var potentialEnergy: Double?
  
  /// The bulk + thermal velocity (in nanometers per picosecond), of each atom.
  public internal(set) var velocities: [SIMD3<Float>]?
  
  internal init() {
    
  }
}

extension MM4ForceField {
  /// Retrieve a frame of the simulation.
  ///
  /// This should be more efficient than using the getters for `forces`,
  /// `positions`, `velocities`, or either of the energies in isolation.
  /// However, the API is less expressive.
  public func state(descriptor: MM4StateDescriptor) -> MM4State {
    if descriptor.energy {
      // Add the energy flag to the OpenMM state data type.
    }
    if descriptor.forces {
      // Add the forces flag to the OpenMM state data type.
    }
    if descriptor.positions {
      // Add the positions flag to the OpenMM state data type.
    }
    if descriptor.velocities {
      // Add the velocities flag to the OpenMM state data type.
    }
   
    let state = MM4State()
    if descriptor.energy {
      // Set the kinetic and potential energy.
    }
    if descriptor.forces {
      // Set the forces.
    }
    if descriptor.positions {
      // Set the positions.
    }
    if descriptor.velocities {
      // Set the velocities.
    }
    return state
  }
}
