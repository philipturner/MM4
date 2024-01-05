//
//  MM4State.swift
//  MM4
//
//  Created by Philip Turner on 9/10/23.
//

import OpenMM

/// A configuration for a frame of a simulation.
public struct MM4StateDescriptor {
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
  
  /// The linear velocity (in nanometers per picosecond), of each atom.
  public internal(set) var velocities: [SIMD3<Float>]?
  
  internal init() {
    
  }
}

extension MM4ForceField {
  /// Retrieve a frame of the simulation.
  ///
  /// This can be efficient than using the getters for `forces`,
  /// `positions`, `velocities`, or either of the energies in isolation.
  /// However, the API is less expressive.
  public func state(descriptor: MM4StateDescriptor) -> MM4State {
    var dataTypes: OpenMM_State.DataType = []
    if descriptor.energy {
      dataTypes = [dataTypes, .energy]
    }
    if descriptor.forces {
      dataTypes = [dataTypes, .forces]
    }
    if descriptor.positions {
      dataTypes = [dataTypes, .positions]
    }
    if descriptor.velocities {
      dataTypes = [dataTypes, .velocities]
    }
    
    let query = context.context.state(types: dataTypes)
    
    // Convert the OpenMM array to a different data type, and map from reordered
    // to original indices.
    func convertArray(_ input: OpenMM_Vec3Array) -> [SIMD3<Float>] {
      // original -> reordered -> original
      system.reorderedIndices.map {
        let index = Int($0)
        return SIMD3<Float>(input[index])
      }
    }
    
    let state = MM4State()
    if descriptor.energy {
      state.kineticEnergy = query.kineticEnergy * MM4ZJPerKJPerMol
      state.potentialEnergy = query.potentialEnergy * MM4ZJPerKJPerMol
    }
    if descriptor.forces {
      state.forces = convertArray(query.forces)
    }
    if descriptor.positions {
      state.positions = convertArray(query.positions)
    }
    if descriptor.velocities {
      state.velocities = convertArray(query.velocities)
    }
    return state
  }
}
