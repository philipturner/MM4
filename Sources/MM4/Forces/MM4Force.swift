//
//  MMForce.swift
//  MM4
//
//  Created by Philip Turner on 10/14/23.
//

import OpenMM

/// Options for customizing the forces that act on a set of atoms.
public struct MM4ForceOptions: OptionSet {
  /// The 32-bit code describing the force options.
  public let rawValue: UInt32
  
  /// Create a set of force options using the specified code.
  @_transparent
  public init(rawValue: UInt32) {
    self.rawValue = rawValue
  }
  
  /// Angle bending with a sextic Taylor expansion.
  public static let bend = MM4ForceOptions(rawValue: 1 << 0)
  
  /// London dispersion, overlap repulsion, and Coulomb forces.
  public static let nonbonded = MM4ForceOptions(rawValue: 1 << 1)
  
  /// Morse bond stretching potential.
  public static let stretch = MM4ForceOptions(rawValue: 1 << 2)
  
  /// Increase in covalent bond length as bond angle shrinks.
  public static let stretchBend = MM4ForceOptions(rawValue: 1 << 3)
}

class MM4Force {
  /// The OpenMM objects containing the forces.
  var forces: [OpenMM_Force]
  
  /// How many times to evaluate this force per timestep.
  var forceGroup: Int
  
  init(forces: [OpenMM_Force?], forceGroup: Int) {
    let activeForces = forces.compactMap { $0 }
    for force in activeForces {
      force.forceGroup = forceGroup
    }
    self.forces = activeForces
    self.forceGroup = forceGroup
  }
  
  required init(system: MM4System, descriptor: MM4ForceFieldDescriptor) {
    fatalError("Not implemented.")
  }
  
  func addForces(to system: OpenMM_System) {
    for force in forces {
      force.transfer()
      system.addForce(force)
    }
  }
}

/// Wraps all the forces owned by a system.
class MM4Forces {
  // Force Group 1
  var electrostatic: MM4ElectrostaticForce
  var external: MM4ExternalForce
  var nonbonded: MM4NonbondedForce
  
  // Force Group 2
  var bend: MM4BendForce
  var stretch: MM4StretchForce
  
  init(system: MM4System, descriptor: MM4ForceFieldDescriptor) {
    // Force Group 1
    self.electrostatic = .init(system: system, descriptor: descriptor)
    self.external = .init(system: system, descriptor: descriptor)
    self.nonbonded = .init(system: system, descriptor: descriptor)
    
    // Force Group 2
    self.bend = .init(system: system, descriptor: descriptor)
    self.stretch = .init(system: system, descriptor: descriptor)
  }
  
  func addForces(to system: OpenMM_System) {
    // Force Group 1
    electrostatic.addForces(to: system)
    external.addForces(to: system)
    nonbonded.addForces(to: system)
    
    // Force Group 2
    bend.addForces(to: system)
    stretch.addForces(to: system)
  }
}
