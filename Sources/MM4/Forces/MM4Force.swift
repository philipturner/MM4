//
//  MMForce.swift
//  MM4
//
//  Created by Philip Turner on 10/14/23.
//

import OpenMM

public enum MM4Force: CaseIterable, Hashable {
  /// Angle bending with a sextic Taylor expansion.
  case bend
  
  /// Coupling of vibrational frequencies of two angles with the same center.
  case bendBend
  
  /// Unchanging potential gradient imposed by something outside the system.
  case external
  
  /// London dispersion, overlap repulsion, and Coulomb forces.
  case nonbonded
  
  /// Morse bond stretching potential.
  case stretch
  
  /// Increase in covalent bond length as bond angle shrinks.
  case stretchBend
  
  /// Coupling of adjacent covalent bonds invoked by electronegative atoms.
  case stretchStretch
  
  /// Torsion and 1-4 nonbonded exception forces.
  case torsion
  
  /// Torsion-bend and bend-torsion-bend for electronegative atoms.
  case torsionBend
  
  /// Decrease in covalent bond length as a torsion reaches the eclipsing
  /// position.
  case torsionStretch
}

class MM4ForceGroup {
  /// The OpenMM objects containing the forces.
  var forces: [OpenMM_Force]
  
  /// Whether each force contains any particles.
  var forcesActive: [Bool]
  
  /// How many times to evaluate this force per timestep.
  var forceGroup: Int
  
  init(forces: [OpenMM_Force], forcesActive: [Bool], forceGroup: Int) {
    for force in forces {
      force.forceGroup = forceGroup
    }
    self.forces = forces
    self.forcesActive = forcesActive
    self.forceGroup = forceGroup
  }
  
  required init(system: MM4System, descriptor: MM4ForceFieldDescriptor) {
    fatalError("Not implemented.")
  }
  
  func addForces(to system: OpenMM_System) {
    for (force, forceActive) in zip(forces, forcesActive) where forceActive {
      system.addForce(force)
    }
  }
}

/// Wraps all the forces owned by a system.
class MM4Forces {
  // Force Group 1
  var electrostatic: MM4ElectrostaticForce
  var electrostaticException: MM4ElectrostaticExceptionForce
  var external: MM4ExternalForce
  var nonbonded: MM4NonbondedForce
  var nonbondedException: MM4NonbondedExceptionForce
  var torsion: MM4TorsionForce
  var torsionExtended: MM4TorsionExtendedForce
  
  // Force Group 2
  var bend: MM4BendForce
  var bendBend: MM4BendBendForce
  var bendExtended: MM4BendExtendedForce
  var stretch: MM4StretchForce
  
  init(system: MM4System, descriptor: MM4ForceFieldDescriptor) {
    // Force Group 1
    self.electrostatic = .init(system: system, descriptor: descriptor)
    self.electrostaticException = .init(system: system, descriptor: descriptor)
    self.external = .init(system: system, descriptor: descriptor)
    self.nonbonded = .init(system: system, descriptor: descriptor)
    self.nonbondedException = .init(system: system, descriptor: descriptor)
    self.torsion = .init(system: system, descriptor: descriptor)
    self.torsionExtended = .init(system: system, descriptor: descriptor)
    
    // Force Group 2
    self.bend = .init(system: system, descriptor: descriptor)
    self.bendBend = .init(system: system, descriptor: descriptor)
    self.bendExtended = .init(system: system, descriptor: descriptor)
    self.stretch = .init(system: system, descriptor: descriptor)
  }
  
  // Make this act explicit instead of performing it in the initializer. This
  // makes the code more declarative; it's easier to debug when one writes
  // functions this way.
  func addForces(to system: OpenMM_System) {
    // Force Group 0
    external.addForces(to: system)
    
    // Force Group 1
    electrostatic.addForces(to: system)
    electrostaticException.addForces(to: system)
    nonbonded.addForces(to: system)
    nonbondedException.addForces(to: system)
    torsion.addForces(to: system)
    torsionExtended.addForces(to: system)
    
    // Force Group 2
    bend.addForces(to: system)
    bendBend.addForces(to: system)
    bendExtended.addForces(to: system)
    stretch.addForces(to: system)
  }
}
