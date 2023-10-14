//
//  MM4ForceField.swift
//  MolecularRenderer
//
//  Created by Philip Turner on 9/10/23.
//

import Foundation

/// A configuration for a force field simulator.
public class MM4ForceFieldDescriptor {
  // TODO: Allow these properties to be changed without making a new object.
  //
  // externalForces
  //   - setParticleParameters, updateParametersInContext
  //   - setIntegrationForceGroups
  // maximumTimeStep
  //   - setStepSize
  // positions
  //   - setState
  // stationaryAtoms
  //   - setParticleParameters, updateParametersInContext
  //   - setIntegrationForceGroups
  // velocities
  //   - setState
  
  /// Optional. The force (in piconewtons) exerted on each atom.
  public var externalForces: [SIMD3<Float>]?
  
  /// Required. The largest time step (in picoseconds) that may be taken during
  /// a simulation.
  ///
  /// The default is approximately 4 femtoseconds.
  public var maximumTimeStep: Double = 0.100 / 23 + 1e-8
  
  /// Required. The set of parameters defining the forcefield.
  public var parameters: MM4Parameters?
  
  /// Required. The position (in nanometers) of each atom's nucleus.
  public var positions: [SIMD3<Float>] = []
  
  /// Required. The one-to-one mapping of atom indices to rigid bodies.
  public var rigidBodies: [[UInt32]] = []
  
  /// Optional. Whether each atom's absolute position should never change.
  ///
  /// This is implemented by setting the particle's mass and velocity to zero.
  public var stationaryAtoms: [Bool]?
  
  /// Required. The temperature (in Kelvin) to initialize thermal velocities at.
  ///
  /// The default is 25 degrees Celsius.
  public var temperature: Double = 298.15
  
  /// Optional. The velocity (in nanometers per picosecond), of each atom at the
  /// start of the simulation.
  ///
  /// These are added to thermal velocities in a way that conserves each rigid
  /// body's overall momentum.
  public var velocities: [SIMD3<Float>]?
  
  // Idea: Set of modes for the forcefield
  // - "Default" - best tradeoff between accuracy and performance.
  // - "Reduced accuracy" - cheaper alternative forcefield that omits torsions
  //   in bulk H/C/Si. Not yet measured whether it's sufficient for replacing
  //   "default", what the artifacts are.
  //   - Measure and rank the stiffness of each force, and their contributions
  //     to a set of specific quantitative properties. No matter how cheap a
  //     force is, there should be a maximally cheap option for situations where
  //     simulation speed vastly outweighs accuracy.
  // - "Energy conserving" - throws an error if HMR is nonzero, uses a much
  //   smaller timestep to minimize energy drift, no mixed precision for now.
  // - "Energy minimizing" - disregards simulation speed, creates checkpoints to
  //   restart the simulation if energy explodes, may use modified forces.
  
  public init() {
    
  }
}

/// A force field simulator.
public class MM4ForceField {
  /// Create a simulator using the specified configuration.
  public init(descriptor: MM4ForceFieldDescriptor) {
    MM4Plugins.global.load()
    
    // Eventually, separate the atoms into two groups of "small" vs "large"
    // atoms, creating different zones of internally contiguous tiles.
  }
  
  /// Simulate the system's evolution for the specified time interval (in
  /// picoseconds).
  public func simulate(time: Double) {
    // If the time doesn't divide evenly into 100 fs, compile a temporary
    // integrator that executes the remainder, potentially with a slightly
    // scaled-down timestep.
  }
}
