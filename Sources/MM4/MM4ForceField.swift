//
//  MM4ForceField.swift
//  MolecularRenderer
//
//  Created by Philip Turner on 9/10/23.
//

import Foundation

/// A configuration for a force field simulator.
public class MM4ForceFieldDescriptor {
  /// Optional. The force (in piconewtons) exerted on each atom.
  public var externalForces: [SIMD3<Float>]?
  
  /// Required. The set of parameters defining the forcefield.
  public var parameters: MM4Parameters?
  
  /// Required. The position (in nanometers) of each atom's nucleus.
  public var positions: [SIMD3<Float>] = []
  
  /// Optional. Whether each atom's absolute position should never change.
  ///
  /// This is implemented by setting the particle's mass and velocity to zero.
  ///
  /// > Warning: Stationary atoms may cause energy measurements to be
  /// nonsensical. This needs to be investigated further.
  public var stationaryAtoms: [Bool]?
  
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
  
  // TODO: Allow these properties to be changed without making a new object.
  // Provide an elegant API for this, but also a batched API for more efficient
  // transactions. The efficient API should use an `MM4ForceFieldDescriptor`,
  // with the requirement that parameters are `nil`.
  //
  // externalForces
  //   - setParticleParameters, updateParametersInContext
  //   - setIntegrationForceGroups, swapping integrators
  // positions
  //   - setState
  // stationaryAtoms
  //   - setParticleParameters, updateParametersInContext
  //   - setIntegrationForceGroups, swapping integrators
  // velocities
  //   - setState
  public func update(descriptor: MM4ForceFieldDescriptor) {
    
  }
  
  /// Simulate the system's evolution for the specified time interval.
  ///
  /// - Parameter time: The time interval, in picoseconds.
  /// - Parameter maximumTimeStep: The largest time step that may be taken
  ///   during the simulation, in picoseconds.
  ///
  /// This is an NVE simulation, using 32-bit single precision and a massive
  /// timestep. Energy is not conserved to the precision of kT. It fluctuates
  /// wildly, but shouldn't systematically drift upward or downward (unless the
  /// timestep exceeds C-H stretching resonance frequency). The default
  /// timestep (4 fs) is so massive that integration error dwarfs numerical
  /// error. Even if the forcefield used 64-bit double or mixed precision, that
  /// would not improve its ability to conserve energy.
  ///
  /// > Note: Accuracy (whether an archer hits the right target) is orthogonal
  ///   to precision (whether the arrows hit the same place every time).
  ///   Conservation of energy only affects the number of significant figures of
  ///   measured energy. It does not make system dynamics more accurate, except
  ///   in the case where energy systematically drifts upward (explodes).
  ///
  /// The integrator uses a multiple time-stepping (MTS) scheme. Cheaper bonded
  /// forces, such as bond-stretch and bond-bend, are only stable at ~2 fs
  /// without constraints. Expensive forces like torsions, nonbonded, and
  /// electrostatic can execute at double the timestep. The value you enter for
  /// `maximumTimeStep` specifies the execution rate of larger forces. Always
  /// assume the C-H stretching forces execute at half the specified timestep.
  /// For example, in the note below, bond stretching forces don't actually
  /// execute at the quoted '2 fs'.
  ///
  /// > Note: To maximize the simulation step, hydrogen mass repetitioning (HMR)
  ///   is enabled by default. This method makes hydrogens heaver and makes
  ///   non-hydrogen atoms lighter, to decrease the C-H or Si-H stretching
  ///   frequency. You can, and should, disable HMR for energy-conserving
  ///   simulations where timestep falls below 2 fs. This can be done by setting
  ///   `hydrogenMassRepartitioning` in the `MM4ParametersDescriptor` to `0`.
  ///   Note that changing the repartitioning will shift the center of mass,
  ///   which is where bulk angular momentum is applied.
  ///
  /// Single precision is an implementation choice to make the simulator more
  /// compatible with different GPU architectures. Most vendors have no or few
  /// FP64 units (unlike CPU, where FP64 is widely supported). Mixed precision
  /// can sometimes run the bulk of computations in FP32, with only the energy
  /// computations in FP64. However, even mixed precision has implementation
  /// issues, being difficult to run on GPUs with zero FP64 units.
  ///
  /// The best method to maximize energy conservation is a timestep of 1 fs.
  /// This requires 4x as much computation for the same trajectory length,
  /// compared to the default (4 fs). It executes expensive forces at 1 fs and
  /// cheap forces at 0.5 fs. Integration error scales O(h^2), so this should be
  /// around 16x more precise. Any smaller timestep, and rounding error from
  /// FP32 will overtake integration error.
  public func simulate(
    time: Double,
    maximumTimeStep: Double = 0.100 / 23 + 1e-8
  ) {
    // If the time doesn't divide evenly into 100 fs, use a temporary
    // integrator for the remainder, potentially with a slightly scaled-down
    // timestep.
    //
    // Dynamically switch between the integrators below, or lazily initialized
    // variants with external forces. There should be multiple OpenMM contexts
    // surrounding the same system.
    // - minimization integrator
    // - single-step custom MTS integrator
    //   - 2, 3, 4-step variants created on demand to encompass a remainder
    // - 5-step custom MTS integrator
    // - 23-step custom MTS integrator
  }
  
  /// Minimize the system's energy using the L-BFGS algorithm.
  ///
  /// This is one of few such algorithms with O(n) computational complexity. It
  /// is a limited-memory version of BFGS, and O(n^2) algorithm. BFGS, in turn,
  /// is an improvement on O(n^3) methods such as Newton's method and the
  /// conjugate gradient method.
  ///
  /// - Parameter tolerance: Accepted uncertainty in potential energy,
  ///   in zeptojoules.
  /// - Parameter maxIterations: Maximum number of force evaluations permitted
  ///   during the minimization. The default value, 0, puts no restrictions on
  ///   the number of evaluations.
  public func minimize(
    tolerance: Double = 10.0 * MM4_ZJPerKJPerMol,
    maxIterations: Int = 0
  ) {
    // Use a different integrator that doesn't incorporate the fusion of
    // multiple timesteps. Also, zero out the system's bulk velocity during the
    // minimization.
  }
  
  /// Create random thermal velocities, while conserving the total momentum of
  /// each rigid body.
  ///
  /// - Parameter temperature: The temperature to randomize thermal velocites
  ///   at, in kelvin.
  /// - Parameter rigidBodies: Atom indices for each rigid body.
  ///
  /// Thermalizing is recommended for any simulation that replicates macroscale
  /// conditions. The default is 298.15 K, but other useful temperatures include
  /// liquid nitrogen (77.00 K) and liquid helium (4.15 K).
  ///
  /// Rigid bodies should have atoms laid out contiguously in memory, in Morton
  /// order. This format ensures spatial locality, which increases performance
  /// of nonbonded forces. Therefore, rigid bodies must be entered as contiguous
  /// ranges of the atom list.
  ///
  /// The set of rigid bodies must cover every atom in the system. No two ranges
  /// may overlap the same atom. If the array of rigid bodies is empty, it
  /// defaults to a range encompassing the entire system. This ensures the
  /// closed system's net momentum stays conserved.
  public func thermalize(
    temperature: Double = 298.15,
    rigidBodies: [Range<Int>] = []
  ) {
    
  }
}
