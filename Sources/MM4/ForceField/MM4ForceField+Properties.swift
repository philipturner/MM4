//
//  MM4ForceField+Properties.swift
//
//
//  Created by Philip Turner on 10/14/23.
//

import OpenMM

// Ergonomic APIs for accessing a force field's state. These delegate to
// a single-property instances of the batched functions, 'update' or 'state'.
// They are not intended to be used as a primary interface for interacting with
// the system (they have high overhead). Rather, the calling application's IR
// should handle most of the data processing. One should only transfer data
// in/out of OpenMM to set up the simulation.

// Idea for adding torques:
// - Can't constrain "anchors" to have constant angular velocity, only constant
//   linear velocity.
// - A simulation can achieve near-constant angular velocity with a large
//   flywheel (which will be provided in the hardware catalog, and documented in
//   MM4).
// - Torques around a computed position:
//   - LocalCoordinatesSite - center of mass of a rigid body, exerts a
//     counterforce on the rigid body as a whole
//   - relative to a set of collinear anchor particles, which must all have the
//     same velocity
//     - check the collinearity constraint every time the user changes particle
//       velocities or modifies the anchors
// - Requires a new type (`MM4Torque`) that wraps a Swift `Quaternion` and an
//   enumeration with an associated value, which specifies the type of origin.

// MARK: - Batched Functions

extension MM4ForceField {
  /// The net varying force (in piconewtons) exerted on each atom.
  ///
  /// > Note: This is a more ergonomic API, but less efficient than the batched
  /// function <doc:MM4ForceField/state(descriptor:)>.
  public var forces: [SIMD3<Float>] {
    get {
      // No need to convert between original and reordered indices.
      var descriptor = MM4StateDescriptor()
      descriptor.forces = true
      return state(descriptor: descriptor).forces!
    }
  }
  
  /// The position (in nanometers) of each atom's nucleus.
  ///
  /// > Note: This is a more ergonomic API, but less efficient than the batched
  /// function <doc:MM4ForceField/state(descriptor:)>.
  public var positions: [SIMD3<Float>] {
    get {
      // No need to convert between original and reordered indices.
      var descriptor = MM4StateDescriptor()
      descriptor.positions = true
      return state(descriptor: descriptor).positions!
    }
    set {
      // reordered -> original -> reordered
      let array = OpenMM_Vec3Array(size: system.parameters.atoms.count)
      for (reordered, original) in system.originalIndices.enumerated() {
        let position = newValue[Int(original)]
        array[reordered] = SIMD3<Double>(position)
      }
      context.context.positions = array
    }
  }
  
  /// The bulk + thermal velocity (in nanometers per picosecond) of each atom.
  ///
  /// > Note: This is a more ergonomic API, but less efficient than the batched
  /// function <doc:MM4ForceField/state(descriptor:)>.
  ///
  /// When thermalizing, the linear and angular momentum over every rigid body
  /// is conserved. Then, the thermal velocities are reinitialized. If you want
  /// more complex motion within the rigid body, fetch the thermalized
  /// velocities. Add the desired bulk velocity component to them, them set the
  /// new velocity values.
  public var velocities: [SIMD3<Float>] {
    get {
      // No need to convert between original and reordered indices.
      var descriptor = MM4StateDescriptor()
      descriptor.velocities = true
      return state(descriptor: descriptor).velocities!
    }
    set {
      // reordered -> original -> reordered
      let array = OpenMM_Vec3Array(size: system.parameters.atoms.count)
      for (reordered, original) in system.originalIndices.enumerated() {
        let velocity = newValue[Int(original)]
        array[reordered] = SIMD3<Double>(velocity)
      }
      context.context.velocities = array
    }
  }
}

// MARK: - Simulation Setup

extension MM4ForceField {
  /// Indices of atoms that should ignore forces exerted on them.
  ///
  /// > Warning: Anchors may cause energy measurements to be
  /// nonsensical. This needs to be investigated further.
  ///
  /// This is implemented by treating particle like it has infinite mass (in
  /// OpenMM, by setting the mass to zero). Either set the particle's velocity
  /// to zero, or give all anchors in same rigid body the same linear velocity.
  /// Otherwise, the divergent paths will eventually rip the rigid body apart.
  /// The library currently doesn't check for adherence to this rule, but may
  /// enforce it in the future.
  public var anchors: [UInt32] {
    get {
      // original -> reordered -> original
      system.reorderedIndices.map {
        let anchor = _anchors[Int($0)]
        return UInt32(truncatingIfNeeded: anchor)
      }
    }
    set {
      // Reset every atom's mass to the one provided by the parameters. Then,
      // selectively set the anchors to zero. This may have more overhead, but
      // not that much more (~overhead of setting positions). It also reduces
      // the chance for bugs in a rarely tested edge case.
      for (index, mass) in system.parameters.atoms.masses.enumerated() {
        let reordered = system.reorderedIndices[Int(index)]
        system.system.setParticleMass(Double(mass), index: Int(reordered))
      }
      
      // reordered -> original -> reordered
      _anchors = system.originalIndices.map {
        let anchor = newValue[Int($0)]
        return Int32(truncatingIfNeeded: anchor)
      }
      _anchors.sort()
      if _anchors.count > 0 {
        for anchorID in 1..<_anchors.count {
          if _anchors[anchorID - 1] == _anchors[anchorID] {
            fatalError("Anchor '\(anchorID)' was entered multiple times.")
          }
        }
      }
      _anchors.forEach { index in
        system.system.setParticleMass(0, index: Int(index))
      }
      
      // TODO: Finish the implementation. Change so that anchors are stored in
      // the pre-reordering format. Perform some extra checks during
      // thermalization to convert between formats.
      fatalError("Account for masses of virtual sites when setting anchors.")
    }
  }
  
  /// The constant force (in piconewtons) exerted on each atom.
  ///
  /// The default value is zero for every atom.
  public var externalForces: [SIMD3<Float>] {
    get {
      _externalForces
    }
    set {
      guard newValue.count == system.parameters.atoms.count else {
        fatalError("Too few atoms.")
      }
      _externalForces = newValue
      
      let force = system.forces.external
      force.updateForces(_externalForces, context: context)
    }
  }
  
  /// Atom indices for each rigid body.
  ///
  /// > Note: This is similar to molecules (`OpenMM_Context.getMolecules()`),
  /// with an additional restriction. The user must enter atoms for each
  /// molecule in one contiguous range of the atom list. Otherwise, the
  /// forcefield cannot initialize. See <doc:MM4ParametersDescriptor/bonds> for
  /// more details.
  ///
  /// Rigid bodies should have atoms laid out contiguously in memory, in Morton
  /// order. This format ensures spatial locality, which increases performance
  /// of nonbonded forces. Therefore, rigid bodies are contiguous ranges of the
  /// atom list.
  ///
  /// The set of rigid bodies must cover every atom in the system. No two ranges
  /// may overlap the same atom. If the array of rigid bodies is unspecified, it
  /// defaults to a range encompassing the entire system. This ensures the
  /// closed system's net momentum stays conserved.
  ///
  /// If a body can flex around a joint, it is no longer a rigid body.
  /// Suppressing macroscopic motion from randomly initialized thermal
  /// velocities becomes non-trivial. For such cases, you must explicitly
  /// thermalize the velocities.
  /// <doc:MM4ForceField/thermalize(temperature:rigidBodies:)> lets you choose
  /// which rigid bodies are thermalized using the MM4 algorithm. Remove
  /// flexible bodies from the list and manually thermalize using a similar
  /// algorithm.
  public var rigidBodies: [Range<UInt32>] {
    get {
      // Create a new array with a different type. This isn't the fastest
      // approach, but the property should rarely be used in most cases. The
      // user should already have a data structure that separates the atoms
      // into rigid bodies during high-level operations.
      system.parameters.rigidBodies.map {
        let lowerBound = UInt32(truncatingIfNeeded: $0.lowerBound)
        let upperBound = UInt32(truncatingIfNeeded: $0.upperBound)
        return lowerBound..<upperBound
      }
    }
  }
}
