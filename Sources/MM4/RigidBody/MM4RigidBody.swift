//
//  MM4RigidBody.swift
//
//
//  Created by Philip Turner on 11/19/23.
//

// TODO: Remove anchors from the forcefield, have the user simply
// mark such atoms as having 0 mass in MM4Parameters. The user defines what to
// do with energy in edge cases. They can compute it manually if they want to.
// Anchors have 0 contribution to energy and their thermal velocity is zero.
//
// Removing anchors from the rigid body requires a little more care. There may
// be some functions that check the edge case where mass is zero.

/// An enclosed group of covalently bonded atoms.
///
/// `MM4RigidBody` is an API to efficiently compute basic properties in rigid
/// body mechanics. Specifically, a non-consecutive memory layout optimized for
/// the vector units on modern CPUs.
///
/// > NOTE: This documentation page is still a draft. It may be inconsistent or
///   difficult to understand.
///
/// The rigid body exists caches the vectorized layout for positions/velocities
/// and efficiently computes properties like moment of inertia.
/// It also efficiently initializes thermal velocities. The user is welcome to
/// create their own data structure wrapping 'MM4Parameters', which avoids the
/// conversion from Vec3 to a swizzled representation. They could also make such
/// data structures only for a tiny fraction of atoms that require more
/// fine-grained access to mutating elements. If such an API were created in
/// MM4, it would have no purpose. There needs to be some kind of unique
/// functionality provided by the code in MM4, such as very optimized methods
/// for computing certain properties + CoW caching. Stuff that would be
/// extremely difficult or tedious to create manually.
///
/// Need performance (API similar to OpenMM, which requires positions and
/// velocities to be set in batches):
/// - MM4RigidBody
/// - API for copying positions/velocities from an atom range to a rigid body
///
/// Need flexibility:
/// - Custom data structure storing an MM4Parameters, along with some positions
///   and velocities.
/// - Can even wrap a rigid body, just with custom or more ergonomic methods to
///   access the positions. Or not; set Vec3 into OpenMM directly if that for
///   some reason has greater performance.
public struct MM4RigidBody {
  /// The rigid body's energy.
  var _energy: MM4RigidBodyEnergy
  
  /// The force field parameters cached for this rigid body.
  public let parameters: MM4Parameters
  
  /// The backing storage object.
  var storage: MM4RigidBodyStorage
  
  /// Create a rigid body using the specified configuration.
  public init(parameters: MM4Parameters) {
    self._energy = MM4RigidBodyEnergy()
    self.parameters = parameters
    self.storage = MM4RigidBodyStorage(
      anchors: [], parameters: parameters)
    
    // Prepare computed properties for access through the public API.
    _ensureReferencesUpdated()
  }
}

// MARK: - Properties

extension MM4RigidBody {
  /// The total mass (in yoctograms) of all atoms, excluding anchors.
  public var mass: Double {
    storage.nonAnchorMass
  }
}
