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

/// A descriptor for a rigid body.
public struct MM4RigidBodyDescriptor {
  /// Optional. Indices of atoms that should be treated as having infinite mass.
  ///
  /// If not specified, the default is an empty set.
  public var anchors: Set<UInt32>?
  
  /// Optional. The constant force (in piconewtons) exerted on each atom.
  ///
  /// If not specified, the default is zero for every atom.
  public var externalForces: [SIMD3<Float>]?
  
  /// Required. The parameters for the rigid body.
  public var parameters: MM4Parameters?
  
  /// Required. The position (in nanometers) of each atom's nucleus.
  public var positions: [SIMD3<Float>]?
  
  /// Optional. The velocity (in nanometers per picosecond) of each atom.
  ///
  /// If not specified, the default is zero for every atom.
  public var velocities: [SIMD3<Float>]?
  
  public init() {
    
  }
}

/// An enclosed group of covalently bonded atoms.
public struct MM4RigidBody {
  /// The rigid body's energy.
  var _energy: MM4RigidBodyEnergy
  
  /// The force field parameters cached for this rigid body.
  public let parameters: MM4Parameters
  
  /// The backing storage object.
  var storage: MM4RigidBodyStorage
  
  /// Create a rigid body using the specified configuration.
  public init(descriptor: MM4RigidBodyDescriptor) {
    // Ensure the required descriptor properties were set.
    guard let descriptorParameters = descriptor.parameters,
          let descriptorPositions = descriptor.positions else {
      fatalError("Descriptor did not have the required properties.")
    }
    self._energy = MM4RigidBodyEnergy()
    self.parameters = descriptorParameters
    
    let atomCount = descriptorParameters.atoms.count
    let anchors = descriptor.anchors ?? []
    if anchors.contains(where: { $0 >= atomCount }) {
      fatalError("An anchor was out of bounds.")
    }
    precondition(anchors.count == 0, "Anchors are not supported yet.")
    
    self.storage = MM4RigidBodyStorage(
      anchors: anchors, parameters: parameters)
    descriptor.externalForces?.withUnsafeBufferPointer {
      // setExternalForces($0)
      _ = $0
    }
    descriptorPositions.withUnsafeBufferPointer {
      setPositions($0)
    }
    descriptor.velocities?.withUnsafeBufferPointer {
      setVelocities($0)
    }
    
    // Prepare computed properties for access through the public API.
    _ensureReferencesUpdated()
  }
}

// MARK: - Properties

extension MM4RigidBody {
  /// Indices of atoms that should be treated as having infinite mass.
  public var anchors: Set<UInt32> {
    // _modify not supported b/c it requires very complex caching logic.
    // Workaround: import a new rigid body initialized with different anchors.
    _read { yield storage.anchors }
  }
  
  /// The number of protons in each atom's nucleus.
  public var atomicNumbers: [UInt8] {
    parameters.atoms.atomicNumbers
  }
  
  /// Pairs of atom indices representing sigma bonds.
  public var bonds: [SIMD2<UInt32>] {
    parameters.bonds.indices
  }
  
  /// The total mass (in yoctograms) of all atoms, excluding anchors.
  public var mass: Double {
    storage.nonAnchorMass
  }
  
  /// The mass (in yoctograms) of each atom.
  ///
  /// This is different than the masses in `MM4Parameters`. `MM4RigidBody`
  /// zeroes out the mass for each anchor, while `MM4Parameters` does not.
  public var masses: [Float] {
    storage.nonAnchorMasses
  }
}
