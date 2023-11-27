//
//  MM4RigidBody.swift
//
//
//  Created by Philip Turner on 11/19/23.
//

/// A descriptor for a rigid body.
public struct MM4RigidBodyDescriptor {
  /// Optional. Indices of atoms that should be treated as having infinite mass.
  public var anchors: Set<UInt32>?
  
  /// Optional. Indices of atoms where external force is applied.
  ///
  /// If not specified, external force is distributed over the entire object.
  public var handles: Set<UInt32>?
  
  /// Required. The parameters for the rigid body.
  public var parameters: MM4Parameters?
  
  /// Required. The position (in nanometers) of each atom's nucleus.
  public var positions: [SIMD3<Float>]?
  
  ///. Optional. The velocity (in nanometers per picosecond) of each atom.
  ///
  /// Velocities attributed to anchors are ignored. They are replaced with a
  /// value determined by the bulk velocities.
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
    var handles: Set<UInt32>
    if let descriptorHandles = descriptor.handles {
      handles = descriptorHandles
      for anchor in anchors {
        if handles.contains(anchor) {
          fatalError("Atom \(anchor) was both an anchor and a handle.")
        }
      }
    } else {
      handles = []
      for atomID in 0..<atomCount {
        handles.insert(UInt32(truncatingIfNeeded: atomID))
      }
    }
    guard !anchors.contains(where: { $0 >= atomCount }),
          !handles.contains(where: { $0 >= atomCount }) else {
      fatalError("An anchor or handle was out of bounds.")
    }
    
    self.storage = MM4RigidBodyStorage(
      anchors: anchors, handles: handles, parameters: parameters)
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

