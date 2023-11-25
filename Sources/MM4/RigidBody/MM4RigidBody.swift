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

/// Wrapper for ergonomically vectorizing computations.
struct MM4RigidBodyAtoms {
  /// A convenient method for accessing the atom count.
  var count: Int
  
  /// A convenient method for accessing the atom vector count.
  var vectorCount: Int
  
  /// A convenient method for addressing the end of the list.
  var vectorMask: SIMDMask<MM4UInt32Vector.MaskStorage>
  
  init(count: Int) {
    self.count = count
    self.vectorCount = count + MM4VectorWidth - 1
    self.vectorCount /= MM4VectorWidth
    
    let lastVectorStart = UInt32((vectorCount - 1) * MM4VectorWidth)
    var lastVector: MM4UInt32Vector = .init(repeating: lastVectorStart)
    for lane in 0..<MM4VectorWidth {
      lastVector[lane] &+= UInt32(lane)
    }
    self.vectorMask = lastVector .>= UInt32(count)
  }
}

/// An enclosed group of covalently bonded atoms.
public struct MM4RigidBody {
  /// A wrapper for ergonomically vectorizing computations.
  var atoms: MM4RigidBodyAtoms
  
  /// The rigid body's energy.
  public var energy: MM4RigidBodyEnergy
  
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
    self.atoms = MM4RigidBodyAtoms(count: descriptorParameters.atoms.count)
    self.energy = MM4RigidBodyEnergy()
    self.parameters = descriptorParameters
    
    self.storage = MM4RigidBodyStorage(atoms: atoms, parameters: parameters)
    descriptorPositions.withUnsafeBufferPointer {
      setPositions($0)
    }
    
    if descriptor.velocities != nil {
      fatalError("Velocity input is not accepted yet.")
    }
    
    // Prepare computed properties for access through the public API.
    ensureReferencesUpdated()
  }
}
