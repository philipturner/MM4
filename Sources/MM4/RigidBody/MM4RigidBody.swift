//
//  MM4RigidBody.swift
//
//
//  Created by Philip Turner on 11/19/23.
//

public struct MM4RigidBodyDescriptor {
  /// Required. The number of protons in each atom's nucleus.
  public var atomicNumbers: [UInt8]?
  
  /// Required. Pairs of atom indices representing sigma bonds.
  public var bonds: [SIMD2<UInt32>]?
  
  /// Required. The amount of mass (in amu) to redistribute from a substituent
  /// atom to each covalently bonded hydrogen.
  ///
  /// The default is 1 amu.
  public var hydrogenMassRepartitioning: Float = 1.0
  
  /// Required. The position (in nanometers) of each atom's nucleus.
  public var positions: [SIMD3<Float>]?
  
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

public struct MM4RigidBody {
  /// A wrapper for ergonomically vectorizing computations.
  var atoms: MM4RigidBodyAtoms
  
  /// The rigid body's energy.
  public var energy: MM4RigidBodyEnergy
  
  /// The force field parameters cached for this rigid body.
  public let parameters: MM4Parameters
  
  /// The backing storage object.
  var storage: MM4RigidBodyStorage
  
  public init(descriptor: MM4RigidBodyDescriptor) throws {
    // Ensure the required descriptor properties were set.
    guard let descriptorAtomicNumbers = descriptor.atomicNumbers,
          let descriptorBonds = descriptor.bonds,
          let descriptorPositions = descriptor.positions else {
      fatalError("Descriptor did not have the required properties.")
    }
    self.atoms = MM4RigidBodyAtoms(count: descriptorAtomicNumbers.count)
    self.energy = MM4RigidBodyEnergy()
    
    var desc = MM4ParametersDescriptor()
    desc.atomicNumbers = descriptorAtomicNumbers
    desc.bonds = descriptorBonds
    desc.hydrogenMassRepartitioning = descriptor.hydrogenMassRepartitioning
    self.parameters = try MM4Parameters(descriptor: desc)
    
    self.storage = MM4RigidBodyStorage(atoms: atoms, parameters: parameters)
    descriptorPositions.withUnsafeBufferPointer {
      setPositions($0)
    }
    
    // Prepare computed properties for access through the public API.
    ensureReferencesUpdated()
  }
}
