//
//  MM4Error.swift
//
//
//  Created by Philip Turner on 10/20/23.
//

import OpenMM

/// Identifying information for an atom.
public struct MM4Address {
  /// The rigid body ID.
  var rigidBodyIndex: Int
  
  /// The position within the rigid body.
  var atomIndex: UInt32
  
  /// The atom's element.
  var atomicNumber: UInt8
}

/// An error thrown by the MM4 library.
public enum MM4Error: Error {
  /// The system's energy exploded.
  ///
  /// Includes the amount of energy drift during the time interval that failed.
  case energyDrift(Double)
  
  // TODO: MM4ForceField should intercept parameter errors, changing the rigid
  // body ID to the correct value.
  
  /// The force field did not have a parameter for a group of atoms.
  ///
  /// Includes the addresses of the atoms.
  case missingParameter([MM4Address])
  
  /// The atom's bond count did not match its valence.
  ///
  /// Includes the address of the atom, and the atoms it was detected as bonding
  /// to.
  case openValenceShell(MM4Address, [MM4Address])
  
  /// Lone atom centers such as methane, silane, germane, and tetrafluoromethane
  /// are not supported.
  ///
  /// Includes the address of the atom center, and the atoms it was detected as
  /// bonding to.
  case unsupportedCenterType(MM4Address, [MM4Address])
  
  /// The force field did not support the ring; it was too small.
  ///
  /// Includes the addresses of the atoms. The number of addresses equals the
  /// ring size.
  case unsupportedRing([MM4Address])
}
