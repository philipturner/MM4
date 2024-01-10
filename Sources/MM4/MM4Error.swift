//
//  MM4Error.swift
//  MM4
//
//  Created by Philip Turner on 10/20/23.
//

import OpenMM

/// Identifying information for an atom.
public struct MM4Address {
  /// The atom's position within the object.
  public var index: UInt32
  
  /// The atom's element.
  public var atomicNumber: UInt8
}

/// An error thrown by the MM4 library.
public enum MM4Error: Error {
  /// The inertia tensor could not be diagonalized.
  ///
  /// This error does not exactly mean the inertia tensor was a defective
  /// matrix. For example, the null matrix is technically diagonalizable. It
  /// just means e.g. the center of mass couldn't be created because the mass
  /// was zero. Therefore, the rigid body couldn't reach the point where it
  /// generates the inertia tensor. The error may also occur if the
  /// characteristic polynomial has repeated roots or the eigensolver fails.
  case defectiveInertiaTensor((SIMD3<Double>, SIMD3<Double>, SIMD3<Double>))
  
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
