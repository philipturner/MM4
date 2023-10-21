//
//  MM4ForceField+Thermalize.swift
//
//
//  Created by Philip Turner on 10/21/23.
//

extension MM4ForceField {
  /// Create random thermal velocities, while conserving the total (bulk)
  /// momentum of each rigid body.
  ///
  /// - Parameter temperature: The temperature to randomize thermal velocites
  ///   at, in kelvin.
  /// - Parameter rigidBodies: Indices of the rigid bodies to thermalize. If not
  ///   specified, it will thermalize the entire system.
  ///
  /// Thermalizing is recommended for any simulation that replicates macroscale
  /// conditions. The default is temperature 298.15 K, but other useful
  /// temperatures include liquid nitrogen (77.00 K) and liquid helium (4.15 K).
  ///
  /// Anchors have no velocity appended to them during thermalization. Before
  /// angular momentum is corrected around the true center of mass, it is first
  /// corrected around a center defined by any anchors. If a rigid body has only
  /// one anchor, that anchor's position is used. Otherwise, the center is the
  /// average position of all anchors in the rigid body.
  public func thermalize(
    temperature: Double = 298.15,
    rigidBodies: [Int]? = nil
  ) {
    // Precursor: test the correctness of the freely available matrix inversion
    // code against Apple's 'simd' library
  }
}
