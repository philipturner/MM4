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
  /// - throws: An error if two anchors have the same position, which interferes
  ///   with the heuristic for conserving angular momentum.
  ///
  /// Thermalizing is recommended for any simulation that replicates macroscale
  /// conditions. The default is temperature 298.15 K, but other useful
  /// temperatures include liquid nitrogen (77.00 K) and liquid helium (4.15 K).
  ///
  /// The velocity of anchors does not change during thermalization. Each rigid
  /// body's bulk velocity is conserved at the average velocity of all the
  /// anchor (where each anchor has the same weight). Angular momentum is
  /// constrained according to the number of anchors present.
  /// - 0 anchors: conserve linear and angular momentum around center of mass
  /// - 1 anchor: conserve linear and angular momentum around anchor
  /// - 2 anchors: conserve linear momentum around average of anchors, constrain
  ///   angular momentum to the axis shared by the two anchors
  /// - 3+ anchors: conserve linear momentum around average of anchors, force
  ///   the object to have zero angular momentum
  /// - 3+ anchors but all are collinear: the same as with 2 anchors. Note that
  ///   this won't be implemented until torques are implemented.
  public func thermalize(
    temperature: Double = 298.15,
    rigidBodies: [Int]? = nil
  ) throws {
    // Notes about angular momentum:
    // https://www.physicsforums.com/threads/how-can-angular-velocity-be-independent-of-the-choice-of-origin.986098/#:~:text=Both%20the%20angular%20momentum%20and,the%20angular%20velocity%20does%20not.
    
    
  }
}
