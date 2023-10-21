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
  /// body's bulk velocity is conserved at the average velocity of all
  /// anchors (where each anchor has the same weight). Angular momentum is
  /// constrained according to the number of anchors present.
  /// - 0 anchors: conserve linear and angular momentum around center of mass.
  /// - 1 anchor: conserve linear and angular momentum around anchor.
  /// - collinear\* anchors: conserve linear momentum around average of
  ///   anchors, constrain angular momentum to the shared axis.
  /// - anchors form plane: conserve linear momentum around average of anchors,
  ///   force angular momentum to zero.
  ///
  /// > \*To be classified as collinear, a group of anchors must share the same
  /// linear velocity, in addition to the same axis (with a tight margin for
  /// floating point error).
  public func thermalize(
    temperature: Double = 298.15,
    rigidBodies: [Int]? = nil
  ) throws {
    // The code currently doesn't recognize the case of 2+ collinear anchors.
    // That will be deferred until the second development round, when constant
    // torques are introduced.
    //
    // Notes about angular momentum:
    // https://www.physicsforums.com/threads/how-can-angular-velocity-be-independent-of-the-choice-of-origin.986098/#:~:text=Both%20the%20angular%20momentum%20and,the%20angular%20velocity%20does%20not.
    
  }
}

// MARK: - Math Utilities for Thermalization

/// Cross-platform implementation of the cross product.
///
/// Source: [Wikipedia](https://en.wikipedia.org/wiki/Cross_product#Computing)
func cross<T: BinaryFloatingPoint & SIMDScalar>(
  _ x: SIMD3<T>, _ y: SIMD3<T>
) -> SIMD3<T> {
  let s1 = x[1] * y[2] - x[2] * y[1]
  let s2 = x[2] * y[0] - x[0] * y[2]
  let s3 = x[0] * y[1] - x[1] * y[0]
  return SIMD3(s1, s2, s3)
}

/// Cross-platform implementation of the 3x3 matrix inverse.
///
/// Source: [Stack Overflow](https://stackoverflow.com/a/18504573)
struct RotationalInertia {
  /// The position to compute rotation around.
  var origin: SIMD3<Float>
  
  /// The accumulator for the rigid body's moment of inertia.
  var columns: (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>)
  
  /// Initialize a moment of inertia with zero mass.
  init(origin: SIMD3<Float>) {
    self.origin = origin
    self.columns = (.zero, .zero, .zero)
  }
}
