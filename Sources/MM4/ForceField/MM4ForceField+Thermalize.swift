//
//  MM4ForceField+Thermalize.swift
//
//
//  Created by Philip Turner on 10/21/23.
//

import OpenMM

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
    
    // Using reordered indices.
    var descriptor = MM4StateDescriptor()
    descriptor.positions = true
    descriptor.velocities = true
    let originalState = state(descriptor: descriptor)
    let positions = originalState.positions!
    let bulkVelocities = originalState.velocities!
    
    // Using reordered indices.
    context.context.setVelocitiesToTemperature(temperature)
    descriptor.positions = false
    let thermalState = state(descriptor: descriptor)
    let thermalVelocities = thermalState.velocities!
    
    // Using reordered indices.
    var anchorStatuses = [Bool](
      repeating: false, count: system.parameters.atoms.count)
    for anchor in _anchors {
      anchorStatuses[Int(anchor)] = true
    }
    
    var activeRigidBodies: [Int]
    if let rigidBodies {
      activeRigidBodies = rigidBodies
    } else {
      activeRigidBodies = []
      for rigidBodyID in system.parameters.rigidBodies.indices {
        activeRigidBodies.append(rigidBodyID)
      }
    }
    
    // Set the system's velocities to these at the end of the function. Note,
    // these are in reordered indices, while the 'velocities' property uses
    // original indices. Do not set these using the 'velocities' property.
    var newVelocities: [SIMD3<Float>] = bulkVelocities
    
    // Don't forget to always map atom IDs to reordered IDs!
    for rigidBodyID in activeRigidBodies {
      // Last anchor ID is reordered.
      let range = system.parameters.rigidBodies[rigidBodyID]
      var anchorCount = 0
      var lastAnchorID: Int32 = -1
      for originalID in range {
        let atomID = system.reorderedIndices[Int(originalID)]
        if anchorStatuses[Int(atomID)] {
          anchorCount += 1
          lastAnchorID = atomID
        }
      }
      
      let atoms = system.parameters.atoms
      var totalMass: Double = .zero
      var centerOfMass: SIMD3<Float>
      if anchorCount == 1 {
        // Reading with a reordered index.
        centerOfMass = positions[Int(lastAnchorID)]
      } else {
        var positionSum: SIMD3<Double> = .zero
        for originalID in range {
          let atomID = system.reorderedIndices[Int(originalID)]
          let mass = Double(atoms.masses[Int(originalID)])
          let position = positions[Int(atomID)]
          totalMass += mass
          positionSum += mass * SIMD3<Double>(position)
        }
        centerOfMass = SIMD3<Float>(positionSum / totalMass)
      }
      
      // Find the original/thermalized linear/angular velocities, while
      // accumulating the moment of inertia.
      var bulkMomentum: SIMD3<Double> = .zero
      var thermalMomentum: SIMD3<Double> = .zero
      var bulkAngularMomentum: SIMD3<Double> = .zero
      var thermalAngularMomentum: SIMD3<Double> = .zero
      var rotationalInertia: RotationalInertia = .init()
      for originalID in range {
        let atomID = system.reorderedIndices[Int(originalID)]
        let mass = Double(atoms.masses[Int(originalID)])
        let position = positions[Int(atomID)]
        let bulkVelocity = bulkVelocities[Int(atomID)]
        let thermalVelocity = thermalVelocities[Int(atomID)]
        bulkMomentum += mass * SIMD3(bulkVelocity)
        thermalMomentum += mass * SIMD3(thermalVelocity)
        
        // From Wikipedia:
        // https://en.wikipedia.org/wiki/Rigid_body_dynamics#Linear_and_angular_momentum
        //
        // L = m * (R - R_cm) cross d/dt (R - R_cm)
        // assume R_cm is stationary
        // L = m * (R - R_cm) cross v
        let relativePosition = position - centerOfMass
        let bulkAngularVelocity = cross(relativePosition, bulkVelocity)
        let thermalAngularVelocity = cross(relativePosition, thermalVelocity)
        bulkAngularMomentum += mass * SIMD3(bulkAngularVelocity)
        thermalAngularMomentum += mass * SIMD3(thermalAngularVelocity)
        rotationalInertia.append(mass: mass, relativePosition: relativePosition)
      }
      
      // If the system has >1 anchor points, suppress the angular velocity.
      if anchorCount > 1 {
        bulkAngularMomentum = .zero
      }
      
      // Matrix:
      // L = I * w
      // (I^{-1}) L = w
      // w = angular velocity
      //
      // Resulting vector:
      // w_x: angular velocity around x-axis (YZ plane)
      // w_y: angular velocity around y-axis (ZX plane)
      // w_z: angular velocity around z-axis (XY plane)
      //
      // Convert into a linear velocity for each particle:
      // v = w cross r
      let inverse = rotationalInertia.inverse
      func project(momentum: SIMD3<Double>) -> SIMD3<Float> {
        var output: SIMD3<Double> = .zero
        output += momentum.x * inverse.0
        output += momentum.y * inverse.1
        output += momentum.z * inverse.2
        return SIMD3<Float>(output)
      }
      
      let bulkVelocity = bulkMomentum / totalMass
      let thermalVelocity = thermalMomentum / totalMass
      let velocityCorrection = -SIMD3<Float>(thermalVelocity - bulkVelocity)
      
      let bulkAngularVelocity = project(momentum: bulkAngularMomentum)
      let thermalAngularVelocity = project(momentum: thermalAngularMomentum)
      let angularVelocityCorrection = -SIMD3<Float>(
        thermalAngularVelocity - bulkAngularVelocity)
      
      // Apply the correction to rotational and angular velocity.
      for originalID in range {
        let atomID = system.reorderedIndices[Int(originalID)]
        let position = positions[Int(atomID)]
        var velocity = thermalVelocities[Int(atomID)]
        velocity += velocityCorrection
        
        let relativePosition = position - centerOfMass
        velocity += cross(angularVelocityCorrection, relativePosition)
        newVelocities[Int(atomID)] = velocity
      }
    }
    
    // Set the system's velocities to the new ones, reverting changes to rigid
    // bodies that shouldn't be thermalized.
    let array = OpenMM_Vec3Array(size: system.parameters.atoms.count)
    for (index, velocity) in newVelocities.enumerated() {
      array[index] = SIMD3<Double>(velocity)
    }
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
  /// The accumulator for the rigid body's moment of inertia.
  var columns: (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>)
  
  /// Initialize a moment of inertia with zero mass.
  init() {
    self.columns = (.zero, .zero, .zero)
  }
  
  /// Add an atom to the accumulator.
  mutating func append(mass: Double, relativePosition: SIMD3<Float>) {
    // From Wikipedia:
    // https://en.wikipedia.org/wiki/Rigid_body_dynamics#Mass_properties
    //
    // I_R = m * (I (S^T S) - S S^T)
    // where S is the column vector R - R_cm
    let positionSquared = relativePosition * relativePosition
    let STS = positionSquared[0] + positionSquared[1] + positionSquared[2]
    var column0 = SIMD3(STS, 0, 0)
    var column1 = SIMD3(0, STS, 0)
    var column2 = SIMD3(0, 0, STS)
    column0 -= relativePosition.x * relativePosition
    column1 -= relativePosition.y * relativePosition
    column2 -= relativePosition.z * relativePosition
    
    // Convert to FP64 before adding to the accumulator. The matrix is
    // symmetric, so it doesn't matter whether you mix up the rows and columns.
    columns.0 += mass * SIMD3<Double>(column0)
    columns.1 += mass * SIMD3<Double>(column1)
    columns.2 += mass * SIMD3<Double>(column2)
  }
  
  // The matrix is symmetric, but not exactly orthonormal. Inversion is not a
  // simple transpose operation.
  var inverse: (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>) {
    // Source: https://stackoverflow.com/a/18504573
    //
    // double det = m(0, 0) * (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) -
    //              m(0, 1) * (m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0)) +
    //              m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));
    let determinant =
    columns.0[0]*(columns.1[1]*columns.2[2]-columns.2[1]*columns.1[2]) -
    columns.0[1]*(columns.1[0]*columns.2[2]-columns.1[2]*columns.2[0]) +
    columns.0[2]*(columns.1[0]*columns.2[1]-columns.1[1]*columns.2[0])
    let invdet = 1 / determinant
    
    // minv(0, 0) = (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) * invdet;
    // minv(0, 1) = (m(0, 2) * m(2, 1) - m(0, 1) * m(2, 2)) * invdet;
    // minv(0, 2) = (m(0, 1) * m(1, 2) - m(0, 2) * m(1, 1)) * invdet;
    let result00 = (columns.1[1]*columns.2[2]-columns.2[1]*columns.1[2])*invdet
    let result01 = (columns.0[2]*columns.2[1]-columns.0[1]*columns.2[2])*invdet
    let result02 = (columns.0[1]*columns.1[2]-columns.0[2]*columns.1[1])*invdet
    
    // minv(1, 0) = (m(1, 2) * m(2, 0) - m(1, 0) * m(2, 2)) * invdet;
    // minv(1, 1) = (m(0, 0) * m(2, 2) - m(0, 2) * m(2, 0)) * invdet;
    // minv(1, 2) = (m(1, 0) * m(0, 2) - m(0, 0) * m(1, 2)) * invdet;
    let result10 = (columns.1[2]*columns.2[0]-columns.1[0]*columns.2[2])*invdet
    let result11 = (columns.0[0]*columns.2[2]-columns.0[2]*columns.2[0])*invdet
    let result12 = (columns.1[0]*columns.0[2]-columns.0[0]*columns.1[2])*invdet
    
    // minv(2, 0) = (m(1, 0) * m(2, 1) - m(2, 0) * m(1, 1)) * invdet;
    // minv(2, 1) = (m(2, 0) * m(0, 1) - m(0, 0) * m(2, 1)) * invdet;
    // minv(2, 2) = (m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1)) * invdet;
    let result20 = (columns.1[0]*columns.2[1]-columns.2[0]*columns.1[1])*invdet
    let result21 = (columns.2[0]*columns.0[1]-columns.0[0]*columns.2[1])*invdet
    let result22 = (columns.0[0]*columns.1[1]-columns.1[0]*columns.0[1])*invdet
    
    let column0 = SIMD3(result00, result10, result20)
    let column1 = SIMD3(result01, result11, result21)
    let column2 = SIMD3(result02, result12, result22)
    return (column0, column1, column2)
  }
}

