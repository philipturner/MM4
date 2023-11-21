//
//  MM4ForceField+Thermalize.swift
//
//
//  Created by Philip Turner on 10/21/23.
//

import OpenMM

/// Cross-platform implementation of the cross product.
///
/// Source: [Wikipedia](https://en.wikipedia.org/wiki/Cross_product#Computing)
fileprivate func cross<T: BinaryFloatingPoint & SIMDScalar>(
  _ x: SIMD3<T>, _ y: SIMD3<T>
) -> SIMD3<T> {
  let s1 = x[1] * y[2] - x[2] * y[1]
  let s2 = x[2] * y[0] - x[0] * y[2]
  let s3 = x[0] * y[1] - x[1] * y[0]
  return SIMD3(s1, s2, s3)
}

extension MM4ForceField {
  // TODO: Remove this function entirely, transferring the functionality to
  // `MM4RigidBody`. The effect can be emulated by changing a rigid body's
  // temperature, then importing it into the simulator.
  // - Keep the new documentation that requires users to specify heat capacity.
  
  // TODO: Copy the OpenMM code for initializing thermal velocities.
  
  // TODO: Rescale velocities after zeroing out the object's momentum, to
  // perfectly match the desired temperature.
  
  /// Create random thermal velocities, while conserving the total (bulk)
  /// momentum of each rigid body.
  ///
  /// - Parameter heatCapacity: The material's heat capacity in J/mol-K at the
  ///   specified temperature.
  /// - Parameter temperature: The temperature to randomize thermal velocites
  ///   at, in kelvin.
  /// - Parameter rigidBodies: Indices of the rigid bodies to thermalize. If not
  ///   specified, it will thermalize the entire system.
  /// - throws: An error if two anchors have the same position, which interferes
  ///   with the heuristic for conserving angular momentum.
  public func thermalize(
    heatCapacity: Double,
    temperature: Double = 298.15,
    rigidBodies: [Int]? = nil
  ) throws {
    // TODO: Replace OpenMM implementation with custom implementation, extract
    // this functionality out into 'MM4RigidBody'.
    fatalError("Incorrect mapping from temperature to energy.")
    
    // The code currently doesn't recognize the case of 2+ collinear anchors.
    // That will be deferred to later, when constant torques are introduced.
    //
    // Notes about angular momentum:
    // https://www.physicsforums.com/threads/how-can-angular-velocity-be-independent-of-the-choice-of-origin.986098/#:~:text=Both%20the%20angular%20momentum%20and,the%20angular%20velocity%20does%20not.
    
    var descriptor = MM4StateDescriptor()
    descriptor.positions = true
    descriptor.velocities = true
    let originalState = state(descriptor: descriptor)
    let positions = originalState.positions!
    let bulkVelocities = originalState.velocities!
    
    context.context.setVelocitiesToTemperature(temperature)
    descriptor.positions = false
    let thermalState = state(descriptor: descriptor)
    let thermalVelocities = thermalState.velocities!
    
    var anchorMask = [Bool](
      repeating: false, count: system.parameters.atoms.count)
    for anchor in _anchors {
      anchorMask[Int(anchor)] = true
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
    
    // Set the system's velocities to these at the end of the function.
    var newVelocities = bulkVelocities
    
    for rigidBodyID in activeRigidBodies {
      let range = system.parameters.rigidBodies[rigidBodyID]
      var anchorCount = 0
      for originalID in range {
        if anchorMask[Int(originalID)] {
          anchorCount += 1
        }
      }
      
      let atoms = system.parameters.atoms
      var anchorMass: Double = .zero
      var nonAnchorMass: Double = .zero
      var centerOfMass: SIMD3<Float>
      do {
        var accumulator: SIMD3<Double> = .zero
        for originalID in range {
          let atomID = Int(originalID)
          let mass = Double(atoms.masses[atomID])
          if (anchorCount == 0) || anchorMask[atomID] {
            anchorMass += mass
            accumulator += mass * SIMD3<Double>(positions[atomID])
          }
          if (anchorCount == 0) || !anchorMask[atomID] {
            nonAnchorMass += mass
          }
        }
        centerOfMass = SIMD3<Float>(accumulator / anchorMass)
      }
      
      // Find the original/thermalized linear/angular velocities, while
      // accumulating the moment of inertia.
      var bulkMomentum: SIMD3<Double> = .zero
      var thermalMomentum: SIMD3<Double> = .zero
      var bulkAngularMomentum: SIMD3<Double> = .zero
      var thermalAngularMomentum: SIMD3<Double> = .zero
      var rotationalInertia: MM4AngularMass = .init()
      for originalID in range {
        let atomID = Int(originalID)
        let mass = Double(atoms.masses[atomID])
        let bulkVelocity = bulkVelocities[atomID]
        let thermalVelocity = thermalVelocities[atomID]
        if (anchorCount == 0) || anchorMask[atomID] {
          bulkMomentum += mass * SIMD3(bulkVelocity)
        }
        if (anchorCount == 0) || !anchorMask[atomID] {
          thermalMomentum += mass * SIMD3(thermalVelocity)
        }
        
        // From Wikipedia:
        // https://en.wikipedia.org/wiki/Rigid_body_dynamics#Linear_and_angular_momentum
        //
        // L = m * (R - R_cm) cross d/dt (R - R_cm)
        // assume R_cm is stationary
        //
        // v = dR / dt
        // L = m * (R - R_cm) cross v
        let relativePosition = positions[atomID] - centerOfMass
        let bulkAngularTerm = cross(relativePosition, bulkVelocity)
        let thermalAngularTerm = cross(relativePosition, thermalVelocity)
        
        // If the system has >1 anchor points, suppress the angular momentum.
        // The singular anchor for exactly 1 anchor point will have zero
        // angular momentum or inertia, as it's defined as the center of mass.
        if (anchorCount <= 1) {
          // Unlike linear velocity, bulk angular velocity can be defined by
          // non-anchor atoms.
          bulkAngularMomentum += mass * SIMD3(bulkAngularTerm)
        }
        if (anchorCount == 0) || !anchorMask[atomID] {
          thermalAngularMomentum += mass * SIMD3(thermalAngularTerm)
          rotationalInertia
            .append(mass: mass, relativePosition: relativePosition)
        }
      }
      
      // Matrix:
      // L = I * w
      // I^{-1} L = w
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
        
        // Provide a fail-safe if the entire rigid body consists of anchors.
        if all(momentum .== 0) {
          return .zero
        } else {
          return SIMD3<Float>(output)
        }
      }
      
      let bulkVelocity = bulkMomentum / anchorMass
      let thermalVelocity = thermalMomentum / nonAnchorMass
      let velocityCorrection = -SIMD3<Float>(thermalVelocity - bulkVelocity)
      
      let bulkAngularVelocity = project(momentum: bulkAngularMomentum)
      let thermalAngularVelocity = project(momentum: thermalAngularMomentum)
      let angularVelocityCorrection = -SIMD3<Float>(
        thermalAngularVelocity - bulkAngularVelocity)
      
      // Apply the correction to rotational and angular velocity.
      for originalID in range {
        let atomID = Int(originalID)
        var velocity = thermalVelocities[atomID]
        velocity += velocityCorrection
        
        let relativePosition = positions[atomID] - centerOfMass
        velocity += cross(angularVelocityCorrection, relativePosition)
        if !anchorMask[atomID] {
          newVelocities[atomID] = velocity
        }
      }
    }
    
    // Set the system's velocities to the new ones, reverting changes to rigid
    // bodies that shouldn't be thermalized.
    self.velocities = newVelocities
  }
}
