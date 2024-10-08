import XCTest
import MM4
import QuaternionModule

// MARK: - Test Execution

final class MM4RigidBodyMomentumTests: XCTestCase {
  
  func testInitialVelocities() throws {
    for (_, descriptor) in MM4RigidBodyTests.descriptors {
      let rigidBody = try MM4RigidBody(descriptor: descriptor)
      for i in rigidBody.masses.indices {
        XCTAssertEqual(
          rigidBody.positions[i], descriptor.positions![i], accuracy: 1e-3)
      }
      XCTAssertEqual(rigidBody.velocities.count, descriptor.positions!.count)
      XCTAssert(rigidBody.velocities.allSatisfy { $0 == .zero })
    }
  }
  
  func testLinearMomentum() throws {
    for (_, descriptor) in MM4RigidBodyTests.descriptors {
      var bulkVelocities: [SIMD3<Double>] = []
      bulkVelocities.append(.zero)
      while true {
        // Add one random velocity with up to 200 m/s for each component.
        let random = SIMD3<Double>.random(in: -0.2...0.2)
        if any(random .!= .zero) {
          bulkVelocities.append(random)
          break
        } else {
          continue
        }
      }
      
      // Test that a pre-existing linear velocity is recognized.
      for bulkVelocity in bulkVelocities {
        var rigidBodyDesc = descriptor
        rigidBodyDesc.velocities = Array(
          repeating: SIMD3(bulkVelocity), count: descriptor.masses!.count)
        let rigidBody = try MM4RigidBody(descriptor: rigidBodyDesc)
        let computedVelocity = rigidBody.linearMomentum / rigidBody.mass
        XCTAssertEqual(computedVelocity, bulkVelocity, accuracy: 1e-5)
      }
      
      // Test that the original linear velocity is zero, but a modified linear
      // velocity is reflected in the 'velocities' property. Ensure the property
      // is all zeroes before linear velocity changes.
      for bulkVelocity in bulkVelocities {
        var rigidBody = try MM4RigidBody(descriptor: descriptor)
        XCTAssertEqual(rigidBody.linearMomentum, .zero, accuracy: 1e-5)
        for i in descriptor.masses!.indices {
          XCTAssertEqual(rigidBody.velocities[i], .zero, accuracy: 1e-5)
        }
        rigidBody.linearMomentum = bulkVelocity * rigidBody.mass
        
        XCTAssertEqual(
          rigidBody.linearMomentum / rigidBody.mass,
          bulkVelocity, accuracy: 1e-5)
        for i in descriptor.masses!.indices {
          XCTAssertEqual(
            rigidBody.velocities[i], SIMD3(bulkVelocity), accuracy: 1e-5)
        }
        XCTAssertEqual(
          rigidBody.linearMomentum / rigidBody.mass,
          bulkVelocity, accuracy: 1e-5)
      }
      
      // Test what happens when you change a rigid body's velocity multiple
      // times.
      var rigidBodyDesc = descriptor
      rigidBodyDesc.velocities = Array(
        repeating: SIMD3(bulkVelocities[1]), count: descriptor.masses!.count)
      var rigidBody = try MM4RigidBody(descriptor: rigidBodyDesc)
      
      XCTAssertEqual(
        rigidBody.linearMomentum / rigidBody.mass,
        (descriptor.masses!.count > 0) ? bulkVelocities[1] : .zero,
        accuracy: 1e-5)
      for i in descriptor.masses!.indices {
        XCTAssertEqual(
          rigidBody.velocities[i], SIMD3(bulkVelocities[1]), accuracy: 1e-5)
      }
      XCTAssertEqual(
        rigidBody.linearMomentum / rigidBody.mass,
        (descriptor.masses!.count > 0) ? bulkVelocities[1] : .zero,
        accuracy: 1e-5)
      
      for _ in 0..<5 {
        // Generate a random 3D vector ranging from -200 m/s to 200 m/s.
        let velocity = SIMD3<Double>.random(in: -0.2...0.2)
        rigidBody.linearMomentum = velocity * rigidBody.mass
        
        XCTAssertEqual(
          rigidBody.linearMomentum,
          (descriptor.masses!.count > 0) ? velocity * rigidBody.mass : .zero,
          accuracy: 1e-5)
        for i in descriptor.masses!.indices {
          XCTAssertEqual(
            rigidBody.velocities[i], SIMD3(velocity), accuracy: 1e-5)
        }
        XCTAssertEqual(
          rigidBody.linearMomentum,
          (descriptor.masses!.count > 0) ? velocity * rigidBody.mass : .zero,
          accuracy: 1e-5)
      }
    }
  }
  
  func testAngularMomentum() throws {
    for descriptor in MM4RigidBodyTests.descriptors {
      try _testAngularMomentum(descriptor.1)
    }
  }
  
  func testForces() throws {
    for (_, descriptor) in MM4RigidBodyTests.descriptors {
      var rigidBody = try MM4RigidBody(descriptor: descriptor)
      XCTAssertNil(rigidBody.forces)
      rigidBody.forces = nil
      XCTAssertNil(rigidBody.forces)
      
      for zeroForces in [false, true] {
        var forces: [SIMD3<Float>] = []
        for _ in rigidBody.masses.indices {
          if zeroForces {
            forces.append(.zero)
          } else {
            forces.append(.random(in: -0.1...0.1))
          }
        }
        
        rigidBody.forces = forces
        XCTAssertNotNil(rigidBody.forces)
        XCTAssertEqual(rigidBody.forces, forces)
        XCTAssertNotNil(rigidBody.netForce)
        XCTAssertNotNil(rigidBody.netTorque)
        
        if zeroForces {
          XCTAssertEqual(rigidBody.netForce, SIMD3<Double>.zero)
          XCTAssertEqual(rigidBody.netTorque, SIMD3<Double>.zero)
          rigidBody.centerOfMass.y += 0.001
        } else {
          XCTAssertNotEqual(rigidBody.netForce, SIMD3<Double>.zero)
          XCTAssertNotEqual(rigidBody.netTorque, SIMD3<Double>.zero)
          
          let rotation = Quaternion<Double>(angle: 0.001, axis: [0, 1, 0])
          rigidBody.rotate(quaternion: rotation)
        }
        XCTAssertNil(rigidBody.forces)
        XCTAssertNil(rigidBody.netForce)
        XCTAssertNil(rigidBody.netTorque)
      }
    }
  }
  
  func testNetForce() throws {
    for (_, descriptor) in MM4RigidBodyTests.descriptors {
      var rigidBody = try MM4RigidBody(descriptor: descriptor)
      
      var forces = [SIMD3<Float>](
        repeating: .zero, count: rigidBody.positions.count)
      forces[1] = [8, 0, 0]
      forces[5] = [0, 9, 0]
      forces[6] = [0, -0.5, -7]
      rigidBody.forces = forces
      XCTAssertEqual(rigidBody.netForce, SIMD3<Double>(8, 8.5, -7))
      
      rigidBody.centerOfMass += SIMD3(repeating: 0.1)
      XCTAssertNil(rigidBody.netForce)
    }
  }
  
  func testNetTorque() throws {
    for (_, descriptor) in MM4RigidBodyTests.descriptors {
      var rigidBody = try MM4RigidBody(descriptor: descriptor)
      
      var forces = [SIMD3<Float>](
        repeating: .zero, count: rigidBody.positions.count)
      var netForce: SIMD3<Double> = .zero
      var netTorque: SIMD3<Double> = .zero
      
      // The first principal axis for the adamantanes is [0, 0, 1]. Create a
      // net torque with a large magnitude along this axis.
      for i in forces.indices {
        let r = rigidBody.positions[i] - SIMD3(rigidBody.centerOfMass)
        var f = SIMD3(-r.y, r.x, 0)
        f /= (f * f).sum().squareRoot()
        let rCrossF = cross(r, f)
        
        netForce += SIMD3(f)
        netTorque += SIMD3(rCrossF)
        forces[i] = f
      }
      rigidBody.forces = forces
      
      // check that torque along the first principal axis = torque.z
      // check that the RMS torque along the remaining two = torque.xy
      let expectedTorque2D = SIMD2(netTorque.x, netTorque.y)
      let actualTorque2D = SIMD2(
        rigidBody.netTorque![1], rigidBody.netTorque![2])
      let expectedTorqueRMS = (
        expectedTorque2D * expectedTorque2D).sum().squareRoot()
      let actualTorqueRMS = (
        actualTorque2D * actualTorque2D).sum().squareRoot()
      XCTAssertEqual(netForce, rigidBody.netForce!, ratioAccuracy: 1e-6)
      XCTAssertEqual(netTorque.z, rigidBody.netTorque![0], accuracy: 1e-6)
      XCTAssertEqual(expectedTorqueRMS, actualTorqueRMS, accuracy: 1e-7)
    }
  }
}

// MARK: - Momentum

// Test that when certain velocities are entered into the object descriptor,
// it automatically recognizes the correct momentum. Test what happens when
// both linear and angular momentum are nonzero.

private func _testAngularMomentum(
  _ descriptor: MM4RigidBodyDescriptor
) throws {
  var bulkVelocities: [SIMD3<Double>] = []
  bulkVelocities.append(.zero)
  while true {
    // Add one random velocity with up to ~0.020 rad/ps angular speed.
    // - This equates to 20 / (2 * .pi) GHz.
    let random = SIMD3<Double>.random(in: -0.02...0.02)
    if any(random .!= .zero) {
      bulkVelocities.append(random)
      break
    } else {
      continue
    }
  }
  
  // Create the center of mass for initializing angular velocities.
  var centerOfMass: SIMD3<Double>
  do {
    let rigidBody = try MM4RigidBody(descriptor: descriptor)
    centerOfMass = rigidBody.centerOfMass
  }
  
  // Helper function for creating expected velocities.
  func createVelocities(
    _ angularVelocity: SIMD3<Double>
  ) -> [SIMD3<Float>] {
    var output: [SIMD3<Float>] = []
    for position in descriptor.positions! {
      let delta = position - SIMD3(centerOfMass)
      let velocity = cross(SIMD3(angularVelocity), delta)
      output.append(velocity)
    }
    return output
  }
  
  // Test that a pre-existing angular velocity is recognized.
  for bulkVelocity in bulkVelocities {
    var rigidBodyDesc = descriptor
    rigidBodyDesc.velocities = createVelocities(bulkVelocity)
    let rigidBody = try MM4RigidBody(descriptor: rigidBodyDesc)
    
    let w = rigidBody.angularMomentum / rigidBody.momentOfInertia
    let Σ = rigidBody.principalAxes
    let Σw = Σ.0 * w[0] + Σ.1 * w[1] + Σ.2 * w[2]
    XCTAssertEqual(
      Σw, (descriptor.masses!.count > 0) ? bulkVelocity : .zero,
      accuracy: 1e-5)
  }
  
  // Test what happens when a linear velocity is also present initially.
  for bulkVelocity in bulkVelocities {
    var originalVelocities = createVelocities(bulkVelocity)
    let linearVelocity = SIMD3<Double>.random(in: -0.2...0.2)
    for i in descriptor.masses!.indices {
      originalVelocities[i] += SIMD3(linearVelocity)
    }
    
    var rigidBodyDesc = descriptor
    rigidBodyDesc.velocities = originalVelocities
    let rigidBody = try MM4RigidBody(descriptor: rigidBodyDesc)
    
    do {
      let w = rigidBody.angularMomentum / rigidBody.momentOfInertia
      let Σ = rigidBody.principalAxes
      let Σw = Σ.0 * w[0] + Σ.1 * w[1] + Σ.2 * w[2]
      XCTAssertEqual(Σw, bulkVelocity, accuracy: 1e-5)
    }
    do {
      let Σ = rigidBody.principalAxes
      let ΣT = (
        SIMD3(Σ.0[0], Σ.1[0], Σ.2[0]),
        SIMD3(Σ.0[1], Σ.1[1], Σ.2[1]),
        SIMD3(Σ.0[2], Σ.1[2], Σ.2[2]))
      let Σw = bulkVelocity
      let w = ΣT.0 * Σw[0] + ΣT.1 * Σw[1] + ΣT.2 * Σw[2]
      let L = w * rigidBody.momentOfInertia
      XCTAssertEqual(rigidBody.angularMomentum, L, accuracy: 1e-5)
    }
    
    XCTAssertEqual(
      rigidBody.linearMomentum / rigidBody.mass,
      (descriptor.masses!.count > 0) ? linearVelocity : .zero,
      accuracy: 1e-5)
    XCTAssertEqual(
      rigidBody.linearMomentum,
      (descriptor.masses!.count > 0) ? linearVelocity * rigidBody.mass : .zero,
      accuracy: 1e-5 * rigidBody.mass)
  }
  
  // Test what happens when a zero velocity becomes a nonzero angular velocity.
  for bulkVelocity in bulkVelocities {
    var rigidBody = try MM4RigidBody(descriptor: descriptor)
    XCTAssert(rigidBody.angularMomentum == .zero)
    for i in descriptor.masses!.indices {
      XCTAssertEqual(rigidBody.velocities[i], .zero, accuracy: 1e-5)
    }
    
    // Fix this.
    let Σ = rigidBody.principalAxes
    let ΣT = (
      SIMD3(Σ.0[0], Σ.1[0], Σ.2[0]),
      SIMD3(Σ.0[1], Σ.1[1], Σ.2[1]),
      SIMD3(Σ.0[2], Σ.1[2], Σ.2[2]))
    let Σw = bulkVelocity
    let w = ΣT.0 * Σw[0] + ΣT.1 * Σw[1] + ΣT.2 * Σw[2]
    rigidBody.angularMomentum = w * rigidBody.momentOfInertia
    
    func createComputedVelocity() -> SIMD3<Double> {
      let w = rigidBody.angularMomentum / rigidBody.momentOfInertia
      let Σ = rigidBody.principalAxes
      let Σw = Σ.0 * w[0] + Σ.1 * w[1] + Σ.2 * w[2]
      return Σw
    }
    
    var computedVelocity = createComputedVelocity()
    XCTAssertEqual(
      computedVelocity,
      (descriptor.masses!.count > 0) ? bulkVelocity : .zero,
      accuracy: 1e-5)
    
    let expectedVelocities = createVelocities(bulkVelocity)
    for i in descriptor.masses!.indices {
      let position = descriptor.positions![i]
      XCTAssertEqual(position, rigidBody.positions[i], accuracy: 1e-3)
      let velocity = expectedVelocities[i]
      XCTAssertEqual(velocity, rigidBody.velocities[i], accuracy: 1e-3)
    }
    
    computedVelocity = createComputedVelocity()
    XCTAssertEqual(
      computedVelocity,
      (descriptor.masses!.count > 0) ? bulkVelocity : .zero,
      accuracy: 1e-5)
  }
  
  // Test a sequence of mutations to linear and/or angular velocity.
  var testOrder: [(Bool, Bool)] = []
  for trialID in 0..<12 {
    let changeLinear = (trialID % 2) >= 1
    let changeAngular = (trialID % 4) >= 2
    testOrder.append((changeLinear, changeAngular))
  }
  testOrder.shuffle()
  
  var rigidBody = try MM4RigidBody(descriptor: descriptor)
  var currentAngularVelocity: SIMD3<Double> = .zero
  var currentLinearVelocity: SIMD3<Double> = .zero
  for (changeLinear, changeAngular) in testOrder {
    if changeLinear {
      currentLinearVelocity = .random(in: -0.2...0.2)
      rigidBody.linearMomentum = currentLinearVelocity * rigidBody.mass
    }
    if changeAngular {
      currentAngularVelocity = .random(in: -0.020...0.020)
      rigidBody.angularMomentum = currentAngularVelocity
      rigidBody.angularMomentum *= rigidBody.momentOfInertia
    }
    
    
    if Bool.random() {
      XCTAssertEqual(
        rigidBody.linearMomentum / rigidBody.mass,
        (descriptor.masses!.count > 0) ? currentLinearVelocity : .zero,
        accuracy: 1e-5)
    } else {
      let expectedVelocity = currentLinearVelocity * rigidBody.mass
      XCTAssertEqual(
        rigidBody.linearMomentum,
        (descriptor.masses!.count > 0) ? expectedVelocity : .zero,
        accuracy: 1e-5 * rigidBody.mass)
    }
    do {
      let w = rigidBody.angularMomentum / rigidBody.momentOfInertia
      XCTAssertEqual(w, currentAngularVelocity, accuracy: 1e-5)
    }
    
    let Σ = {
      let f64 = rigidBody.principalAxes
      let f32 = (
        SIMD3<Float>(f64.0),
        SIMD3<Float>(f64.1),
        SIMD3<Float>(f64.2))
      return f32
    }()
    let ΣT = (
      SIMD3(Σ.0[0], Σ.1[0], Σ.2[0]),
      SIMD3(Σ.0[1], Σ.1[1], Σ.2[1]),
      SIMD3(Σ.0[2], Σ.1[2], Σ.2[2]))
    
    for i in descriptor.masses!.indices {
      let position = descriptor.positions![i]
      var r = position - SIMD3(centerOfMass)
      r = ΣT.0 * r[0] + ΣT.1 * r[1] + ΣT.2 * r[2]
      
      var v = cross(SIMD3(currentAngularVelocity), r)
      v = Σ.0 * v[0] + Σ.1 * v[1] + Σ.2 * v[2]
      v += SIMD3(currentLinearVelocity)
      
      let actualVelocity = rigidBody.velocities[i]
      XCTAssertEqual(v, actualVelocity, accuracy: 1e-5)
    }
  }
}

// MARK: - Utilities

func cross(_ a: SIMD3<Float>, _ b: SIMD3<Float>) -> SIMD3<Float> {
  var c: SIMD3<Float> = .zero
  c.x = a.y * b.z - a.z * b.y
  c.y = a.z * b.x - a.x * b.z
  c.z = a.x * b.y - a.y * b.x
  return c
}
