import XCTest
import MM4

// MARK: - Test Execution

final class MM4RigidBodyVelocityTests: XCTestCase {
  #if false
  func testRigidBodyVelocity() throws {
    let descriptors = try MM4RigidBodyTests.createDescriptors()
    for descriptor in descriptors {
      testLinearMomentum(descriptor)
      testAngularMomentum(descriptor)
    }
  }
  #endif
}

// MARK: - Velocity

// Test that when certain velocities are entered into the object descriptor,
// it automatically recognizes the correct momentum. Test what happens when
// both linear and angular momentum are nonzero.

private func testLinearMomentum(_ descriptor: MM4RigidBodyDescriptor) {
  let parameters = descriptor.parameters!
  
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
      repeating: SIMD3(bulkVelocity), count: parameters.atoms.count)
    let rigidBody = MM4RigidBody(descriptor: rigidBodyDesc)
    let computedVelocity = rigidBody.linearMomentum / rigidBody.mass
    let expectedVelocity = (parameters.atoms.count > 0) ? bulkVelocity : .zero
    XCTAssertEqual(computedVelocity, expectedVelocity, accuracy: 1e-5)
  }
  
  // Test that the original linear velocity is zero, but a modified linear
  // velocity is reflected in the 'velocities' property. Ensure the property is
  // all zeroes before linear velocity changes.
  for bulkVelocity in bulkVelocities {
    var rigidBody = MM4RigidBody(descriptor: descriptor)
    XCTAssertEqual(rigidBody.linearMomentum, .zero, accuracy: 1e-5)
    for i in parameters.atoms.indices {
      XCTAssertEqual(rigidBody.velocities[i], .zero, accuracy: 1e-5)
    }
    rigidBody.linearMomentum = bulkVelocity * rigidBody.mass
    
    let expectedVelocity = (parameters.atoms.count > 0) ? bulkVelocity : .zero
    XCTAssertEqual(
      rigidBody.linearMomentum / rigidBody.mass,
      expectedVelocity, accuracy: 1e-5)
    for i in parameters.atoms.indices {
      XCTAssertEqual(
        rigidBody.velocities[i], SIMD3(bulkVelocity), accuracy: 1e-5)
    }
    XCTAssertEqual(
      rigidBody.linearMomentum / rigidBody.mass,
      expectedVelocity, accuracy: 1e-5)
  }
  
  // Test what happens when you change a rigid body's velocity multiple times.
  var rigidBodyDesc = descriptor
  rigidBodyDesc.velocities = Array(
    repeating: SIMD3(bulkVelocities[1]), count: parameters.atoms.count)
  var rigidBody = MM4RigidBody(descriptor: descriptor)
  
  XCTAssertEqual(
    rigidBody.linearMomentum / rigidBody.mass,
    (parameters.atoms.count > 0) ? bulkVelocities[1] : .zero,
    accuracy: 1e-5)
  for i in parameters.atoms.indices {
    XCTAssertEqual(
      rigidBody.velocities[i], SIMD3(bulkVelocities[1]), accuracy: 1e-5)
  }
  XCTAssertEqual(
    rigidBody.linearMomentum / rigidBody.mass,
    (parameters.atoms.count > 0) ? bulkVelocities[1] : .zero,
    accuracy: 1e-5)
  
  for _ in 0..<5 {
    // Generate a random 3D vector ranging from -200 m/s to 200 m/s.
    let velocity = SIMD3<Double>.random(in: -0.2...0.2)
    rigidBody.linearMomentum = velocity * rigidBody.mass
    
    XCTAssertEqual(
      rigidBody.linearMomentum,
      (parameters.atoms.count > 0) ? velocity * rigidBody.mass : .zero,
      accuracy: 1e-5)
    for i in parameters.atoms.indices {
      XCTAssertEqual(rigidBody.velocities[i], SIMD3(velocity), accuracy: 1e-5)
    }
    XCTAssertEqual(
      rigidBody.linearMomentum,
      (parameters.atoms.count > 0) ? velocity * rigidBody.mass : .zero,
      accuracy: 1e-5)
  }
}

// TODO: This function will likely fail because the angular velocity is w.r.t.
// the local reference frame (diagonalized). Modify the test accordingly.

private func testAngularMomentum(_ descriptor: MM4RigidBodyDescriptor) {
  let parameters = descriptor.parameters!
  
  func cross(_ a: SIMD3<Float>, _ b: SIMD3<Float>) -> SIMD3<Float> {
    var c: SIMD3<Float> = .zero
    c.x = a.y * b.z - a.z * b.y
    c.y = a.z * b.x - a.x * b.z
    c.z = a.x * b.y - a.y * b.x
    return c
  }
  
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
    let rigidBody = MM4RigidBody(descriptor: descriptor)
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
    let rigidBody = MM4RigidBody(descriptor: rigidBodyDesc)
    
    let computedVelocity = rigidBody.angularMomentum
    XCTAssertEqual(
      computedVelocity,
      (parameters.atoms.count > 0) ? bulkVelocity : .zero,
      accuracy: 1e-5)
  }
  
  // Test what happens when a linear velocity is also present initially.
  for bulkVelocity in bulkVelocities {
    var originalVelocities = createVelocities(bulkVelocity)
    let linearVelocity = SIMD3<Double>.random(in: -0.2...0.2)
    for i in parameters.atoms.indices {
      originalVelocities[i] += SIMD3(linearVelocity)
    }
    
    var rigidBodyDesc = descriptor
    rigidBodyDesc.velocities = originalVelocities
    let rigidBody = MM4RigidBody(descriptor: rigidBodyDesc)
    
    // TODO: Test both methods for angular velocity. Check that the reported
    // angular momentum matches what you expect, and shape the angular momentum
    // into an angular velocity equalling 'bulkVelocity'.
    let computedVelocity = rigidBody.angularMomentum
    XCTAssertEqual(
      computedVelocity,
      (parameters.atoms.count > 0) ? bulkVelocity : .zero,
      accuracy: 1e-5)
    XCTAssertEqual(
      rigidBody.linearMomentum / rigidBody.mass,
      (parameters.atoms.count > 0) ? linearVelocity : .zero,
      accuracy: 1e-5)
    XCTAssertEqual(
      rigidBody.linearMomentum,
      (parameters.atoms.count > 0) ? linearVelocity * rigidBody.mass : .zero,
      accuracy: 1e-5)
  }
  
  // Test what happens when a zero velocity becomes a nonzero angular velocity.
  for bulkVelocity in bulkVelocities {
    var rigidBody = MM4RigidBody(descriptor: descriptor)
    XCTAssert(rigidBody.angularMomentum == .zero)
    for i in parameters.atoms.indices {
      XCTAssertEqual(rigidBody.velocities[i], .zero, accuracy: 1e-5)
    }
    
    // Fix this.
    rigidBody.angularMomentum = bulkVelocity
    
    var computedVelocity = rigidBody.angularMomentum
    XCTAssertEqual(
      computedVelocity,
      (parameters.atoms.count > 0) ? bulkVelocity : .zero,
      accuracy: 1e-5)
    
    let expectedVelocities = createVelocities(bulkVelocity)
    for i in parameters.atoms.indices {
      let position = descriptor.positions![i]
      XCTAssertEqual(position, rigidBody.positions[i], accuracy: 1e-3)
      let velocity = expectedVelocities[i]
      XCTAssertEqual(velocity, rigidBody.velocities[i], accuracy: 1e-3)
    }
    
    computedVelocity = rigidBody.angularMomentum
    XCTAssertEqual(
      computedVelocity,
      (parameters.atoms.count > 0) ? bulkVelocity : .zero,
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
  
  var rigidBody = MM4RigidBody(descriptor: descriptor)
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
    }
    
    let computedVelocity = rigidBody.angularMomentum
    if Bool.random() {
      XCTAssertEqual(
        rigidBody.linearMomentum / rigidBody.mass,
        (parameters.atoms.count > 0) ? currentLinearVelocity : .zero,
        accuracy: 1e-5)
    } else {
      let expectedVelocity = currentLinearVelocity * rigidBody.mass
      XCTAssertEqual(
        rigidBody.linearMomentum,
        (parameters.atoms.count > 0) ? expectedVelocity : .zero,
        accuracy: 1e-5)
    }
    if Bool.random() {
      XCTAssertEqual(
        computedVelocity,
        (parameters.atoms.count > 0) ? currentAngularVelocity : .zero,
        accuracy: 1e-5)
    } else {
      // TODO: Two different methods of checking angular motion.
    }
    
    // Fix this.
    for i in parameters.atoms.indices {
      let position = descriptor.positions![i]
      let delta = position - SIMD3(centerOfMass)
      var velocity: SIMD3<Float> = .zero
      velocity += cross(SIMD3(currentAngularVelocity), delta)
      velocity += SIMD3(currentLinearVelocity)
      
      let actualVelocity = rigidBody.velocities[i]
      XCTAssertEqual(velocity, actualVelocity, accuracy: 1e-5)
    }
  }
}
