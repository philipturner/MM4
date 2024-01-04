import XCTest
import MM4
import Numerics

// Test the correctness of functionality that initializes and mutates
// rigid body velocities.
//
// Edge cases to test when adding support for anchors:
// - Set an anchor to one linear velocity, when the rest of the object has an
//   entirely different linear velocity. Ensure the rest of the atoms' velocity
//   does not appear in the bulk linear velocity.

// MARK: - Test Execution

final class MM4RigidBodyVelocityTests: XCTestCase {
  func testRigidBodyVelocity() throws {
    let references = try MM4RigidBodyTests.createRigidBodyReferences()
    for reference in references {
      testLinearVelocity(reference)
      testAngularVelocity(reference)
    }
  }
}

// MARK: - Velocity

// Test that when certain velocities are entered into the object descriptor, it
// automatically recognizes the correct velocity. Then, test that a stationary
// object with its velocity mutated shows the expected atom velocities.
// Also, test what happens when both linear and angular velocity are nonzero.

private func testLinearVelocity(_ reference: MM4RigidBody) {
  let parameters = reference.parameters
  
  var bulkVelocities: [SIMD3<Float>] = []
  bulkVelocities.append(.zero)
  while true {
    // Add one random velocity with up to 200 m/s for each component.
    let random = SIMD3<Float>.random(in: -0.2...0.2)
    if any(random .!= .zero) {
      bulkVelocities.append(random)
      break
    } else {
      continue
    }
  }
  
  // Test that a pre-existing linear velocity is recognized.
  for bulkVelocity in bulkVelocities {
    let originalVelocities = [SIMD3<Float>](
      repeating: bulkVelocity, count: parameters.atoms.count)
    
    var rigidBody = MM4RigidBody(parameters: parameters)
    rigidBody.setPositions(reference.positions)
    rigidBody.setVelocities(originalVelocities)
    let computedVelocity = rigidBody.linearVelocity
    let expectedVelocity = (parameters.atoms.count > 0) ? bulkVelocity : .zero
    XCTAssertEqual(computedVelocity, expectedVelocity, accuracy: 1e-5)
  }
  
  // Test that the original linear velocity is zero, but a modified linear
  // velocity is reflected in the 'velocities' property. Ensure the property is
  // all zeroes before linear velocity changes.
  for bulkVelocity in bulkVelocities {
    var rigidBody = MM4RigidBody(parameters: parameters)
    rigidBody.setPositions(reference.positions)
    XCTAssertEqual(rigidBody.linearVelocity, .zero, accuracy: 1e-5)
    for i in parameters.atoms.indices {
      XCTAssertEqual(rigidBody.velocities[i], .zero, accuracy: 1e-5)
    }
    rigidBody.linearVelocity = bulkVelocity
    
    let expectedVelocity = (parameters.atoms.count > 0) ? bulkVelocity : .zero
    XCTAssertEqual(rigidBody.linearVelocity, expectedVelocity, accuracy: 1e-5)
    for i in parameters.atoms.indices {
      XCTAssertEqual(rigidBody.velocities[i], bulkVelocity, accuracy: 1e-5)
    }
    XCTAssertEqual(rigidBody.linearVelocity, expectedVelocity, accuracy: 1e-5)
  }
  
  // Test what happens when you change a rigid body's velocity multiple times.
  var rigidBody = MM4RigidBody(parameters: parameters)
  rigidBody.setPositions(reference.positions)
  rigidBody.setVelocities(
    Array(repeating: bulkVelocities[1], count: parameters.atoms.count))
  
  XCTAssertEqual(
    rigidBody.linearVelocity,
    (parameters.atoms.count > 0) ? bulkVelocities[1] : .zero,
    accuracy: 1e-5)
  for i in parameters.atoms.indices {
    XCTAssertEqual(rigidBody.velocities[i], bulkVelocities[1], accuracy: 1e-5)
  }
  XCTAssertEqual(
    rigidBody.linearVelocity,
    (parameters.atoms.count > 0) ? bulkVelocities[1] : .zero,
    accuracy: 1e-5)
  
  for _ in 0..<5 {
    // Generate a random 3D vector ranging from -200 m/s to 200 m/s.
    let velocity = SIMD3<Float>.random(in: -0.2...0.2)
    rigidBody.linearVelocity = velocity
    
    XCTAssertEqual(
      rigidBody.linearVelocity,
      (parameters.atoms.count > 0) ? velocity : .zero,
      accuracy: 1e-5)
    for i in parameters.atoms.indices {
      XCTAssertEqual(rigidBody.velocities[i], velocity, accuracy: 1e-5)
    }
    XCTAssertEqual(
      rigidBody.linearVelocity,
      (parameters.atoms.count > 0) ? velocity : .zero,
      accuracy: 1e-5)
  }
}

private func testAngularVelocity(_ reference: MM4RigidBody) {
  let parameters = reference.parameters
  
  func cross(_ a: SIMD3<Float>, _ b: SIMD3<Float>) -> SIMD3<Float> {
    var c: SIMD3<Float> = .zero
    c.x = a.y * b.z - a.z * b.y
    c.y = a.z * b.x - a.x * b.z
    c.z = a.x * b.y - a.y * b.x
    return c
  }
  
  var bulkVelocities: [SIMD3<Float>] = []
  bulkVelocities.append(.zero)
  while true {
    // Add one random velocity with up to ~0.020 rad/ps angular speed.
    // - This equates to 20 / (2 * .pi) GHz.
    let random = SIMD3<Float>.random(in: -0.02...0.02)
    if any(random .!= .zero) {
      bulkVelocities.append(random)
      break
    } else {
      continue
    }
  }
  
  // Create the center of mass for initializing angular velocities.
  var centerOfMass: SIMD3<Float>
  do {
    var rigidBody = MM4RigidBody(parameters: reference.parameters)
    rigidBody.setPositions(reference.positions)
    centerOfMass = rigidBody.centerOfMass
  }
  
  // Helper function for creating expected velocities.
  func createVelocities(
    _ angularVelocity: SIMD3<Float>
  ) -> [SIMD3<Float>] {
    var output: [SIMD3<Float>] = []
    for position in reference.positions {
      let delta = position - centerOfMass
      let velocity = cross(angularVelocity, delta)
      output.append(velocity)
    }
    return output
  }
  
  // Helper function for extracting a quaternion's rotation vector.
  func createRotationVector(
    _ quaternion: Quaternion<Float>
  ) -> SIMD3<Float> {
    if quaternion.angle.isNaN {
      return .zero
    } else {
      return quaternion.rotationVector
    }
  }
  
  // Test that a pre-existing angular velocity is recognized.
  for bulkVelocity in bulkVelocities {
    let originalVelocities = createVelocities(bulkVelocity)
    
    var rigidBody = MM4RigidBody(parameters: reference.parameters)
    rigidBody.setPositions(reference.positions)
    rigidBody.setVelocities(originalVelocities)
    
    let quaternion = rigidBody.angularVelocity
    let computedVelocity = createRotationVector(quaternion)
    XCTAssertEqual(
      computedVelocity,
      (parameters.atoms.count > 0) ? bulkVelocity : .zero,
      accuracy: 1e-5)
  }
  
  // Test what happens when a linear velocity is also present initially.
  for bulkVelocity in bulkVelocities {
    var originalVelocities = createVelocities(bulkVelocity)
    let linearVelocity = SIMD3<Float>.random(in: -0.2...0.2)
    for i in parameters.atoms.indices {
      originalVelocities[i] += linearVelocity
    }
    
    var rigidBody = MM4RigidBody(parameters: reference.parameters)
    rigidBody.setPositions(reference.positions)
    rigidBody.setVelocities(originalVelocities)
    
    let quaternion = rigidBody.angularVelocity
    let computedVelocity = createRotationVector(quaternion)
    XCTAssertEqual(
      computedVelocity,
      (parameters.atoms.count > 0) ? bulkVelocity : .zero,
      accuracy: 1e-5)
    XCTAssertEqual(
      rigidBody.linearVelocity,
      (parameters.atoms.count > 0) ? linearVelocity : .zero,
      accuracy: 1e-5)
  }
  
  // Test what happens when a zero velocity becomes a nonzero angular velocity.
  for bulkVelocity in bulkVelocities {
    var rigidBody = MM4RigidBody(parameters: reference.parameters)
    rigidBody.setPositions(reference.positions)
    XCTAssert(rigidBody.angularVelocity.angle.isNaN)
    for i in parameters.atoms.indices {
      XCTAssertEqual(rigidBody.velocities[i], .zero, accuracy: 1e-5)
    }
    
    let angle = (bulkVelocity * bulkVelocity).sum().squareRoot()
    let axis = (angle == 0) ? SIMD3<Float>(1, 0, 0) : bulkVelocity / angle
    rigidBody.angularVelocity = Quaternion(angle: angle, axis: axis)
    
    var quaternion = rigidBody.angularVelocity
    var computedVelocity = createRotationVector(quaternion)
    XCTAssertEqual(
      computedVelocity,
      (parameters.atoms.count > 0) ? bulkVelocity : .zero,
      accuracy: 1e-5)
    
    let expectedVelocities = createVelocities(bulkVelocity)
    for i in parameters.atoms.indices {
      let position = reference.positions[i]
      XCTAssertEqual(position, rigidBody.positions[i], accuracy: 1e-3)
      let velocity = expectedVelocities[i]
      XCTAssertEqual(velocity, rigidBody.velocities[i], accuracy: 1e-3)
    }
    
    quaternion = rigidBody.angularVelocity
    computedVelocity = createRotationVector(quaternion)
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
  
  var rigidBody = MM4RigidBody(parameters: reference.parameters)
  rigidBody.setPositions(reference.positions)
  var currentAngularVelocity: SIMD3<Float> = .zero
  var currentLinearVelocity: SIMD3<Float> = .zero
  for (changeLinear, changeAngular) in testOrder {
    if changeLinear {
      currentLinearVelocity = .random(in: -0.2...0.2)
      rigidBody.linearVelocity = currentLinearVelocity
    }
    if changeAngular {
      currentAngularVelocity = .random(in: -0.020...0.020)
      
      let velocity = currentAngularVelocity
      let angle = (velocity * velocity).sum().squareRoot()
      let axis = (angle == 0) ? SIMD3<Float>(1, 0, 0) : velocity / angle
      rigidBody.angularVelocity = Quaternion(angle: angle, axis: axis)
    }
    
    let quaternion = rigidBody.angularVelocity
    let computedVelocity = createRotationVector(quaternion)
    XCTAssertEqual(
      rigidBody.linearVelocity,
      (parameters.atoms.count > 0) ? currentLinearVelocity : .zero,
      accuracy: 1e-5)
    XCTAssertEqual(
      computedVelocity,
      (parameters.atoms.count > 0) ? currentAngularVelocity : .zero,
      accuracy: 1e-5)
    
    for i in parameters.atoms.indices {
      let position = reference.positions[i]
      let delta = position - centerOfMass
      var velocity: SIMD3<Float> = .zero
      velocity += cross(currentAngularVelocity, delta)
      velocity += currentLinearVelocity
      
      let actualVelocity = rigidBody.velocities[i]
      XCTAssertEqual(velocity, actualVelocity, accuracy: 1e-5)
    }
  }
}


