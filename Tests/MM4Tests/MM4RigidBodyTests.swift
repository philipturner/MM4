import XCTest
import MM4
import Numerics

// Test the basic functionality of the rigid body API, and bulk properties
// derived from rigid body mechanics.

// MARK: - Test Execution

final class MM4RigidBodyTests: XCTestCase {
  func testAdamantane() throws {
    try testAdamantaneVariant(atomCode: .alkaneCarbon)
  }
  
  func testSilaAdamantane() throws {
    try testAdamantaneVariant(atomCode: .silicon)
  }
  
  func testGold() throws {
    let positions = GoldTests.createGoldPositions()
    
    // Although it is good practice to make the gold atoms all be anchors, that
    // would be cross-coverage for this test case.
    var rigidBodyDesc = MM4RigidBodyDescriptor()
    rigidBodyDesc.parameters = try GoldTests.createGoldParameters()
    rigidBodyDesc.positions = positions
    
    try testRigidBody(descriptor: rigidBodyDesc)
  }
}

private func testAdamantaneVariant(atomCode: MM4AtomCode) throws {
  let adamantane = Adamantane(atomCode: atomCode)
  
  var paramsDesc = MM4ParametersDescriptor()
  paramsDesc.atomicNumbers = adamantane.atomicNumbers
  paramsDesc.bonds = adamantane.bonds
  
  var rigidBodyDesc = MM4RigidBodyDescriptor()
  rigidBodyDesc.parameters = try MM4Parameters(descriptor: paramsDesc)
  rigidBodyDesc.positions = adamantane.positions
  
  try testRigidBody(descriptor: rigidBodyDesc)
}

private func testRigidBody(descriptor: MM4RigidBodyDescriptor) throws {
  // Test whether the initial velocities are zero. Also, test whether the
  // positions and velocities are different after the rigid body is changed with
  // 'setPositions' and 'setVelocities'.
  var rigidBody = MM4RigidBody(descriptor: descriptor)
  XCTAssertEqual(rigidBody.positions, descriptor.positions!)
  XCTAssertEqual(rigidBody.velocities.count, descriptor.positions!.count)
  XCTAssert(rigidBody.velocities.allSatisfy { $0 == .zero })
  testInertia(rigidBody)
  
  let originalCenterOfMass = rigidBody.centerOfMass
  
  // Ensure all four permutations of changing positions and velocities result
  // in the expected state.
  let zeroPositions = [SIMD3<Float>](
    repeating: .zero, count: rigidBody.parameters.atoms.count)
  let zeroVelocities = [SIMD3<Float>](
    repeating: .zero, count: rigidBody.parameters.atoms.count)
  var newPositions = zeroPositions
  var newVelocities = zeroVelocities
  for i in rigidBody.parameters.atoms.indices {
    newPositions[i] = .random(in: -2...2)
    newVelocities[i] = .random(in: -0.1...0.1)
  }
  
  // Pseudo-randomly switch between the two APIs for setting data.
  var setMethodCounter = 0
  func setPositions(_ array: [SIMD3<Float>]) {
    setMethodCounter += 1
    if setMethodCounter % 3 == 0 {
      rigidBody.setPositions(array)
    } else {
      array.withUnsafeBufferPointer {
        rigidBody.setPositions($0)
      }
    }
  }
  func setVelocities(_ array: [SIMD3<Float>]) {
    setMethodCounter += 1
    if setMethodCounter % 3 == 0 {
      rigidBody.setVelocities(array)
    } else {
      array.withUnsafeBufferPointer {
        rigidBody.setVelocities($0)
      }
    }
  }
  
  for changePositions in [false, true] {
    for changeVelocities in [false, true] {
      // Reset the state to zero, to avoid interference from previous trials.
      setPositions(zeroPositions)
      setVelocities(zeroVelocities)
      XCTAssert(rigidBody.positions.allSatisfy { $0 == .zero })
      XCTAssert(rigidBody.velocities.allSatisfy { $0 == .zero })
      
      // Change the object based on the current permutation of state changes.
      if changePositions {
        setPositions(newPositions)
      }
      if changeVelocities {
        setVelocities(newVelocities)
      }
      
      // Test both the position/velocity properties and low-level getters into
      // copies of the arrays.
      if changePositions {
        XCTAssertEqual(newPositions, rigidBody.positions)
      }
      if changeVelocities {
        XCTAssertEqual(newVelocities, rigidBody.velocities)
      }
      let expectedPositions = changePositions ? newPositions : zeroPositions
      let expectedVelocities = changeVelocities ? newVelocities : zeroVelocities
      var outputPositions = zeroPositions
      var outputVelocities = zeroVelocities
      outputPositions.withUnsafeMutableBufferPointer {
        rigidBody.getPositions($0)
      }
      outputVelocities.withUnsafeMutableBufferPointer {
        rigidBody.getVelocities($0)
      }
      XCTAssertEqual(outputPositions, expectedPositions)
      XCTAssertEqual(outputVelocities, expectedVelocities)
    }
  }
  
  // Repeat the inertia tests after changing the positions and velocities. If
  // random initialization doesn't change the CoM, the tests fail to cover the
  // functionality.
  XCTAssert(
    rigidBody.centerOfMass != originalCenterOfMass,
    "tests did not cover functionality")
  testInertia(rigidBody)
  
  // Run the bulk velocity tests on similar rigid bodies created with the same
  // descriptor.
  testLinearVelocity(descriptor)
  testAngularVelocity(descriptor)
  testCoW(descriptor)
}

// MARK: - Inertia

// This function does not validate whether the inertia is correct after
// modification. Another function should mutate some inertia properties and
// determine whether they behave correctly.
private func testInertia(_ rigidBody: MM4RigidBody) {
  // Assert that all quantities are equal within 1e-3 tolerance. This is an
  // internally consistent unit system where the fundamental quantities are
  // roughly 1. Therefore, absolute errors should be on the same order as what
  // you expect for relative errors when using a multiplicative factor. Using
  // absolute errors enables use of more ergonomic testing functions.
  
  let mass = deriveMass(rigidBody)
  XCTAssertEqual(mass, rigidBody.mass, accuracy: 1e-3, "mass")
  
  // Floating-point error starts to make the moment of inertia tests flaky for
  // gold, as the absolute value of the numbers is so large. Therefore, we
  // change the accuracy with a variable declared outside 'compareVectors'.
  var accuracy: Float = 1e-3
  func compareVectors(_ x: SIMD3<Float>, _ y: SIMD3<Float>, _ message: String) {
    XCTAssertEqual(x.x, y.x, accuracy: accuracy, message)
    XCTAssertEqual(x.y, y.y, accuracy: accuracy, message)
    XCTAssertEqual(x.z, y.z, accuracy: accuracy, message)
  }
  
  let centerOfMass = deriveCenterOfMass(rigidBody)
  compareVectors(centerOfMass, rigidBody.centerOfMass, "center of mass")
  
  if rigidBody.parameters.atoms.codes.contains(.gold) {
    accuracy = 1e-2
  }
  let momentOfInertia = deriveMomentOfInertia(rigidBody)
  compareVectors(
    momentOfInertia.0, rigidBody.momentOfInertia.0, "moment of inertia")
  compareVectors(
    momentOfInertia.1, rigidBody.momentOfInertia.1, "moment of inertia")
  compareVectors(
    momentOfInertia.2, rigidBody.momentOfInertia.2, "moment of inertia")
}

private func deriveMass(_ rigidBody: MM4RigidBody) -> Double {
  var output: Double = .zero
  for i in rigidBody.parameters.atoms.indices {
    output += Double(rigidBody.parameters.atoms.masses[i])
  }
  return output
}

private func deriveCenterOfMass(
  _ rigidBody: MM4RigidBody
) -> SIMD3<Float> {
  var output: SIMD3<Double> = .zero
  for i in rigidBody.parameters.atoms.indices {
    let position = rigidBody.positions[i]
    let mass = rigidBody.parameters.atoms.masses[i]
    output += SIMD3<Double>(mass * position)
  }
  output /= deriveMass(rigidBody)
  return SIMD3<Float>(output)
}

private func deriveMomentOfInertia(
  _ rigidBody: MM4RigidBody
) -> (SIMD3<Float>, SIMD3<Float>, SIMD3<Float>) {
  let centerOfMass = deriveCenterOfMass(rigidBody)
  var output: (
    SIMD3<Double>, SIMD3<Double>, SIMD3<Double>
  ) = (.zero, .zero, .zero)
  
  for i in rigidBody.parameters.atoms.indices {
    let delta = rigidBody.positions[i] - centerOfMass
    let mass = rigidBody.parameters.atoms.masses[i]
    let STS = (delta * delta).sum()
    output.0[0] += Double(mass * STS)
    output.1[1] += Double(mass * STS)
    output.2[2] += Double(mass * STS)
    
    output.0 -= SIMD3<Double>(mass * delta * delta.x)
    output.1 -= SIMD3<Double>(mass * delta * delta.y)
    output.2 -= SIMD3<Double>(mass * delta * delta.z)
  }
  return (
    SIMD3<Float>(output.0),
    SIMD3<Float>(output.1),
    SIMD3<Float>(output.2)
  )
}

private func invertMomentOfInertia(
  _ momentOfInertia: (SIMD3<Float>, SIMD3<Float>, SIMD3<Float>)
) -> (SIMD3<Float>, SIMD3<Float>, SIMD3<Float>) {
  let col = (
    SIMD3<Double>(momentOfInertia.0),
    SIMD3<Double>(momentOfInertia.1),
    SIMD3<Double>(momentOfInertia.2)
  )
  let determinant =
  col.0[0] * (col.1[1] * col.2[2] - col.2[1] * col.1[2]) -
  col.0[1] * (col.1[0] * col.2[2] - col.1[2] * col.2[0]) +
  col.0[2] * (col.1[0] * col.2[1] - col.1[1] * col.2[0])
  let invdet = 1 / determinant
  
  let result00 = (col.1[1] * col.2[2] - col.2[1] * col.1[2]) * invdet
  let result01 = (col.0[2] * col.2[1] - col.0[1] * col.2[2]) * invdet
  let result02 = (col.0[1] * col.1[2] - col.0[2] * col.1[1]) * invdet
  
  let result10 = (col.1[2] * col.2[0] - col.1[0] * col.2[2]) * invdet
  let result11 = (col.0[0] * col.2[2] - col.0[2] * col.2[0]) * invdet
  let result12 = (col.1[0] * col.0[2] - col.0[0] * col.1[2]) * invdet
  
  let result20 = (col.1[0] * col.2[1] - col.2[0] * col.1[1]) * invdet
  let result21 = (col.2[0] * col.0[1] - col.0[0] * col.2[1]) * invdet
  let result22 = (col.0[0] * col.1[1] - col.1[0] * col.0[1]) * invdet
  
  let column0 = SIMD3(result00, result10, result20)
  let column1 = SIMD3(result01, result11, result21)
  let column2 = SIMD3(result02, result12, result22)
  return (
    SIMD3<Float>(column0),
    SIMD3<Float>(column1),
    SIMD3<Float>(column2)
  )
}

// MARK: - Velocity

// Test that when certain velocities are entered into the object descriptor, it
// automatically recognizes the correct velocity. Then, test that a stationary
// object with its velocity mutated shows the expected atom velocities.
// Also, test what happens when both linear and angular velocity are nonzero.

private func testLinearVelocity(_ descriptor: MM4RigidBodyDescriptor) {
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
    let parameters = descriptor.parameters!
    let originalVelocities = [SIMD3<Float>](
      repeating: bulkVelocity, count: parameters.atoms.count)
    
    var desc = descriptor
    desc.velocities = originalVelocities
    let rigidBody = MM4RigidBody(descriptor: desc)
    let computedVelocity = rigidBody.linearVelocity
    XCTAssertEqual(computedVelocity, bulkVelocity, accuracy: 1e-5)
  }
  
  // Test that the original linear velocity is zero, but a modified linear
  // velocity is reflected in the 'velocities' property. Ensure the property is
  // all zeroes before linear velocity changes.
  for bulkVelocity in bulkVelocities {
    var rigidBody = MM4RigidBody(descriptor: descriptor)
    XCTAssertEqual(rigidBody.linearVelocity, .zero, accuracy: 1e-5)
    for i in descriptor.parameters!.atoms.indices {
      XCTAssertEqual(rigidBody.velocities[i], .zero, accuracy: 1e-5)
    }
    
    rigidBody.linearVelocity = bulkVelocity
    XCTAssertEqual(rigidBody.linearVelocity, bulkVelocity, accuracy: 1e-5)
    for i in descriptor.parameters!.atoms.indices {
      XCTAssertEqual(rigidBody.velocities[i], bulkVelocity, accuracy: 1e-5)
    }
    XCTAssertEqual(rigidBody.linearVelocity, bulkVelocity, accuracy: 1e-5)
  }
  
  // Test what happens when you change a rigid body's velocity multiple times.
  var desc = descriptor
  let parameters = descriptor.parameters!
  desc.velocities = Array(
    repeating: bulkVelocities[1], count: parameters.atoms.count)
  
  var rigidBody = MM4RigidBody(descriptor: desc)
  XCTAssertEqual(rigidBody.linearVelocity, bulkVelocities[1], accuracy: 1e-5)
  for i in parameters.atoms.indices {
    XCTAssertEqual(rigidBody.velocities[i], bulkVelocities[1], accuracy: 1e-5)
  }
  XCTAssertEqual(rigidBody.linearVelocity, bulkVelocities[1], accuracy: 1e-5)
  
  for _ in 0..<5 {
    // Generate a random 3D vector ranging from -200 m/s to 200 m/s.
    let velocity = SIMD3<Float>.random(in: -0.2...0.2)
    rigidBody.linearVelocity = velocity
    XCTAssertEqual(rigidBody.linearVelocity, velocity, accuracy: 1e-5)
    for i in parameters.atoms.indices {
      XCTAssertEqual(rigidBody.velocities[i], velocity, accuracy: 1e-5)
    }
    XCTAssertEqual(rigidBody.linearVelocity, velocity, accuracy: 1e-5)
  }
}

private func testAngularVelocity(_ descriptor: MM4RigidBodyDescriptor) {
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
    let rigidBody = MM4RigidBody(descriptor: descriptor)
    centerOfMass = rigidBody.centerOfMass
  }
  
  // Helper function for creating expected velocities.
  func createVelocities(
    _ angularVelocity: SIMD3<Float>
  ) -> [SIMD3<Float>] {
    var output: [SIMD3<Float>] = []
    for position in descriptor.positions! {
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
    
    var desc = descriptor
    desc.velocities = originalVelocities
    let rigidBody = MM4RigidBody(descriptor: desc)
    let quaternion = rigidBody.angularVelocity
    let computedVelocity = createRotationVector(quaternion)
    XCTAssertEqual(computedVelocity, bulkVelocity, accuracy: 1e-5)
  }
  
  // Test what happens when a linear velocity is also present initially.
  for bulkVelocity in bulkVelocities {
    let parameters = descriptor.parameters!
    var originalVelocities = createVelocities(bulkVelocity)
    let linearVelocity = SIMD3<Float>.random(in: -0.2...0.2)
    for i in parameters.atoms.indices {
      originalVelocities[i] += linearVelocity
    }
    
    var desc = descriptor
    desc.velocities = originalVelocities
    let rigidBody = MM4RigidBody(descriptor: desc)
    let quaternion = rigidBody.angularVelocity
    let computedVelocity = createRotationVector(quaternion)
    XCTAssertEqual(computedVelocity, bulkVelocity, accuracy: 1e-5)
    XCTAssertEqual(rigidBody.linearVelocity, linearVelocity, accuracy: 1e-5)
  }
  
  // Test what happens when a zero velocity becomes a nonzero angular velocity.
  for bulkVelocity in bulkVelocities {
    var rigidBody = MM4RigidBody(descriptor: descriptor)
    XCTAssert(rigidBody.angularVelocity.angle.isNaN)
    for i in descriptor.parameters!.atoms.indices {
      XCTAssertEqual(rigidBody.velocities[i], .zero, accuracy: 1e-5)
    }
    
    let angle = (bulkVelocity * bulkVelocity).sum().squareRoot()
    let axis = (angle == 0) ? SIMD3<Float>(1, 0, 0) : bulkVelocity / angle
    rigidBody.angularVelocity = Quaternion(angle: angle, axis: axis)
    
    var quaternion = rigidBody.angularVelocity
    var computedVelocity = createRotationVector(quaternion)
    XCTAssertEqual(computedVelocity, bulkVelocity, accuracy: 1e-5)
    
    let expectedVelocities = createVelocities(bulkVelocity)
    for i in descriptor.parameters!.atoms.indices {
      let position = descriptor.positions![i]
      XCTAssertEqual(position, rigidBody.positions[i], accuracy: 1e-3)
      let velocity = expectedVelocities[i]
      XCTAssertEqual(velocity, rigidBody.velocities[i], accuracy: 1e-3)
    }
    
    quaternion = rigidBody.angularVelocity
    computedVelocity = createRotationVector(quaternion)
    XCTAssertEqual(computedVelocity, bulkVelocity, accuracy: 1e-5)
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
      rigidBody.linearVelocity, currentLinearVelocity, accuracy: 1e-5)
    XCTAssertEqual(computedVelocity, currentAngularVelocity, accuracy: 1e-5)
    
    for i in descriptor.parameters!.atoms.indices {
      let position = descriptor.positions![i]
      let delta = position - centerOfMass
      var velocity: SIMD3<Float> = .zero
      velocity += cross(currentAngularVelocity, delta)
      velocity += currentLinearVelocity
      
      let actualVelocity = rigidBody.velocities[i]
      XCTAssertEqual(velocity, actualVelocity, accuracy: 1e-5)
    }
  }
}

// MARK: - Copy on Write

private func testCoW(_ descriptor: MM4RigidBodyDescriptor) {
  var rigidBody1 = MM4RigidBody(descriptor: descriptor)
  var rigidBody2 = rigidBody1
  
  let centerOfMass1 = rigidBody1.centerOfMass
  XCTAssertEqual(centerOfMass1, rigidBody1.centerOfMass, accuracy: 1e-3)
  XCTAssertEqual(centerOfMass1, rigidBody2.centerOfMass, accuracy: 1e-3)
  XCTAssertEqual(rigidBody1.centerOfMass, rigidBody2.centerOfMass)
  
  let delta2 = SIMD3<Float>.random(in: -2...2)
  rigidBody2.centerOfMass += delta2
  let centerOfMass2 = centerOfMass1 + delta2
  XCTAssertEqual(centerOfMass1, rigidBody1.centerOfMass, accuracy: 1e-3)
  XCTAssertEqual(centerOfMass2, rigidBody2.centerOfMass, accuracy: 1e-3)
  XCTAssertNotEqual(centerOfMass1, centerOfMass2)
  XCTAssertNotEqual(rigidBody1.centerOfMass, rigidBody2.centerOfMass)
  
  let rigidBody3 = rigidBody1
  let velocity1 = SIMD3<Float>.random(in: -0.2...0.2)
  rigidBody1.linearVelocity = velocity1
  XCTAssertEqual(velocity1, rigidBody1.linearVelocity, accuracy: 1e-3)
  XCTAssertEqual(.zero, rigidBody2.linearVelocity, accuracy: 1e-3)
  XCTAssertEqual(.zero, rigidBody3.linearVelocity, accuracy: 1e-3)
  XCTAssertNotEqual(rigidBody1.linearVelocity, rigidBody2.linearVelocity)
  XCTAssertNotEqual(rigidBody1.linearVelocity, rigidBody3.linearVelocity)
  XCTAssertEqual(rigidBody2.linearVelocity, rigidBody3.linearVelocity)
  
  XCTAssertEqual(rigidBody1.centerOfMass, rigidBody3.centerOfMass)
  XCTAssertNotEqual(rigidBody2.centerOfMass, rigidBody3.centerOfMass)
}
