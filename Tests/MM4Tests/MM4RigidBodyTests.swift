import XCTest
import MM4

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
    
    var rigidBodyDesc = MM4RigidBodyDescriptor()
    rigidBodyDesc.parameters = try GoldTests.createGoldParameters()
    rigidBodyDesc.positions = positions
    
    let rigidBody = MM4RigidBody(descriptor: rigidBodyDesc)
    try testRigidBody(rigidBody, positions: positions)
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
  
  let rigidBody = MM4RigidBody(descriptor: rigidBodyDesc)
  try testRigidBody(rigidBody, positions: adamantane.positions)
}

private func testRigidBody(
  _ rigidBody: MM4RigidBody, positions: [SIMD3<Float>]
) throws {
  XCTAssertEqual(rigidBody.positions, positions)
  XCTAssertEqual(rigidBody.velocities.count, positions.count)
  XCTAssert(rigidBody.velocities.allSatisfy { $0 == .zero })
  
  testInertia(rigidBody)
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
  
  func compareVectors(_ x: SIMD3<Float>, _ y: SIMD3<Float>, _ message: String) {
    XCTAssertEqual(x.x, y.x, accuracy: 1e-3, message)
    XCTAssertEqual(x.y, y.y, accuracy: 1e-3, message)
    XCTAssertEqual(x.z, y.z, accuracy: 1e-3, message)
  }
  
  let centerOfMass = deriveCenterOfMass(rigidBody)
  compareVectors(centerOfMass, rigidBody.centerOfMass, "center of mass")
  
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
  fatalError("Not implemented.")
}

// MARK: - Velocity

// Test that when certain velocities are entered into the object descriptor, it
// automatically recognizes the correct velocity. Then, test that a stationary
// object with its velocity mutated shows the expected atom velocities.
