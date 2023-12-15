import XCTest
import MM4

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
  XCTAssertEqual(mass, rigidBody.mass, accuracy: 1e-3)
  
  // TODO: continue
  // https://gist.github.com/philipturner/b07137b478bac692eceaa079f4993f32
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
  // derive the mass
  fatalError("Not implemented.")
}

private func deriveMomentOfInertia(
  _ rigidBody: MM4RigidBody
) -> (SIMD3<Float>, SIMD3<Float>, SIMD3<Float>) {
  // derive the center of mass
  fatalError("Not implemented.")
}

private func invertMomentOfInertia(
  _ momentOfInertia: (SIMD3<Float>, SIMD3<Float>, SIMD3<Float>)
) -> (SIMD3<Float>, SIMD3<Float>, SIMD3<Float>) {
  fatalError("Not implemented.")
}

// MARK: - Velocity

// MARK: - Mutation
