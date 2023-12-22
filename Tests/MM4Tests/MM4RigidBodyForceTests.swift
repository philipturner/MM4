import XCTest
import MM4

// Test all functionality related to forced motions - anchors, handles, and
// external forces. This includes every edge case that calculates inertia and
// velocity differently.
//
// Edge cases:
// - Set an anchor to one linear velocity, when the rest of the object has an
//   entirely different linear velocity. Ensure the rest of the atoms' velocity
//   does not appear in the bulk linear velocity.

// MARK: - Test Execution

final class MM4RigidBodyForceTests: XCTestCase {
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
    
    testRigidBodyForce(descriptor: rigidBodyDesc)
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
  
  testRigidBodyForce(descriptor: rigidBodyDesc)
}

private func testRigidBodyForce(descriptor: MM4RigidBodyDescriptor) {
  
}
