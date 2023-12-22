import XCTest
import MM4

// Test all functionality related to forced motions - anchors, handles, and
// external forces. This includes every edge case that calculates inertia and
// velocity differently.
//
// NOTE: This part of the test suite will be incomplete until MM4ForceField is
// fully functioning. We need to define exactly how external forces are
// imported/exported to the rigid body.
//
// Edge cases:
// - Sort the atoms by farthest distance from the CoM, then choose some
//   interesting ones for the anchors.
// - Set an anchor to one linear velocity, when the rest of the object has an
//   entirely different linear velocity. Ensure the rest of the atoms' velocity
//   does not appear in the bulk linear velocity.
// - Test combinations where anchors and handles collide. You may need to add
//   some MM4 errors for invalid anchors and/or invalid handles.
// - Test what happens when every atom is an anchor.

// MARK: - Test Execution

final class MM4RigidBodyForceTests: XCTestCase {
  func testRigidBodyForce() throws {
    let descriptors = try MM4RigidBodyTests.createRigidBodyDescriptors()
    for descriptor in descriptors {
      _ = descriptor
    }
  }
}
