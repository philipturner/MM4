import XCTest
import MM4

// Test the basic functionality of the rigid body API.

// Edge cases to test when adding support for anchors:
// - Sort the atoms by farthest distance from the CoM, then choose some
//   interesting ones for the anchors. Test that the inertia conforms to the
//   expected rules regarding anchors.
// - Test what happens when every atom is an anchor.
// - Ensure there is a fatal error when an anchor has an external force.
//   This will not be permanently part of the test suite, but the check should
//   still be performed.

// MARK: - Test Execution

final class MM4RigidBodyTests: XCTestCase {
  func testRigidBody() throws {
    let references = try MM4RigidBodyTests.createRigidBodyReferences()
    for reference in references {
      // Test whether the initial velocities are zero. Also, test whether the
      // positions and velocities are different after the rigid body is changed
      // with 'setPositions' and 'setVelocities'.
      var rigidBody = MM4RigidBody(parameters: reference.parameters)
      rigidBody.setPositions(reference.positions)
      XCTAssertEqual(rigidBody.positions, reference.positions)
      XCTAssertEqual(rigidBody.velocities.count, reference.positions.count)
      XCTAssert(rigidBody.velocities.allSatisfy { $0 == .zero })
      testInertia(rigidBody)
      
      let originalCenterOfMass = rigidBody.centerOfMass
      
      // Ensure all four permutations of changing positions and velocities
      // result in the expected state.
      let zeroExternalForces = [SIMD3<Float>](
        repeating: .zero, count: rigidBody.parameters.atoms.count)
      let zeroPositions = [SIMD3<Float>](
        repeating: .zero, count: rigidBody.parameters.atoms.count)
      let zeroVelocities = [SIMD3<Float>](
        repeating: .zero, count: rigidBody.parameters.atoms.count)
      var newExternalForces = zeroExternalForces
      var newPositions = zeroPositions
      var newVelocities = zeroVelocities
      for i in rigidBody.parameters.atoms.indices {
        newExternalForces[i] = .random(in: -0.1...0.1)
        newPositions[i] = .random(in: -2...2)
        newVelocities[i] = .random(in: -0.1...0.1)
      }
      
      // Pseudo-randomly switch between the two APIs for setting data.
      var setMethodCounter = 0
      func setExternalForces(_ array: [SIMD3<Float>]) {
        setMethodCounter += 1
        if setMethodCounter % 5 == 0 {
          rigidBody.externalForces = array
        } else {
          rigidBody.externalForces = array
        }
      }
      func setPositions(_ array: [SIMD3<Float>]) {
        setMethodCounter += 1
        if setMethodCounter % 5 == 0 {
          rigidBody.setPositions(array)
        } else {
          array.withUnsafeBufferPointer {
            rigidBody.setPositions($0)
          }
        }
      }
      func setVelocities(_ array: [SIMD3<Float>]) {
        setMethodCounter += 1
        if setMethodCounter % 5 == 0 {
          rigidBody.setVelocities(array)
        } else {
          array.withUnsafeBufferPointer {
            rigidBody.setVelocities($0)
          }
        }
      }
      
      for changeExternalForces in [false, true] {
        for changePositions in [false, true] {
          for changeVelocities in [false, true] {
            // Reset the state to zero, to avoid interference from previous
            // trials.
            setExternalForces(zeroExternalForces)
            setPositions(zeroPositions)
            setVelocities(zeroVelocities)
            XCTAssert(rigidBody.externalForces.allSatisfy { $0 == .zero })
            XCTAssert(rigidBody.positions.allSatisfy { $0 == .zero })
            XCTAssert(rigidBody.velocities.allSatisfy { $0 == .zero })
            
            // Change the object based on the current permutation of state
            // changes.
            if changeExternalForces {
              setExternalForces(newExternalForces)
            }
            if changePositions {
              setPositions(newPositions)
            }
            if changeVelocities {
              setVelocities(newVelocities)
            }
            
            // Test both the position/velocity properties and low-level getters
            // into copies of the arrays.
            if changeExternalForces {
              XCTAssertEqual(newExternalForces, rigidBody.externalForces)
            }
            if changePositions {
              XCTAssertEqual(newPositions, rigidBody.positions)
            }
            if changeVelocities {
              XCTAssertEqual(newVelocities, rigidBody.velocities)
            }
            let expectedExternalForces =
            changeExternalForces ? newExternalForces : zeroExternalForces
            let expectedPositions =
            changePositions ? newPositions : zeroPositions
            let expectedVelocities =
            changeVelocities ? newVelocities : zeroVelocities
            
            var outputExternalForces = zeroExternalForces
            var outputPositions = zeroPositions
            var outputVelocities = zeroVelocities
            outputExternalForces = rigidBody.externalForces
            outputPositions.withUnsafeMutableBufferPointer {
              rigidBody.getPositions($0)
            }
            outputVelocities.withUnsafeMutableBufferPointer {
              rigidBody.getVelocities($0)
            }
            XCTAssertEqual(outputExternalForces, expectedExternalForces)
            XCTAssertEqual(outputPositions, expectedPositions)
            XCTAssertEqual(outputVelocities, expectedVelocities)
          }
        }
      }
      
      // Repeat the inertia tests after changing the positions and velocities.
      // If random initialization doesn't change the CoM, the tests fail to
      // cover the functionality.
      if rigidBody.parameters.atoms.count > 0 {
        XCTAssert(
          rigidBody.centerOfMass != originalCenterOfMass,
          "tests did not cover functionality")
      }
      testInertia(rigidBody)
      testCoW(reference)
    }
  }
}

extension MM4RigidBodyTests {
  static func createRigidBodyReferences() throws -> [MM4RigidBody] {
    var output: [MM4RigidBody] = []
    for atomCode in [MM4AtomCode.alkaneCarbon, .silicon] {
      let adamantane = Adamantane(atomCode: atomCode)
      
      var paramsDesc = MM4ParametersDescriptor()
      paramsDesc.atomicNumbers = adamantane.atomicNumbers
      paramsDesc.bonds = adamantane.bonds
      let params = try MM4Parameters(descriptor: paramsDesc)
      
      var rigidBody = MM4RigidBody(parameters: params)
      rigidBody.setPositions(adamantane.positions)
      output.append(rigidBody)
    }
    
    do {
      var paramsDesc = MM4ParametersDescriptor()
      paramsDesc.atomicNumbers = []
      paramsDesc.bonds = []
      let params = try MM4Parameters(descriptor: paramsDesc)
      
      var rigidBody = MM4RigidBody(parameters: params)
      rigidBody.setPositions(Array<SIMD3<Float>>())
      output.append(rigidBody)
    }
    return output
  }
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
  XCTAssertEqual(Float(mass), rigidBody.mass, accuracy: 1e-3)
  
  let centerOfMass = deriveCenterOfMass(rigidBody)
  XCTAssertEqual(centerOfMass, rigidBody.centerOfMass, accuracy: 1e-3)
  
  let momentOfInertia = deriveMomentOfInertia(rigidBody)
  XCTAssertEqual(momentOfInertia.0, rigidBody.momentOfInertia.0, accuracy: 1e-3)
  XCTAssertEqual(momentOfInertia.1, rigidBody.momentOfInertia.1, accuracy: 1e-3)
  XCTAssertEqual(momentOfInertia.2, rigidBody.momentOfInertia.2, accuracy: 1e-3)
}

// MARK: - Copy on Write

private func testCoW(_ reference: MM4RigidBody) {
  var rigidBody1 = MM4RigidBody(parameters: reference.parameters)
  rigidBody1.setPositions(reference.positions)
  var rigidBody2 = rigidBody1
  let parameters = rigidBody1.parameters
  
  let centerOfMass1 = rigidBody1.centerOfMass
  XCTAssertEqual(centerOfMass1, rigidBody1.centerOfMass, accuracy: 1e-3)
  XCTAssertEqual(centerOfMass1, rigidBody2.centerOfMass, accuracy: 1e-3)
  XCTAssertEqual(rigidBody1.centerOfMass, rigidBody2.centerOfMass)
  
  let delta2 = SIMD3<Float>.random(in: -2...2)
  rigidBody2.centerOfMass += delta2
  let centerOfMass2 = centerOfMass1 + delta2
  XCTAssertEqual(centerOfMass1, rigidBody1.centerOfMass, accuracy: 1e-3)
  XCTAssertEqual(
    (parameters.atoms.count > 0) ? centerOfMass2 : .zero,
    rigidBody2.centerOfMass,
    accuracy: 1e-3)
  XCTAssertNotEqual(centerOfMass1, centerOfMass2)
  if parameters.atoms.count > 0 {
    XCTAssertNotEqual(rigidBody1.centerOfMass, rigidBody2.centerOfMass)
  }
  
  let rigidBody3 = rigidBody1
  let velocity1 = SIMD3<Float>.random(in: -0.2...0.2)
  rigidBody1.linearVelocity = velocity1
  XCTAssertEqual(
    (parameters.atoms.count > 0) ? velocity1 : .zero,
    rigidBody1.linearVelocity,
    accuracy: 1e-3)
  XCTAssertEqual(.zero, rigidBody2.linearVelocity, accuracy: 1e-3)
  XCTAssertEqual(.zero, rigidBody3.linearVelocity, accuracy: 1e-3)
  if parameters.atoms.count > 0 {
    XCTAssertNotEqual(rigidBody1.linearVelocity, rigidBody2.linearVelocity)
    XCTAssertNotEqual(rigidBody1.linearVelocity, rigidBody3.linearVelocity)
  }
  XCTAssertEqual(rigidBody2.linearVelocity, rigidBody3.linearVelocity)
  
  XCTAssertEqual(rigidBody1.centerOfMass, rigidBody3.centerOfMass)
  if parameters.atoms.count > 0 {
    XCTAssertNotEqual(rigidBody2.centerOfMass, rigidBody3.centerOfMass)
  }
  for i in parameters.atoms.indices {
    let original = reference.positions[i]
    let modified = original + delta2
    XCTAssertEqual(rigidBody1.positions[i], original, accuracy: 1e-3)
    XCTAssertEqual(rigidBody2.positions[i], modified, accuracy: 1e-3)
    XCTAssertEqual(rigidBody3.positions[i], original, accuracy: 1e-3)
  }
  
  var rigidBody4 = rigidBody3
  let velocity4 = SIMD3<Float>.random(in: -0.2...0.2)
  let velocities4 = Array(repeating: velocity4, count: parameters.atoms.count)
  rigidBody4.setVelocities(velocities4)
  for i in parameters.atoms.indices {
    XCTAssertEqual(rigidBody3.velocities[i], .zero, accuracy: 1e-5)
    XCTAssertEqual(rigidBody4.velocities[i], velocity4, accuracy: 1e-5)
  }
  XCTAssertEqual(.zero, rigidBody3.linearVelocity, accuracy: 1e-5)
  XCTAssertEqual(
    (parameters.atoms.count > 0) ? velocity4 : .zero,
    rigidBody4.linearVelocity,
    accuracy: 1e-5)
}
