import XCTest
import MM4

// TODO: Most equality tests will likely fail. You need to assert they're equal
// to within a specific tolerance. But get the test suite compiling before you
// add more complexity.
//
// In addition, the tests for the empty rigid body will likely fail, because the
// inertia tensor cannot be diagonalized.

// MARK: - Test Execution

final class MM4RigidBodyTests: XCTestCase {
  #if false
  func testRigidBody() throws {
    let descriptors = try MM4RigidBodyTests.createDescriptors()
    for descriptor in descriptors {
      testCoW(descriptor)
      
      // Test whether the initial velocities are zero.
      var rigidBody = MM4RigidBody(descriptor: descriptor)
      XCTAssertEqual(rigidBody.positions, descriptor.positions!)
      XCTAssertEqual(rigidBody.velocities.count, descriptor.positions!.count)
      XCTAssert(rigidBody.velocities.allSatisfy { $0 == .zero })
      testInertia(rigidBody)
      
      rigidBody.forces = nil
      XCTAssertNil(rigidBody.forces)
      
      for zeroForces in [false, true] {
        var forces: [SIMD3<Float>] = []
        for _ in rigidBody.parameters.atoms.indices {
          if zeroForces {
            forces.append(.zero)
          } else {
            forces.append(.random(in: -0.1...0.1))
          }
        }
        rigidBody.forces = forces
        XCTAssertEqual(rigidBody.forces, forces)
      }
    }
  }
  #endif
}

extension MM4RigidBodyTests {
  static func createDescriptors() throws -> [MM4RigidBodyDescriptor] {
    var output: [MM4RigidBodyDescriptor] = []
    for atomCode in [MM4AtomCode.alkaneCarbon, .silicon] {
      let adamantane = Adamantane(atomCode: atomCode)
      
      var paramsDesc = MM4ParametersDescriptor()
      paramsDesc.atomicNumbers = adamantane.atomicNumbers
      paramsDesc.bonds = adamantane.bonds
      let params = try MM4Parameters(descriptor: paramsDesc)
      
      var rigidBodyDesc = MM4RigidBodyDescriptor()
      rigidBodyDesc.parameters = params
      rigidBodyDesc.positions = adamantane.positions
      output.append(rigidBodyDesc)
    }
    
    do {
      var paramsDesc = MM4ParametersDescriptor()
      paramsDesc.atomicNumbers = []
      paramsDesc.bonds = []
      let params = try MM4Parameters(descriptor: paramsDesc)
      
      var rigidBodyDesc = MM4RigidBodyDescriptor()
      rigidBodyDesc.parameters = params
      rigidBodyDesc.positions = []
      output.append(rigidBodyDesc)
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
  XCTAssertEqual(mass, rigidBody.mass, accuracy: 1e-3)
  
  let centerOfMass = deriveCenterOfMass(rigidBody)
  XCTAssertEqual(centerOfMass, rigidBody.centerOfMass, accuracy: 1e-3)
  
  let I = deriveInertiaTensor(rigidBody)
  let Σ = rigidBody.principalAxes
  let Λ = rigidBody.momentOfInertia
  let ΛΣT = (
    Λ[0] * SIMD3(Σ.0[0], Σ.1[0], Σ.2[0]),
    Λ[1] * SIMD3(Σ.0[1], Σ.1[1], Σ.2[1]),
    Λ[2] * SIMD3(Σ.0[2], Σ.1[2], Σ.2[2]))
  
  func gemv(
    matrix: (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>),
    vector: SIMD3<Double>
  ) -> SIMD3<Double> {
    matrix.0 * vector[0] + matrix.1 * vector[1] + matrix.2 * vector[2]
  }
  XCTAssertEqual(I.0, gemv(matrix: Σ, vector: ΛΣT.0), accuracy: 1e-3)
  XCTAssertEqual(I.1, gemv(matrix: Σ, vector: ΛΣT.1), accuracy: 1e-3)
  XCTAssertEqual(I.2, gemv(matrix: Σ, vector: ΛΣT.2), accuracy: 1e-3)
}

// MARK: - Copy on Write

private func testCoW(_ descriptor: MM4RigidBodyDescriptor) {
  var rigidBody1 = MM4RigidBody(descriptor: descriptor)
  var rigidBody2 = rigidBody1
  let parameters = rigidBody1.parameters
  
  let centerOfMass1 = rigidBody1.centerOfMass
  XCTAssertEqual(centerOfMass1, rigidBody1.centerOfMass, accuracy: 1e-3)
  XCTAssertEqual(centerOfMass1, rigidBody2.centerOfMass, accuracy: 1e-3)
  XCTAssertEqual(rigidBody1.centerOfMass, rigidBody2.centerOfMass)
  
  let delta2 = SIMD3<Double>.random(in: -2...2)
  do {
    rigidBody2.centerOfMass += delta2
    let centerOfMass2 = centerOfMass1 + delta2
    XCTAssertEqual(centerOfMass1, rigidBody1.centerOfMass, accuracy: 1e-3)
    XCTAssertEqual(
      (parameters.atoms.count > 0) ? centerOfMass2 : .zero,
      rigidBody2.centerOfMass, accuracy: 1e-3)
    XCTAssertNotEqual(centerOfMass1, centerOfMass2)
    if parameters.atoms.count > 0 {
      XCTAssertNotEqual(rigidBody1.centerOfMass, rigidBody2.centerOfMass)
    }
  }
    
  let rigidBody3 = rigidBody1
  do {
    let momentum1 = SIMD3<Double>.random(in: -0.2...0.2)
    rigidBody1.linearMomentum = momentum1
    XCTAssertEqual(
      (parameters.atoms.count > 0) ? momentum1 : .zero,
      rigidBody1.linearMomentum, accuracy: 1e-3)
    
    XCTAssertEqual(.zero, rigidBody2.linearMomentum, accuracy: 1e-3)
    XCTAssertEqual(.zero, rigidBody3.linearMomentum, accuracy: 1e-3)
    if parameters.atoms.count > 0 {
      XCTAssertNotEqual(rigidBody1.linearMomentum, rigidBody2.linearMomentum)
      XCTAssertNotEqual(rigidBody1.linearMomentum, rigidBody3.linearMomentum)
    }
    XCTAssertEqual(rigidBody2.linearMomentum, rigidBody3.linearMomentum)
  }
  
  do {
    XCTAssertEqual(rigidBody1.centerOfMass, rigidBody3.centerOfMass)
    if parameters.atoms.count > 0 {
      XCTAssertNotEqual(rigidBody2.centerOfMass, rigidBody3.centerOfMass)
    }
    for i in parameters.atoms.indices {
      let original = descriptor.positions![i]
      let modified = original + SIMD3<Float>(delta2)
      XCTAssertEqual(rigidBody1.positions[i], original, accuracy: 1e-3)
      XCTAssertEqual(rigidBody2.positions[i], modified, accuracy: 1e-3)
      XCTAssertEqual(rigidBody3.positions[i], original, accuracy: 1e-3)
    }
  }
  
  var rigidBody4 = rigidBody3
  do {
    let velocity4 = SIMD3<Double>.random(in: -0.2...0.2)
    rigidBody4.linearMomentum = velocity4 * rigidBody4.mass
    for i in parameters.atoms.indices {
      XCTAssertEqual(rigidBody3.velocities[i], .zero, accuracy: 1e-5)
      XCTAssertEqual(rigidBody4.velocities[i], SIMD3(velocity4), accuracy: 1e-5)
    }
    XCTAssertEqual(.zero, rigidBody3.linearMomentum, accuracy: 1e-5)
    XCTAssertEqual(
      (parameters.atoms.count > 0) ? velocity4 : .zero,
      rigidBody4.linearMomentum / rigidBody4.mass,
      accuracy: 1e-5)
  }
}

// MARK: - Utilities

func XCTAssertEqual<T: SIMD>(
  _ lhs: T,
  _ rhs: T,
  accuracy: T.Scalar,
  _ message: String = "",
  file: StaticString = #filePath,
  line: UInt = #line
)
where T.Scalar : FloatingPoint {
  for lane in 0..<T.scalarCount {
    XCTAssertEqual(
      lhs[lane], rhs[lane], accuracy: accuracy,
      "v[\(lane)] did not match for '\(message)'.",
      file: file, line: line)
  }
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
) -> SIMD3<Double> {
  var output: SIMD3<Double> = .zero
  for i in rigidBody.parameters.atoms.indices {
    let position = rigidBody.positions[i]
    let mass = rigidBody.parameters.atoms.masses[i]
    output += SIMD3<Double>(mass * position)
  }
  
  let mass = deriveMass(rigidBody)
  if mass == 0 {
    return .zero
  } else {
    output /= mass
    return SIMD3<Double>(output)
  }
}

private func deriveInertiaTensor(
  _ rigidBody: MM4RigidBody
) -> (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>) {
  let centerOfMass = deriveCenterOfMass(rigidBody)
  var output: (
    SIMD3<Double>, SIMD3<Double>, SIMD3<Double>
  ) = (.zero, .zero, .zero)
  
  for i in rigidBody.parameters.atoms.indices {
    let delta = rigidBody.positions[i] - SIMD3<Float>(centerOfMass)
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
    SIMD3<Double>(output.0),
    SIMD3<Double>(output.1),
    SIMD3<Double>(output.2)
  )
}

private func invertInertiaTensor(
  _ momentOfInertia: (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>)
) -> (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>) {
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
    SIMD3<Double>(column0),
    SIMD3<Double>(column1),
    SIMD3<Double>(column2)
  )
}
