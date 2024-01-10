import XCTest
import MM4

// MARK: - Test Execution

final class MM4RigidBodyTests: XCTestCase {
  
  // TODO: Test the rotate() function.
  //
  // After adding the test for rotate(), this should be ready to merge into the
  // main branch.
  
  func testInertia() throws {
    for descriptor in MM4RigidBodyTests.descriptors {
      // Assert that all quantities are equal within 1e-3 tolerance. This is an
      // internally consistent unit system where the fundamental quantities are
      // roughly 1. Therefore, absolute errors should be on the same order as
      // what you expect for relative errors when using a multiplicative factor.
      // Using absolute errors enables use of more ergonomic testing functions.
      
      let rigidBody = try MM4RigidBody(descriptor: descriptor)
      let mass = deriveMass(descriptor)
      XCTAssertEqual(mass, rigidBody.mass, accuracy: 1e-3)
      
      let centerOfMass = deriveCenterOfMass(descriptor)
      XCTAssertEqual(centerOfMass, rigidBody.centerOfMass, accuracy: 1e-3)
      
      let I = deriveInertiaTensor(descriptor)
      let Σ = rigidBody.principalAxes
      let Λ = rigidBody.momentOfInertia
      let ΣT = (
        SIMD3(Σ.0[0], Σ.1[0], Σ.2[0]),
        SIMD3(Σ.0[1], Σ.1[1], Σ.2[1]),
        SIMD3(Σ.0[2], Σ.1[2], Σ.2[2]))
      let ΛΣT = (Λ * ΣT.0, Λ * ΣT.1, Λ * ΣT.2)
      
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
  }
  
  func testEmptyRigidBody() throws {
    for descriptor in MM4RigidBodyTests.descriptors {
      XCTAssertNotEqual(descriptor.parameters!.atoms.count, 0)
      XCTAssertNoThrow(try MM4RigidBody(descriptor: descriptor))
    }
    do {
      let descriptor = MM4RigidBodyTests.emptyDescriptor
      XCTAssertEqual(descriptor.parameters!.atoms.count, 0)
      XCTAssertThrowsError(try MM4RigidBody(descriptor: descriptor))
    }
  }
  
  func testCoW() throws {
    for descriptor in MM4RigidBodyTests.descriptors {
      var rigidBody1 = try MM4RigidBody(descriptor: descriptor)
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
          XCTAssertNotEqual(
            rigidBody1.linearMomentum, rigidBody2.linearMomentum)
          XCTAssertNotEqual(
            rigidBody1.linearMomentum, rigidBody3.linearMomentum)
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
          XCTAssertEqual(
            rigidBody4.velocities[i], SIMD3(velocity4), accuracy: 1e-5)
        }
        XCTAssertEqual(.zero, rigidBody3.linearMomentum, accuracy: 1e-5)
        XCTAssertEqual(
          (parameters.atoms.count > 0) ? velocity4 : .zero,
          rigidBody4.linearMomentum / rigidBody4.mass,
          accuracy: 1e-5)
      }
    }
  }
}

// MARK: - Descriptors

extension MM4RigidBodyTests {
  static var descriptors: [MM4RigidBodyDescriptor] {
    var output: [MM4RigidBodyDescriptor] = []
    for atomCode in [MM4AtomCode.alkaneCarbon, .silicon] {
      let adamantane = Adamantane(atomCode: atomCode)
      
      var paramsDesc = MM4ParametersDescriptor()
      paramsDesc.atomicNumbers = adamantane.atomicNumbers
      paramsDesc.bonds = adamantane.bonds
      let params = try! MM4Parameters(descriptor: paramsDesc)
      
      var rigidBodyDesc = MM4RigidBodyDescriptor()
      rigidBodyDesc.parameters = params
      rigidBodyDesc.positions = adamantane.positions
      output.append(rigidBodyDesc)
    }
    return output
  }
  
  static var emptyDescriptor: MM4RigidBodyDescriptor {
    var paramsDesc = MM4ParametersDescriptor()
    paramsDesc.atomicNumbers = []
    paramsDesc.bonds = []
    let params = try! MM4Parameters(descriptor: paramsDesc)
    
    var rigidBodyDesc = MM4RigidBodyDescriptor()
    rigidBodyDesc.parameters = params
    rigidBodyDesc.positions = []
    return rigidBodyDesc
  }
}

// MARK: - Utilities

private func deriveMass(_ descriptor: MM4RigidBodyDescriptor) -> Double {
  var output: Double = .zero
  for i in descriptor.parameters!.atoms.indices {
    output += Double(descriptor.parameters!.atoms.masses[i])
  }
  return output
}

private func deriveCenterOfMass(
  _ descriptor: MM4RigidBodyDescriptor
) -> SIMD3<Double> {
  var output: SIMD3<Double> = .zero
  for i in descriptor.parameters!.atoms.indices {
    let position = descriptor.positions![i]
    let mass = descriptor.parameters!.atoms.masses[i]
    output += SIMD3<Double>(mass * position)
  }
  
  let mass = deriveMass(descriptor)
  if mass == 0 {
    return .zero
  } else {
    output /= mass
    return SIMD3<Double>(output)
  }
}

private func deriveInertiaTensor(
  _ descriptor: MM4RigidBodyDescriptor
) -> (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>) {
  let centerOfMass = deriveCenterOfMass(descriptor)
  var output: (
    SIMD3<Double>, SIMD3<Double>, SIMD3<Double>
  ) = (.zero, .zero, .zero)
  
  for i in descriptor.parameters!.atoms.indices {
    let delta = descriptor.positions![i] - SIMD3<Float>(centerOfMass)
    let mass = descriptor.parameters!.atoms.masses[i]
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

private func invertMatrix(
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
