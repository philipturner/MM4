import XCTest
#if DEBUG
@testable import MM4
#else
import MM4
#endif
import Numerics

final class QuaternionTests: XCTestCase {
  #if DEBUG
  func testQuaternionToVector() throws {
    // Use a known rotation - 90° around the Z axis - to test the validity of
    // the formula for generating rotational velocity. Atoms' individual
    // velocities should come out at right angles to the w vector.
    let positions: [SIMD3<Float>] = [
      SIMD3<Float>(1, 0, 0),
      SIMD3<Float>(0, 1, 0),
      SIMD3<Float>(0, 1, 1),
      SIMD3<Float>(0, -1, 0),
      SIMD3<Float>(Float(0.5).squareRoot(), -Float(0.5).squareRoot(), 0)
    ]
    let expectedVelocities: [SIMD3<Float>] = [
      SIMD3<Float>(0, 1, 0),
      SIMD3<Float>(-1, 0, 0),
      SIMD3<Float>(-1, 0, 0),
      SIMD3<Float>(1, 0, 0),
      SIMD3<Float>(Float(0.5).squareRoot(), Float(0.5).squareRoot(), 0),
    ]
    
    // Test the scalarized method for computing angular velocity.
    let angularVelocity: Quaternion<Float> =
      .init(angle: 1.0, axis: [0, 0, 1])
    let w = quaternionToVector(angularVelocity)
    for i in positions.indices {
      let r = positions[i]
      let vExp = expectedVelocities[i]
      let v = SIMD3<Float>(
        w.y * r.z - w.z * r.y,
        w.z * r.x - w.x * r.z,
        w.x * r.y - w.y * r.x)
      
      XCTAssertEqual(v[0], vExp[0], accuracy: 1e-3, "v.x did not match.")
      XCTAssertEqual(v[1], vExp[1], accuracy: 1e-3, "v.y did not match.")
      XCTAssertEqual(v[2], vExp[2], accuracy: 1e-3, "v.z did not match.")
    }
    
    // Test the vectorized method for computing angular velocity.
    var vectorCount = positions.count
    vectorCount += MM4VectorWidth - 1
    vectorCount /= MM4VectorWidth
    for vID in 0..<vectorCount {
      var rX: MM4FloatVector = .zero
      var rY: MM4FloatVector = .zero
      var rZ: MM4FloatVector = .zero
      var vExpX: MM4FloatVector = .zero
      var vExpY: MM4FloatVector = .zero
      var vExpZ: MM4FloatVector = .zero
      
      let vStart = vID * MM4VectorWidth
      let vEnd = min(vStart + MM4VectorWidth, positions.count)
      for lane in 0..<vEnd - vStart {
        rX[lane] = positions[vStart + lane].x
        rY[lane] = positions[vStart + lane].y
        rZ[lane] = positions[vStart + lane].z
        vExpX[lane] = expectedVelocities[vStart + lane].x
        vExpY[lane] = expectedVelocities[vStart + lane].y
        vExpZ[lane] = expectedVelocities[vStart + lane].z
      }
      
      let vX = w.y * rZ - w.z * rY
      let vY = w.z * rX - w.x * rZ
      let vZ = w.x * rY - w.y * rX
      
      XCTAssertEqual(vX, vExpX, accuracy: 1e-3)
      XCTAssertEqual(vY, vExpY, accuracy: 1e-3)
      XCTAssertEqual(vZ, vExpZ, accuracy: 1e-3)
    }
  }
  
  func testQuaternionAct() throws {
    // Test three different methods for computing rotations. They must all
    // simultaneously agree for the test to pass:
    // - quaternion.act(on:)
    // - 3x3 rotation matrix
    // - vQuaternionPrepare/vQuaternionAct
    let quaternionAmount: Int = 10
    let positionsAmount: Int = 10
    var quaternions: [Quaternion<Float>] = []
    var positions: [SIMD3<Float>] = []
    
    // Assert the zero-valued edge cases behave as expected.
    quaternions.append(Quaternion(angle: .zero, axis: [1, 0, 0]))
    positions.append(SIMD3<Float>.zero)
    
    // Randomly initialize the quaternions and vectors.
    while quaternions.count < quaternionAmount {
      let angle = Float.random(in: -4...4)
      var axis = SIMD3<Float>.random(in: -1...1)
      let length = sqrt(axis.x * axis.x + axis.y * axis.y + axis.z * axis.z)
      axis /= length
      
      if angle.magnitude < 1e-3 || length < 1e-3 || length > 1 {
        continue
      }
      let quaternion = Quaternion<Float>(angle: angle, axis: axis)
      precondition(quaternion.length != 0)
      quaternions.append(quaternion)
    }
    while positions.count < positionsAmount {
      let radius: Float = 5
      let pos = SIMD3<Float>.random(in: -radius...radius)
      let length = sqrt(pos.x * pos.x + pos.y * pos.y + pos.z * pos.z)
      
      if length < 1e-3 || length > radius {
        continue
      }
      positions.append(pos)
    }
    
    // Test N quaternions against M vectors; O(N * M).
    for quaternion in quaternions {
      var rotatedPositions: [SIMD3<Float>] = []
      rotatedPositions.append(.zero)
      for position in positions[1...] {
        let rotated = quaternion.act(on: position)
        rotatedPositions.append(rotated)
      }
      
      // Assert the expectation value for the zero-valued edge case.
      XCTAssertEqual(rotatedPositions[0], .zero, accuracy: 1e-3, "0")
      
      // Test the rotation matrix method.
      let (x, y, z) = (
        quaternion.act(on: [1, 0, 0]),
        quaternion.act(on: [0, 1, 0]),
        quaternion.act(on: [0, 0, 1]))
      for positionID in positions.indices {
        let position = positions[positionID]
        let actual = x * position.x + y * position.y + z * position.z
        let expected = rotatedPositions[positionID]
        XCTAssertEqual(actual, expected, accuracy: 1e-3)
      }
      
      // Test the vQuaternionPrepare/vQuaternionAct method.
      let prepared = vQuaternionPrepare(quaternion)
      var vectorCount = positions.count
      vectorCount += MM4VectorWidth - 1
      vectorCount /= MM4VectorWidth
      for vID in 0..<vectorCount {
        var rX: MM4FloatVector = .zero
        var rY: MM4FloatVector = .zero
        var rZ: MM4FloatVector = .zero
        var expX: MM4FloatVector = .zero
        var expY: MM4FloatVector = .zero
        var expZ: MM4FloatVector = .zero
        
        let vStart = vID * MM4VectorWidth
        let vEnd = min(vStart + MM4VectorWidth, positions.count)
        for lane in 0..<vEnd - vStart {
          rX[lane] = positions[vStart + lane].x
          rY[lane] = positions[vStart + lane].y
          rZ[lane] = positions[vStart + lane].z
          expX[lane] = rotatedPositions[vStart + lane].x
          expY[lane] = rotatedPositions[vStart + lane].y
          expZ[lane] = rotatedPositions[vStart + lane].z
        }
        
        vQuaternionAct(prepared, &rX, &rY, &rZ)
        XCTAssertEqual(rX, expX, accuracy: 1e-3)
        XCTAssertEqual(rY, expY, accuracy: 1e-3)
        XCTAssertEqual(rZ, expZ, accuracy: 1e-3)
      }
    }
  }
  #endif
}
