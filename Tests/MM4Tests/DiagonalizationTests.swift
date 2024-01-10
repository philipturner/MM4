import XCTest
#if DEBUG
@testable import MM4
#else
import MM4
#endif

// Test the correctness of math functions for diagonalizing the inertia tensor.
// There are 6 different test cases in the original source file. They should
// be copied over into the unit tests here.
// - 3 tests for factoring a cubic polynomial
// - 3 tests for diagonalizing a matrix
// - test that a rigid body for NCFPart successfully initializes

final class DiagonalizationTests: XCTestCase {
  #if DEBUG
  func testFactorCubicPolynomial() throws {
    do {
      let coefficients: SIMD4<Double> = [1, -6, 11, -6]
      let (root0, root1, root2) = factorCubicPolynomial(
        coefficients: coefficients)
      XCTAssertEqual(root0, 1.0)
      XCTAssertEqual(root1, 3.0)
      XCTAssertEqual(root2, 2.0)
    }
    
    do {
      let coefficients: SIMD4<Double> = [1, 1, 7, 2]
      let (root0, root1, root2) = factorCubicPolynomial(
        coefficients: coefficients)
      XCTAssertNotNil(root0)
      XCTAssertNil(root1)
      XCTAssertNil(root2)
      if let root0 {
        XCTAssertEqual(root0, -0.2944532632827759, accuracy: 1e-8)
      }
    }
    
    do {
      let coefficients: SIMD4<Double> = [
        -1.0, 220458.03662109375, -12580728532.801353, 47331747116183.29
      ]
      let (root0, root1, root2) = factorCubicPolynomial(
        coefficients: coefficients)
      XCTAssertNotNil(root0)
      XCTAssertNotNil(root1)
      XCTAssertNotNil(root2)
      if let root0 {
        XCTAssertEqual(root0, 109983.7421875, accuracy: 1e-3)
      }
      if let root1 {
        XCTAssertEqual(root1, 4043.4931640625, accuracy: 1e-3)
      }
      if let root2 {
        XCTAssertEqual(root2, 106430.8046875, accuracy: 1e-2)
      }
    }
  }
  
  func testDiagonalize() throws {
    
  }
  #endif
}
