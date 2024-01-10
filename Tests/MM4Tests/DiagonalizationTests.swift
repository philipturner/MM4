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

#if DEBUG
final class DiagonalizationTests: XCTestCase {
  
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
        XCTAssertEqual(root0, 109983.7421875, accuracy: root0 * 1e-5)
      }
      if let root1 {
        XCTAssertEqual(root1, 4043.4931640625, accuracy: root1 * 1e-5)
      }
      if let root2 {
        XCTAssertEqual(root2, 106430.8046875, accuracy: root2 * 1e-5)
      }
    }
  }
  
  func testInvert() throws {
    do {
      let matrix = (
        SIMD3(5.973220001159007, -782.010009765625, -0.0),
        SIMD3(-782.010009765625, 102381.33528054804, -0.0),
        SIMD3(-0.0, -0.0, 105940.24934304804))
      let expectedInverse = (
        SIMD3(15549.384063716987, 118.76944122867943, 0.0),
        SIMD3(118.76944122867943, 0.9071955512260657, 0.0),
        SIMD3(0.0, 0.0, 9.439283050510099e-06))
      let inverse = invert(matrix: matrix)!
      XCTAssertEqual(inverse.0, expectedInverse.0, ratioAccuracy: 1e-8)
      XCTAssertEqual(inverse.1, expectedInverse.1, ratioAccuracy: 1e-8)
      XCTAssertEqual(inverse.2, expectedInverse.2, ratioAccuracy: 1e-8)
    }
    
    do {
      let matrix = (
        SIMD3(5.973163166685026, -782.010498046875, -7.963180541992188e-05),
        SIMD3(-782.010498046875, 102381.33497957293, -6.258487701416016e-06),
        SIMD3(-7.963180541992188e-05, -6.258487701416016e-06, 105940.24904207293))
      let expectedInverse = (
        SIMD3(205266640498.44598, 1567870430.6856866, 154.3848142804133),
        SIMD3(1567870430.6856866, 11975729.136752428, 1.1792241772427414),
        SIMD3(154.3848142804133, 1.1792241772427414, 9.552917198195528e-06))
      let inverse = invert(matrix: matrix)!
      XCTAssertEqual(inverse.0, expectedInverse.0, ratioAccuracy: 1e-8)
      XCTAssertEqual(inverse.1, expectedInverse.1, ratioAccuracy: 1e-8)
      XCTAssertEqual(inverse.2, expectedInverse.2, ratioAccuracy: 1e-8)
    }
  }
  
  func testDiagonalize() throws {
    func checkEigenPairs(
      _ actual: (
        eigenValues: SIMD3<Double>,
        eigenVectors: (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>)
      )?,
      _ expected: (
        eigenValues: SIMD3<Double>,
        eigenVectors: (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>)
      ),
      file: StaticString = #filePath,
      line: UInt = #line
    ) {
      guard let actual else {
        XCTAssert(false, "Matrix failed to diagonalize", file: file, line: line)
        return
      }
      
      XCTAssertEqual(
        actual.eigenValues, expected.eigenValues, ratioAccuracy: 1e-5,
        file: file, line: line)
      XCTAssertEqual(
        actual.eigenVectors.0, expected.eigenVectors.0, accuracy: 1e-4,
        file: file, line: line)
      XCTAssertEqual(
        actual.eigenVectors.1, expected.eigenVectors.1, accuracy: 1e-4,
        file: file, line: line)
      XCTAssertEqual(
        actual.eigenVectors.2, expected.eigenVectors.2, accuracy: 1e-4,
        file: file, line: line)
    }
    
    do {
      let matrix = (
        SIMD3(4049.46630859375, -782.010498046875, -7.963180541992188e-05),
        SIMD3(-782.010498046875, 106424.828125, -6.258487701416016e-06),
        SIMD3(-7.963180541992188e-05, -6.258487701416016e-06, 109983.7421875))
      let expectedEigenValues = SIMD3(
        109983.74, 106430.805, 4043.4932)
      let expectedEigenVectors = (
        SIMD3(0, 0, 1.0),
        SIMD3(-0.0076379897, 0.99997085, 0),
        SIMD3(-0.99997085, -0.0076379897, 0))
      
      let eigenPairs = diagonalize(matrix: matrix)
      checkEigenPairs(
        eigenPairs, (expectedEigenValues, expectedEigenVectors))
    }
    
    do {
      let matrix = (
        SIMD3(4049.4675941467285, -782.0110626220703, -2.2411346435546875e-05),
        SIMD3(-782.0110626220703, 106424.83778762817, 4.26173210144043e-06),
        SIMD3(-2.2411346435546875e-05, 4.26173210144043e-06, 109983.75480651855))
      let expectedEigenValues = SIMD3(
        109983.7548065182, 106430.81095893068, 4043.494422844558)
      let expectedEigenVectors = (
        SIMD3(0, 0, 1.0),
        SIMD3(-0.0076379897, 0.99997085, 0),
        SIMD3(-0.99997085, -0.0076379897, 0))
      
      let eigenPairs = diagonalize(matrix: matrix)
      checkEigenPairs(
        eigenPairs, (expectedEigenValues, expectedEigenVectors))
    }
    
    do {
      let matrix = (
        SIMD3(4049.468017578125, -770.1312255859375, -135.79461669921875),
        SIMD3(-770.1312255859375, 106532.140625, -608.6104736328125),
        SIMD3(-135.79461669921875, -608.6104736328125, 109876.4296875))
      let expectedEigenValues = SIMD3(
        109983.74, 106430.8, 4043.4949)
      let expectedEigenVectors = (
        SIMD3(0, -0.17364809, 0.9848078),
        SIMD3(-0.0076380027, 0.98477906, 0.17364302),
        SIMD3(-0.99997085, -0.0075219646, -0.0013263224))
      
      let eigenPairs = diagonalize(matrix: matrix)
      checkEigenPairs(
        eigenPairs, (expectedEigenValues, expectedEigenVectors))
    }
    
    // This test case required an additional trial to find 1 eigenvector.
    do {
      let matrix = (
        SIMD3(82663.453125, 19038.958984375, 41616.65234375),
        SIMD3(19038.958984375, 91756.046875, -30655.078125),
        SIMD3(41616.65234375, -30655.078125, 46038.53515625))
      let expectedEigenValues = SIMD3(
        109983.74, 106430.8, 4043.494)
      let expectedEigenVectors = (
        SIMD3(0.77459675, -0.19999957, 0.6),
        SIMD3(0.3834677, 0.90293205, -0.19407801),
        SIMD3(-0.5029437, 0.38041282, 0.7761016))
      
      let eigenPairs = diagonalize(matrix: matrix)
      checkEigenPairs(
        eigenPairs, (expectedEigenValues, expectedEigenVectors))
    }
    
    // This test case required an additional trial to find 1 eigenvector.
    do {
      let matrix = (
        SIMD3(27748.875, 26787.69921875, -34814.97265625),
        SIMD3(26787.69921875, 97504.8984375, 10952.2744140625),
        SIMD3(-34814.97265625, 10952.2744140625, 95204.203125))
      let expectedEigenValues = SIMD3(
        109983.72, 106430.77, 4043.49)
      let expectedEigenVectors = (
        SIMD3(-0.41388798, -0.09375615, 0.9054869),
        SIMD3(0.23297852, 0.9506457, 0.20492388),
        SIMD3(-0.88001007, 0.29577452, -0.37161765))
      
      let eigenPairs = diagonalize(matrix: matrix)
      checkEigenPairs(
        eigenPairs, (expectedEigenValues, expectedEigenVectors))
    }
    
    // This test case required an additional trial to find 2 eigenvectors.
    do {
      let matrix = (
        SIMD3(93399.5, 17936.953125, 34045.8984375),
        SIMD3(17936.953125, 87377.6328125, -35676.8359375),
        SIMD3(34045.8984375, -35676.8359375, 39680.83984375))
      let expectedEigenValues = SIMD3(
        109983.72, 106430.766, 4043.4954)
      let expectedEigenVectors = (
        SIMD3(0.9168911, 0.1343581, 0.37584394),
        SIMD3(0.05351307, 0.89175814, -0.44933698),
        SIMD3(-0.39553395, 0.43210563, 0.8104552))
      
      let eigenPairs = diagonalize(matrix: matrix)
      checkEigenPairs(
        eigenPairs, (expectedEigenValues, expectedEigenVectors))
    }
  }
}
#endif

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

func XCTAssertEqual<T: SIMD>(
  _ lhs: T,
  _ rhs: T,
  ratioAccuracy: T.Scalar,
  _ message: String = "",
  file: StaticString = #filePath,
  line: UInt = #line
)
where T.Scalar : FloatingPoint {
  for lane in 0..<T.scalarCount {
    let maxValue = max(lhs[lane].magnitude, rhs[lane].magnitude)
    XCTAssertEqual(
      lhs[lane], rhs[lane],
      accuracy: maxValue * ratioAccuracy,
      "v[\(lane)] did not match for '\(message)'.",
      file: file, line: line)
  }
}
