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
    
    do {
      let coefficients: SIMD4<Double> = [
        -1.0, 681.8423500061035, -154951.57960266608, 11736402.461415779
      ]
      let (root0, root1, root2) = factorCubicPolynomial(
        coefficients: coefficients)
      XCTAssertNotNil(root0)
      XCTAssertNotNil(root1)
      XCTAssertNotNil(root2)
      if let root0 {
        XCTAssertEqual(root0, 229.73597, accuracy: root0 * 1e-5)
      }
      if let root1 {
        XCTAssertEqual(root1, 222.37041, accuracy: root1 * 1e-5)
      }
      if let root2 {
        XCTAssertEqual(root2, 229.73597, accuracy: root2 * 1e-5)
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
      _ matrix: (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>),
      _ expected: (
        eigenValues: SIMD3<Double>,
        eigenVectors: (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>)
      ),
      file: StaticString = #filePath,
      line: UInt = #line
    ) {
      let (Λ, Σ, failureReason) = diagonalize(matrix: matrix)
      guard let Λ, let Σ else {
        XCTAssert(
          false, "Matrix failed to diagonalize: \(failureReason!)",
          file: file, line: line)
        return
      }
      
      XCTAssertEqual(
        Λ, expected.eigenValues, ratioAccuracy: 1e-5,
        file: file, line: line)
      XCTAssertEqual(
        Σ.0, expected.eigenVectors.0, accuracy: 1e-4,
        file: file, line: line)
      XCTAssertEqual(
        Σ.1, expected.eigenVectors.1, accuracy: 1e-4,
        file: file, line: line)
      XCTAssertEqual(
        Σ.2, expected.eigenVectors.2, accuracy: 1e-4,
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
      checkEigenPairs(
        matrix, (expectedEigenValues, expectedEigenVectors))
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
      checkEigenPairs(
        matrix, (expectedEigenValues, expectedEigenVectors))
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
      checkEigenPairs(
        matrix, (expectedEigenValues, expectedEigenVectors))
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
      checkEigenPairs(
        matrix, (expectedEigenValues, expectedEigenVectors))
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
      checkEigenPairs(
        matrix, (expectedEigenValues, expectedEigenVectors))
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
      checkEigenPairs(
        matrix, (expectedEigenValues, expectedEigenVectors))
    }
    
    // This test case caused convergence issues.
    do {
      let matrix = (
        SIMD3(35702.36472773552, -1.3455748558044434e-05, 7.459503173828125),
        SIMD3(-1.3455748558044434e-05, 82421.14924621582, -5.0455331802368164e-05),
        SIMD3(7.459503173828125, -5.0455331802368164e-05, 47511.179409742355))
      let (Λ, Σ, failureReason) = diagonalize(matrix: matrix)
      XCTAssertNotNil(Λ)
      XCTAssertNotNil(Σ, failureReason ?? "")
    }
    
    // This test case caused convergence issues.
    do {
      let matrix = (
        SIMD3(35702.347987532616, 0.007739812135696411, 7.47467041015625),
        SIMD3(0.007739812135696411, 82421.11811828613, 0.007063537836074829),
        SIMD3(7.47467041015625, 0.007063537836074829, 47511.16352403164))
      let (Λ, Σ, failureReason) = diagonalize(matrix: matrix)
      XCTAssertNotNil(Λ)
      XCTAssertNotNil(Σ, failureReason ?? "")
    }
    
    // This test case caused convergence issues.
    do {
      let matrix = (
        SIMD3(229.3068447113037, -9.5367431640625e-07, 3.2782554626464844e-07),
        SIMD3(-9.5367431640625e-07, 229.3068389892578, -9.685754776000977e-07),
        SIMD3(3.2782554626464844e-07, -9.685754776000977e-07, 222.27340507507324))
      let (Λ, Σ, failureReason) = diagonalize(matrix: matrix)
      XCTAssertNotNil(Λ)
      XCTAssertNotNil(Σ, failureReason ?? "")
    }
    
    // This test case caused convergence issues.
    do {
      let matrix = (
        SIMD3(229.30638122558594, 0.00022840499877929688, 0.00021526217460632324),
        SIMD3(0.00022840499877929688, 229.30638122558594, 0.00021557509899139404),
        SIMD3(0.00021526217460632324, 0.00021557509899139404, 222.27295684814453))
      let (Λ, Σ, failureReason) = diagonalize(matrix: matrix)
      XCTAssertNotNil(Λ)
      XCTAssertNotNil(Σ, failureReason ?? "")
    }
    
    // This test case has known repeated roots.
    do {
      let matrix = (
        SIMD3(1.0, 0, 0),
        SIMD3(0, 1.0, 0),
        SIMD3(0, 0, 1.0))
      let (Λ, Σ, failureReason) = diagonalize(matrix: matrix)
      XCTAssertNotNil(Λ)
      XCTAssertNotNil(Σ, failureReason ?? "")
    }
    
    // This test case has known repeated roots.
    do {
      let matrix = (
        SIMD3(2.0, 0, 0),
        SIMD3(0, 2.0, 0),
        SIMD3(0, 0, 1.0))
      let (Λ, Σ, failureReason) = diagonalize(matrix: matrix)
      XCTAssertNotNil(Λ)
      XCTAssertNotNil(Σ, failureReason ?? "")
    }
    
    // This test case had a tricky characteristic polynomial.
    do {
      let matrix = (
        SIMD3(229.7359676361084, -1.7881393432617188e-07, 2.3543834686279297e-06),
        SIMD3(-1.7881393432617188e-07, 229.73597145080566, -2.1904706954956055e-06),
        SIMD3(2.3543834686279297e-06, -2.1904706954956055e-06, 222.37041091918945))
      let (Λ, Σ, failureReason) = diagonalize(matrix: matrix)
      XCTAssertNotNil(Λ)
      XCTAssertNotNil(Σ, failureReason ?? "")
    }
    
    // This test case caused convergence issues.
    do {
      let matrix = (
        SIMD3(230.0067958831787, 2.384185791015625e-07, 2.562999725341797e-06),
        SIMD3(2.384185791015625e-07, 230.00679206848145, -1.2218952178955078e-06),
        SIMD3(2.562999725341797e-06, -1.2218952178955078e-06, 222.01723098754883))
      let (Λ, Σ, failureReason) = diagonalize(matrix: matrix)
      XCTAssertNotNil(Λ)
      XCTAssertNotNil(Σ, failureReason ?? "")
    }
    
    // This test case caused convergence issues.
    do {
      let matrix = (
        SIMD3(230.44046783447266, -7.748603820800781e-07, -1.9371509552001953e-06),
        SIMD3(-7.748603820800781e-07, 230.4404697418213, -2.6226043701171875e-06),
        SIMD3(-1.9371509552001953e-06, -2.6226043701171875e-06, 223.87350273132324))
      let (Λ, Σ, failureReason) = diagonalize(matrix: matrix)
      XCTAssertNotNil(Λ)
      XCTAssertNotNil(Σ, failureReason ?? "")
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
