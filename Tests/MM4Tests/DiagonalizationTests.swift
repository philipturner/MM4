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
  
  #endif
}
