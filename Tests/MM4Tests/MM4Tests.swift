import XCTest
import MM4

final class MM4Tests: XCTestCase {
  func testExample() throws {
    // XCTest Documentation
    // https://developer.apple.com/documentation/xctest
    
    // Defining Test Cases and Test Methods
    // https://developer.apple.com/documentation/xctest/defining_test_cases_and_test_methods
  }
}

final class MM4ParametersTests: XCTestCase {
  // TODO: Double check the new parameters for O and B-N.
  
  // TODO: Use cyclosilane ab initio structures as a test case.
  
  // TODO: Test that the torsion-stretch parameter actually increases system
  // energy by the expected amount. We can't rule out an error in unit
  // conversions with the custom formula.
  
  // TODO: Look the most unusual edge cases in the parameters, ensure they show
  // up correctly.
  
  // TODO: Some profiling tests. See whether multicore CPU could help for
  // setting up quick energy minimizations of large nanosystems.
}
