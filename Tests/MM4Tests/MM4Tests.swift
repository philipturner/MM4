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

// TODO: Pick some examples experimental and computed forcefield values from the
// MM3/MM4 research papers. Measure how well the current implementation
// reproduces them. This should be a litmus test of how good the forcefield is.
//
// Then, this should expand to measuring bulk properties of materials we want to
// handle. Diamond, lonsdaleite, silicon carbide, silicon, germanium carbide,
// and germanium should all be checked for correctness. We may need to tweak the
// parameters for germanium to get correct results.

final class MM4ParametersTests: XCTestCase {
  // TODO: Measure material properties of solid germanium to check that the
  // Ge-Ge-Ge equilibrium angle is correct.
  
  // TODO: Measure the match of perfluoroethane to the MM4 forcefield
  // calculations in the paper. This should help me reverse-engineer MM4's
  // mechanism for 1,4 electrostatics exceptions. One of the following:
  // - Zero out the contributions from dipoles on torsions
  // - Run them at half-strength
  // - Run them at full-strength
  // - Include 1,3 interactions as well
  //
  // This may also determine whether computing charge-charge or classical
  // dipole-dipole is measurably more accurate. The initial implementation will
  // use dipoles because that's the most likely algorithm used.
  
  // TODO: Use cyclosilane ab initio structures as a test case.
}

// TODO: In another file testing basic functionality, see what happens when
// 'MM4RigidBody' and/or 'MM4ForceField' has zero atoms. There should be no
// crashes or errors.
//
// TODO: Test mass of hydrogen atoms with the new rule. Does it agree with
// what you expect from the old rule?
//
// TODO: Test whether the code enters an infinite loop while generating random
// velocities.
