import XCTest
import MM4

final class MM4Tests: XCTestCase {
  func testExample() throws {
    
  }
}

// MARK: - Tests for Carbon-Only Molecules

// Every part of the API should have defined behavior in every edge case. This
// includes places where the atom count is zero, every atom is an anchor, etc.
// Each edge case should have a unit test.
//


// MARK: - Tests for Non-Carbon Elements

// Pick some examples experimental and computed forcefield values from the
// MM3/MM4 research papers. Measure how well the current implementation
// reproduces them. This should be a litmus test of how good the forcefield is.
// Once can use ab initio structures from GFN2-xTB as an independent evaluation
// of the correct structure. Examining how GFN calculates vibrational
// frequencies could also enable tests of stiffness.
//
// Then, this should expand to measuring bulk properties of materials we want to
// handle. Diamond, lonsdaleite, silicon carbide, silicon, germanium carbide,
// and germanium should all be checked for correctness. We may need to tweak the
// parameters for germanium to get correct results.

// Measure material properties of solid germanium to check that the
// Ge-Ge-Ge equilibrium angle is correct.
//
// Measure the match of perfluoroethane to the MM4 forcefield
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
//
// Use cyclosilane ab initio structures as a test case.
