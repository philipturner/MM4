import XCTest
import MM4

// Test the correctness of APIs that deal with rigid body energy. This includes
// simultaneous adherence to the Boltzmann distribution and conservation of
// momentum. In addition, test that the heat capacity heuristic accesses the
// lookup table correctly.
//
// Test whether the code enters an infinite loop while generating random
// velocities. This is a hypothesized bug in the untested implementation.
// Plot a histogram of the random velocities and ensure they match one from
// OpenMM at 1.5 kT. Maybe even fall back to OpenMM's implementation requiring
// 1.5 kT for now.
//
// Edge cases to test when adding support for anchors:
// - Ensure the energy of anchor atoms is separated from the total energy.
