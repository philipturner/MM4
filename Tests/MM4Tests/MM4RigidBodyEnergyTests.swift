import XCTest
import MM4

// Test the correctness of APIs that deal with rigid body energy. This includes
// simultaneous adherence to the Boltzmann distribution and conservation of
// momentum. In addition, test that the heat capacity heuristic accesses the
// lookup table correctly.
//
// Test the energy of anchor atoms is separated from the total energy.
