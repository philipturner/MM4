import XCTest
import MM4

// Test the correctness of APIs that deal with rigid body energy. This includes
// simultaneous adherence to the Boltzmann distribution and conservation of
// momentum. In addition, test that the heat capacity heuristic accesses the
// lookup table correctly.
//
// This file should also handle the edge case where an object has thermal
// velocities, and is then rotated. The thermal velocities should rotate with
// the object. A similar test exists in the file 'MM4RigidBodyTests.swift', but
// thermal velocities are out of scope there.
