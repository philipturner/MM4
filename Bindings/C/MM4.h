//
//  MM4.h
//
//
//  Created by Philip Turner on 10/16/23.
//

#ifndef MM4_h
#define MM4_h

#include <cstdint>

// MARK: - Types

__attribute__((aligned(16)))
struct MM4Float3 {
  float x;
  float y;
  float z;
};

__attribute__((aligned(16)))
struct MM4Range {
  long lowerBound;
  long upperBound;
}

typedef void MM4Error;
typedef void MM4ForceField;
typedef void MM4ForceFieldDescriptor;
typedef void MM4ForceFieldUpdateDescriptor;

#define MM4_ARRAY(expr) expr, int64_t* size

// The parameters must be arranged in a different way than the Swift API. Swift
// optionals and tuples are not compatible with C. Create a data layout for
// bridging to C. Request that the user enter a pointer to write the data to. If
// the value is 'nil', return 'false' from the function to indicate failure.

// All arrays for getters/setters should have length equal to the number of
// atoms in the system. The exception is rigid bodies, which should equal the
// number of atoms. If an array is null, the function will return the needed
// number of elements in the 'size' argument.

// Instances of 'Bool' in Swift are replaced with 'uint8_t' in C.

// Arrays have an efficient method of getting/setting. Instead of entering the
// entire at once, you can access the data at a per-element granularity. If an
// index is out of bounds, it returns something invalid (getter) or false
// (setter). Most arrays have a pre-determined element count, but those that
// don't have a separate...
// - Haven't decided how to proceed with this.
// - Note: the APIs are supposed to be batched. In Swift, there is a warning
//   about how efficient fetching even a contiguous array can be. However, an
//   element-by-element API may be extremely useful, even for Swift.
// - Consider how to implement this internally, e.g. by caching and
//   automatically flushing changes from a dirty array. This may be more work
//   and debugging, but a more ergonomic API overall.
// - Add this API after the first round of debugging/testing. It is bound to
//   introduce additional bugs, although that doesn't mean it's a bad idea from
//   a software engineering prospective. A minimum viable product with
//   sub-optimal performance for some use cases is okay.
// - TODO: Add this per-element API and finish the C bindings, once development
//   on MM4 resumes, to add the remaining optimizations.

// MARK: - Functions

// ForceField/MM4ForceField.swift
MM4ForceFieldDescriptor* MM4ForceFieldDescriptor_init();
void MM4ForceFieldDescriptor_deinit(MM4ForceFieldDescriptor* target);
MM4Parameters* MM4ForceFieldDescriptor_getParameters();
void MM4ForceFieldDescriptor_setParameters(MM4ForceFieldParameters* parameters);
// getter for parameters
// getter for positions
// setter for parameters
// setter for positions

MM4ForceField* MM4ForceField_init(MM4ForceFieldDescriptor* descriptor);
void MM4ForceField_deinit(MM4ForceField* target);

// ForceField/MM4ForceField+Actions.swift
void MM4ForceField_simulate(MM4ForceField* target, double time, double maximumTimeStep, MM4Error** error);
void MM4ForceField_minimize(MM4ForceField* target, double tolerance, int64_t maxIterations, MM4Error** error);
void MM4ForceField_thermalize(MM4ForceField* target, double temperature, MM4_ARRAY(const int64_t* rigidBodies));

// ForceField/MM4ForceField+Properties.swift
void MM4ForceField_getExternalForces(MM4ForceField* target, MM4_ARRAY(MM4Float3* externalForces));
void MM4ForceField_getForces(MM4ForceField* target, MM4_ARRAY(MM4Float3* forces));
double MM4ForceField_getKineticEnergy(MM4ForceField* target);
void MM4ForceField_getPositions(MM4ForceField* target, MM4_ARRAY(MM4Float3* positions));
double MM4ForceField_getPotentialEnergy(MM4ForceField* target);
void MM4ForceField_getRigidBodies(MM4ForceField* target, MM4_ARRAY(MM4Range* rigidBodies));
void MM4ForceField_getStationaryAtoms(MM4ForceField* target, MM4_ARRAY(uint8_t* stationaryAtoms));
void MM4ForceField_getVelocities(MM4ForceField* target, MM4_ARRAY(MM4Float3* velocities));

void MM4ForceField_setExternalForces(MM4ForceField* target, MM4_ARRAY(const MM4Float3* externalForces));
void MM4ForceField_setPositions(MM4ForceField* target, MM4_ARRAY(const MM4Float3* positions));
void MM4ForceField_setRigidBodies(MM4ForceField* target, MM4_ARRAY(const MM4Range* rigidBodies));
void MM4ForceField_setStationaryAtoms(MM4ForceField* target, MM4_ARRAY(const uint8_t* stationaryAtoms));
void MM4ForceField_setVelocities(MM4ForceField* target, MM4_ARRAY(const MM4Float3* velocities));

// ForceField/MM4ForceField+Update.swift
MM4ForceFieldUpdateDescriptor* MM4ForceFieldUpdateDescriptor_init();
void MM4ForceFieldUpdateDescriptor_deinit(MM4ForceFieldUpdateDescriptor* target);

// MM4Error.swift
MM4Error* MM4Error_init(const char* description);
void MM4Error_deinit(MM4Error* target);
const char* MM4Error_description(MM4Error* target);

#undef MM4_ARRAY

#endif /* MM4_h */
