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

// The C API will be drafted and implemented right now. The Python API will have
// to wait until the second near-term round of development. During this round,
// the forcefield will be debugged and Swift test cases will be created. The
// development round will begin after sufficient progress has been made on the
// geometry compiler.

__attribute__((aligned(16)))
struct MM4Float3 {
  float x;
  float y;
  float z;
};

__attribute__((aligned(8)))
struct MM4UInt2 {
  uint32_t x;
  uint32_t y;
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
typedef void MM4ParametersDescriptor;
typedef void MM4Parameters;
typedef void MM4StateDescriptor;
typedef void MM4State;

#define MM4_ARRAY(expr) expr, int64_t* size

// All arrays for getters/setters should have length equal to the number of
// atoms in the system. The exception is rigid bodies, which should equal the
// number of atoms. If an array is null, the function will return the needed
// number of elements in the 'size' argument.
//
// Unless explicitly stated in this header, you do not need to deallocate an
// object returned from a function.
//
// Instances of 'Bool' in Swift are replaced with 'uint8_t' in C.

// MARK: - CoreObjects/MM4State.swift

MM4StateDescriptor* MM4StateDescriptor_init();
void MM4StateDescriptor_deinit(MM4StateDescriptor* target);

uint8_t MM4StateDescriptor_getEnergy(MM4StateDescriptor* target);
uint8_t MM4StateDescriptor_getForces(MM4StateDescriptor* target);
uint8_t MM4StateDescriptor_getPositions(MM4StateDescriptor* target);
uint8_t MM4StateDescriptor_getVelocities(MM4StateDescriptor* target);

void MM4StateDescriptor_setEnergy(MM4StateDescriptor* target, uint8_t energy);
void MM4StateDescriptor_setForces(MM4StateDescriptor* target, uint8_t forces);
void MM4StateDescriptor_setPositions(MM4StateDescriptor* target, uint8_t positions);
void MM4StateDescriptor_setVelocities(MM4StateDescriptor* target, uint8_t velocities);

void MM4State_destroy(MM4State* target);
void MM4State_getForces(MM4State* target, MM4_ARRAY(MM4Float3* forces));
double MM4State_getKineticEnergy(MM4State* target);
void MM4State_getPositions(MM4State* target, MM4_ARRAY(MM4Float3* positions));
double MM4State_getPotentialEnergy(MM4State* target);
void MM4State_getVelocities(MM4State* target, MM4_ARRAY(MM4Float3* velocities));

/// You own the object created by this function. You must destroy it.
MM4State* MM4ForceField_state(MM4StateDescriptor* descriptor);

// MARK: - ForceField/MM4ForceField.swift

MM4ForceFieldDescriptor* MM4ForceFieldDescriptor_init();
void MM4ForceFieldDescriptor_deinit(MM4ForceFieldDescriptor* target);
MM4Parameters* MM4ForceFieldDescriptor_getParameters(MM4ForceFieldDescriptor* target);
void MM4ForceFieldDescriptor_getPositions(MM4ForceFieldDescriptor* target, MM4_ARRAY(MM4Float3* positions));
void MM4ForceFieldDescriptor_setParameters(MM4ForceFieldDescriptor* target, MM4ForceFieldParameters* parameters);
void MM4ForceFieldDescriptor_setPositions(MM4ForceFieldDescriptor* target, MM4_ARRAY(const MM4Float3* positions));

MM4ForceField* MM4ForceField_init(MM4ForceFieldDescriptor* descriptor);
void MM4ForceField_deinit(MM4ForceField* target);

// MARK: - ForceField/MM4ForceField+Actions.swift

void MM4ForceField_simulate(MM4ForceField* target, double time, double maximumTimeStep, MM4Error** error);
void MM4ForceField_minimize(MM4ForceField* target, double tolerance, int64_t maxIterations, MM4Error** error);
void MM4ForceField_thermalize(MM4ForceField* target, double temperature, MM4_ARRAY(const int64_t* rigidBodies));

// MARK: - ForceField/MM4ForceField+Properties.swift

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
void MM4ForceField_setStationaryAtoms(MM4ForceField* target, MM4_ARRAY(const uint8_t* stationaryAtoms));
void MM4ForceField_setVelocities(MM4ForceField* target, MM4_ARRAY(const MM4Float3* velocities));

// MARK: - ForceField/MM4ForceField+Update.swift

MM4ForceFieldUpdateDescriptor* MM4ForceFieldUpdateDescriptor_init();
void MM4ForceFieldUpdateDescriptor_deinit(MM4ForceFieldUpdateDescriptor* target);
void MM4ForceFieldUpdateDescriptor_getPositions(MM4ForceFieldDescriptor* target, MM4_ARRAY(MM4Float3* positions));
void MM4ForceFieldUpdateDescriptor_getVelocities(MM4ForceFieldDescriptor* target, MM4_ARRAY(MM4Float3* velocities));
void MM4ForceFieldUpdateDescriptor_setPositions(MM4ForceFieldDescriptor* target, MM4_ARRAY(const MM4Float3* positions));
void MM4ForceFieldUpdateDescriptor_setVelocities(MM4ForceFieldDescriptor* target, MM4_ARRAY(const MM4Float3* velocities));

void MM4ForceField_update(MM4ForceFieldUpdateDescriptor* descriptor);

// MARK: - MM4Error.swift

MM4Error* MM4Error_init(const char* description);
void MM4Error_deinit(MM4Error* target);
const char* MM4Error_description(MM4Error* target);

// MARK: - Parameters/MM4Parameters.swift

// The parameters must be arranged in a different way than the Swift API. Swift
// optionals and tuples are not compatible with C. Create a data layout for
// bridging to C. Request that the user enter a pointer to write the data to. If
// the value is 'nil', return 'false' from the function to indicate failure.
//
// Therefore, as of now, you cannot inspect the individual parameters in an
// MM4Parameters object without using the Swift API.

MM4ParametersDescriptor* MM4ParametersDescriptor_init();
void MM4ParametersDescriptor_destroy(MM4ParametersDescriptor* target);

void MM4ParametersDescriptor_getAtomicNumbers(MM4ParametersDescriptor* target, MM4_ARRAY(uint8_t* atomicNumbers));
void MM4ParametersDescriptor_getBonds(MM4ParametersDescriptor* target, MM4_ARRAY(MM4UInt2* bonds));
void MM4ParametersDescriptor_getBondOrders(MM4ParametersDescriptor* target, MM4_ARRAY(float* bondOrders));
double MM4ParametersDescriptor_getHydrogenMassRepartitioning(MM4ParametersDescriptor* target);

void MM4ParametersDescriptor_setAtomicNumbers(MM4ParametersDescriptor* target, MM4_ARRAY(const uint8_t* atomicNumbers));
void MM4ParametersDescriptor_setBonds(MM4ParametersDescriptor* target, MM4_ARRAY(const MM4UInt2* bonds));
void MM4ParametersDescriptor_setBondOrders(MM4ParametersDescriptor* target, MM4_ARRAY(const float* bondOrders));
void MM4ParametersDescriptor_setHydrogenMassRepartitioning(MM4ParametersDescriptor* target, double hydrogenMassRepartitioning);

#undef MM4_ARRAY

#endif /* MM4_h */
