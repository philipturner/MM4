//
//  MM4ForceField+Properties.swift
//
//
//  Created by Philip Turner on 10/14/23.
//

// Ergonomic APIs for accessing a force field's state. These delegate to
// a single-property instances of the batched functions, 'update' or 'state'.
//
// Use the following OpenMM functions to update the system.
// externalForces
//   - setParticleParameters, updateParametersInContext
//   - setIntegrationForceGroups, swapping integrators
// positions
//   - setState
// stationaryAtoms
//   - setParticleParameters, updateParametersInContext
//   - setIntegrationForceGroups, swapping integrators
// velocities
//   - setState
