//
//  MM4ForceField+Update.swift
//
//
//  Created by Philip Turner on 11/21/23.
//

import OpenMM

struct MM4UpdateRecord {
  var externalForces: Bool = false
  var positions: Bool = false
  var velocities: Bool = false
  
  func active() -> Bool {
    externalForces || positions || velocities
  }
  
  mutating func erase() {
    externalForces = false
    positions = false
    velocities = false
  }
}

extension MM4ForceField {
  func ensurePositionsAndVelocitiesCached() {
    if cachedState.positions == nil ||
        cachedState.velocities == nil {
      if updateRecord.active() {
        fatalError(
          "Fetched new positions or velocities while previous updates were not flushed.")
      }
      var descriptor = MM4StateDescriptor()
      descriptor.positions = true
      descriptor.velocities = true
      
      let state = self.state(descriptor: descriptor)
      cachedState.positions = state.positions!
      cachedState.velocities = state.velocities!
    }
  }
  
  func ensureForcesAndEnergyCached() {
    if updateRecord.active() {
      // Forces and energies must be computed using the most recent system
      // state.
      flushUpdateRecord()
      invalidateForcesAndEnergy()
    }
    
    if cachedState.forces == nil ||
        cachedState.kineticEnergy == nil ||
        cachedState.potentialEnergy == nil {
      var descriptor = MM4StateDescriptor()
      descriptor.forces = true
      descriptor.energy = true
      
      let state = self.state(descriptor: descriptor)
      cachedState.forces = state.forces!
      cachedState.kineticEnergy = state.kineticEnergy!
      cachedState.potentialEnergy = state.potentialEnergy!
    }
  }
  
  func flushUpdateRecord() {
    if updateRecord.positions || updateRecord.velocities {
      guard let positions = cachedState.positions,
            let velocities = cachedState.velocities else {
        fatalError("Positions or velocities not fetched before update.")
      }
      
      let arrayP = OpenMM_Vec3Array(size: system.reorderedIndices.count)
      let arrayV = OpenMM_Vec3Array(size: system.reorderedIndices.count)
      for (original, reordered) in system.reorderedIndices.enumerated() {
        let position = positions[Int(original)]
        let velocity = velocities[Int(original)]
        arrayP[Int(reordered)] = SIMD3<Double>(position)
        arrayV[Int(reordered)] = SIMD3<Double>(velocity)
      }
      context.context.positions = arrayP
      context.context.velocities = arrayV
    }
    
    if updateRecord.externalForces {
      guard _externalForces.count == system.parameters.atoms.count else {
        fatalError("Number of external forces does not match atom count.")
      }
      
      // Do not enforce the restriction that anchors must have nonzero external
      // force here. The check is not performed here for the bulk linear
      // velocity either. Rather, the error should appear when exporting to a
      // rigid body.
      let force = system.forces.external
      force.updateForces(_externalForces, system: system)
      force.updateParametersInContext(context)
    }
    
    updateRecord.erase()
  }
  
  func invalidatePositionsAndVelocities() {
    cachedState.positions = nil
    cachedState.velocities = nil
  }
  
  func invalidateForcesAndEnergy() {
    cachedState.forces = nil
    cachedState.kineticEnergy = nil
    cachedState.potentialEnergy = nil
  }
}
