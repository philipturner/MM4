//
//  MM4ForceField+Update.swift
//  MM4
//
//  Created by Philip Turner on 11/21/23.
//

import OpenMM

struct MM4UpdateRecord {
  var positions: Bool = false
  var velocities: Bool = false
  
  func active() -> Bool {
    positions || velocities
  }
  
  mutating func erase() {
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
      
      let particleCount =
      system.parameters.atoms.count + system.virtualSiteCount
      let arrayP = OpenMM_Vec3Array(size: particleCount)
      let arrayV = OpenMM_Vec3Array(size: particleCount)
      for (original, reordered) in system.reorderedIndices.enumerated() {
        let position = positions[Int(original)]
        let velocity = velocities[Int(original)]
        arrayP[Int(reordered)] = SIMD3<Double>(position)
        arrayV[Int(reordered)] = SIMD3<Double>(velocity)
      }
      context.context.positions = arrayP
      context.context.velocities = arrayV
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
