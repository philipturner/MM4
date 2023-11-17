//
//  MM4System+VirtualSites.swift
//
//
//  Created by Philip Turner on 11/17/23.
//

import OpenMM

extension MM4System {
  func createReorderedIndices(parameters: MM4Parameters) {
    let atomicNumbers = parameters.atoms.atomicNumbers
    var virtualSiteCount = 0
    for originalID in 0..<atomicNumbers.count {
      if atomicNumbers[originalID] == 1 {
        virtualSiteCount += 1
      }
    }
    self.virtualSiteIndices.reserveCapacity(virtualSiteCount)
    
    let particleCount = virtualSiteCount + atomicNumbers.count
    self.reorderedIndices = Array(repeating: -1, count: atomicNumbers.count)
    self.originalIndices = Array(repeating: -1, count: particleCount)
    
    for originalID in 0..<atomicNumbers.count {
      if atomicNumbers[originalID] == 1 {
        let reorderedID = virtualSiteIndices.count
        reorderedIndices[originalID] = Int32(truncatingIfNeeded: reorderedID)
        originalIndices[reorderedID] = Int32(truncatingIfNeeded: originalID)
        
        let virtualSiteID = virtualSiteCount + originalID
        virtualSiteIndices.append(Int32(truncatingIfNeeded: originalID))
        originalIndices[virtualSiteID] = Int32(truncatingIfNeeded: originalID)
      } else {
        let reorderedID = virtualSiteCount + originalID
        reorderedIndices[originalID] = Int32(truncatingIfNeeded: reorderedID)
        originalIndices[reorderedID] = Int32(truncatingIfNeeded: originalID)
      }
    }
  }
  
  func createMasses(parameters: MM4Parameters) {
    let masses = parameters.atoms.masses
    for reorderedID in reorderedIndices {
      var mass: Double = 0
      if reorderedID > -1 {
        let originalID = originalIndices[Int(reorderedID)]
        mass = Double(masses[Int(originalID)])
      }
      system.addParticle(mass: mass)
    }
    
    fatalError("Change this function to account for virtual sites.")
  }
  
  func createVirtualSites(parameters: MM4Parameters) {
    fatalError("Not implemented.")
  }
}
