//
//  MM4System+VirtualSites.swift
//  MM4
//
//  Created by Philip Turner on 11/17/23.
//

import OpenMM

extension MM4System {
  func createReorderedIndices() {
    let atomicNumbers = parameters.atoms.atomicNumbers
    for originalID in 0..<atomicNumbers.count {
      if atomicNumbers[originalID] == 1 {
        virtualSiteCount += 1
      }
    }
    
    let particleCount = virtualSiteCount + atomicNumbers.count
    self.reorderedIndices = Array(repeating: .max, count: atomicNumbers.count)
    self.originalIndices = Array(repeating: .max, count: particleCount)
    
    var virtualSitePointer = 0
    for originalID in 0..<atomicNumbers.count {
      if atomicNumbers[originalID] == 1 {
        let reorderedID = virtualSitePointer
        virtualSitePointer += 1
        reorderedIndices[originalID] = UInt32(truncatingIfNeeded: reorderedID)
        originalIndices[reorderedID] = UInt32(truncatingIfNeeded: originalID)
        
        let virtualSiteID = virtualSiteCount + originalID
        originalIndices[virtualSiteID] = UInt32(truncatingIfNeeded: originalID)
      } else {
        let reorderedID = virtualSiteCount + originalID
        reorderedIndices[originalID] = UInt32(truncatingIfNeeded: reorderedID)
        originalIndices[reorderedID] = UInt32(truncatingIfNeeded: originalID)
      }
    }
    guard virtualSitePointer == virtualSiteCount else {
      fatalError("This should never happen.")
    }
  }
  
  func createMasses() {
    for (index, originalID) in originalIndices.enumerated() {
      var mass = parameters.atoms.masses[Int(originalID)]
      let atomicNumber = parameters.atoms.atomicNumbers[Int(originalID)]
      if index >= virtualSiteCount {
        if atomicNumber == 1 {
          mass = 0
        }
      } else {
        guard atomicNumber == 1 else {
          fatalError("This should never happen.")
        }
      }
      system.addParticle(mass: Double(mass))
    }
  }
  
  func createVirtualSites() {
    for index in 0..<virtualSiteCount {
      let originalID = originalIndices[index]
      let map = parameters.atomsToBondsMap[Int(originalID)]
      guard map[0] != -1, map[1] == -1, map[2] == -1, map[3] == -1 else {
        fatalError("Invalid virtual site.")
      }
      
      let bondID = map[0]
      let otherID = parameters.other(atomID: originalID, bondID: bondID)
      
      let otherParameters = parameters.atoms.parameters[Int(otherID)]
      let reductionFactor = Double(otherParameters.hydrogenReductionFactor)
      let weights = SIMD2(1 - reductionFactor, reductionFactor)
      
      let reordered = self.reorder(SIMD2(
        UInt32(truncatingIfNeeded: otherID), originalID))
      let virtualSite = OpenMM_TwoParticleAverageSite(
        particles: reordered, weights: weights)
      let virtualSiteID = virtualSiteCount + Int(originalID)
      system.setVirtualSite(virtualSite, index: virtualSiteID)
    }
  }
  
  func createVirtualSiteMask() {
    virtualSiteMask.reserveCapacity(reorderedIndices.count)
    for reorderedID in reorderedIndices.indices {
      if reorderedID >= virtualSiteCount {
        let originalID = Int(originalIndices[reorderedID])
        let atomicNumber = parameters.atoms.atomicNumbers[originalID]
        if atomicNumber == 1 {
          virtualSiteMask.append(false)
          continue
        }
      }
      virtualSiteMask.append(true)
    }
  }
}

extension MM4System {
  func createExceptions(force: OpenMM_CustomNonbondedForce) {
    for bond in parameters.bonds.indices {
      let reordered = self.virtualSiteReorder(bond)
      force.addExclusion(particles: reordered)
    }
    for exception in parameters.nonbondedExceptions13 {
      let reordered = self.virtualSiteReorder(exception)
      force.addExclusion(particles: reordered)
    }
    
    let group = OpenMM_IntSet()
    for atomID in parameters.atoms.indices {
      let reordered = self.virtualSiteReorder(atomID)
      group.insert(reordered)
    }
    force.addInteractionGroup(set1: group, set2: group)
  }
  
  @inline(__always)
  func virtualSiteReorder(_ index: Int) -> Int {
    return virtualSiteCount + index
  }
  
  @inline(__always)
  func virtualSiteReorder(_ indices: SIMD2<UInt32>) -> SIMD2<Int> {
    var output: SIMD2<UInt32> = indices
    let virtualSiteCount = UInt32(truncatingIfNeeded: virtualSiteCount)
    output &+= SIMD2(repeating: virtualSiteCount)
    return SIMD2<Int>(truncatingIfNeeded: output)
  }
  
  @inline(__always)
  func virtualSiteReorder(_ indices: SIMD4<UInt32>) -> SIMD4<Int> {
    var output: SIMD4<UInt32> = indices
    let virtualSiteCount = UInt32(truncatingIfNeeded: virtualSiteCount)
    output &+= SIMD4(repeating: virtualSiteCount)
    return SIMD4<Int>(truncatingIfNeeded: output)
  }
}
