//
//  MM4System+VirtualSites.swift
//
//
//  Created by Philip Turner on 11/17/23.
//

import OpenMM

extension MM4System {
  func createReorderedIndices() {
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
  
  func createMasses() {
    for (index, reorderedID) in reorderedIndices.enumerated() {
      let originalID = originalIndices[Int(reorderedID)]
      let atomicNumber = parameters.atoms.atomicNumbers[Int(originalID)]
      var mass = parameters.atoms.masses[Int(originalID)]
      
      if index >= virtualSiteIndices.count {
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
    fatalError("Not implemented.")
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
    return virtualSiteIndices.count + index
  }
  
  @inline(__always)
  func virtualSiteReorder(_ indices: SIMD2<Int32>) -> SIMD2<Int> {
    var output: SIMD2<Int32> = indices
    let virtualSiteCount = Int32(truncatingIfNeeded: virtualSiteIndices.count)
    output &+= SIMD2(repeating: virtualSiteCount)
    return SIMD2<Int>(truncatingIfNeeded: output)
  }
  
  @inline(__always)
  func virtualSiteReorder(_ indices: SIMD4<Int32>) -> SIMD4<Int> {
    var output: SIMD4<Int32> = indices
    let virtualSiteCount = Int32(truncatingIfNeeded: virtualSiteIndices.count)
    output &+= SIMD4(repeating: virtualSiteCount)
    return SIMD4<Int>(truncatingIfNeeded: output)
  }
}
