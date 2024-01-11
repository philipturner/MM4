//
//  MM4System+VirtualSites.swift
//  MM4
//
//  Created by Philip Turner on 11/17/23.
//

import OpenMM

extension MM4System {
  func createReorderedIndices() {
    for atomID in parameters.atoms.indices {
      let reorderedID = atomID + virtualSiteCount
      reorderedIndices.append(UInt32(truncatingIfNeeded: reorderedID))
      
      if parameters.atoms.atomicNumbers[atomID] == 1 {
        virtualSiteCount += 1
      }
    }
  }
  
  func createMasses() {
    for atomID in parameters.atoms.indices {
      let mass = parameters.atoms.masses[atomID]
      system.addParticle(mass: Double(mass))
      
      let atomicNumber = parameters.atoms.atomicNumbers[atomID]
      if atomicNumber == 1 {
        system.addParticle(mass: 0)
      }
    }
  }
  
  func createVirtualSites() {
    for hydrogenID in parameters.atoms.indices {
      guard parameters.atoms.atomicNumbers[hydrogenID] == 1 else {
        continue
      }
      let map = parameters.atomsToBondsMap[Int(hydrogenID)]
      guard map[0] != -1, map[1] == -1, map[2] == -1, map[3] == -1 else {
        fatalError("Invalid virtual site.")
      }
      
      let bondID = map[0]
      let bond = parameters.bonds.indices[Int(bondID)]
      let otherID = (bond[0] == hydrogenID) ? bond[1] : bond[0]
      
      let otherParameters = parameters.atoms.parameters[Int(otherID)]
      let reductionFactor = Double(otherParameters.hydrogenReductionFactor)
      let weights = SIMD2(1 - reductionFactor, reductionFactor)
      
      let reordered = self.reorder(SIMD2(
        UInt32(truncatingIfNeeded: otherID),
        UInt32(truncatingIfNeeded: hydrogenID)))
      let virtualSite = OpenMM_TwoParticleAverageSite(
        particles: reordered, weights: weights)
      system.setVirtualSite(virtualSite, index: Int(reordered[1] &+ 1))
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
  }
  
  @_transparent
  func virtualSiteReorder(_ index: UInt32) -> Int {
    var reorderedID = reorderedIndices[Int(index)]
    if parameters.atoms.atomicNumbers[Int(index)] == 1 {
      reorderedID &+= 1
    }
    return Int(reorderedID)
  }
  
  @_transparent
  func virtualSiteReorder(_ indices: SIMD2<UInt32>) -> SIMD2<Int> {
    let reordered0 = virtualSiteReorder(indices[0])
    let reordered1 = virtualSiteReorder(indices[1])
    return SIMD2(reordered0, reordered1)
  }
  
  @_transparent
  func virtualSiteReorder(_ indices: SIMD4<UInt32>) -> SIMD4<Int> {
    let reordered0 = virtualSiteReorder(indices[0])
    let reordered1 = virtualSiteReorder(indices[1])
    let reordered2 = virtualSiteReorder(indices[2])
    let reordered3 = virtualSiteReorder(indices[3])
    return SIMD4(reordered0, reordered1, reordered2, reordered3)
  }
}
