//
//  MM4Parameters+Topology.swift
//  MM4
//
//  Created by Philip Turner on 10/7/23.
//

import Atomics
import Dispatch

/// Parameters for a group of 3 to 5 atoms.
public struct MM4Rings {
  /// Groups of atom indices that form a ring.
  ///
  /// The unused vector lanes are set to `UInt32.max`.
  public var indices: [SIMD8<UInt32>] = []
  
  /// Map from a group of atoms to a ring index.
  public var map: [SIMD8<UInt32>: UInt32] = [:]
  
  /// The number of atoms in the ring.
  public var ringTypes: [UInt8] = []
  
  mutating func append(contentsOf other: Self, atomOffset: UInt32) {
    let ringOffset = UInt32(self.indices.count)
    self.indices += other.indices.map {
      let modified = $0 &+ atomOffset
      return $0.replacing(with: modified, where: $0 .< UInt32.max)
    }
    for var key in other.map.keys {
      let value = other.map[key].unsafelyUnwrapped
      let modifiedKey = key &+ atomOffset
      key.replace(with: modifiedKey, where: key .< UInt32.max)
      self.map[key] = value &+ ringOffset
    }
    self.ringTypes += other.ringTypes
  }
}

extension MM4Parameters {
  /// - throws: `.openValenceShell`
  mutating func createAtomsToBondsMap() throws {
    atomsToBondsMap = Array(
      repeating: SIMD4(repeating: -1), count: atoms.count)
    
    for bondID in 0..<bonds.indices.count {
      // Initialize the bonds map, so it can be used later when generating the
      // topology of angles/torsions/rings.
      let bond = bonds.indices[bondID]
      bonds.map[bond] = .max
      
      for j in 0..<2 {
        let atomID = Int(bond[j])
        var map = atomsToBondsMap[atomID]
        var succeeded = false
        for k in 0..<4 {
          if map[k] == -1 {
            map[k] = Int32(bondID)
            succeeded = true
            break
          }
        }
        
        if !succeeded {
          var neighbors: [MM4Address] = []
          for lane in 0..<4 where map[lane] != -1 {
            let bond = bonds.indices[Int(map[lane])]
            let otherID = (bond[0] == bond[j]) ? bond[1] : bond[0]
            neighbors.append(createAddress(otherID))
          }
          neighbors.append(createAddress(bond[1 - j]))
          
          // There is an error because we exceeded the bond capacity (4).
          let address = createAddress(bond[j])
          throw MM4Error.openValenceShell(address, neighbors)
        }
        atomsToBondsMap[atomID] = map
      }
    }
  }
  
  /// - throws: Nothing
  mutating func createAtomsToAtomsMap() throws {
    atomsToAtomsMap.reserveCapacity(atoms.count)
    for atomID in 0..<atoms.count {
      let bondsMap = atomsToBondsMap[atomID]
      var atomsMap = SIMD4<Int32>(repeating: -1)
      
      for lane in 0..<4 where bondsMap[lane] != -1 {
        let bond = bonds.indices[Int(bondsMap[lane])]
        let otherID = (bond[0] == atomID) ? bond[1] : bond[0]
        atomsMap[lane] = Int32(truncatingIfNeeded: otherID)
      }
      
      var duplicates: Int32 = .zero
      for lane in 0..<4 where atomsMap[lane] != -1 {
        let zero = SIMD4<Int32>(repeating: 0)
        let one = SIMD4<Int32>(repeating: 1)
        
        var matchMask = zero
        matchMask.replace(with: one, where: atomsMap .== atomsMap[lane])
        matchMask[lane] = 0
        duplicates &+= matchMask.wrappedSum()
      }
      
      // If the bonds are all unique, we can employ some interesting tricks
      // during the topology search. We can pre-allocate a fixed amount of
      // memory for each atom's list of angles and torsions.
      if duplicates > 0 {
        fatalError("The same bond was entered twice.")
      }
      
      atomsToAtomsMap.append(atomsMap)
    }
  }
  
  /// - throws: `.unsupportedRing`
  mutating func createTopology(forces: MM4ForceOptions) throws {
    // Map from atoms to connected atoms that can be efficiently traversed.
    var vAtomsToAtomsMap: UnsafeMutablePointer<SIMD4<Int32>>
    vAtomsToAtomsMap = .allocate(capacity: atoms.count + 1)
    vAtomsToAtomsMap += 1
    vAtomsToAtomsMap[-1] = SIMD4(repeating: -1)
    vAtomsToAtomsMap.initialize(from: atomsToAtomsMap, count: atoms.count)
    defer {
      (vAtomsToAtomsMap - 1).deallocate()
    }
    
    let includeTorsions =
    forces.contains(.torsion) ||
    forces.contains(.torsionBend) ||
    forces.contains(.torsionStretch)
    
    // Use a conservative metric to determine whether angles are included. If
    // torsions are included but angles aren't, there's still some torsion
    // cross-terms that depend on equilibrium angle.
    let includeAngles =
    forces.contains(.bend) ||
    forces.contains(.bendBend) ||
    forces.contains(.stretchBend) ||
    forces.contains(.stretchStretch) ||
    includeTorsions
    
    let angleAtomics: UnsafeMutablePointer<UInt16.AtomicRepresentation> =
      .allocate(capacity: atoms.count)
    let angleCounts = UnsafeMutablePointer<UInt16>(
      OpaquePointer(angleAtomics))
    let angleBuckets: UnsafeMutablePointer<SIMD3<UInt32>> =
      .allocate(capacity: 6 * atoms.count)
    
    let torsionAtomics: UnsafeMutablePointer<UInt16.AtomicRepresentation> =
      .allocate(capacity: atoms.count)
    let torsionCounts =  UnsafeMutablePointer<UInt16>(
      OpaquePointer(torsionAtomics))
    let torsionBuckets: UnsafeMutablePointer<SIMD4<UInt32>> =
      .allocate(capacity: 36 * atoms.count)
    
    angleCounts.initialize(repeating: .zero, count: atoms.count)
    torsionCounts.initialize(repeating: .zero, count: atoms.count)
    defer {
      angleAtomics.deallocate()
      angleBuckets.deallocate()
      torsionAtomics.deallocate()
      torsionBuckets.deallocate()
    }
    
    @_transparent
    func wrap(_ index: Int) -> Int {
      (index + 5) % 5
    }
    
    // Scope the rings map into a local dictionary per-thread. Merge the
    // partial results on single-core after the loop is over.
    let taskSize = 128
    let taskCount = (atoms.count + taskSize - 1) / taskSize
    var localRingsMaps = [[SIMD8<UInt32>: Bool]](
      repeating: [:], count: taskCount)
    var localErrors = [MM4Error?](repeating: nil, count: taskCount)
    
    DispatchQueue.concurrentPerform(iterations: taskCount) { z in
      do {
        try execute(taskID: z)
      } catch let error {
        localErrors[z] = error as? MM4Error
      }
    }
    for error in localErrors {
      if let error {
        throw error
      }
    }
    var ringsMap: [SIMD8<UInt32>: Bool] = [:]
    for localRingsMap in localRingsMaps {
      for key in localRingsMap.keys {
        ringsMap[key] = true
      }
    }
    
    func execute(taskID: Int) throws {
      var ringsMap: [SIMD8<UInt32>: Bool] = [:]
      defer {
        localRingsMaps[taskID] = ringsMap
      }
      
      let atomStart = taskID * taskSize
      let atomEnd = min(atomStart + taskSize, atoms.count)
      for atom1 in UInt32(atomStart)..<UInt32(atomEnd) {
        let map1 = vAtomsToAtomsMap[Int(atom1)]
        var ringType: UInt8 = 6
        
        for lane2 in 0..<4 where map1[lane2] != -1 {
          let atom2 = UInt32(truncatingIfNeeded: map1[lane2])
          let map2 = vAtomsToAtomsMap[Int(atom2)]
          
          for lane3 in 0..<4 where map2[lane3] != -1 {
            let atom3 = UInt32(truncatingIfNeeded: map2[lane3])
            if atom1 == atom3 { continue }
            if includeAngles, atom1 < atom3 {
              let angle = SIMD3(atom1, atom2, atom3)
              let atomID = Int(atom2)
              let atomic = UnsafeAtomic<UInt16>(
                at: angleAtomics.advanced(by: atomID))
              
              let count = atomic
                .loadThenWrappingIncrement(ordering: .relaxed)
              angleCounts[atomID] = count &+ 1
              angleBuckets[6 &* atomID &+ Int(count)] = angle
            }
            let map3 = vAtomsToAtomsMap[Int(atom3)]
            
            for lane4 in 0..<4 where map3[lane4] != -1 {
              let atom4 = UInt32(truncatingIfNeeded: map3[lane4])
              if atom2 == atom4 {
                continue
              } else if atom1 == atom4 {
                ringType = min(3, ringType)
                continue
              }
              
              let map4 = vAtomsToAtomsMap[Int(atom4)]
              var maskA: SIMD4<Int32> = .init(repeating: 6)
              var maskB: SIMD4<Int32> = .init(repeating: 6)
              var maskC: SIMD4<Int32> = .init(repeating: 6)
              var maskD: SIMD4<Int32> = .init(repeating: 6)
              let mapA = vAtomsToAtomsMap[Int(map4[0])]
              let mapB = vAtomsToAtomsMap[Int(map4[1])]
              let mapC = vAtomsToAtomsMap[Int(map4[2])]
              let mapD = vAtomsToAtomsMap[Int(map4[3])]
              
              let fives = SIMD4<Int32>(repeating: 5)
              let atom = Int32(truncatingIfNeeded: atom1)
              maskA.replace(with: fives, where: atom .== mapA)
              maskB.replace(with: fives, where: atom .== mapB)
              maskC.replace(with: fives, where: atom .== mapC)
              maskD.replace(with: fives, where: atom .== mapD)
              
              maskA.replace(with: maskB, where: maskB .< maskA)
              maskC.replace(with: maskD, where: maskD .< maskC)
              maskA.replace(with: maskC, where: maskC .< maskA)
              maskA.replace(with: .init(repeating: 4), where: atom .== map4)
              ringType = min(ringType, UInt8(truncatingIfNeeded: maskA.min()))
              
              if atom2 < atom3 {
                if includeTorsions {
                  let torsion = SIMD4(atom1, atom2, atom3, atom4)
                  let atomID = Int(atom2)
                  let atomic = UnsafeAtomic<UInt16>(
                    at: torsionAtomics.advanced(by: atomID))
                  
                  let count = atomic
                    .loadThenWrappingIncrement(ordering: .relaxed)
                  torsionCounts[atomID] = count &+ 1
                  torsionBuckets[36 &* atomID &+ Int(count)] = torsion
                }
                
                var match1: SIMD4<Int32> = .zero
                var match2: SIMD4<Int32> = .zero
                var match3: SIMD4<Int32> = .zero
                var match4: SIMD4<Int32> = .zero
                match1.replace(with: .one, where: map1[0] .== map4)
                match2.replace(with: .one, where: map1[1] .== map4)
                match3.replace(with: .one, where: map1[2] .== map4)
                match4.replace(with: .one, where: map1[3] .== map4)
                
                match1 &+= match2
                match3 &+= match4
                match1 &+= match3
                match1.replace(with: SIMD4.zero, where: map4 .== -1)
                
                // Atom indices within the angles/torsions are already sorted.
                // Here, we initialize rings with correct ordering as well.
                for lane in 0..<4 where match1[lane] > 0 {
                  let array = SIMD8<UInt32>(
                    atom1, atom2, atom3, atom4,
                    UInt32(truncatingIfNeeded: map4[lane]), 0, 0, 0)
                  
                  var minIndex = 0
                  for lane in 0..<5 {
                    if array[lane] < array[minIndex] {
                      minIndex = lane
                    }
                  }
                  
                  let prev = array[wrap(minIndex &- 1)]
                  let next = array[wrap(minIndex &+ 1)]
                  let increment = (next > prev) ? +1 : -1
                  
                  var output: SIMD8<UInt32> = .init(repeating: .max)
                  for lane in 0..<5 {
                    let index = wrap(minIndex &+ lane &* increment)
                    output[lane] = array[index]
                  }
                  ringsMap[output] = true
                }
              }
            }
            
            if ringType < 5 {
              var addresses: [MM4Address] = []
              addresses.append(createAddress(atom1))
              addresses.append(createAddress(atom2))
              addresses.append(createAddress(atom3))
              
              for lane4 in 0..<4 where map3[lane4] != -1 {
                let atom4 = UInt32(map3[lane4])
                if atom2 == atom4 {
                  continue
                } else if atom1 == atom4 {
                  throw MM4Error.unsupportedRing(addresses)
                }
                
                let map4 = vAtomsToAtomsMap[Int(atom4)]
                for lane5 in 0..<4 where atom1 == map4[lane5] {
                  addresses.append(createAddress(atom4))
                  throw MM4Error.unsupportedRing(addresses)
                }
              }
              fatalError(
                "Ring type was less than 5, but faulting atom group was not found.")
            }
          }
        }
      }
    }
    
    for atomID in atoms.indices {
      var angleBucket = UnsafeMutableBufferPointer(
        start: angleBuckets.advanced(by: 6 &* atomID),
        count: Int(angleCounts[atomID]))
      angleBucket.sort(by: { x, y in
        if x[0] != y[0] { return x[0] < y[0] }
        if x[2] != y[2] { return x[2] < y[2] }
        return true
      })
      angles.indices += angleBucket
      
      var torsionBucket = UnsafeMutableBufferPointer(
        start: torsionBuckets.advanced(by: 36 &* atomID),
        count: Int(torsionCounts[atomID]))
      torsionBucket.sort(by: { x, y in
        if x[2] != y[2] { return x[2] < y[2] }
        if x[0] != y[0] { return x[0] < y[0] }
        if x[3] != y[3] { return x[3] < y[3] }
        return true
      })
      torsions.indices += torsionBucket
    }
    
    rings.indices = ringsMap.keys.map { $0 }
    rings.indices.sort(by: compareRing)
    
    atoms.ringTypes = .init(repeating: 6, count: atoms.count)
    bonds.ringTypes = .init(repeating: 6, count: bonds.indices.count)
    angles.ringTypes = .init(repeating: 6, count: angles.indices.count)
    torsions.ringTypes = .init(repeating: 6, count: torsions.indices.count)
    rings.ringTypes = .init(repeating: 6, count: rings.indices.count)
    
    guard bonds.indices.count < Int32.max,
          angles.indices.count < Int32.max,
          torsions.indices.count < Int32.max,
          rings.indices.count < Int32.max else {
      fatalError("Too many bonds, angles, torsions, or rings.")
    }
    for (index, angle) in bonds.indices.enumerated() {
      bonds.map[angle] = UInt32(truncatingIfNeeded: index)
    }
    for (index, angle) in angles.indices.enumerated() {
      angles.map[angle] = UInt32(truncatingIfNeeded: index)
    }
    for (index, torsion) in torsions.indices.enumerated() {
      torsions.map[torsion] = UInt32(truncatingIfNeeded: index)
    }
    for (index, ring) in rings.indices.enumerated() {
      rings.map[ring] = UInt32(truncatingIfNeeded: index)
    }
    
    for ringID in rings.indices.indices {
      let ring = rings.indices[ringID]
      for lane in 0..<5 {
        let atomID = ring[lane]
        let unsortedBond = SIMD2(atomID, ring[wrap(lane &+ 1)])
        let unsortedAngle = SIMD3(unsortedBond, ring[wrap(lane &+ 2)])
        let unsortedTorsion = SIMD4(unsortedAngle, ring[wrap(lane &+ 3)])
        let bond = sortBond(unsortedBond)
        let angle = sortAngle(unsortedAngle)
        let torsion = sortTorsion(unsortedTorsion)
        
        guard atomID < .max, let bondID = bonds.map[bond] else {
          fatalError("Invalid atom or bond in ring.")
        }
        atoms.ringTypes[Int(atomID)] = 5
        bonds.ringTypes[Int(bondID)] = 5
        
        if includeAngles, let angleID = angles.map[angle]  {
          angles.ringTypes[Int(angleID)] = 5
        }
        if includeTorsions, let torsionID = torsions.map[torsion] {
          torsions.ringTypes[Int(torsionID)] = 5
        }
      }
      rings.ringTypes[ringID] = 5
    }
  }
  
  /// - throws: `.unsupportedCenterType`
  mutating func createCenterTypes() throws {
    let permittedAtomicNumbers: [UInt8] = [6, 7, 8, 14, 15, 16, 32]
    let blacklistedAtomicNumbers: [UInt8] = [1, 9]
    for atomID in atoms.indices {
      let atomicNumber = atoms.atomicNumbers[atomID]
      guard permittedAtomicNumbers.contains(atomicNumber) else {
        precondition(
          blacklistedAtomicNumbers.contains(atomicNumber),
          "Atomic number \(atomicNumber) not recognized.")
        atoms.centerTypes.append(nil)
        continue
      }
      
      let map = atomsToAtomsMap[atomID]
      var otherElements: SIMD4<UInt8> = .zero
      for lane in 0..<4 where map[lane] != -1 {
        if map[lane] == -1 {
          otherElements[lane] = 1
        } else {
          let otherID = map[lane]
          otherElements[lane] = atoms.atomicNumbers[Int(otherID)]
        }
      }
      
      let halogenMask =
      (otherElements .== 1) .|
      (otherElements .== 9) .|
      (otherElements .== 17) .|
      (otherElements .== 35) .|
      (otherElements .== 53)
      
      if all(halogenMask) {
        let address = createAddress(atomID)
        let neighbors = createAddresses(map)
        throw MM4Error.unsupportedCenterType(address, neighbors)
      }
      
      // In MM4, fluorine is treated like carbon when determining carbon types.
      // Allinger notes this may be a weakness of the forcefield. This idea has
      // been extended to encapsulate all non-hydrogen atoms.
      var matchMask: SIMD4<UInt8> = .zero
      matchMask.replace(with: .one, where: otherElements .!= 1)
      
      var carbonType: MM4CenterType
      switch matchMask.wrappedSum() {
      case 4:
        carbonType = .quaternary
      case 3:
        carbonType = .tertiary
      case 2:
        carbonType = .secondary
      case 1:
        carbonType = .primary
      default:
        fatalError("This should never happen.")
      }
      atoms.centerTypes.append(carbonType)
    }
  }
}
