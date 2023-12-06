//
//  MM4Parameters+Topology.swift
//
//
//  Created by Philip Turner on 10/7/23.
//

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
      $0 &+ atomOffset
    }
    for key in other.map.keys {
      let value = other.map[key].unsafelyUnwrapped
      self.map[key &+ atomOffset] = value &+ ringOffset
    }
    self.ringTypes += other.ringTypes
  }
  
  mutating func reserveCapacity(_ minimumCapacity: Int) {
    indices.reserveCapacity(minimumCapacity)
    map.reserveCapacity(minimumCapacity)
    ringTypes.reserveCapacity(minimumCapacity)
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
            let bondID = Int(map[lane])
            let bond = bonds.indices[Int(map[lane])]
            let otherID = other(atomID: bond[j], bondID: bondID)
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
        let otherID = other(atomID: atomID, bondID: bondsMap[lane])
        atomsMap[lane] = Int32(truncatingIfNeeded: otherID)
      }
      atomsToAtomsMap.append(atomsMap)
    }
  }
  
  /// - throws: `.unsupportedRing`
  mutating func createTopology() throws {
    // Map from atoms to connected atoms that can be efficiently traversed.
    var vAtomsToAtomsMap: UnsafeMutablePointer<SIMD4<Int32>>
    vAtomsToAtomsMap = .allocate(capacity: atoms.count + 1)
    vAtomsToAtomsMap += 1
    vAtomsToAtomsMap[-1] = SIMD4(repeating: -1)
    vAtomsToAtomsMap.initialize(from: atomsToAtomsMap, count: atoms.count)
    defer {
      (vAtomsToAtomsMap - 1).deallocate()
    }
    
    for atom1 in 0..<UInt32(atoms.count) {
      let map1 = vAtomsToAtomsMap[Int(atom1)]
      var ringType: UInt8 = 6
      
      for lane2 in 0..<4 where map1[lane2] != -1 {
        let atom2 = UInt32(truncatingIfNeeded: map1[lane2])
        let map2 = vAtomsToAtomsMap[Int(atom2)]
        
        for lane3 in 0..<4 where map2[lane3] != -1 {
          let atom3 = UInt32(truncatingIfNeeded: map2[lane3])
          if atom1 == atom3 { continue }
          if atom1 < atom3 {
            angles.map[SIMD3(atom1, atom2, atom3)] = .max
          }
          let map3 = vAtomsToAtomsMap[Int(atom3)]
          
          for lane4 in 0..<4 where map3[lane4] != -1 {
            let atom4 = UInt32(truncatingIfNeeded: map3[lane4])
            if atom2 == atom4 {
              continue
            } else if atom1 == atom4 {
              ringType = min(3, ringType)
              continue
            } else if atom2 < atom3 {
              torsions.map[SIMD4(atom1, atom2, atom3, atom4)] = .max
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
          }
          
          if ringType != 6 {
            print("Found a ring at atom \(atom1)")
          }
          
          if _slowPath(ringType < 5) {
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
            fatalError("Ring type was less than 5, but faulting atom group was not found.")
          }
        }
      }
    }
    angles.indices = angles.map.keys.map { $0 }
    torsions.indices = torsions.map.keys.map { $0 }
    
    func wrap(_ index: Int) -> Int {
      (index + 5) % 5
    }
    var ringsMap: [SIMD8<UInt32>: Bool] = [:]
    for torsion in torsions.indices {
      // Mask out the -1 indices, then check whether any atoms from the first
      // atom's map match the fourth atom's map.
      let map1 = vAtomsToAtomsMap[Int(torsion[0])]
      let map4 = vAtomsToAtomsMap[Int(torsion[1])]
      
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
      
      for lane in 0..<4 where match1[lane] > 0 {
        var array: [UInt32] = []
        array.reserveCapacity(5)
        for i in 0..<4 {
          array.append(torsion[i])
        }
        array.append(UInt32(truncatingIfNeeded: map4[lane]))
        
        // Create a sorted list of atom indices.
        let minIndex = (0..<5).min(by: { array[$0] < array[$1] })!
        let prev = array[wrap(minIndex - 1)]
        let next = array[wrap(minIndex + 1)]
        let increment = (next > prev) ? +1 : -1
        
        var output: SIMD8<UInt32> = .init(repeating: .max)
        for lane in 0..<5 {
          let index = wrap(minIndex + lane * increment)
          output[lane] = array[index]
        }
        ringsMap[output] = true
      }
    }
    
    rings.indices = ringsMap.keys.map { $0 }
    atoms.ringTypes = .init(repeating: 6, count: atoms.count)
    bonds.ringTypes = .init(repeating: 6, count: bonds.indices.count)
    angles.ringTypes = .init(repeating: 6, count: angles.indices.count)
    torsions.ringTypes = .init(repeating: 6, count: torsions.indices.count)
    
    guard bonds.indices.count < Int32.max,
          angles.indices.count < Int32.max,
          torsions.indices.count < Int32.max,
          rings.indices.count < Int32.max else {
      fatalError("Too many bonds, angles, torsions, or rings.")
    }
    for (index, angle) in bonds.indices.enumerated() {
      bonds.map[angle]! = UInt32(truncatingIfNeeded: index)
    }
    for (index, angle) in angles.indices.enumerated() {
      angles.map[angle]! = UInt32(truncatingIfNeeded: index)
    }
    for (index, torsion) in torsions.indices.enumerated() {
      torsions.map[torsion]! = UInt32(truncatingIfNeeded: index)
    }
    for (index, ring) in rings.indices.enumerated() {
      rings.map[ring]! = UInt32(truncatingIfNeeded: index)
    }
    
    for ring in rings.indices {
      for lane in 0..<5 {
        let atomID = ring[lane]
        let bond = sortBond(SIMD2(atomID, ring[wrap(lane + 1)]))
        let angle = sortAngle(SIMD3(bond, ring[wrap(lane + 2)]))
        let torsion = sortTorsion(SIMD4(angle, ring[wrap(lane + 3)]))
        
        guard atomID < .max,
              let bondID = bonds.map[bond],
              let angleID = angles.map[angle],
              let torsionID = torsions.map[torsion] else {
          fatalError("Invalid atom, bond, angle, or torsion in ring.")
        }
        atoms.ringTypes[Int(atomID)] = 5
        bonds.ringTypes[Int(bondID)] = 5
        angles.ringTypes[Int(angleID)] = 5
        torsions.ringTypes[Int(torsionID)] = 5
      }
    }
  }
  
  /// - throws: `.unsupportedCenterType`
  mutating func createCenterTypes() throws {
    let permittedAtomicNumbers: [UInt8] = [6, 7, 8, 14, 15, 16, 32]
    for atomID in atoms.indices {
      let atomicNumber = atoms.atomicNumbers[atomID]
      guard permittedAtomicNumbers.contains(atomicNumber) else {
        precondition(
          atomicNumber == 1 || atomicNumber == 9,
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
