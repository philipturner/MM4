//
//  MM4Parameters+Topology.swift
//
//
//  Created by Philip Turner on 10/7/23.
//

// MARK: - Functions for generating the topology and assigning ring types.

/// Parameters for a group of 5 or less atoms.
public struct MM4Rings {
  /// Groups of atom indices that form a ring.
  ///
  /// The unused vector lanes are set to `UInt32.max`.
  public internal(set) var indices: [SIMD8<UInt32>] = []
  
  /// Map from a group of atoms to a ring index.
  public internal(set) var map: [SIMD8<UInt32>: UInt32] = [:]
  
  /// The number of atoms in the ring.
  public internal(set) var ringTypes: [UInt8] = []
}

extension MM4Parameters {
  func createTopology() throws {
    // Traverse the bond topology.
    for atom1 in 0..<UInt32(atoms.count) {
      let map1 = atomsToAtomsMap[Int(atom1)]
      var ringType: UInt8 = 6
      
      for lane2 in 0..<4 where map1[lane2] != -1 {
        let atom2 = UInt32(truncatingIfNeeded: map1[lane2])
        let map2 = atomsToAtomsMap[Int(atom2)]
        
        for lane3 in 0..<4 where map2[lane3] != -1 {
          let atom3 = UInt32(truncatingIfNeeded: map2[lane3])
          if atom1 == atom3 { continue }
          if atom1 < atom3 {
            angles.map[SIMD3(atom1, atom2, atom3)] = .max
          }
          let map3 = atomsToAtomsMap[Int(atom3)]
          
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
            
            let map4 = atomsToAtomsMap[Int(atom4)]
            var maskA: SIMD4<Int32> = .init(repeating: 6)
            var maskB: SIMD4<Int32> = .init(repeating: 6)
            var maskC: SIMD4<Int32> = .init(repeating: 6)
            var maskD: SIMD4<Int32> = .init(repeating: 6)
            let mapA = atomsToAtomsMap[Int(map4[0])]
            let mapB = atomsToAtomsMap[Int(map4[1])]
            let mapC = atomsToAtomsMap[Int(map4[2])]
            let mapD = atomsToAtomsMap[Int(map4[3])]
            
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
          
          if _slowPath(ringType < 5) {
            func makeAddress(_ atomID: UInt32) -> MM4Address {
              MM4Address(
                rigidBodyIndex: 0,
                atomIndex: atomID,
                atomicNumber: atoms.atomicNumbers[Int(atomID)])
            }
            var addresses: [MM4Address] = []
            addresses.append(makeAddress(atom1))
            addresses.append(makeAddress(atom2))
            addresses.append(makeAddress(atom3))
            
            for lane4 in 0..<4 where map3[lane4] != -1 {
              let atom4 = UInt32(truncatingIfNeeded: map3[lane4])
              if atom2 == atom4 {
                continue
              } else if atom1 == atom4 {
                throw MM4Error.unsupportedRing(addresses, 3)
              }
              
              let map4 = atomsToAtomsMap[Int(atom4)]
              for lane5 in 0..<4 where atom1 == map4[lane5] {
                addresses.append(makeAddress(atom4))
                throw MM4Error.unsupportedRing(addresses, 4)
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
      let map1 = atomsToAtomsMap[Int(torsion[0])]
      let map4 = atomsToAtomsMap[Int(torsion[1])]
      
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
}

