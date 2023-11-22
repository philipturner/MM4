//
//  MM4Parameters+Utilities.swift
//
//
//  Created by Philip Turner on 10/13/23.
//

extension MM4Parameters {
  @inline(__always)
  func other<T: FixedWidthInteger, U: FixedWidthInteger>(
    atomID: T, bondID: U
  ) -> Int32 {
    let bond = bondsToAtomsMap[Int(bondID)]
    return (bond[0] == atomID) ? bond[1] : bond[0]
  }
  
  func createAddress<T: FixedWidthInteger>(_ atomID: T) -> MM4Address {
    MM4Address(
      rigidBodyIndex: 0,
      atomIndex: UInt32(atomID),
      atomicNumber: atoms.atomicNumbers[Int(atomID)])
  }
  
  func createAddresses(_ map: SIMD4<Int32>) -> [MM4Address] {
    var output: [MM4Address] = []
    for lane in 0..<4 where map[lane] != -1 {
      output.append(createAddress(map[lane]))
    }
    return output
  }
  
  /// Overloads the function `createAtomCodes()`, but is a more elegant API,
  /// that follows the de facto naming convention established for MM4 and
  /// related projects.
  @inline(__always)
  func createAtomCodes<T: SIMD, U: SIMD>(group: T, zero: U) -> U
  where T.Scalar == UInt32, U.Scalar == UInt8 {
    var codes = zero
    for lane in 0..<T.scalarCount {
      let atomID = group[lane]
      codes[lane] = atoms.codes[Int(atomID)].rawValue
    }
    return codes
  }
  
  /// Fetch atomic numbers from an atom-to-atom map, some of which may be -1.
  /// Return invalid atoms as zero.
  @inline(__always)
  func createAtomicNumbers(map: SIMD4<Int32>) -> SIMD4<UInt8> {
    var atomicNumbers: SIMD4<UInt8> = .zero
    for lane in 0..<4 where map[lane] != -1 {
      let atomID = map[lane]
      atomicNumbers[lane] = atoms.atomicNumbers[Int(atomID)]
    }
    return atomicNumbers
  }
  
  static let excluding5RingElements: [UInt8] = [11, 19, 25, 31]
  
  @inline(__always)
  func with5RingsRemoved<T: SIMD>(_ codes: () -> T) -> T
  where T.Scalar == UInt8 {
    var output = codes()
    for element in Self.excluding5RingElements {
      if any(output .== element) {
        output.replace(with: .init(repeating: 1), where: output .== 123)
      }
    }
    return output
  }
  
  static let nonCarbonElements: [UInt8] = [6, 8, 11, 15, 19, 25, 31]
  
  @inline(__always)
  func containsTwoNonCarbons<T: SIMD>(_ codes: T) -> Bool
  where T.Scalar == UInt8 {
    var with5RingsRemoved = codes
    with5RingsRemoved.replace(with: .init(repeating: 1), where: codes .== 123)
    
    var nonCarbonElementCount: Int = 0
    for element in Self.nonCarbonElements {
      if any(with5RingsRemoved .== element) {
        nonCarbonElementCount += 1
      }
    }
    return nonCarbonElementCount > 1
  }
  
  /// `codes` could also contain atom indices, for sorting the angle while
  /// generating the bond topology.
  @inline(__always)
  func sortBond<T>(_ codes: SIMD2<T>) -> SIMD2<T>
  where T: FixedWidthInteger {
    if codes[0] > codes[1] {
      return SIMD2(codes[1], codes[0])
    } else {
      return codes
    }
  }
  
  /// `codes` could also contain atom indices, for sorting the angle while
  /// generating the bond topology.
  @inline(__always)
  func sortAngle<T>(_ codes: SIMD3<T>) -> SIMD3<T>
  where T: FixedWidthInteger {
    if codes[0] > codes[2] {
      return SIMD3(codes[2], codes[1], codes[0])
    } else {
      return codes
    }
  }
  
  /// `codes` could also contain atom indices, for sorting the torsion while
  /// generating the bond topology.
  @inline(__always)
  func sortTorsion<T>(_ codes: SIMD4<T>) -> SIMD4<T>
  where T: FixedWidthInteger {
    var reorder = false
    if codes[1] > codes[2] {
      reorder = true
    } else if codes[1] == codes[2] {
      if codes[0] > codes[3] {
        reorder = true
      }
    }
    
    if reorder {
      return SIMD4(codes[3], codes[2], codes[1], codes[0])
    } else {
      return codes
    }
  }
}

