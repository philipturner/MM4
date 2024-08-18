//
//  MM4Parameters+Utilities.swift
//  MM4
//
//  Created by Philip Turner on 10/13/23.
//

// MARK: - Locating Atoms

extension MM4Parameters {
  func createAddress<T: FixedWidthInteger>(_ atomID: T) -> MM4Address {
    MM4Address(
      index: UInt32(atomID),
      atomicNumber: atoms.atomicNumbers[Int(atomID)])
  }
  
  func createAddresses(_ map: SIMD4<Int32>) -> [MM4Address] {
    var output: [MM4Address] = []
    for lane in 0..<4 where map[lane] != -1 {
      output.append(createAddress(map[lane]))
    }
    return output
  }
}

// MARK: - Handling Different Elements

extension MM4Parameters {
  /// Overloads the function `createAtomCodes()`, but is a more elegant API,
  /// that follows the de facto naming convention established for MM4 and
  /// related projects.
  @_transparent
  @_specialize(where T == SIMD2<UInt32>, U == SIMD2<UInt8>)
  @_specialize(where T == SIMD3<UInt32>, U == SIMD3<UInt8>)
  @_specialize(where T == SIMD4<UInt32>, U == SIMD4<UInt8>)
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
  @_transparent
  func createAtomicNumbers(map: SIMD4<Int32>) -> SIMD4<UInt8> {
    var atomicNumbers: SIMD4<UInt8> = .zero
    for lane in 0..<4 where map[lane] != -1 {
      let atomID = map[lane]
      atomicNumbers[lane] = atoms.atomicNumbers[Int(atomID)]
    }
    return atomicNumbers
  }
  
  static let excluding5RingElements: [UInt8] = [11, 19, 25, 31]
  
  @_transparent
  @_specialize(where T == SIMD2<UInt8>)
  @_specialize(where T == SIMD3<UInt8>)
  @_specialize(where T == SIMD4<UInt8>)
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
  
  @_transparent
  @_specialize(where T == SIMD2<UInt8>)
  @_specialize(where T == SIMD3<UInt8>)
  @_specialize(where T == SIMD4<UInt8>)
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
}

// MARK: - Sorting Within Bonds

extension MM4Parameters {
  /// `codes` could also contain atom indices, for sorting the bond while
  /// generating the bond topology.
  @_transparent
  @_specialize(where T == UInt8)
  @_specialize(where T == UInt32)
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
  @_transparent
  @_specialize(where T == UInt8)
  @_specialize(where T == UInt32)
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
  @_transparent
  @_specialize(where T == UInt8)
  @_specialize(where T == UInt32)
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

// MARK: - Sorting Between Bonds

extension MM4Parameters {
  @_transparent
  @_specialize(where T == UInt8)
  @_specialize(where T == UInt32)
  func compareAngle<T>(_ x: SIMD3<T>, _ y: SIMD3<T>) -> Bool
  where T: FixedWidthInteger {
    if x[1] != y[1] { return x[1] < y[1] }
    if x[0] != y[0] { return x[0] < y[0] }
    if x[2] != y[2] { return x[2] < y[2] }
    return true
  }
  
  @_transparent
  @_specialize(where T == UInt8)
  @_specialize(where T == UInt32)
  func compareTorsion<T>(_ x: SIMD4<T>, _ y: SIMD4<T>) -> Bool
  where T: FixedWidthInteger {
    if x[1] != y[1] { return x[1] < y[1] }
    if x[2] != y[2] { return x[2] < y[2] }
    if x[0] != y[0] { return x[0] < y[0] }
    if x[3] != y[3] { return x[3] < y[3] }
    return true
  }
  
  @_transparent
  @_specialize(where T == UInt8)
  @_specialize(where T == UInt32)
  func compareRing<T>(_ x: SIMD8<T>, _ y: SIMD8<T>) -> Bool
  where T: FixedWidthInteger {
    if x[0] != y[0] { return x[0] < y[0] }
    if x[1] != y[1] { return x[1] < y[1] }
    if x[2] != y[2] { return x[2] < y[2] }
    if x[3] != y[3] { return x[3] < y[3] }
    if x[4] != y[4] { return x[4] < y[4] }
    return true
  }
}
