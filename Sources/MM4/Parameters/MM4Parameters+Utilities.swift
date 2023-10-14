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
  
  /// Overloads the function `createAtomCodes()`, but is a more elegant API,
  /// that follows the de facto naming convention established for MM4 and
  /// related projects.
  @inline(__always)
  func createAtomCodes<T: SIMD, U: SIMD>(group: T, zero: U) -> U
  where T.Scalar == Int32, U.Scalar == UInt8 {
    var codes = zero
    for lane in 0..<T.scalarCount {
      let atomID = group[lane]
      codes[lane] = atoms.codes[Int(atomID)].rawValue
    }
    return codes
  }
  
  private static let excluding5RingElements: [UInt8] = [11, 19, 25]
  
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
  
  private static let nonCarbonElements: [UInt8] = [8, 11, 15, 19, 25]
  
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
