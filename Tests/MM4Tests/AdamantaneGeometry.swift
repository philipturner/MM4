//
//  AdamantaneGeometry.swift
//
//
//  Created by Philip Turner on 12/13/23.
//

import MM4

// An adamantane-like cage with one sidewall carbon replaced with a sigma bond,
// forming two 5-membered rings.
struct Adamantane {
  var atomicNumbers: [UInt8]
  var bonds: [SIMD2<UInt32>]
  var positions: [SIMD3<Float>]
  
  init(atomCode: MM4AtomCode) {
    var rawAtoms: [SIMD4<Float>]
    var rawBonds: [SIMD2<UInt32>]
    
    switch atomCode {
    case .alkaneCarbon:
      rawAtoms = carbonAtoms
      rawBonds = carbonBonds
    case .silicon:
      rawAtoms = siliconAtoms
      rawBonds = siliconBonds
    default:
      fatalError("Unrecognized atom code for adamantane.")
    }
    
    self.atomicNumbers = rawAtoms.map {
      UInt8(exactly: $0.w)!
    }
    self.bonds = rawBonds
    self.positions = rawAtoms.map {
      SIMD3($0.x, $0.y, $0.z)
    }
  }
}

private let carbonAtoms: [SIMD4<Float>] = [
  SIMD4<Float>(0.000, 0.000, -0.178, 6),
  SIMD4<Float>(0.063, 0.063, -0.241, 1),
  SIMD4<Float>(-0.063, -0.063, -0.241, 1),
  SIMD4<Float>(-0.089, -0.089, 0.089, 6),
  SIMD4<Float>(-0.152, -0.152, 0.152, 1),
  SIMD4<Float>(0.000, -0.178, 0.000, 6),
  SIMD4<Float>(0.063, -0.241, 0.063, 1),
  SIMD4<Float>(-0.063, -0.241, -0.063, 1),
  SIMD4<Float>(0.089, -0.089, -0.089, 6),
  SIMD4<Float>(0.152, -0.152, -0.152, 1),
  SIMD4<Float>(-0.178, 0.000, 0.000, 6),
  SIMD4<Float>(-0.241, -0.063, -0.063, 1),
  SIMD4<Float>(-0.241, 0.063, 0.063, 1),
  SIMD4<Float>(-0.089, 0.089, -0.089, 6),
  SIMD4<Float>(-0.152, 0.152, -0.152, 1),
  SIMD4<Float>(0.000, 0.178, 0.000, 6),
  SIMD4<Float>(0.063, 0.241, -0.063, 1),
  SIMD4<Float>(-0.063, 0.241, 0.063, 1),
  SIMD4<Float>(0.178, 0.000, 0.000, 6),
  SIMD4<Float>(0.241, -0.063, 0.063, 1),
  SIMD4<Float>(0.241, 0.063, -0.063, 1),
  SIMD4<Float>(0.089, 0.089, 0.089, 6),
  SIMD4<Float>(0.152, 0.152, 0.152, 1),
]

private let siliconAtoms: [SIMD4<Float>] = [
  SIMD4<Float>(0.000, 0.000, -0.272, 14),
  SIMD4<Float>(0.085, 0.085, -0.357, 1),
  SIMD4<Float>(-0.085, -0.085, -0.357, 1),
  SIMD4<Float>(0.000, -0.272, 0.000, 14),
  SIMD4<Float>(-0.085, -0.357, -0.085, 1),
  SIMD4<Float>(0.085, -0.357, 0.085, 1),
  SIMD4<Float>(-0.272, 0.000, 0.000, 14),
  SIMD4<Float>(-0.357, 0.085, 0.085, 1),
  SIMD4<Float>(-0.357, -0.085, -0.085, 1),
  SIMD4<Float>(0.136, -0.136, -0.136, 14),
  SIMD4<Float>(0.221, -0.221, -0.221, 1),
  SIMD4<Float>(0.272, 0.000, 0.000, 14),
  SIMD4<Float>(0.357, -0.085, 0.085, 1),
  SIMD4<Float>(0.357, 0.085, -0.085, 1),
  SIMD4<Float>(-0.136, 0.136, -0.136, 14),
  SIMD4<Float>(-0.221, 0.221, -0.221, 1),
  SIMD4<Float>(0.000, 0.272, 0.000, 14),
  SIMD4<Float>(0.085, 0.357, -0.085, 1),
  SIMD4<Float>(-0.085, 0.357, 0.085, 1),
  SIMD4<Float>(-0.136, -0.136, 0.136, 14),
  SIMD4<Float>(-0.221, -0.221, 0.221, 1),
  SIMD4<Float>(0.136, 0.136, 0.136, 14),
  SIMD4<Float>(0.221, 0.221, 0.221, 1),
]

private let carbonBonds: [SIMD2<UInt32>] = [
  SIMD2<UInt32>(0, 1),
  SIMD2<UInt32>(0, 2),
  SIMD2<UInt32>(3, 4),
  SIMD2<UInt32>(3, 5),
  SIMD2<UInt32>(5, 6),
  SIMD2<UInt32>(5, 7),
  SIMD2<UInt32>(0, 8),
  SIMD2<UInt32>(5, 8),
  SIMD2<UInt32>(8, 9),
  SIMD2<UInt32>(3, 10),
  SIMD2<UInt32>(10, 11),
  SIMD2<UInt32>(10, 12),
  SIMD2<UInt32>(0, 13),
  SIMD2<UInt32>(10, 13),
  SIMD2<UInt32>(13, 14),
  SIMD2<UInt32>(13, 15),
  SIMD2<UInt32>(15, 16),
  SIMD2<UInt32>(15, 17),
  SIMD2<UInt32>(8, 18),
  SIMD2<UInt32>(18, 19),
  SIMD2<UInt32>(18, 20),
  SIMD2<UInt32>(15, 21),
  SIMD2<UInt32>(18, 21),
  SIMD2<UInt32>(21, 22),
  SIMD2<UInt32>(3, 21),
]

private let siliconBonds: [SIMD2<UInt32>] = [
  SIMD2<UInt32>(0, 1),
  SIMD2<UInt32>(0, 2),
  SIMD2<UInt32>(3, 4),
  SIMD2<UInt32>(3, 5),
  SIMD2<UInt32>(6, 7),
  SIMD2<UInt32>(6, 8),
  SIMD2<UInt32>(0, 9),
  SIMD2<UInt32>(3, 9),
  SIMD2<UInt32>(9, 10),
  SIMD2<UInt32>(9, 11),
  SIMD2<UInt32>(11, 12),
  SIMD2<UInt32>(11, 13),
  SIMD2<UInt32>(0, 14),
  SIMD2<UInt32>(6, 14),
  SIMD2<UInt32>(14, 15),
  SIMD2<UInt32>(14, 16),
  SIMD2<UInt32>(16, 17),
  SIMD2<UInt32>(16, 18),
  SIMD2<UInt32>(3, 19),
  SIMD2<UInt32>(6, 19),
  SIMD2<UInt32>(19, 20),
  SIMD2<UInt32>(11, 21),
  SIMD2<UInt32>(16, 21),
  SIMD2<UInt32>(21, 22),
  SIMD2<UInt32>(19, 21),
]

// TODO: - Copy parameters from the image, ensure they are the same.
