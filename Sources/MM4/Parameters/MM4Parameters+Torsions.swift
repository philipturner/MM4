//
//  MM4Parameters+Torsions.swift
//  MM4
//
//  Created by Philip Turner on 10/7/23.
//

/// Parameters for a group of 4 atoms.
public struct MM4Torsions {
  /// Each value corresponds to the torsion at the same array index.
  public var extendedParameters: [MM4TorsionExtendedParameters?] = []
  
  /// Groups of atom indices that form a torsion.
  public var indices: [SIMD4<UInt32>] = []
  
  /// Map from a group of atoms to a torsion index.
  public var map: [SIMD4<UInt32>: UInt32] = [:]
  
  /// Each value corresponds to the torsion at the same array index.
  public var parameters: [MM4TorsionParameters] = []
  
  /// The smallest ring each torsion is involved in.
  public var ringTypes: [UInt8] = []
  
  mutating func append(contentsOf other: Self, atomOffset: UInt32) {
    let torsionOffset = UInt32(self.indices.count)
    self.extendedParameters += other.extendedParameters
    self.indices += other.indices.map {
      $0 &+ atomOffset
    }
    for key in other.map.keys {
      let value = other.map[key].unsafelyUnwrapped
      self.map[key &+ atomOffset] = value &+ torsionOffset
    }
    self.parameters += other.parameters
    self.ringTypes += other.ringTypes
  }
}

/// Parameters for a torsion between three nonpolar bonds.
///
/// V1 term:
/// - zeroed out for X-C-C-H
/// - present for C-C-C-C
/// - present for 5-membered rings
/// - present for X-C-C-F
///
/// V3 term:
/// - present for X-C-C-H
/// - present for C-C-C-C
/// - present for 5-membered rings
/// - present for X-C-C-F
///
/// Vn term:
/// - 6 for some cases of X-C-C-H
/// - 2 for some cases of C-C-C-C
/// - zeroed out for 5-membered rings
/// - 2 for X-C-C-F
///
/// 1-term torsion-stretch:
/// - present for X-C-C-H
/// - present for C-C-C-C
/// - present for 5-membered rings
/// - zeroed out for X-C-C-F, to prioritize conciseness over performance
public struct MM4TorsionParameters {
  /// Units: kilocalorie / mole
  ///
  /// > WARNING: Convert kcal/mol to kJ/mol.
  public var V1: Float
  
  /// Units: kilocalorie / mole
  ///
  /// > WARNING: Convert kcal/mol to kJ/mol.
  public var Vn: Float
  
  /// Units: kilocalorie / mole
  ///
  /// > WARNING: Convert kcal/mol to kJ/mol.
  public var V3: Float
  
  /// The factor to multiply the angle with, inside the cosine term for Vn.
  ///
  /// The value of `n` is most often 2. It must be an even integer.
  public var n: Float
  
  /// Units: kilocalorie / angstrom \* mole
  ///
  /// > WARNING: Convert kcal/mol to kJ/mol.
  public var Kts3: Float
}

/// Parameters for torsion forces unique to highly electronegative elements.
///
/// The parameters include V4, V6, 3-term torsion-stretch, and torsion-bend. In
/// addition, parameters for the bend-torsion-bend force. This force is omitted
/// from C-H torsions for efficiency.
public struct MM4TorsionExtendedParameters {
  /// Units: kilocalorie / mole
  ///
  /// > WARNING: Convert kcal/mol to kJ/mol.
  public var V4: Float
  
  /// Units: kilocalorie / mole
  ///
  /// > WARNING: Convert kcal/mol to kJ/mol.
  public var V6: Float
  
  /// The V1-like term contributing to torsion-stretch stiffness.
  ///
  /// > WARNING: Convert kcal/mol to kJ/mol.
  public var Kts1: (left: Float, central: Float, right: Float)
  
  /// The V2-like term contributing to torsion-stretch stiffness.
  ///
  /// > WARNING: Convert kcal/mol to kJ/mol.
  public var Kts2: (left: Float, central: Float, right: Float)
  
  /// The V3-like term contributing to torsion-stretch stiffness.
  ///
  /// > WARNING: Convert kcal/mol to kJ/mol.
  public var Kts3: (left: Float, central: Float, right: Float)
  
  /// The V1-like term contributing to torsion-bend stiffness.
  ///
  /// > WARNING: Convert radians to degrees.
  ///
  /// > WARNING: Convert kcal/mol to kJ/mol.
  public var Ktb1: (left: Float, right: Float)
  
  /// The V2-like term contributing to torsion-bend stiffness.
  ///
  /// > WARNING: Convert radians to degrees.
  ///
  /// > WARNING: Convert kcal/mol to kJ/mol.
  public var Ktb2: (left: Float, right: Float)
  
  /// The V3-like term contributing to torsion-bend stiffness.
  ///
  /// > WARNING: Convert radians to degrees.
  ///
  /// > WARNING: Convert kcal/mol to kJ/mol.
  public var Ktb3: (left: Float, right: Float)
  
  /// Bend-torsion-bend constant.
  ///
  /// > WARNING: Convert radians to degrees.
  ///
  /// > WARNING: Convert kcal/mol to kJ/mol.
  public var Kbtb: Float
}

extension MM4Parameters {
  /// - throws: `.missingParameter`
  mutating func createTorsionParameters(forces: MM4ForceOptions) throws {
    for torsionID in torsions.indices.indices {
      let torsion = torsions.indices[torsionID]
      let ringType = torsions.ringTypes[torsionID]
      let codes = with5RingsRemoved {
        createAtomCodes(group: torsion, zero: SIMD4<UInt8>.zero)
      }
      
      func createTorsionError() -> MM4Error {
        let map = SIMD4<Int32>(truncatingIfNeeded: torsion)
        let addresses = createAddresses(map)
        return MM4Error.missingParameter(addresses)
      }
      if containsTwoNonCarbons(codes) {
        throw createTorsionError()
      }
      
      /// "Note that unless they are explicitly in the table, parameters
      /// involving five-membered ring atom types (122 and 123) are assigned
      /// parameters which involve the regular atom types (types 2 and 1). This
      /// is true for all parameters."
      ///
      /// - Parameter closure: The function to try once with 5-ring carbons,
      ///   then again with all 5-ring carbons set to 6-ring. Returns whether
      ///   the first attempt succeeded.
      /// - Returns: Whether the sorted codes are different from the original
      ///   ones.
      @_transparent
      @discardableResult
      func with5RingAttempt(_ closure: (SIMD4<UInt8>) -> Bool) throws -> Bool {
        var sortedCodes = sortTorsion(codes)
        if closure(sortedCodes) {
          return any(sortedCodes .!= codes)
        }
        
        let newCodes = codes.replacing(with: .one, where: codes .== 123)
        sortedCodes = sortTorsion(newCodes)
        if closure(sortedCodes) {
          return any(sortedCodes .!= newCodes)
        } else {
          throw createTorsionError()
        }
      }
      
      // MARK: - Torsion, Bend-Torsion-Bend
      
      var V1: Float = 0.000
      var Vn: Float = 0.000
      var V3: Float = 0.000
      var n: Float = 2
      
      var V4: Float?
      var V6: Float?
      var Kbtb: Float?
      
      try with5RingAttempt { codes in
        switch (codes[0], codes[1], codes[2], codes[3]) {
          // Carbon
        case (1, 1, 1, 1):   (V1, Vn, V3) = (0.239, 0.024, 0.637)
        case (1, 1, 1, 5):             V3 = 0.290
        case (5, 1, 1, 5):    (V3, Vn, n) = (0.260, 0.008, 6)
        case (1, 123, 123, 1):   (V1, V3) = (0.160, 0.550)
        case (1, 123, 123, 123): (V1, V3) = (0.160, 0.550)
        case (5, 123, 123, 5):         V3 = 0.300
        case (5, 123, 123, 123):       V3 = 0.290
        case (5, 1, 123, 123):         V3 = 0.306
        case (1, 123, 123, 5):         V3 = 0.306
        case (5, 1, 123, 5):           V3 = 0.260
        case (123, 123, 123, 123):
          if ringType == 5 {     (V1, V3) = (-0.150, 0.160) }
          else {                 (V1, V3) = (-0.120, 0.550) }
          
          // Nitrogen
        case (1, 1, 1, 8):     (V1, Vn, V3) = (1.139, 1.348, 1.582)
          /**/               (V4, V6, Kbtb) = (-0.140, 0.172, -0.05)
        case (5, 1, 1, 8):       (V3, Kbtb) = (0.455, -0.05)
        case (8, 1, 1, 8):     (V1, Vn, V3) = (2.545, -2.520, 3.033)
        case (1, 1, 8, 1):     (V1, Vn, V3) = (1.193, -0.337, 0.870)
          /**/                     (V4, V6) = (0.228, -0.028)
        case (5, 1, 8, 1):     (V1, Vn, V3) = (0.072, -0.512, 0.562)
          /**/                         Kbtb = 0.05
        case (5, 1, 8, 123): (Vn, V3, Kbtb) = (-0.450, 0.170, -0.12)
        case (1, 8, 123, 5):       (Vn, V3) = (0.550, 0.100)
        case (1, 8, 123, 123):     (Vn, V3) = (-0.520, 0.180)
        case (123, 8, 123, 5): (V1, Vn, V3) = (0.072, -0.012, 0.597)
        case (5, 123, 123, 8):   (V3, Kbtb) = (0.374, -0.09)
        case (123, 8, 123, 123):
          if ringType == 5 {   (V1, Vn, V3) = (1.150, -0.040, 0.860) }
          else { return false }
        case (8, 123, 123, 123):
          if ringType == 5 {             V3 = 0.699 }
          else {                   (Vn, V3) = (-0.850, 0.200) }
          
          // Oxygen
        case (1, 1, 1, 6):       (V1, Vn, V3) = (-0.333, 0.037, 0.552)
        case (5, 1, 1, 6):       (V1, Vn, V3) = (-0.593, 0.554, 0.474)
          /**/                     (V6, Kbtb) = (0.070, -0.100)
        case (6, 1, 1, 6):       (V1, Vn, V3) = (-0.917, -0.631, 0.641)
          /**/                             V6 = -0.100
        case (1, 1, 6, 1):       (V1, Vn, V3) = (1.900, -0.500, 1.250)
        case (5, 1, 6, 1):     (V3, V6, Kbtb) = (0.730, 0.028, -0.050)
        case (6, 1, 6, 1):       (V1, Vn, V3) = (-0.350, -0.900, -0.020)
          /**/                     (V4, Kbtb) = (0.503, -0.150)
        case (123, 6, 123, 5): return false
        case (123, 6, 123, 6):
          if ringType == 5 {     (V1, Vn, V3) = (0.350, -1.900, 1.550) }
          else {                 (V1, Vn, V3) = (-0.350, -0.900, -1.550)
            /**/                         Kbtb = 0.150 }
        case (5, 123, 123, 6): (V1, V3, Kbtb) = (0.200, 0.480, -0.100)
        case (6, 123, 123, 6):
          if ringType == 5 {     (V1, Vn, V3) = (0.846, -2.314, 1.653)
            /**/                     (V4, V6) = (-0.101, 0.350) }
          else {                 (V1, Vn, V3) = (-0.217, -0.851, 0.441)
            /**/               (V4, V6, Kbtb) = (-0.101, 0.350, 0.080) }
        case (6, 123, 123, 123):
          if ringType == 5 {         (Vn, V3) = (0.500, 0.930) }
          else {                 (V1, Vn, V3) = (-0.128, -0.645, 0.934) }
        case (123, 6, 123, 123):
          if ringType == 5 {         (V1, V3) = (-0.900, 0.460) }
          else { return false }
          
          // Fluorine
        case (1, 1, 1, 11):  (V1, Vn, V3) = (-0.360, 0.380, 0.978)
          /**/             (V4, V6, Kbtb) = (0.240, 0.010, -0.06)
        case (5, 1, 1, 11):  (V1, Vn, V3) = (-0.460, 1.190, 0.420)
          /**/                       Kbtb = 0.06
        case (11, 1, 1, 11): (V1, Vn, V3) = (-1.350, 0.305, 0.355)
          /**/                       Kbtb = -0.06
          
          // Silicon
        case (5, 1, 1, 19):   V3 = 0.200
        case (19, 1, 1, 19):  V3 = 0.167
        case (1, 1, 19, 5):   V3 = 0.295
        case ( 5, 1, 19, 1):  V3 = 0.195
        case ( 5, 1, 19, 5):  V3 = 0.177
        case (19, 1, 19, 1):  V3 = 0.100
        case (19, 1, 19, 5):  V3 = 0.167
        case ( 1, 1, 19, 19): V3 = 0.300
        case ( 5, 1, 19, 19): V3 = 0.270
        case (1, 19, 19, 5):  V3 = 0.127
        case (1, 19, 19, 1):  V3 = 0.107
        case (1, 19, 19, 19): V3 = 0.350
        case (5, 19, 19, 5):  V3 = 0.132
        case (5, 19, 19, 19): V3 = 0.070
        case (1, 1, 1, 19):
          if ringType == 5 {  V3 = 0.850 }
          else {        (Vn, V3) = (0.050, 0.240) }
        case (1, 1, 19, 1):
          if ringType == 5 {  Vn = 0.800 }
          else {              V3 = 0.167 }
        case (19, 19, 19, 19):
          if ringType == 5 {  V3 = 0.175 }
          else {              V3 = 0.125 }
          
          // Phosphorus
        case (5, 1, 25, 1):     (V3, V6) = (0.300, -0.050)
        case (5, 1, 1, 25):     (Vn, V3) = (0.200, 0.360)
        case (1, 1, 25, 1): (V1, Vn, V3) = (0.800, -0.400, 0.100)
          
          // Sulfur
          //
          // Grabbing the 15-1-15-1 torsion parameters from MM3.
        case (15, 1, 15, 1):      (Vn, V3) = (-0.900, 0.300)
        case (5, 1, 15, 1):     (V3, Kbtb) = (0.540, 0.080)
        case (1, 1, 15, 1): (V1, V3, Kbtb) = (0.410, 0.600, 0.004)
        case (15, 1, 1, 15):  (V1, Vn, V3) = (0.461, 0.144, 1.511)
          /**/                        Kbtb = 0.130
        case (5, 1, 1, 15):     (V3, Kbtb) = (0.460, 0.050)
        case (1, 1, 1, 15):   (V1, Vn, V3) = (0.420, 0.100, 0.200)
          /**/                        Kbtb = 0.090
        case (5, 123, 123, 15): (V3, Kbtb) = (0.200, 0.020)
        case (123, 15, 123, 1):   (V1, V3) = (0.100, 0.200)
        case (5, 1, 123, 15):     (Vn, V3) = (0.330, 0.200)
          /**/                        Kbtb = 0.050
        case (123, 15, 123, 5): (V3, Kbtb) = (0.450, 0.020)
        case (15, 123, 123, 123):
          if ringType == 5 {
            (V1, Vn, V3, Kbtb) = (0.040, 0.200, 0.300, 0.100)
          } else {
            (V1, Vn, V3, Kbtb) = (0.520, 0.080, 0.250, 0.100)
          }
        case (123, 15, 123, 123):
          if ringType == 5 {  (V1, Vn, V3) = (0.440, 0.300, 0.500) }
          else { return false }
          
          // Germanium
          //
          // There are two values for 1-1-1-31 in the MM3 paper. It hints that
          // one is from the preliminary MM3(1996), but doesn't explicitly
          // state which. The Tinker implementation suggests the one without
          // the "b" footnote.
          //
          // There are no germanium parameters for the following torsions in
          // MM3. As with angles, it looks like silicon has the same ballpark
          // value for torsions - neither consistently greater or smaller.
          // Note that crystolecules are supposedly very insensitive to
          // torsions. I wouldn't call these gold standard, just borderline
          // okay for trying out germanium in diamondoid nanomachines.
          // - 31-1-1-31
          // - 31-1-31-1
          // - 31-1-31-5
          // - 1-1-31-31
          // - 5-1-31-31
          // - 1-31-31-5
          // - 1-31-31-31
          // - 5-31-31-31
        case (5, 1, 1, 31):        V3 = 0.185
        case (31, 1, 1, 31):       V3 = 0.167
        case (1, 1, 31, 5):        V3 = 0.172
        case (5, 1, 31, 1):        V3 = 0.127
        case (5, 1, 31, 5):        V3 = 0.132
        case (31, 1, 31, 1):       V3 = 0.100
        case (31, 1, 31, 5):       V3 = 0.167
        case (1, 1, 31, 31):       V3 = 0.300
        case (5, 1, 31, 31):       V3 = 0.270
        case (1, 31, 31, 5):       V3 = 0.127
        case (1, 31, 31, 1):       V3 = 0.112
        case (1, 31, 31, 31):      V3 = 0.350
        case (5, 31, 31, 5):       V3 = 0.165
        case (5, 31, 31, 31):      V3 = 0.070
        case (31, 31, 31, 31):     V3 = 0.112
        case (1, 1, 1, 31):
          if ringType == 5 {       V3 = 0.520 }
          else {             (V1, V3) = (-0.200, 0.112) }
        case (1, 1, 31, 1):
          if ringType == 5 { (V1, V3) = (-0.200, 0.100) }
          else {         (V1, Vn, V3) = (-0.200, 0.085, 0.112) }
          
        default:
          return false
        }
        return true
      }
      
      if !forces.contains(.torsionBend) {
        Kbtb = nil
      }
      
      // MARK: - Torsion-Bend
      
      var Ktb_l: SIMD3<Float>?
      var Ktb_r: SIMD3<Float>?
      
      if forces.contains(.torsionBend) {
        let swappedCodesForTorsionBend = try with5RingAttempt { codes in
          switch (codes[0], codes[1], codes[2], codes[3]) {
            // Carbon
          case (1, 1, 1, 1): fallthrough
          case (1, 1, 1, 5): fallthrough
          case (5, 1, 1, 5): fallthrough
          case (1, 123, 123, 1): fallthrough
          case (1, 123, 123, 123): fallthrough
          case (5, 123, 123, 5): fallthrough
          case (5, 123, 123, 123): fallthrough
          case (5, 1, 123, 123): fallthrough
          case (1, 123, 123, 5): fallthrough
          case (5, 1, 123, 5): fallthrough
          case (123, 123, 123, 123): break
            
            // Nitrogen
          case (1, 1, 1, 8):       Ktb_l = SIMD3(0.007, 0.000, -0.005)
            /**/                   Ktb_r = SIMD3(0.010, -0.004, -0.006)
          case (5, 1, 1, 8):       Ktb_r = SIMD3(-0.003, -0.003, 0.000)
          case (8, 1, 1, 8):       Ktb_l = SIMD3(-0.008, 0.008, -0.018)
            /**/                   Ktb_r = SIMD3(-0.008, 0.008, -0.018)
          case (1, 1, 8, 1):       Ktb_l = SIMD3(-0.028, -0.020, 0.007)
            /**/                   Ktb_r = SIMD3(0.003, -0.004, -0.009)
          case (5, 1, 8, 1):       Ktb_l = SIMD3(-0.010, -0.025, 0.000)
            /**/                   Ktb_r = SIMD3(-0.003, -0.003, 0.000)
          case (5, 1, 8, 123):     Ktb_l = SIMD3(-0.003, -0.003, 0.000)
          case (1, 8, 123, 5):     Ktb_r = SIMD3(-0.015, -0.028, 0.000)
          case (1, 8, 123, 123):   Ktb_r = SIMD3(-0.003, -0.003, 0.000)
          case (123, 8, 123, 5):   Ktb_r = SIMD3(-0.015, -0.028, 0.000)
          case (5, 123, 123, 8):   Ktb_r = SIMD3(-0.003, -0.003, 0.000)
          case (123, 8, 123, 123):
            if ringType == 5 {     Ktb_r = SIMD3(-0.008, -0.004, 0.000) }
            else { return false }
          case (8, 123, 123, 123): Ktb_l = SIMD3(-0.008, -0.004, 0.000)
            
            // Oxygen
          case (1, 1, 1, 6):       Ktb_l = SIMD3(-0.003, -0.003, 0.000)
            /**/                   Ktb_r = SIMD3(0.005, -0.002, 0.001)
          case (5, 1, 1, 6):       Ktb_l = SIMD3(0.006, -0.006, 0.000)
            /**/                   Ktb_r = SIMD3(0.008, -0.010, 0.000)
          case (6, 1, 1, 6):       Ktb_l = SIMD3(-0.004, -0.010, 0.000)
            /**/                   Ktb_r = SIMD3(-0.004, -0.010, 0.000)
          case (1, 1, 6, 1):       Ktb_l = SIMD3(-0.043, -0.013, -0.002)
            /**/                   Ktb_r = SIMD3(-0.015, 0.017, -0.014)
          case (5, 1, 6, 1):       Ktb_l = SIMD3(-0.018, -0.008, 0.000)
            /**/                   Ktb_r = SIMD3(-0.005, 0.006, 0.000)
          case (6, 1, 6, 1):       Ktb_l = SIMD3(0.030, -0.040, 0.009)
            /**/                   Ktb_r = SIMD3(-0.003, -0.001, 0.000)
          case (123, 6, 123, 5):   Ktb_l = SIMD3(-0.005, 0.006, 0.000)
            /**/                   Ktb_r = SIMD3(-0.018, -0.008, 0.000)
          case (123, 6, 123, 6):
            if ringType == 5 {     Ktb_l = SIMD3(0.000, 0.000, 0.001) }
            else {                 Ktb_l = SIMD3(0.030, -0.001, 0.000)
              /**/                 Ktb_r = SIMD3(0.030, -0.040, 0.009) }
          case (5, 123, 123, 6):   Ktb_l = SIMD3(0.008, -0.010, 0.000)
            /**/                   Ktb_r = SIMD3(0.006, -0.006, 0.000)
          case (6, 123, 123, 6):   Ktb_l = SIMD3(-0.005, -0.014, 0.000)
            /**/                 Ktb_r = SIMD3(-0.005, -0.014, 0.000)
          case (6, 123, 123, 123): Ktb_l = SIMD3(0.005, -0.002, 0.001)
            /**/                 Ktb_r = SIMD3(-0.003, -0.003, -0.000)
          case (123, 6, 123, 123):
            if ringType == 5 {     Ktb_l = SIMD3(-0.010, 0.017, 0.000)
              /**/                 Ktb_r = SIMD3(-0.043, -0.013, -0.002) }
            else { return false }
            
            // Fluorine
          case (1, 1, 1, 11):  Ktb_l = SIMD3(0.000, -0.012, -0.009)
            /**/               Ktb_r = SIMD3(0.005, 0.004, 0.003)
          case (5, 1, 1, 11):  Ktb_l = SIMD3(0.002, -0.022, 0.000)
            /**/               Ktb_r = SIMD3(0.000, 0.000, -0.001)
          case (11, 1, 1, 11): Ktb_l = SIMD3(0.000, -0.015, -0.003)
            /**/               Ktb_r = SIMD3(0.000, -0.015, -0.003)
            
            // Silicon
          case (19, 1, 1, 19): fallthrough
          case ( 5, 1, 19, 1): fallthrough
          case ( 5, 1, 19, 5): fallthrough
          case (19, 1, 19, 1): fallthrough
          case (19, 1, 19, 5): fallthrough
          case ( 1, 1, 19, 19): fallthrough
          case ( 5, 1, 19, 19): fallthrough
          case (1, 19, 19, 5): fallthrough
          case (1, 19, 19, 1): fallthrough
          case (1, 19, 19, 19): fallthrough
          case (5, 19, 19, 5): fallthrough
          case (5, 19, 19, 19): fallthrough
          case (1, 1, 1, 19): fallthrough
          case (1, 1, 19, 1): fallthrough
          case (19, 19, 19, 19): break
            
            // Phosphorus
            //
            // WARNING: Entering phosphorus in the reverse order that it appears
            // in the research paper. This may cause mistakes from confusing two
            // quantities. Place Ktb_r before Ktb_l in this section. Not all
            // cases have reversed order, but pay careful attention to the ones
            // that do.
          case (5, 1, 25, 1): Ktb_r = SIMD3(0.000, 0.000, -0.001)
            /**/              Ktb_l = SIMD3(-0.001, -0.004, -0.001)
          case (5, 1, 1, 25): Ktb_r = SIMD3(0.005, 0.000, 0.000)
          case (1, 1, 25, 1): Ktb_r = SIMD3(-0.002, -0.003, -0.001)
            /**/              Ktb_l = SIMD3(-0.015, 0.007, -0.006)
            
            // Sulfur
            //
            // The 15-1-1-15 torsion is symmetric, so setting only one of the
            // angles to 0.003 seems suspicious. I will set both angles.
          case (15, 1, 15, 1): break
          case (5, 1, 15, 1):  Ktb_l = SIMD3(0.006, 0.020, 0.000)
          case (1, 1, 15, 1):  Ktb_l = SIMD3(0.030, 0.023, 0.000)
          case (15, 1, 1, 15): Ktb_l = SIMD3(0.000, 0.003, 0.000)
            /**/               Ktb_r = SIMD3(0.000, 0.003, 0.000)
          case (5, 1, 1, 15): fallthrough
          case (1, 1, 1, 15): fallthrough
          case (5, 123, 123, 15): fallthrough
          case (123, 15, 123, 1): fallthrough
          case (5, 1, 123, 15): fallthrough
          case (123, 15, 123, 5): fallthrough
          case (15, 123, 123, 123): fallthrough
          case (123, 15, 123, 123): break
            
          default:
            return false
          }
          return true
        }
        
        if swappedCodesForTorsionBend {
          swap(&Ktb_l, &Ktb_r)
        }
      }
      
      // MARK: - Torsion-Stretch
      
      // The formula from the MM4 alkene paper was ambiguous, specifying "-k":
      //   -k * Δl * Kts * (1 + cos(3ω))
      // The formula from the MM3 original paper was:
      //   11.995 * (Kts/2) * Δl * (1 + cos(3ω))
      // After running several parameters through Desmos, and comparing similar
      // ones (https://www.desmos.com/calculator/p5wqbw7tku), I think I have
      // identified a typo. "-k" was supposed to be "K_s^-1" or "1/K_s". This
      // would result in:
      // - C-C having ~33% less TS stiffness in MM4 than in MM3
      // - C-Csp2 (MM4) having ~17% less stiffness than Si-Csp2 (MM3)
      // - C-S (MM4) having nearly identical stiffness to C-Si (MM3)
      // - Central TS for H-C-C-F (MM4) having 42% less peak stiffness than C-S
      // - Central TS for C-C-C-F (MM4) having 2x larger peaks than C-S, and a
      //   new trough with ~2x the magnitude of the C-S peak
      // - Central TS for F-C-C-F (MM4) having 1% less peak stiffness than C-S,
      //   but a new trough with ~3.7x the magnitude of the C-S peak
      //
      // New formula:
      //   Δl * (Kts / Ks) * (1 + cos(3ω))
      
      // Update: During an experiment with lonsdaleite, the TS caused very
      // incorrect behavior with the old MM4 implementation (~10 too much
      // stiffness). It didn't show up with cubic diamond, where the torsion
      // angle was 180 degrees.
      //
      // While double-checking torsion stretch, I found a simpler explanation.
      //
      // Bond         MM3     MM4   MM4/MM3 Ratio
      // 1-1          0.059   0.66  11.18644068
      // 123-123(5)   0.059   0.84  14.23728814
      // 56-56(4)     0.059   0.94  15.93220339
      // 1-2          0.27    2.159 7.996296296
      // 1-6          0.1     1.559 15.59
      // 1-8          0.061   0.58  9.508196721
      // 1-15         0.17    1.559 9.170588235
      // 1-25         0.104   1.247 11.99038462
      //
      // The average of all eight (MM4/MM3 Ratio) values is 11.951. This is very
      // close to 11.995, and not a coincidence. The actual formula is shown
      // below. However, I will leave the previous (incorrect) reasoning there
      // for reference. "-k" could have been "0.5". The "-" key on the keyboard
      // is close to "0". lowercase "k" is also sort-of close to "." and also
      // not activated by "shift".
      //
      // Furthermore, those two keys are both offset by ~1 key from the correct
      // one. The third key could have simply been tapped out-of-bounds of the
      // keyboard. Another extrapolation: it activated a distracting Fn key that
      // interrupted the author. This suggests a very plausible environment that
      // would facilitate the hypothesized mistake.
      //
      // (Kts/2) * Δl * (1 + cos(3ω))
      
      var Kts_l: SIMD3<Float>?
      var Kts_c: SIMD3<Float>?
      var Kts_r: SIMD3<Float>?
      
      if forces.contains(.torsionStretch) {
        let swappedCodesForTorsionStretch = try with5RingAttempt { codes in
          let middle = SIMD2(codes[1], codes[2])
          var middle6Ring = middle.replacing(with: .one, where: middle .== 123)
          middle6Ring = sortBond(middle6Ring)
          
          if all(middle .== 123) {
            Kts_c = SIMD3(0.000, 0.000, 0.840)
          } else {
            // For silicon, multiply the MM3 torsion-stretch constant by 11.995.
            switch (middle6Ring[0], middle6Ring[1]) {
            case (1, 1): Kts_c = SIMD3(0.000, 0.000, 0.640)
            case (1, 19): Kts_c = SIMD3(0.000, 0.000, 0.036 * 11.995)
            case (19, 19): Kts_c = SIMD3(0.000, 0.000, 0.012 * 11.995)
            case (1, 15): Kts_c = SIMD3(0.000, 0.000, 1.559)
            default: break
            }
          }
          
          switch (codes[0], codes[1], codes[2], codes[3]) {
            // Carbon
          case (1, 1, 1, 1): fallthrough
          case (1, 1, 1, 5): fallthrough
          case (5, 1, 1, 5): fallthrough
          case (1, 123, 123, 1): fallthrough
          case (1, 123, 123, 123): fallthrough
          case (5, 123, 123, 5): fallthrough
          case (5, 123, 123, 123): fallthrough
          case (5, 1, 123, 123): fallthrough
          case (1, 123, 123, 5): fallthrough
          case (5, 1, 123, 5): fallthrough
          case (123, 123, 123, 123): break
            
            // Deviating from the convention of making everything tabular, for
            // nitrogen, oxygen, and fluorine. This makes it easier to read
            // and/or less effort to format than tabular form.
            
            // Nitrogen
          case (1, 1, 1, 8):
            Kts_l = SIMD3(0.000, 0.000, 1.300)
            Kts_c = SIMD3(3.000, -3.100, 2.860)
            Kts_r = SIMD3(-0.600, 1.000, 1.500)
          case (5, 1, 1, 8): break
          case (8, 1, 1, 8):
            Kts_l = SIMD3(0.000, 5.000, 0.000)
            Kts_c = SIMD3(4.000, -5.560, 2.660)
            Kts_r = SIMD3(0.000, 5.000, 0.000)
          case (1, 1, 8, 1):
            Kts_l = SIMD3(5.000, 5.000, 0.000)
            Kts_c = SIMD3(1.000, 0.000, 4.580)
          case (5, 1, 8, 1):
            Kts_l = SIMD3(3.100, 8.200, 0.000)
            Kts_c = SIMD3(0.000, -1.100, 0.580)
          case (5, 1, 8, 123): return false
          case (1, 8, 123, 5): return false
          case (1, 8, 123, 123): return false
          case (123, 8, 123, 5):
            Kts_c = SIMD3(0.000, 0.000, 0.580)
            Kts_r = SIMD3(0.000, 8.038, 0.000)
          case (5, 123, 123, 8): return false
          case (123, 8, 123, 123): return false
          case (8, 123, 123, 123):
            Kts_l = SIMD3(1.250, 1.200, 0.000)
            Kts_c = SIMD3(4.000, -3.000, 2.660)
            Kts_r = SIMD3(1.500, 0.000, 2.500)
            
            // Oxygen
          case (1, 1, 1, 6):
            Kts_l = SIMD3(1.000, -1.000, 4.500)
            Kts_c = SIMD3(0.000, 0.000, 0.660)
            Kts_r = SIMD3(-0.500, 2.000, 0.700)
          case (5, 1, 1, 6):
            Kts_c = SIMD3(0.000, 0.000, 0.660)
          case (6, 1, 1, 6):
            Kts_l = SIMD3(-5.000, 6.000, 0.000)
            Kts_c = SIMD3(0.000, -7.000, 3.800)
            Kts_r = SIMD3(-5.000, 6.000, 0.000)
          case (1, 1, 6, 1):
            Kts_l = SIMD3(-1.500, 3.500, 4.500)
            Kts_c = SIMD3(0.000, 0.000, 1.559)
            Kts_r = SIMD3(0.000, 2.000, 0.000)
          case (5, 1, 6, 1):
            Kts_l = SIMD3(5.997, 7.798, 0.000)
            Kts_c = SIMD3(0.000, 0.000, 1.559)
          case (6, 1, 6, 1):
            Kts_l = SIMD3(5.300, 9.650, 0.000)
            Kts_c = SIMD3(-7.000, -6.520, 1.559)
            Kts_r = SIMD3(0.000, 3.900, 0.000)
          case (123, 6, 123, 6):
            Kts_l = SIMD3(0.000, 14.000, 0.000)
            Kts_c = SIMD3(0.000, -6.920, -2.200)
            Kts_r = SIMD3(5.300, 9.650, 0.000)
          case (123, 6, 123, 123):
            Kts_l = SIMD3(0.000, 9.000, -5.000)
            Kts_c = SIMD3(0.000, 0.000, 1.559)
          case (6, 123, 123, 6):
            Kts_l = SIMD3(5.000, 0.000, 0.000)
            Kts_c = SIMD3(12.000, -7.000, 3.800)
            Kts_r = SIMD3(5.000, 0.000, 0.000)
          case (6, 123, 123, 123):
            Kts_l = SIMD3(-0.500, 2.000, 0.700)
            Kts_c = SIMD3(0.000, 0.000, 5.000)
            Kts_r = SIMD3(1.000, -1.000, 4.500)
            
            // Fluorine
          case (1, 1, 1, 11):
            Kts_c = SIMD3(5.300, -4.800, 4.500)
          case (5, 1, 1, 11):
            Kts_c = SIMD3(0.000, -1.650, 2.400)
            Kts_r = SIMD3(0.000, 0.000, 0.550)
          case (11, 1, 1, 11):
            Kts_l = SIMD3(-5.550, 4.500, 0.000)
            Kts_c = SIMD3(4.800, -5.000, 1.559)
            Kts_r = SIMD3(-5.550, 4.500, 0.000)
            
            // Silicon
          case (19, 1, 1, 19): fallthrough
          case ( 5, 1, 19, 1): fallthrough
          case ( 5, 1, 19, 5): fallthrough
          case (19, 1, 19, 1): fallthrough
          case (19, 1, 19, 5): fallthrough
          case ( 1, 1, 19, 19): fallthrough
          case ( 5, 1, 19, 19): fallthrough
          case (1, 19, 19, 5): fallthrough
          case (1, 19, 19, 1): fallthrough
          case (1, 19, 19, 19): fallthrough
          case (5, 19, 19, 5): fallthrough
          case (5, 19, 19, 19): fallthrough
          case (1, 1, 1, 19): fallthrough
          case (1, 1, 19, 1): fallthrough
          case (19, 19, 19, 19): break
            
            // Phosphorus
          case (5, 1, 25, 1): fallthrough
          case (5, 1, 1, 25): fallthrough
          case (1, 1, 25, 1): Kts_c = SIMD3(0.000, 0.000, 0.000)
            
            // Sulfur
          case (15, 1, 15, 1): break
          case (5, 1, 15, 1):  Kts_l = SIMD3(0.000, 0.900, 0.000)
          case (1, 1, 15, 1):  Kts_l = SIMD3(0.000, 1.439, 0.000)
          case (15, 1, 1, 15): Kts_l = SIMD3(1.919, 1.919, 0.000)
            /**/               Kts_r = SIMD3(1.919, 1.919, 0.000)
          case (5, 1, 1, 15):  break
          case (1, 1, 1, 15):  Kts_l = SIMD3(4.798, 4.798, 0.000)
          case (5, 123, 123, 15): return false
          case (123, 15, 123, 1): return false
          case (5, 1, 123, 15): return false
          case (123, 15, 123, 5): return false
          case (15, 123, 123, 123): return false
          case (123, 15, 123, 123): return false
          
          default:
            return false
          }
          return true
        }
        
        if swappedCodesForTorsionStretch {
          swap(&Kts_l, &Kts_r)
        }
      } else {
        Kts_c = .zero
      }
      guard let Kts_c else {
        // The central torsion-stretch parameter was missing.
        throw createTorsionError()
      }
      
      // If all parameters besides Kts_c[2] are zero, send this to the cheaper
      // non-extended torsion force.
      var Kts3: Float
      if V4 == nil, V6 == nil, Kbtb == nil,
         Ktb_l == nil, Ktb_r == nil,
         Kts_l == nil, Kts_r == nil,
         Kts_c[0] == 0, Kts_c[1] == 0 {
        Kts3 = Kts_c[2]
        
        // There are no extended parameters for this torsion.
        torsions.extendedParameters.append(nil)
      } else {
        Kts3 = 0.000
        
        let Ks_l = Kts_l ?? .zero
        let Ks_c = Kts_c
        let Ks_r = Kts_r ?? .zero
        let Kb_l = Ktb_l ?? .zero
        let Kb_r = Ktb_r ?? .zero
        torsions.extendedParameters.append(
          MM4TorsionExtendedParameters(
            V4: V4 ?? 0, V6: V6 ?? 0,
            Kts1: (Ks_l[0], Ks_c[0], Ks_r[0]),
            Kts2: (Ks_l[1], Ks_c[1], Ks_r[1]),
            Kts3: (Ks_l[2], Ks_c[2], Ks_r[2]),
            Ktb1: (Kb_l[0], Kb_r[0]),
            Ktb2: (Kb_l[1], Kb_r[1]),
            Ktb3: (Kb_l[2], Kb_r[2]),
            Kbtb: Kbtb ?? 0))
      }
      torsions.parameters.append(
        MM4TorsionParameters(
          V1: V1, Vn: Vn, V3: V3, n: n, Kts3: Kts3))
    }
  }
}
