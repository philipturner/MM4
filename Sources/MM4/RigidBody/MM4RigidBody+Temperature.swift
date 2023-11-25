//
//  MM4RigidBody+Temperature.swift
//
//
//  Created by Philip Turner on 11/20/23.
//

extension MM4RigidBody {
  /// Estimate of the heat capacity in kT.
  ///
  /// This may not be the most appropriate number for characterizing thermal
  /// properties. Molecular dynamics does not simulate certain quantum effects,
  /// such as freezing of higher energy vibrational modes. Freezing is the
  /// primary reason for diamond's exceptionally low heat capacity. Perform
  /// simulations at both 3 kT and the estimated heat capacity (0.75-2.5 kT).
  /// Report whether the system functions correctly and efficiently in both
  /// sets of conditions.
  ///
  /// Heat capacity is derived from data for C, SiC, and Si. The object is
  /// matched to one of these materials based on its elemental composition.
  /// - Elements with Z=6 to Z=8 are treated like carbon.
  /// - Elements with Z=14 to Z=32 are treated like silicon.
  /// - Heat capacity of octane (0.87 kT) is close to diamond (0.74 kT) at 298 K
  ///   ([Gang et al., 1998](https://doi.org/10.1016/S0301-0104(97)00369-8)).
  ///   Therefore, hydrogens and halogens likely have the same thermodynamic
  ///   characteristics as whatever bulk they are attached to. These atoms are
  ///   omitted from the enthalpy derivation.
  /// - The elemental composition is mapped to a spectrum: 100% carbon to
  ///   100% silicon. Moissanite falls at the halfway point. The result is
  ///   interpolated between the two closest materials.
  public func heatCapacity(temperature: Double) -> Double {
    var vNumCarbons: MM4UInt32Vector = .zero
    var vNumSilicons: MM4UInt32Vector = .zero
    atomicNumbers.withUnsafeBufferPointer { buffer in
      let rawBaseAddress = OpaquePointer(buffer.baseAddress)
      let vBaseAddress = UnsafeRawPointer(rawBaseAddress)!
        .assumingMemoryBound(to: MM4UInt8Vector.self)
      
      for vID in 0..<storage.atoms.vectorCount {
        let atomicNumber = MM4UInt32Vector(
          truncatingIfNeeded: vBaseAddress[vID])
        
        let isHalogen =
        (atomicNumber .== 1) .|
        (atomicNumber .== 9) .|
        (atomicNumber .== 17) .|
        (atomicNumber .== 35) .|
        (atomicNumber .== 53)
        
        var isCarbon =
        (atomicNumber .>= 6) .&
        (atomicNumber .<= 8) .&
        (.!isHalogen)
        
        var isSilicon =
        (atomicNumber .>= 14) .&
        (atomicNumber .<= 32) .&
        (.!isHalogen)
        
        if vID == storage.atoms.vectorCount - 1 {
          isCarbon = isCarbon .& storage.atoms.vectorMask
          isSilicon = isSilicon .& storage.atoms.vectorMask
        }
        vNumCarbons.replace(with: vNumCarbons &+ 1, where: isCarbon)
        vNumSilicons.replace(with: vNumSilicons &+ 1, where: isSilicon)
      }
    }
    
    let numCarbons = vNumCarbons.wrappedSum()
    let numSilicons = vNumSilicons.wrappedSum()
    let numTotal = numCarbons + numSilicons
    guard numTotal > 0 else {
      fatalError("Zero non-hydrogen or non-halogen atoms.")
    }
    let closenessC = Double(numCarbons) / Double(numTotal)
    let closenessSi = Double(numSilicons) / Double(numTotal)
    
    func lookupCapacity(keys: [Double], values: [Double]) -> Double {
      let key = binarySearch(keys, element: temperature)
      let value =
      values[key.lowerIndex] * key.lowerWeight +
      values[key.upperIndex] * key.upperWeight
      return value
    }
    var capacity: Double = 0
    
    if closenessC < 0.5 {
      capacity += 2 * (closenessC - 0.5) * lookupCapacity(
        keys: diamondLookupTemperatures, values: diamondLookupCapacities)
      capacity += 2 * (1 - closenessC) * lookupCapacity(
        keys: moissaniteLookupTemperatures, values: moissaniteLookupCapacities)
    } else {
      capacity += 2 * (closenessSi - 0.5) * lookupCapacity(
        keys: siliconLookupTemperatures, values: siliconLookupCapacities)
      capacity += 2 * (1 - closenessSi) * lookupCapacity(
        keys: moissaniteLookupTemperatures, values: moissaniteLookupCapacities)
    }
    return capacity
  }
  
  /// Set the thermal kinetic energy to match a given temperature, assuming
  /// positions are energy-minimized at 0 K.
  ///
  /// - Parameter temperature: The temperature to match the thermal energy to,
  ///   in kelvin.
  /// - Parameter heatCapacity: The partitioning of overall thermal energy in
  ///   thermodynamic units per atom (kT or R).
  ///
  /// Some of the energy will be lost to thermal potential energy during a
  /// simulation. This information can technically be recovered from the atoms'
  /// positions. Typical use cases minimize the system at 0 K, then initialize
  /// the simulator at room temperature. It is not anticipated that users will
  /// extract temperature (e.g. local temperature differentials) from the
  /// simulation.
  ///
  /// > WARNING: There is no trivial method to translate between thermal energy
  /// and temperature. Therefore, you must find a heat capacity lookup table
  /// from an external source. Diamond has
  /// [significantly different](https://physics.stackexchange.com/a/583043) heat
  /// capacity characteristics than other solids. In 1957, C. V. Raman devised a
  /// [theoretical function](http://dspace.rri.res.in/bitstream/2289/1763/1/1957%20Proc%20Indian%20Acad%20Sci%20A%20V46%20p323-332.pdf)
  /// to map temperature to heat capacity for diamond. Experimental measurements
  /// matched the prediction with around 1% margin of error.
  ///
  /// Heat capacity in kT/atom equals the value in J/mol-K divided by 8.314. For
  /// reference, here are some common heat capacities:
  /// - By the equipartition theorem, ideal gases are 1.5 kT.
  /// - Most crystalline solids approach 3.0 kT at high temperatures.
  /// - Diamond: 0.74 kT at 298 K, 1.62 kT at 500 K ([Raman, 1957](http://dspace.rri.res.in/bitstream/2289/1763/1/1957%20Proc%20Indian%20Acad%20Sci%20A%20V46%20p323-332.pdf)).
  /// - Moissanite: 1.62 kT at 298 K, 2.31 kT at 500 K ([Chekhovskoy, 1971](https://doi.org/10.1016/S0021-9614(71)80045-9)).
  /// - Silicon: 2.41 kT at 298 K, 2.84 kT at 500 K ([Desai, 1985](https://srd.nist.gov/JPCRD/jpcrd298.pdf)).
  ///
  /// ![Material Heat Capacities](MaterialHeatCapacities)
  public mutating func setThermalKineticEnergy(
    temperature: Double,
    heatCapacity: Double
  ) {
    // E = thermal energy
    // C = heat capacity
    // N = number of atoms
    // k = Boltzmann constant
    // T = temperature
    //
    // E = C N kT
    let N = Double(storage.atoms.count)
    let kT = MM4BoltzInZJPerK * temperature
    energy.kinetic.thermal = heatCapacity * N * kT
  }
}

// MARK: - Utilities

fileprivate typealias BinaryReturn = (
  lowerIndex: Int, lowerWeight: Double,
  upperIndex: Int, upperWeight: Double
)

/*
 // https://en.wikipedia.org/wiki/Binary_search_algorithm
 function binary_search(A, n, T) is
     L := 0
     R := n − 1
     while L ≤ R do
         m := floor((L + R) / 2)
         if A[m] < T then
             L := m + 1
         else if A[m] > T then
             R := m − 1
         else:
             return m
     return unsuccessful
 */
fileprivate func binarySearch(
  _ array: [Double],
  element: Double
) -> BinaryReturn {
  func equalLower(_ lower: Int) -> BinaryReturn {
    return (lower, 1.000, lower + 1, 0.000)
  }
  func equalUpper(_ upper: Int) -> BinaryReturn {
    return (upper - 1, 0.000, upper, 1.000)
  }
  if element < array.first! {
    return equalLower(0)
  } else if element > array.last! {
    return equalUpper(array.count - 1)
  } else if element.isNaN {
    fatalError("Invalid search element: \(element)")
  }
  
  var L = 0
  var R = array.count - 1
  while L <= R {
    let m = (L + R) / 2
    if array[m] < element {
      L = m + 1
    } else if array[m] > element {
      R = m - 1
    } else {
      break
    }
  }
  
  if L <= R {
    let m = (L + R) / 2
    guard array[m] == element else {
      fatalError("Middle of binary search was invalid.")
    }
    if m == 0 {
      return equalLower(0)
    } else {
      return equalUpper(m)
    }
  } else {
    swap(&L, &R)
    guard L < R else {
      fatalError("Swapped L was not less than swapped R.")
    }
    guard L + 1 == R else {
      fatalError("L and R were not separated by 1.")
    }
    
    let distanceL = element - array[L]
    let distanceR = array[R] - element
    let distanceLR = distanceL + distanceR
    guard distanceL > 0, distanceR > 0, distanceLR > 0 else {
      fatalError("Distances were either zero or negative.")
    }
    
    let closenessL = 1 - distanceL / distanceLR
    let closenessR = 1 - distanceR / distanceLR
    return (L, closenessL, R, closenessR)
  }
}

/// Units: kelvin.
fileprivate let diamondLookupTemperatures: [Double] = [
  0,
  1,
  2,
  3,
  4,
  5,
  6,
  7,
  8,
  9,
  10,
  11,
  12,
  13,
  14,
  15,
  20,
  25,
  30,
  40,
  50,
  60,
  70,
  80,
  90,
  100,
  120,
  140,
  160,
  175,
  200,
  225,
  250,
  275,
  300,
  350,
  400,
  450,
  500,
  600,
  700,
  800,
  900,
  1000,
  1100,
]

/// Units: kelvin.
fileprivate let moissaniteLookupTemperatures: [Double] = [
  0,
  1,
  2,
  3,
  4,
  5,
  6,
  7,
  8,
  9,
  10,
  11,
  12,
  13,
  14,
  15,
  20,
  30,
  40,
  50,
  60,
  70,
  80,
  90,
  100,
  110,
  120,
  130,
  140,
  150,
  160,
  170,
  180,
  190,
  200,
  210,
  220,
  230,
  240,
  250,
  260,
  270,
  273.15,
  280,
  290,
  298.15,
  300,
  400,
  500,
  600,
  700,
  800,
  900,
  1000,
  1100,
  1200,
  1300,
  1400,
  1500,
  1600,
  1700,
  1800,
  1900,
  2000,
  2100,
  2200,
  2300,
  2400,
  2500,
  2600,
  2700,
  2800,
  2900,
]

/// Units: kelvin.
fileprivate let siliconLookupTemperatures: [Double] = [
  0,
  1,
  2,
  3,
  4,
  5,
  6,
  7,
  8,
  9,
  10,
  15,
  20,
  25,
  30,
  35,
  40,
  45,
  50,
  60,
  70,
  75,
  80,
  90,
  100,
  110,
  120,
  125,
  130,
  140,
  150,
  160,
  170,
  175,
  180,
  190,
  200,
  210,
  220,
  225,
  230,
  240,
  250,
  260,
  270,
  273.15,
  280,
  290,
  298.15,
  300,
  350,
  400,
  450,
  500,
  550,
  600,
  650,
  700,
  750,
  800,
  850,
  900,
  950,
  1000,
  1100,
  1200,
  1300,
  1400,
  1500,
  1600,
  1687,
]

/// Units: kT at the corresponding temperature.
fileprivate let diamondLookupCapacities: [Double] = [
  0,
  0.00000002087545327,
  0.0000001670036262,
  0.0000005636372384,
  0.00000133602901,
  0.000002609431659,
  0.000004509097907,
  0.000007160280473,
  0.00001068823208,
  0.00001521820544,
  0.00002087545327,
  0.00002778522831,
  0.00003607278326,
  0.00004586337084,
  0.00005728224378,
  0.0000704546548,
  0.0001660716863,
  0.0003220784219,
  0.0005535722877,
  0.001308443589,
  0.002566562425,
  0.004428578302,
  0.007297089247,
  0.01142371903,
  0.01741236469,
  0.02596757277,
  0.05223709406,
  0.09229559779,
  0.1458411354,
  0.1937503007,
  0.2868510945,
  0.3919795045,
  0.5016371422,
  0.6268451287,
  0.7500401251,
  0.9939138802,
  1.224401251,
  1.432242483,
  1.615927833,
  1.913850373,
  2.135279288,
  2.301854222,
  2.425653115,
  2.521773394,
  2.595247534,
]

/// Units: kT at the corresponding temperature.
fileprivate let moissaniteLookupCapacities: [Double] = [
  0,
  0.000000283076738,
  0.000002264613904,
  0.000007643071927,
  0.00001811691123,
  0.00003538459225,
  0.00006114457541,
  0.00009709532115,
  0.0001449352899,
  0.000206362942,
  0.000283076738,
  0.0003767751383,
  0.0004891566033,
  0.0006219195935,
  0.0007767625692,
  0.0009553839909,
  0.002264613904,
  0.007297089247,
  0.0176136637,
  0.0344724561,
  0.05887996151,
  0.09184267501,
  0.1401544383,
  0.1962665384,
  0.2556497474,
  0.3220784219,
  0.3917782054,
  0.4639942266,
  0.5357070002,
  0.6089295165,
  0.6821520327,
  0.7581424104,
  0.8318681742,
  0.9045874429,
  0.9763002165,
  1.048516238,
  1.120732259,
  1.189425547,
  1.256860717,
  1.321276401,
  1.382420977,
  1.441300938,
  1.459417849,
  1.505213375,
  1.568370941,
  1.616682704,
  1.627250902,
  2.054004811,
  2.310157806,
  2.482016839,
  2.606822228,
  2.702439259,
  2.780442627,
  2.845864806,
  2.903738273,
  2.951546789,
  2.994322829,
  3.032066394,
  3.064777484,
  3.094972336,
  3.12265095,
  3.145297089,
  3.165426991,
  3.183040654,
  3.200654318,
  3.215751744,
  3.225816695,
  3.238397883,
  3.248462834,
  3.258527784,
  3.266076497,
  3.27362521,
  3.278657686,
]

/// Units: kT at the corresponding temperature.
fileprivate let siliconLookupCapacities: [Double] = [
  0,
  0.0000009309598268,
  0.000007447678614,
  0.00002516237671,
  0.0000593817657,
  0.0001162496993,
  0.0002007457301,
  0.0003188597546,
  0.0004763050277,
  0.0006780129901,
  0.0009309598268,
  0.003673322107,
  0.01137358672,
  0.0286865528,
  0.05785422179,
  0.09742602839,
  0.1487851816,
  0.205677171,
  0.2652152995,
  0.3885013231,
  0.5113062305,
  0.5731296608,
  0.6354342074,
  0.7569160452,
  0.8667308155,
  0.9965119076,
  1.112581188,
  1.16923262,
  1.225162377,
  1.335939379,
  1.437936012,
  1.536685109,
  1.630502766,
  1.674404619,
  1.718306471,
  1.80298292,
  1.881404859,
  1.951286986,
  2.015997113,
  2.04919413,
  2.079745008,
  2.138320904,
  2.191363964,
  2.240197258,
  2.287106086,
  2.301299014,
  2.331128217,
  2.374067837,
  2.406422901,
  2.413519365,
  2.564830407,
  2.677171037,
  2.76641809,
  2.83714217,
  2.892350253,
  2.937214337,
  2.976665865,
  3.012990137,
  3.047149387,
  3.080105846,
  3.111017561,
  3.140485927,
  3.168631224,
  3.19569401,
  3.246090931,
  3.290834737,
  3.332571566,
  3.373225884,
  3.412557133,
  3.448881405,
  3.479672841,
]

