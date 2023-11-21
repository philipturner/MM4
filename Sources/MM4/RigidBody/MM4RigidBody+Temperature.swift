//
//  MM4RigidBody+Temperature.swift
//  
//
//  Created by Philip Turner on 11/20/23.
//

extension MM4RigidBody {
  /// Set the thermal energy to match a given temperature.
  ///
  /// - Parameter enthalpy: The material's thermodynamic enthalpy at the specified
  ///   temperature, in kT per atom.
  /// - Parameter temperature: The temperature to match the thermal energy to,
  ///   in kelvin.
  ///
  /// Typical use case: minimize at zero kelvin, initialize the simulator at
  /// room temperature, and never attempt to extract temperature from the
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
  /// Enthalpy in kT equals the number of J/mol divided by 8.314. For
  /// reference, here are some common enthalpies:
  /// - By the equipartition theorem, ideal gases are 1.5 kT.
  /// - Most crystalline solids approach 3.0 kT at high temperatures.
  /// - Diamond: 0.21 kT at 298 K, 0.61 kT at 500 K ([Raman, 1957](http://dspace.rri.res.in/bitstream/2289/1763/1/1957%20Proc%20Indian%20Acad%20Sci%20A%20V46%20p323-332.pdf)).
  /// - Moissanite: 0.66 kT at 298 K, 1.21 kT at 500 K ([Chekhovskoy, 1971](https://doi.org/10.1016/S0021-9614(71)80045-9)).
  /// - Silicon: 0.65 kT at 298 K, 0.92 kT at 500 K ([Desai, 1985](https://srd.nist.gov/JPCRD/jpcrd298.pdf)).
  ///
  /// Take the integral of heat capacity over the entire temperature
  /// range. The result is the absolute deviation from zero-point energy, or the
  /// enthalpy. This process be quite laborious, so it has already been done for
  /// common diamondoid materials: diamond (C), moissanite (SiC), and silicon
  /// (Si). The chart shows enthalpies over the range from absolute zero to
  /// thermal decomposition temperature.
  ///
  /// ![Material Enthalpies](MaterialEnthalpies)
  ///
  /// If enthalpy is not specified, it will be derived from data for C, SiC, and
  /// Si.
  /// - Hydrogen and halogens are omitted (see explanation below).
  /// - Elements with Z=6 to Z=8 are treated like carbon.
  /// - Elements with Z=14 to Z=32 are treated like silicon.
  /// - The elemental composition is mapped to a spectrum: 100% carbon to
  ///   100% silicon. The enthalpy is interpolated between the two closest
  ///   materials.
  ///
  /// Diamond has a moderately low heat capacity at 298 K, roughly 0.8 kT.
  /// When integrated over temperature, that translates to an extremely low
  /// enthalpy, just 0.2 kT. This reflects that diamond is an excellent heat
  /// conductor. It expels thermal energy that would pile up in other materials.
  ///
  /// Heat capacity of octane is 0.87 kT/atom at 298 K, close to diamond
  /// ([Gang et al., 1998](https://sci-hub.live/https://doi.org/10.1016/S0301-0104(97)00369-8)).
  /// Enthalpy likely matches diamond as well. Therefore, hydrogens and
  /// halogens likely have the same thermodynamic characteristics as whatever
  /// bulk they are attached to. These atoms are omitted from the enthalpy
  /// derivation.
  public mutating func setTemperature(
    _ temperature: Double,
    enthalpy: Double? = nil
  ) {
    // TODO: Derive enthalpy from elemental composition.
    
    // H is enthalpy
    // E is thermal energy
    // N is number of atoms
    // E = H N kT
    let kT = MM4BoltzInZJPerK * temperature
    thermalEnergy = enthalpy! * Double(atomCount) * kT
  }
  
  /// Units: kelvin.
  fileprivate static let diamondLookupTemperatures: [Double] = [
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
  fileprivate static let moissaniteLookupTemperatures: [Double] = [
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
  fileprivate static let siliconLookupTemperatures: [Double] = [
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
  fileprivate static let diamondLookupEnthalpies: [Double] = [
    0,
    0.00000001043738984,
    0.0000000521869492,
    0.0000001565608476,
    0.0000003548712546,
    0.0000006784303397,
    0.000001158550272,
    0.000001826543222,
    0.000002713721359,
    0.000003851396851,
    0.00000527088187,
    0.000007003488583,
    0.000009080529162,
    0.00001153331577,
    0.00001439316059,
    0.00001769137578,
    0.00004283337047,
    0.00008308013208,
    0.0001422019813,
    0.0003393959603,
    0.0006590048661,
    0.001132080306,
    0.001807875205,
    0.002751903567,
    0.004048089463,
    0.005812207402,
    0.01136035145,
    0.02006087468,
    0.03243633093,
    0.04420952307,
    0.06871995066,
    0.09879599475,
    0.1335957859,
    0.172743708,
    0.2157167667,
    0.3094642095,
    0.4094214054,
    0.5115166943,
    0.6127686228,
    0.8047792115,
    0.9790249674,
    1.133958743,
    1.270594149,
    1.390898077,
    1.497037153,
  ]
  
  /// Units: kT at the corresponding temperature.
  fileprivate static let moissaniteLookupEnthalpies: [Double] = [
    0,
    0.00000006919430319,
    0.0000005535544255,
    0.000001868246186,
    0.000004428435404,
    0.000008649287898,
    0.00001494596949,
    0.00002373364599,
    0.00003542748323,
    0.00005044264702,
    0.00006919430319,
    0.00009209761754,
    0.0001195677559,
    0.0001520198841,
    0.0001898691679,
    0.0002335307733,
    0.0005535544255,
    0.00183679423,
    0.00435295071,
    0.008504608901,
    0.01480338729,
    0.02329241998,
    0.03459715159,
    0.04948441076,
    0.067181378,
    0.08715051079,
    0.1096624866,
    0.1341304954,
    0.1603151128,
    0.1877052734,
    0.2163894572,
    0.2459912982,
    0.2764976398,
    0.3076332369,
    0.3393037013,
    0.3713128062,
    0.403843115,
    0.4364984501,
    0.4692631834,
    0.5020235408,
    0.5347800271,
    0.5672534997,
    0.5774770262,
    0.599564144,
    0.6318155684,
    0.6580911698,
    0.6640136949,
    0.9630588925,
    1.208761573,
    1.407370191,
    1.570081643,
    1.705954093,
    1.821138145,
    1.920582241,
    2.007206646,
    2.084006604,
    2.152281542,
    2.213858251,
    2.269573145,
    2.320210794,
    2.366667183,
    2.409219829,
    2.448485113,
    2.484704524,
    2.518433002,
    2.54978148,
    2.578950993,
    2.606109074,
    2.631698385,
    2.655512838,
    2.678029213,
    2.699117001,
    2.71901075,
  ]
  
  /// Units: kT at the corresponding temperature.
  fileprivate static let siliconLookupEnthalpies: [Double] = [
    0,
    0.0000002327324468,
    0.000001163662234,
    0.000003493191747,
    0.000007903732211,
    0.00001510427566,
    0.00002579461309,
    0.00004066640831,
    0.00006043120491,
    0.00008577998134,
    0.0001174250057,
    0.0004619611177,
    0.001286872295,
    0.003032440181,
    0.006132782741,
    0.0108022152,
    0.01714579034,
    0.02508656128,
    0.03434983706,
    0.05586217792,
    0.08001681391,
    0.09275570797,
    0.1058416773,
    0.1327566389,
    0.1600708368,
    0.1878640016,
    0.2161466899,
    0.2303182241,
    0.2444821168,
    0.2727515926,
    0.3007979179,
    0.3284750033,
    0.3557277337,
    0.3691698038,
    0.3824748205,
    0.4086758269,
    0.4342953967,
    0.4592404751,
    0.4834472272,
    0.4952876226,
    0.5069596629,
    0.5297731125,
    0.5518776396,
    0.5732613672,
    0.5939475875,
    0.6003261982,
    0.6139709882,
    0.6333603009,
    0.648715153,
    0.6521317115,
    0.7373209279,
    0.8088803821,
    0.8704178545,
    0.9234966379,
    0.9698769419,
    1.010431552,
    1.046504859,
    1.078627589,
    1.107830406,
    1.134209763,
    1.158617667,
    1.181048598,
    1.201941313,
    1.221346132,
    1.256722345,
    1.288207109,
    1.316467152,
    1.342193485,
    1.365812667,
    1.387644618,
    1.40545128,
  ]
}
