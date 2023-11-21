//
//  MM4RigidBody+Temperature.swift
//  
//
//  Created by Philip Turner on 11/20/23.
//

extension MM4RigidBody {
  /// Set the thermal energy to match a given temperature.
  ///
  /// - Parameter enthalpy: The material's heat capacity at the specified
  ///   temperature, in kT/atom.
  /// - Parameter temperature: The temperature to randomize thermal velocites
  ///   at, in kelvin.
  ///
  /// > WARNING: There is no trivial method to translate thermal energy into
  /// temperature. Therefore, you must look up the heat capacity and cite the
  /// source. Diamond has
  /// [significantly different](https://physics.stackexchange.com/a/583043) heat
  /// capacity characteristics than other solids. In 1957, C. V. Raman devised a
  /// [theoretical function](http://dspace.rri.res.in/bitstream/2289/1763/1/1957%20Proc%20Indian%20Acad%20Sci%20A%20V46%20p323-332.pdf)
  /// to map temperature to heat capacity. Experimental measurements matched the
  /// prediction with a very small margin of error.
  ///
  /// Heat capacity in kT equals the number of J/mol-K divided by 8.314. For
  /// reference, here are some common heat capacities:
  /// - Most gases are 1.5 kT at room temperature.
  /// - Most crystalline solids approach 3.0 kT at high temperatures.
  /// - Raman calculated the heat capacity of diamond as ~0.889 kT at 298 K.
  /// - Silicon carbide has ~1.25x the heat capacity of diamond at 500 K, and
  ///   ~1.40x the heat capacity at 1000 K. Elemental silicon is even smaller.
  ///
  /// Next, take the integral of heat capacity over the entire temperature
  /// range. The result is the absolute deviation from zero-point vibrational
  /// energy, which is assumed to be zero. This process be quite laborious, so
  /// it has already been some for common diamondoid materials: diamond,
  /// moissanite, and silicon. If enthalpy is not specified, it will be
  /// interpolated from these three materials' curves, based on the rigid body's
  /// elemental composition.
  ///
  /// ![Material Enthalpies](MaterialEnthalpies)
  public func setTemperature(_ temperature: Double, enthalpy: Double? = nil) {
    // C is heat capacity
    // E is thermal energy
    // N is number of atoms
    // E = C N kT
    let kT = MM4BoltzInZJPerK * temperature
    
    // 8.314462618 J/mol-K
    
    fatalError("Not implemented.")
  }
  
  public func getTemperature(enthalpy: Double? = nil) -> Double {
    // C is heat capacity
    // E is thermal energy
    // N is number of atoms
    // E = C N kT
    // T = E / C N k
    fatalError("Not implemented.")
  }
}
