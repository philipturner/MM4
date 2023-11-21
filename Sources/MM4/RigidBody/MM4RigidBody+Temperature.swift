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
  public func setTemperature(_ temperature: Double, enthalpy: Double? = nil) {
    // C is heat capacity
    // E is thermal energy
    // N is number of atoms
    // E = C N kT
    let kT = MM4BoltzInZJPerK * temperature
    
    // 8.314462618 J/mol-K
    
    fatalError("Not implemented.")
  }
}
