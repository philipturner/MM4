//
//  MM4RigidBody+Temperature.swift
//  
//
//  Created by Philip Turner on 11/20/23.
//

extension MM4RigidBody {
  /// Set the thermal energy to match a given temperature.
  ///
  /// - Parameter heatCapacity: The material's heat capacity in kT at the
  ///   specified temperature.
  /// - Parameter temperature: The temperature to randomize thermal velocites
  ///   at, in kelvin.
  ///
  /// Heat capacity in kT equals the number of J/mol-K divided by 8.314. For
  /// reference, here are some common heat capacities:
  /// - Most gases are 1.5 kT at room temperature.
  /// - Most crystalline solids approach 3.0 kT at high temperatures.
  /// - Raman calculated the heat capacity of diamond as ~0.889 kT at 298 K.
  ///
  /// > WARNING: There is no trivial method to translate thermal energy into
  /// temperature. Therefore, you must look up the heat capacity and cite the
  /// source. Diamond has
  /// [significantly different](https://physics.stackexchange.com/a/583043) heat
  /// capacity characteristics than other solids. In 1957, C. V. Raman devised a
  /// [theoretical function](http://dspace.rri.res.in/bitstream/2289/1763/1/1957%20Proc%20Indian%20Acad%20Sci%20A%20V46%20p323-332.pdf)
  /// to map temperature to heat capacity. Experimental measurements matched the
  /// prediction with a very small margin of error.
  public func setTemperature(_ temperature: Double, heatCapacity: Double) {
    fatalError("Not implemented.")
  }
  
  public func getTemperature(heatCapacity: Double) -> Double {
    fatalError("Not implemented.")
  }
}
