//
//  MM4Error.swift
//
//
//  Created by Philip Turner on 10/20/23.
//

import OpenMM

/// An error returned by the MM4 library.
public enum MM4Error: Error {
  /// The system's energy exploded.
  ///
  /// Includes the amount of energy drift during the time interval that failed.
  case energyDrift(Double)
}
