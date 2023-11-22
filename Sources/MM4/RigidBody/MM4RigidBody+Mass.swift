//
//  MM4RigidBody+Mass.swift
//
//
//  Created by Philip Turner on 11/20/23.
//

/// Wrapper class to bypass the issue of mutating a `struct`. A mutation would
/// occur when lazily recomputing the center of mass.
class MM4CenterOfMass {
  var value: SIMD3<Double>? = nil
}

extension MM4RigidBody {
  mutating func createMasses() {
    let capacity = atomVectorCount * MM4VectorWidth
    self.masses = Array(unsafeUninitializedCapacity: capacity) {
      self.mass = 0
      for i in 0..<atomCount {
        let atomicNumber = atomicNumbers[i]
        let mass = MM4MassParameters.global.mass(atomicNumber: atomicNumber)
        $0.baseAddress.unsafelyUnwrapped[i] = mass
        self.mass += Double(mass)
      }
      
      // Pad with zeroes to remove the need for bounds checking in certain
      // vectorized computations.
      for i in atomCount..<atomVectorCount * MM4VectorWidth {
        $0.baseAddress.unsafelyUnwrapped[i] = 0
      }
      $1 = atomCount
    }
    
    self.masses.reserveCapacity(atomVectorCount * MM4VectorWidth)
    self.masses = atomicNumbers.map { atomicNumber in
      MM4MassParameters.global.mass(atomicNumber: atomicNumber)
    }
    
    for bond in bonds {
      let atomicNumber1 = atomicNumbers[Int(bond[0])]
      let atomicNumber2 = atomicNumbers[Int(bond[1])]
      guard atomicNumber1 == 1 || atomicNumber2 == 1 else {
        continue
      }
      if atomicNumber1 == atomicNumber2 {
        fatalError("Hydrogen cannot be bonded to another hydrogen.")
      }
      
      let hydrogen = (atomicNumber1 == 1) ? bond[0] : bond[1]
      let nonHydrogen = (atomicNumber1 == 1) ? bond[1] : bond[0]
      masses[Int(hydrogen)] += hydrogenMassRepartitioning
      masses[Int(nonHydrogen)] -= hydrogenMassRepartitioning
    }
  }
  
  func ensureCenterOfMassCached() {
    guard _centerOfMass.value == nil else {
      return
    }
    fatalError("Not implemented.")
  }
  
  // implementation of center of mass
  
  // public API that shifts all the positions by the offset, and changes the
  // center of mass accordingly (accepting slight inaccuracy from floating-point
  // roundoff error, positions drifting away from cached center of mass)
}
