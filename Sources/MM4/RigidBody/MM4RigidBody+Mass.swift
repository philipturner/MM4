//
//  MM4RigidBody+Mass.swift
//
//
//  Created by Philip Turner on 11/20/23.
//

extension MM4RigidBody {
  mutating func createMasses() {
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
  
  // implementation of center of mass
  
  // public API that shifts all the positions by the offset, and changes the
  // center of mass accordingly (accepting slight inaccuracy from floating-point
  // roundoff error, positions drifting away from cached center of mass)
}
