//
//  MM4Parameters+Bonds.swift
//  MM4
//
//  Created by Philip Turner on 10/7/23.
//

import Atomics
import Dispatch

/// Parameters for a group of 2 atoms.
public struct MM4Bonds {
  /// Each value corresponds to the bond at the same array index.
  public var extendedParameters: [MM4BondExtendedParameters?] = []
  
  /// Groups of atom indices that form a bond.
  public var indices: [SIMD2<UInt32>] = []
  
  /// Map from a group of atoms to a bond index.
  public var map: [SIMD2<UInt32>: UInt32] = [:]
  
  /// Each value corresponds to the bond at the same array index.
  public var parameters: [MM4BondParameters] = []
  
  /// The smallest ring each bond is involved in.
  public var ringTypes: [UInt8] = []
  
  mutating func append(contentsOf other: Self, atomOffset: UInt32) {
    let bondOffset = UInt32(self.indices.count)
    self.extendedParameters += other.extendedParameters
    self.indices += other.indices.map {
      $0 &+ atomOffset
    }
    for key in other.map.keys {
      let value = other.map[key].unsafelyUnwrapped
      self.map[key &+ atomOffset] = value &+ bondOffset
    }
    self.parameters += other.parameters
    self.ringTypes += other.ringTypes
  }
}

/// Morse stretching parameters for a covalent bond.
public struct MM4BondParameters {
  /// Units: millidyne \* angstrom
  ///
  /// The parameter's name originates from its description in
  /// Nanosystems 3.3.3(a).
  ///
  /// > WARNING: Convert aJ to kJ/mol.
  public var potentialWellDepth: Float
  
  /// Units: millidyne / angstrom
  ///
  /// > WARNING: Convert nanometers to angstroms.
  public var stretchingStiffness: Float
  
  /// Units: angstrom
  ///
  /// > WARNING: Convert nanometers to angstroms.
  public var equilibriumLength: Float
}

/// Parameters for covalent bonds that create partial charges.
public struct MM4BondExtendedParameters {
  /// Units: debye
  public var dipoleMoment: Float
}

extension MM4Parameters {
  /// - throws: `.missingParameter`
  mutating func createBondParameters(forces: MM4ForceOptions) throws {
    for bondID in bonds.indices.indices {
      let bond = bonds.indices[bondID]
      let ringType = bonds.ringTypes[bondID]
      let codes = with5RingsRemoved {
        createAtomCodes(group: bond, zero: SIMD2<UInt8>.zero)
      }
      let sortedCodes = sortBond(codes)
      var sortedIDs = bond
      if any(sortedCodes .!= codes) {
        sortedIDs = SIMD2(bond[1], bond[0])
      }
      
      var potentialWellDepth: Float
      var stretchingStiffness: Float
      var equilibriumLength: Float
      var dipoleMoment: Float?
      
      switch (sortedCodes[0], sortedCodes[1]) {
        // Carbon
      case (1, 1):
        potentialWellDepth = 1.130
        stretchingStiffness = 4.5500
        equilibriumLength = 1.5270
      case (1, 5):
        potentialWellDepth = 0.854
        equilibriumLength = 1.1120
        
        let centerType = atoms.centerTypes[Int(sortedIDs[0])]
        switch centerType {
        case .tertiary:
          stretchingStiffness = 4.7400
        case .secondary:
          stretchingStiffness = 4.6700
        case .primary:
          stretchingStiffness = 4.7400
        default:
          fatalError("This should never happen.")
        }
      case (5, 123):
        potentialWellDepth = 0.854
        equilibriumLength = 1.1120
        
        let centerType = atoms.centerTypes[Int(sortedIDs[1])]
        switch centerType {
        case .tertiary:
          stretchingStiffness = 4.7000
        case .secondary:
          stretchingStiffness = 4.6400
        default:
          fatalError("This should never happen.")
        }
      case (1, 123):
        potentialWellDepth = 1.130
        stretchingStiffness = 4.5600
        equilibriumLength = 1.5270
      case (123, 123):
        potentialWellDepth = 1.130
        stretchingStiffness = (ringType == 5) ? 4.9900 : 4.5600
        equilibriumLength = (ringType == 5) ? 1.5290 : 1.5270
        
        // Nitrogen
      case (1, 8):
        potentialWellDepth = 1.140
        stretchingStiffness = 5.20
        equilibriumLength = 1.4585
        dipoleMoment = (codes[1] == 8) ? +0.64 : -0.64
      case (8, 123):
        potentialWellDepth = 1.140
        stretchingStiffness = 4.90
        equilibriumLength = (ringType == 5) ? 1.4640 : 1.4520
        dipoleMoment = (codes[1] == 123) ? -1.40 : +1.40
        
        // Oxygen
      case (1, 6):
        // The oxygen well depth seems too low, lower than neighbors nitrogen
        // and fluorine. In general, all the Morse parameters here are lower
        // than in Nanosystems, even though they're in the same units (aJ).
        // Looking through the other oxygen atom types in the Morse parameters
        // chart, the pattern remains for sp2 carbon (where O also follows the
        // pattern of being close to H). Either this is an artifact of the
        // CCSD(T) (highly accurate) simulator, or an actual quantum effect that
        // cannot be explained rationally.
        potentialWellDepth = 0.851
        stretchingStiffness = 4.90
        equilibriumLength = 1.4190
        dipoleMoment = (codes[1] == 6) ? +1.160 : -1.160
      case (6, 123):
        potentialWellDepth = 0.851
        stretchingStiffness = 4.90
        equilibriumLength = (ringType == 5) ? 1.4096 : 1.4199
        dipoleMoment = (codes[1] == 123) ? -1.160 : +1.160
        
        // Fluorine
      case (1, 11):
        potentialWellDepth = 0.989
        stretchingStiffness = 6.10
        equilibriumLength = 1.3859
        dipoleMoment = (codes[1] == 11) ? +1.82 : -1.82
        
        // Silicon
      case (1, 19):
        potentialWellDepth = 0.812
        stretchingStiffness = (ringType == 5) ? 2.85 : 3.05
        equilibriumLength = (ringType == 5) ? 1.884 : 1.876
        
        var dipoleMagnitude: Float = (ringType == 5) ? 0.55 : 0.70
        dipoleMagnitude *= (codes[1] == 19) ? -1 : +1
        dipoleMoment = dipoleMagnitude
      case (5, 19):
        potentialWellDepth = 0.777
        stretchingStiffness = 2.65
        equilibriumLength = 1.483
      case (19, 19):
        potentialWellDepth = 0.672
        stretchingStiffness = 1.65
        equilibriumLength = (ringType == 5) ? 2.336 : 2.322
        
        // Phosphorus
      case (1, 25):
        potentialWellDepth = 0.702
        stretchingStiffness = 2.9273
        equilibriumLength = 1.8514
        dipoleMoment = (codes[1] == 25) ? +0.9254 : -0.9254
        
        // Sulfur
      case (1, 15):
        potentialWellDepth = 0.651
        stretchingStiffness = 2.92
        equilibriumLength = 1.814
        dipoleMoment = (codes[1] == 15) ? +0.70 : -0.70
      case (15, 123):
        potentialWellDepth = 0.651
        stretchingStiffness = (ringType == 5) ? 3.20 : 2.92
        equilibriumLength = (ringType == 5) ? 1.821 : 1.814
        dipoleMoment = (codes[1] == 15) ? +0.70 : -0.70
        
        // Germanium
      case (1, 31):
        potentialWellDepth = 0.744
        stretchingStiffness = (ringType == 5) ? 2.95 : 2.72
        equilibriumLength = (ringType == 5) ? 1.944 : 1.949
        
        var dipoleMagnitude: Float = (ringType == 5) ? 0.495 : 0.635
        dipoleMagnitude *= (codes[1] == 31) ? -1 : +1
        dipoleMoment = dipoleMagnitude
      case (5, 31):
        potentialWellDepth = 0.689
        stretchingStiffness = 2.55
        equilibriumLength = 1.529
      case (31, 31):
        potentialWellDepth = 0.542
        stretchingStiffness = 1.45
        equilibriumLength = 2.404
        
      default:
        var addresses: [MM4Address] = []
        for lane in 0..<2 {
          addresses.append(createAddress(bond[lane]))
        }
        throw MM4Error.missingParameter(addresses)
      }
      
      if !forces.contains(.stretch) {
        potentialWellDepth = 0
        stretchingStiffness = 0
      }
      if !forces.contains(.nonbonded) {
        dipoleMoment = nil
      }
      
      bonds.parameters.append(
        MM4BondParameters(
          potentialWellDepth: potentialWellDepth,
          stretchingStiffness: stretchingStiffness,
          equilibriumLength: equilibriumLength))
      
      if let dipoleMoment {
        bonds.extendedParameters.append(
          MM4BondExtendedParameters(dipoleMoment: dipoleMoment))
      } else {
        bonds.extendedParameters.append(nil)
      }
    }
  }
  
  /// - Parameter dipoleMoment: Original dipole moment parameter in elementary
  /// charge-angstroms.
  /// - Parameter bondID: Usage of 32-bit integers for `bondID` reflects that
  /// bond IDs are most often stored in compact 32-bit form, not extended
  /// 64-bit. This is not something exposed in a public API, so the choice is
  /// permissible.
  func projectDipole(_ dipoleMoment: Float, bondID: UInt32) -> SIMD2<Float> {
    // Units: angstrom
    let length = bonds.parameters[Int(bondID)].equilibriumLength
    
    // Units: elementary charge * angstrom
    //
    // Dipole moment is positive if it points from a positive charge to a
    // negative charge. For example, a bond from C -> F should have a positive
    // dipole moment. A bond from F -> C should have a negative moment.
    let partialCharge = dipoleMoment * Float(MM4EAngstromPerDebye) / length
    
    // Dipole points from positive to negative (+->)
    return SIMD2(+partialCharge, -partialCharge)
  }
  
  private func electrostaticEffect(sign: Float) -> [Float] {
    func correction(
      atomID: Int32, endID: Int32, bondID: Int32
    ) -> (
      correction: Float, bohlmann: Float?, decay: Float, beta: Float
    )? {
      // Try once with 5-ring carbons. If that doesn't work, try again with
      // 6-ring carbons.
      for attemptID in 0..<2 {
        let bond = bonds.indices[Int(bondID)]
        let otherID = (bond[0] == endID) ? bond[1] : bond[0]
        var codeActing = atoms.codes[Int(atomID)].rawValue
        var codeEnd = atoms.codes[Int(endID)].rawValue
        var codeOther = atoms.codes[Int(otherID)].rawValue
        if attemptID == 1 {
          codeActing = (codeActing == 123) ? 1 : codeActing
          codeEnd = (codeEnd == 123) ? 1 : codeEnd
          codeOther = (codeOther == 123) ? 1 : codeOther
        }
        
        let bondCodes = (min(codeEnd, codeOther), max(codeEnd, codeOther))
        let presentCodes = SIMD3(codeActing, codeEnd, codeOther)
        var nonCarbonElementCount: Int = 0
        for nonCarbonElement in MM4Parameters.nonCarbonElements {
          if any(presentCodes .== nonCarbonElement) {
            nonCarbonElementCount += 1
          }
        }
        if nonCarbonElementCount > 1 {
          // No parameters for electronegativity correction between two
          // non-carbon elements. This should never happen because torsions are
          // assigned before electronegativity corrections.
          fatalError("Encountered two non-carbon elements while computing electrostatic effect.")
        }
        if nonCarbonElementCount == 0 {
          continue
        }
        
        switch (bondCodes.0, bondCodes.1, codeEnd, codeActing) {
          // Nitrogen
        case (1, 1, 1, 8):       return (-0.0195, nil, 0.62, 0.20)
        case (1, 5, 1, 8):       return (-0.0118, nil, 0.62, 0.20)
        case (1, 8, 8, 1):       return (-0.0015, nil, 0.62, 0.20)
        case (1, 8, 8, 123):     return (-0.0030, nil, 0.62, 0.40)
        case (5, 123, 123, 8):   return (-0.0100, nil, 0.62, 0.20)
        case (8, 123, 8, 1):     return (-0.0200, nil, 0.62, 0.20)
        case (8, 123, 8, 123):   return (0.0000, nil, 0.62, 0.20)
        case (123, 123, 123, 8): return (-0.0140, nil, 0.62, 0.20)
          
          // Oxygen
          //
          // The supplementary information for the oxygen paper has "4%" for the
          // secondary 1-1-1-6 correction. I have never seen anything this low -
          // the lowest was 0.20. Given that the 123-123-123 parameter is 40%,
          // 4% is probably a typo. The correct number should be 40%.
        case (1, 1, 1, 6):       return (-0.0095, nil, 0.62, 0.40)
        case (1, 5, 1, 6):       return (-0.0034, nil, 0.62, 0.40)
        case (1, 6, 1, 6):       return (-0.0217, nil, 0.62, 0.20)
        case (1, 6, 6, 21):      return (0.0146, nil, 0.62, 0.20)
        case (5, 123, 123, 6):   return (-0.0034, nil, 0.62, 0.40)
        case (123, 123, 123, 6): return (-0.0097, nil, 0.62, 0.40)
          
          // Fluorine
        case (1, 1, 1, 11):
          var sum: Float = 0.00
          var decay: Float = 1.00
          var count: Int = 0
          let neighbors = atomsToAtomsMap[Int(endID)]
          for lane in 0..<4 where neighbors[lane] != -1 {
            let atomID = neighbors[lane]
            let otherElement = atoms.atomicNumbers[Int(atomID)]
            if otherElement == 6 {
              count += 1
              sum += decay
              decay *= 0.38
            }
          }
          precondition(
            count > 0,
            "Carbon with a C-C bond didn't have any carbon neighbors.")
          
          let units = sum / Float(count)
          return (-0.0193 * units, nil, 0.38, 0.05)
        case (1, 5, 1, 11):  return (-0.0052, 0.0011, 0.55, 0.30)
        case (1, 11, 1, 1):  return (0.0127, nil, 0.62, 0.40)
        case (1, 11, 1, 11): return (-0.0268, -0.0028, 0.33, 0.67)
          
          // Silicon
        case (1, 1, 1, 19):    return (0.009, nil, 0.62, 0.40)
        case (5, 19, 19, 19):  return (0.003, nil, 0.62, 0.40)
        case (19, 19, 19, 19): return (-0.002, nil, 0.62, 0.40)
        case (19, 19, 19, 5):  return (0.004, nil, 0.62, 0.40)
        case (1, 19, 19, 19):  return (0.009, nil, 0.62, 0.40)
        case (1, 19, 1, 19):   return (-0.004, nil, 0.62, 0.40)
          
          // Phosphorus
        case (1, 25, 25, 1): return (-0.0036, nil, 0.62, 0.40)
        case (1, 5, 1, 25):  return (-0.0015, nil, 0.62, 0.40)
        case (1, 25, 1, 1):  return (0.0005, nil, 0.62, 0.40)
          
          // Sulfur
        case (1, 1, 1, 15):       return (-0.0010, nil, 0.62, 0.40)
        case (1, 5, 1, 15):       return (-0.015, nil, 0.62, 0.40)
        case (5, 123, 123, 15):   return (-0.022, nil, 0.62, 0.40)
        case (123, 123, 123, 15): return (-0.015, nil, 0.62, 0.40)
          
          // Germanium
        case (5, 31, 31, 31): return (0.006, nil, 0.62, 0.40)
        case (5, 31, 31, 1):  return (0.008, nil, 0.62, 0.40)
          
        default: break
        }
      }
      return nil
    }
    
    let bondCapacity = bonds.indices.count
    let primaryNeighborContributions: UnsafeMutablePointer<SIMD2<Float>> =
      .allocate(capacity: 64 * bondCapacity)
    let secondaryNeighborContributions: UnsafeMutablePointer<Float> =
      .allocate(capacity: 64 * bondCapacity)
    let bohlmannEffectContributions: UnsafeMutablePointer<Float> =
      .allocate(capacity: 64 * bondCapacity)
    
    typealias AtomicPointer = UnsafeMutablePointer<UInt16.AtomicRepresentation>
    let primaryNeighborAtomics: AtomicPointer =
      .allocate(capacity: bondCapacity)
    let secondaryNeighborAtomics: AtomicPointer =
      .allocate(capacity: bondCapacity)
    let bohlmannEffectAtomics: AtomicPointer =
      .allocate(capacity: bondCapacity)
    let primaryNeighborCounts = UnsafeMutablePointer<UInt16>(
      OpaquePointer(primaryNeighborAtomics))
    let secondaryNeighborCounts = UnsafeMutablePointer<UInt16>(
      OpaquePointer(secondaryNeighborAtomics))
    let bohlmannEffectCounts = UnsafeMutablePointer<UInt16>(
      OpaquePointer(bohlmannEffectAtomics))
    
    primaryNeighborCounts.initialize(repeating: .zero, count: bondCapacity)
    secondaryNeighborCounts.initialize(repeating: .zero, count: bondCapacity)
    bohlmannEffectCounts.initialize(repeating: .zero, count: bondCapacity)
    defer {
      primaryNeighborAtomics.deallocate()
      secondaryNeighborAtomics.deallocate()
      bohlmannEffectAtomics.deallocate()
      primaryNeighborContributions.deallocate()
      secondaryNeighborContributions.deallocate()
      bohlmannEffectContributions.deallocate()
    }
    
    let taskSize = 128
    let taskCount = (atoms.count + taskSize - 1) / taskSize
    DispatchQueue.concurrentPerform(iterations: taskCount) { z in
      execute(taskID: z)
    }
    
    func execute(taskID: Int) {
      let atomStart = taskID * taskSize
      let atomEnd = min(atomStart + taskSize, atoms.count)
      for atom0 in Int32(atomStart)..<Int32(atomEnd) {
        let zeroLevelAtoms = atomsToAtomsMap[Int(atom0)]
        
        for lane in 0..<4 where zeroLevelAtoms[lane] != -1 {
          let atom1 = zeroLevelAtoms[lane]
          let primaryLevelAtoms = atomsToAtomsMap[Int(atom1)]
          let primaryLevelBonds = atomsToBondsMap[Int(atom1)]
          
          for lane in 0..<4 where primaryLevelAtoms[lane] != -1 {
            let atom2 = primaryLevelAtoms[lane]
            if atom0 == atom2 { continue }
            
            let bondID = primaryLevelBonds[lane]
            let corr = correction(
              atomID: atom0, endID: atom1, bondID: bondID)
            if let corr, corr.correction * sign > 0 {
              // Atomically accumulate the primary neighbor contribution.
              let contribution = SIMD2(corr.correction, corr.decay)
              let atomic = UnsafeAtomic<UInt16>(
                at: primaryNeighborAtomics.advanced(by: Int(bondID)))
              
              let slotID = atomic
                .loadThenWrappingIncrement(ordering: .relaxed)
              if slotID >= 64 {
                fatalError("Exceeded slot capacity for primary neighbors.")
              }
              let address = 64 &* Int(bondID) &+ Int(slotID)
              primaryNeighborContributions[address] = contribution
            }
            if let bohlmann = corr?.bohlmann {
              // Atomically accumulate the Bohlmann Effect contribution.
              let contribution = bohlmann
              let atomic = UnsafeAtomic<UInt16>(
                at: bohlmannEffectAtomics.advanced(by: Int(bondID)))
              
              let slotID = atomic
                .loadThenWrappingIncrement(ordering: .relaxed)
              if slotID >= 64 {
                fatalError("Exceeded slot capacity for Bohlmann Effect.")
              }
              let address = 64 &* Int(bondID) &+ Int(slotID)
              bohlmannEffectContributions[address] = contribution
            }
            
            let secondaryLevelAtoms = atomsToAtomsMap[Int(atom2)]
            let secondaryLevelBonds = atomsToBondsMap[Int(atom2)]
            for lane in 0..<4 where secondaryLevelAtoms[lane] != -1 {
              let atom3 = secondaryLevelAtoms[lane]
              if atom1 == atom3 { continue }
              if atom0 == atom3 { fatalError("This should never happen.") }
              
              let bondID = secondaryLevelBonds[lane]
              let corr = correction(
                atomID: atom0, endID: atom2, bondID: bondID)
              if let corr, corr.correction * sign > 0 {
                // Atomically accumulate the secondary neighbor contribution.
                let contribution = corr.correction * corr.beta
                let atomic = UnsafeAtomic<UInt16>(
                  at: secondaryNeighborAtomics.advanced(by: Int(bondID)))
                
                let slotID = atomic
                  .loadThenWrappingIncrement(ordering: .relaxed)
                if slotID >= 64 {
                  fatalError("Exceeded slot capacity for secondary neighbors.")
                }
                let address = 64 &* Int(bondID) &+ Int(slotID)
                secondaryNeighborContributions[address] = contribution
              }
            }
          }
        }
      }
    }
    
    return bonds.indices.indices.map { bondID in
      var correction: Float = 0
      
      struct PrimaryCorrection: Equatable {
        var correction: Float
        var decay: Float
      }
      
      let neighborCount = Int(primaryNeighborCounts[bondID])
      var neighbors: [PrimaryCorrection] = []
      neighbors.reserveCapacity(neighborCount)
      for slotID in 0..<neighborCount {
        let address = 64 &* bondID &+ slotID
        let contribution = primaryNeighborContributions[address]
        neighbors.append(
          PrimaryCorrection(
            correction: contribution[0], decay: contribution[1]))
      }
      neighbors.sort(by: { $0.correction.magnitude > $1.correction.magnitude })
      
      var decay: Float = 1
      for (neighborID, neighbor) in neighbors.enumerated() {
        if neighborID > 0 {
          decay *= neighbor.decay
        }
        correction += neighbor.correction * decay
      }
      for slotID in 0..<Int(secondaryNeighborCounts[bondID]) {
        let address = 64 &* bondID &+ slotID
        let contribution = secondaryNeighborContributions[address]
        correction += contribution
      }
      for slotID in 0..<Int(bohlmannEffectCounts[bondID]) {
        let address = 64 &* bondID &+ slotID
        let contribution = bohlmannEffectContributions[address]
        correction += contribution
      }
      return correction
    }
  }
  
  mutating func createElectronegativityEffectCorrections() {
    // Add electronegativity corrections to bond length.
    var electronegativeCorrections: [Float] = []
    var electropositiveCorrections: [Float] = []
    DispatchQueue.concurrentPerform(iterations: 2) { z in
      if z == 0 {
        electronegativeCorrections = electrostaticEffect(sign: -1)
      } else if z == 1 {
        electropositiveCorrections = electrostaticEffect(sign: +1)
      }
    }
    
    for i in bonds.indices.indices {
      // We are not adding electronegativity effects to bond stiffness, due to
      // the results being questionable. It would be quite time-consuming to
      // implement and test whether the results were correct, using ab initio
      // calculations performed by a different person with the appropriate
      // expertise. In the future, this may be something to consider. For now,
      // it creates a mostly conservative estimate of machine stiffness.
      //
      // There are also electronegativity effect corrections to equilibrium
      // bond angle. The best approach would be to omit them, but consider
      // which structures would be affected by the term. The greatest example
      // may be hydrofluorocarbon storage tapes.
      var correction: Float = 0
      correction += electronegativeCorrections[i]
      correction += electropositiveCorrections[i]
      bonds.parameters[i].equilibriumLength += correction
    }
  }
  
  mutating func createPartialCharges() {
    // Compute partial charges.
    for (bondID, parameters) in bonds.extendedParameters.enumerated() {
      guard let parameters else { continue }
      let atomCharges = projectDipole(
        parameters.dipoleMoment, bondID: UInt32(truncatingIfNeeded: bondID))
      
      let bond = bonds.indices[bondID]
      for lane in 0..<2 {
        let atomID = bond[lane]
        let partialCharge = atomCharges[lane]
        atoms.parameters[Int(atomID)].charge += partialCharge
      }
    }
  }
  
  mutating func removeBondLengths(forces: MM4ForceOptions) {
    // If the user omitted '.stretch' from the forces, zero out the bond
    // lengths. Do this after the partial charges have been initialized.
    //
    // If and only if they omitted '.nonbonded', zero out the dipole moment.
    // The latter was already done when bond parameters were initialized.
    // To get all the information about how to project a bond dipole onto
    // charges, the user must request parameters for both '.stretch' and
    // '.nonbonded'.
    if !forces.contains(.stretch) {
      for bondID in bonds.indices.indices {
        bonds.parameters[bondID].equilibriumLength = .zero
      }
    }
  }
}
