//
//  MM4TestsReference.swift
//
//
//  Created by Philip Turner on 12/15/23.
//

// Precursor to the MM4RigidBody test suite, copied into the source tree for
// reference.

#if false
// Save as GitHub Gist instead of polluting the molecular-renderer source tree.

import Foundation
import HDL
import MM4
import MolecularRenderer
import Numerics

// Make sure to save this with the other code in the GitHub gist.
func adamantaneLattice() -> Lattice<Cubic> {
  let lattice = Lattice<Cubic> { h, k, l in
    Bounds { 4 * h + 4 * k + 4 * l }
    Material { .elemental(.silicon) }
    
    Volume {
      Origin { 2 * h + 2 * k + 2 * l }
      Origin { 0.25 * (h + k - l) }
      
      // Remove the front plane.
      Convex {
        Origin { 0.25 * (h + k + l) }
        Plane { h + k + l }
      }
      
      func triangleCut(sign: Float) {
        Convex {
          Origin { 0.25 * sign * (h - k - l) }
          Plane { sign * (h - k / 2 - l / 2) }
        }
        Convex {
          Origin { 0.25 * sign * (k - l - h) }
          Plane { sign * (k - l / 2 - h / 2) }
        }
        Convex {
          Origin { 0.25 * sign * (l - h - k) }
          Plane { sign * (l - h / 2 - k / 2) }
        }
      }
      
      // Remove three sides forming a triangle.
      triangleCut(sign: +1)
      
      // Remove their opposites.
      triangleCut(sign: -1)
      
      // Remove the back plane.
      Convex {
        Origin { -0.25 * (h + k + l) }
        Plane { -(h + k + l) }
      }
      
      Replace { .empty }
    }
  }
  return lattice
}

func adamantaneDiamondoid() -> Diamondoid {
  let lattice = adamantaneLattice()
  
  let latticeAtoms = lattice.entities.map(MRAtom.init)
  var diamondoid = Diamondoid(atoms: latticeAtoms)
  diamondoid.translate(offset: -diamondoid.createCenterOfMass())
  
  // Remove a sidewall carbon, creating two 5-membered rings.
  do {
    #if true
    // Detect the sidewall carbon farthest in Z.
    var maxZValue: Float = -.greatestFiniteMagnitude
    var maxZIndex: Int = -1
    for (index, atom) in diamondoid.atoms.enumerated() {
      if atom.element == 1 {
        continue
      }
      if atom.origin.z > maxZValue {
        maxZValue = atom.origin.z
        maxZIndex = index
      }
    }
    var removedAtoms = [maxZIndex]
    
    // Detect all hydrogens farther in Z than the removed sidewall.
    for (index, atom) in diamondoid.atoms.enumerated() {
      if atom.element != 1 {
        continue
      }
      if atom.origin.z > maxZValue {
        removedAtoms.append(index)
      }
    }
    
    // Create a new bond between the atoms that are about to become free
    // radicals.
    var neighbors: [Int] = []
    for var bond in diamondoid.bonds {
      guard Int(bond[0]) == maxZIndex ||
              Int(bond[1]) == maxZIndex else {
        continue
      }
      if Int(bond[0]) == maxZIndex {
        bond = SIMD2(bond[1], bond[0])
      }
      
      let atom = diamondoid.atoms[Int(bond[0])]
      if atom.element == 1 {
        continue
      }
      neighbors.append(Int(bond[0]))
    }
    guard neighbors.count == 2 else {
      fatalError("Unrecognized number of neighbors.")
    }
    diamondoid.bonds.append(SIMD2(
      Int32(neighbors[0]),
      Int32(neighbors[1])))
    
    // Remove all bonds containing the removed sidewall.
    diamondoid.bonds.removeAll(where: {
      Int($0[0]) == maxZIndex ||
      Int($0[1]) == maxZIndex
    })
    
    // Remove the atoms one at a time, fixing the bonds with a simple
    // O(n^2) method.
    removedAtoms.sort()
    for atomID in removedAtoms.reversed() {
      for bondID in diamondoid.bonds.indices {
        var bond = diamondoid.bonds[bondID]
        if any(bond .== Int32(atomID)) {
          fatalError("A bond remained that contained a removed atom.")
        }
        let shifted = bond &- 1
        bond.replace(with: shifted, where: bond .>= Int32(atomID))
        diamondoid.bonds[bondID] = bond
      }
      diamondoid.atoms.remove(at: atomID)
    }
    #endif
  }
  
  #if false
  var paramsDesc = MM4ParametersDescriptor()
  paramsDesc.atomicNumbers = diamondoid.atoms.map { $0.element }
  paramsDesc.bonds = diamondoid.bonds.map {
    SIMD2<UInt32>(truncatingIfNeeded: $0)
  }
  
  let params = try! MM4Parameters(descriptor: paramsDesc)
  print("atomic numbers (Z):", params.atoms.atomicNumbers)
  print("atomic parameters (r, eps, Hred):", params.atoms.parameters.map {
    ($0.radius.default, $0.epsilon.default, $0.hydrogenReductionFactor)
  })
  print("atom ringTypes:", params.atoms.ringTypes)
  print("rings:", params.rings.indices)
  print()
  print("bond ringTypes:", params.bonds.ringTypes)
  print("bond base parameters (ks, l):", params.bonds.parameters.map { ($0.stretchingStiffness, $0.equilibriumLength) })
  print("bond extended parameters (complex cross-terms):", params.bonds.extendedParameters)
  print()
  print("angle ringTypes:", params.angles.ringTypes)
  print("angle base parameters (kθ, θ, kθθ):",params.angles.parameters.map { ($0.bendingStiffness, $0.equilibriumAngle, $0.bendBendStiffness) })
  print("angle extended parameters (complex cross-terms):", params.angles.extendedParameters)
  print()
  print("torsion ringTypes:", params.torsions.ringTypes)
  print("torsion base parameters (V1, V2, V3, Kts):", params.torsions.parameters.map {
    ($0.V1, $0.Vn, $0.V3, $0.Kts3)
  })
  print("torsion extended parameters (complex cross-terms):", params.torsions.extendedParameters)
  #endif
  
  print("atoms:")
  for atom in diamondoid.atoms {
    var output: String = "SIMD4<Float>("
    for lane in 0..<3 {
      output += String(format: "%.3f", atom.origin[lane])
      output += ", "
    }
    output += "\(atom.element)"
    output += "),"
    print(output)
  }
  print()
  print("bonds:")
  for bond in diamondoid.bonds {
    var output: String = "SIMD2<UInt32>("
    output += "\(bond[0])"
    output += ", "
    output += "\(bond[1])"
    output += "),"
    print(output)
  }
  print()
  
  return diamondoid
}

func example() {
  var diamondoid = adamantaneDiamondoid()
  self.atomProvider = ArrayAtomProvider(diamondoid.atoms)
  
  var paramsDesc = MM4ParametersDescriptor()
  paramsDesc.atomicNumbers = diamondoid.atoms.map { $0.element }
  paramsDesc.bonds = diamondoid.bonds.map {
    SIMD2<UInt32>(truncatingIfNeeded: $0)
  }
  let params = try! MM4Parameters(descriptor: paramsDesc)
  
  var rigidBodyDesc = MM4RigidBodyDescriptor()
  rigidBodyDesc.positions = diamondoid.atoms.map { $0.origin }
  rigidBodyDesc.parameters = params
  rigidBodyDesc.velocities = diamondoid.createVelocities()
  var rigidBody = MM4RigidBody(descriptor: rigidBodyDesc)
  
  let angularVelocity = Quaternion<Float>(angle: 1, axis: [
    0.70710678, 0.70710678, 0
  ])
  diamondoid.angularVelocity = angularVelocity
  rigidBody.angularVelocity = angularVelocity
  do {
    var differences: [SIMD3<Float>] = []
    for (lhs, rhs) in zip(diamondoid.createVelocities(), rigidBody.velocities) {
      differences.append(lhs - rhs)
    }
  }
  rigidBody.velocities.withUnsafeBufferPointer {
    rigidBody.setVelocities($0)
  }
  print("expected angular velocity:", angularVelocity)
  print("expected angular velocity:", angularVelocity.angle * angularVelocity.axis)
  
  // rigid body angular velocity: (0.87758255, 0.33900505, 0.33900505, 0.0)
  
  let repartitionedMasses = rigidBody.parameters.atoms.masses.map(Double.init)
  let statePositions = diamondoid.atoms.map { $0.origin }.map(SIMD3<Double>.init)
  let stateVelocities = diamondoid.createVelocities().map(SIMD3<Double>.init)
  
  var totalMass: Double = 0
  var totalMomentum: SIMD3<Double> = .zero
  var centerOfMass: SIMD3<Double> = .zero
  
  // Conserve mtoomentum after the velocities are randomly initialized.
  for atomID in diamondoid.atoms.indices {
    let mass = repartitionedMasses[atomID]
    let position = statePositions[atomID]
    let velocity = stateVelocities[atomID]
    totalMass += Double(mass)
    totalMomentum += mass * velocity
    centerOfMass += mass * position
  }
  centerOfMass /= totalMass
  
  // Should be SIMD3<Float>(1.3476824e-07, 1.3006077e-07, -0.03191334)
  //           SIMD3<Float>(1.3132467e-07, 1.277512e-07, -0.03190822)
  print("           center of mass:", centerOfMass)
  print("rigid body center of mass:", rigidBody.centerOfMass)
  
  // Conserve angular momentum along the three cardinal axes, therefore
  // conserving angular momentum along any possible axis.
  var totalAngularMomentum: SIMD3<Double> = .zero
  var totalMomentOfInertia: cross_platform_double3x3 = .init(diagonal: .zero)
  for atomID in diamondoid.atoms.indices {
    let mass = repartitionedMasses[atomID]
    let delta = statePositions[atomID] - centerOfMass
    let velocity = stateVelocities[atomID]
    
    // From Wikipedia:
    // https://en.wikipedia.org/wiki/Rigid_body_dynamics#Mass_properties
    //
    // I_R = m * (I (S^T S) - S S^T)
    // where S is the column vector R - R_cm
    let STS = cross_platform_dot(delta, delta)
    var momentOfInertia = cross_platform_double3x3(diagonal: .init(repeating: STS))
    momentOfInertia -= cross_platform_double3x3(rows: [
      SIMD3(delta.x * delta.x, delta.x * delta.y, delta.x * delta.z),
      SIMD3(delta.y * delta.x, delta.y * delta.y, delta.y * delta.z),
      SIMD3(delta.z * delta.x, delta.z * delta.y, delta.z * delta.z),
    ])
    momentOfInertia *= mass
    totalMomentOfInertia += momentOfInertia
    
    // From Wikipedia:
    // https://en.wikipedia.org/wiki/Rigid_body_dynamics#Linear_and_angular_momentum
    //
    // L = m * (R - R_cm) cross d/dt (R - R_cm)
    // assume R_cm is stationary
    // L = m * (R - R_cm) cross v
    let angularMomentum = mass * cross_platform_cross(delta, velocity)
    totalAngularMomentum += angularMomentum
  }
  let totalAngularVelocity = totalMomentOfInertia
    .inverse * totalAngularMomentum
  print("     total angular velocity:", totalAngularVelocity)
  print("rigid body angular velocity:", rigidBody.angularVelocity)
  
  // Next, query the linear velocity, see what happens when you mutate:
  // - center of mass
  // - linear velocity
  //
  // Inertia | Position | Velocity |
  // mass    | centerOfMass | linearVelocity |
  // momentOfInertia | rotate() | angularVelocity |
  
  // Next, call rotate(), which likely will face bugs in the custom quaternion
  // act implementation.
  
  
  // Before the code that computes moment of inertia above, rotate both objects. Confirm that both objects' moment of inertia become a different value, and the rigid body's moment is un-cached upon mutation.
  //
  // Should be different than: (SIMD3<Float>(19.959436, -0.04870774, -2.4214387e-07), SIMD3<Float>(-0.04870774, 19.959436, -2.9802322e-07), SIMD3<Float>(-2.4214387e-07, -2.9802322e-07, 24.405872))
  // And whatever the new value is, now that you have changed the position
  // and orientation:
  print("     total moment of inertia:", totalMomentOfInertia)
  print("rigid body moment of inertia:", rigidBody.momentOfInertia)
  
  // Next, see what happens when you add anchors, handles, and external forces.
  
  // Next, see what happens when you read/write energy, including the heuristic for heat capacity.
  
  // TODO: Switch over to using Swift package tests for MM4; copy the silaadamantane coordinates and bonds to Swift source code.
}
#endif
