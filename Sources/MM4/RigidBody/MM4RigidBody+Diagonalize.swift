//
//  MM4RigidBodyStorage+Utilities.swift
//  MM4
//
//  Created by Philip Turner on 11/25/23.
//

import ComplexModule

// Source: https://en.wikipedia.org/wiki/Cubic_equation#General_cubic_formula
func factorCubicPolynomial(coefficients: SIMD4<Double>) -> (
  Double?, Double?, Double?
) {
  let a = coefficients[0]
  let b = coefficients[1]
  let c = coefficients[2]
  let d = coefficients[3]
  
  let Δ0 = b * b - 3 * a * c
  let Δ1 = 2 * b * b * b - 9 * a * b * c + 27 * a * a * d
  
  // The square root term may be negative, producing an imaginary number.
  let squareRootTerm = Δ1 * Δ1 - 4 * Δ0 * Δ0 * Δ0
  let squareRootMagnitude = squareRootTerm.magnitude.squareRoot()
  var cubeRootTerm = Complex(Δ1)
  if squareRootTerm < 0 {
    cubeRootTerm.imaginary += squareRootMagnitude
  } else {
    cubeRootTerm.real += squareRootMagnitude
  }
  cubeRootTerm /= 2
  
  // The cube root is a complex number.
  var (length, phase) = cubeRootTerm.polar
  length = Double.root(length, 3)
  phase /= 3
  
  let cbrtTermSIMD = SIMD2(cubeRootTerm.real, cubeRootTerm.imaginary)
  let cbrtTermLengthSq = (cbrtTermSIMD * cbrtTermSIMD).sum()
  if cbrtTermLengthSq.magnitude < .leastNormalMagnitude {
    // Handle the case of a 3-fold repeated real root. We don't yet handle the
    // case of a 2-fold repeated real root.
    if Δ0.magnitude < .leastNormalMagnitude {
      let repeatedRoot = b / Double(-3 * a)
      return (repeatedRoot, repeatedRoot, repeatedRoot)
    }
  }
  
  func x(k: Int) -> Double? {
    let phaseShift = Double(k) * (2 * .pi / 3)
    let cubeRoot = Complex(length: length, phase: phase + phaseShift)
    let numerator = Complex(b) + cubeRoot + Complex(Δ0) / cubeRoot
    
    let imaginaryProportion = numerator.imaginary / numerator.real
    guard imaginaryProportion.magnitude < 1e-4 else {
      return nil
    }
    return numerator.real / Double(-3 * a)
  }
  return (x(k: 0), x(k: 1), x(k: 2))
}

@_transparent
func gemv(
  matrix: (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>),
  vector: SIMD3<Double>
) -> SIMD3<Double> {
  matrix.0 * vector[0] + matrix.1 * vector[1] + matrix.2 * vector[2]
}

@_transparent
func determinant(
  matrix: (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>)
) -> Double {
  matrix.0[0] * (matrix.1[1] * matrix.2[2] - matrix.2[1] * matrix.1[2]) -
  matrix.0[1] * (matrix.1[0] * matrix.2[2] - matrix.1[2] * matrix.2[0]) +
  matrix.0[2] * (matrix.1[0] * matrix.2[1] - matrix.1[1] * matrix.2[0])
}

@_transparent
func cross(
  leftVector x: SIMD3<Double>,
  rightVector y: SIMD3<Double>
) -> SIMD3<Double> {
  SIMD3(
    x[1] * y[2] - x[2] * y[1],
    x[2] * y[0] - x[0] * y[2],
    x[0] * y[1] - x[1] * y[0])
}

@_transparent
func normalize(
  vector: SIMD3<Double>
) -> SIMD3<Double>? {
  let length = (vector * vector).sum().squareRoot()
  guard length >= .leastNormalMagnitude else {
    return nil
  }
  return vector / length
}

// Source: https://stackoverflow.com/a/18504573
func invert(
  matrix: (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>)
) -> (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>)? {
  let determinant = determinant(matrix: matrix)
  guard determinant.magnitude > .leastNormalMagnitude else {
    return nil
  }
  
  let entry00 = matrix.1[1] * matrix.2[2] - matrix.2[1] * matrix.1[2]
  let entry10 = matrix.0[2] * matrix.2[1] - matrix.0[1] * matrix.2[2]
  let entry20 = matrix.0[1] * matrix.1[2] - matrix.0[2] * matrix.1[1]
  
  let entry01 = matrix.1[2] * matrix.2[0] - matrix.1[0] * matrix.2[2]
  let entry11 = matrix.0[0] * matrix.2[2] - matrix.0[2] * matrix.2[0]
  let entry21 = matrix.1[0] * matrix.0[2] - matrix.0[0] * matrix.1[2]
  
  let entry02 = matrix.1[0] * matrix.2[1] - matrix.2[0] * matrix.1[1]
  let entry12 = matrix.2[0] * matrix.0[1] - matrix.0[0] * matrix.2[1]
  let entry22 = matrix.0[0] * matrix.1[1] - matrix.1[0] * matrix.0[1]
  
  let column0 = SIMD3(entry00, entry10, entry20) / determinant
  let column1 = SIMD3(entry01, entry11, entry21) / determinant
  let column2 = SIMD3(entry02, entry12, entry22) / determinant
  return (column0, column1, column2)
}

func diagonalize(
  matrix: (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>)
) -> (
  eigenValues: SIMD3<Double>?,
  eigenVectors: (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>)?,
  failureReason: String?
) {
  // Characteristic polynomial:
  // [-1, tr(A), -0.5 * [tr(A)^2 - tr(A^2)], det(A)]
  var characteristicPolynomial: SIMD4<Double> = .zero
  do {
    let matrix2 = (
      gemv(matrix: matrix, vector: matrix.0),
      gemv(matrix: matrix, vector: matrix.1),
      gemv(matrix: matrix, vector: matrix.2))
    let trace = matrix.0[0] + matrix.1[1] + matrix.2[2]
    let trace2 = matrix2.0[0] + matrix2.1[1] + matrix2.2[2]
    
    characteristicPolynomial[0] = -1
    characteristicPolynomial[1] = trace
    characteristicPolynomial[2] = -0.5 * (trace * trace - trace2)
    characteristicPolynomial[3] = determinant(matrix: matrix)
  }
  let (root0, root1, root2) = factorCubicPolynomial(
    coefficients: characteristicPolynomial)
  guard let root0,
        let root1,
        let root2 else {
    return (nil, nil, """
      Failed to factor characteristic polynomial:
      a = \(characteristicPolynomial[0])
      b = \(characteristicPolynomial[1])
      c = \(characteristicPolynomial[2])
      d = \(characteristicPolynomial[3])
      root0 = \(root0 == nil ? String(describing: root0) : "nil")
      root1 = \(root1 == nil ? String(describing: root1) : "nil")
      root2 = \(root2 == nil ? String(describing: root2) : "nil")
      """)
  }
  
  // Sort the eigenvalues in descending order.
  var eigenValuesArray = [root0, root1, root2]
  eigenValuesArray.sort(by: { $0 > $1 })
  let eigenValues = SIMD3(eigenValuesArray[0],
                          eigenValuesArray[1],
                          eigenValuesArray[2])
  
  // Find the eigenvector corresponding to each eigenvalue. Start by
  // preconditioning the eigensolver with direct matrix inversion, if possible.
  var candidateEigenVectors: [SIMD3<Double>] = []
  for i in 0..<3 {
    let eigenValue = eigenValues[i]
    var B = (SIMD3<Double>(SIMD3<Float>(matrix.0)),
             SIMD3<Double>(SIMD3<Float>(matrix.1)),
             SIMD3<Double>(SIMD3<Float>(matrix.2)))
    B.0[0] -= eigenValue
    B.1[1] -= eigenValue
    B.2[2] -= eigenValue
    
    // Make 2 attempts to invert the B matrix. The second attempt slightly
    // modifies the diagonal, bringing it closer to something invertible.
    for _ in 0..<2 {
      if let inverseB = invert(matrix: B),
         let column = normalize(vector: inverseB.2) {
        candidateEigenVectors.append(column)
        break
      } else {
        B.0[0] *= 1 + 1e-10
        B.1[1] *= 1 + 1e-10
        B.2[2] *= 1 + 1e-10
      }
    }
  }
  
  // If we cannot find all eigenvectors through direct matrix inversion, start
  // with the slow path that rotates the cardinal axes into the eigenbasis.
  var x: SIMD3<Double>
  var y: SIMD3<Double>
  var z: SIMD3<Double>
  if candidateEigenVectors.count == 3 {
    x = candidateEigenVectors[0]
    y = candidateEigenVectors[1]
    z = candidateEigenVectors[2]
  } else {
    x = SIMD3(0.0, 0.0, 1.0)
    y = SIMD3(0.0, 1.0, 0.0)
    z = SIMD3(1.0, 0.0, 0.0)
  }
  
  // Execute something like conjugate gradient descent.
  for trialID in 1...160 {
    // Perform an Arnoldi iteration on each candidate vector.
    @_transparent
    func createEigenPair(_ vector: SIMD3<Double>) -> SIMD4<Double> {
      let Av = gemv(matrix: matrix, vector: vector)
      let λ = (Av * Av).sum().squareRoot()
      return SIMD4(Av, λ)
    }
    var eigenPair0 = createEigenPair(x)
    var eigenPair1 = createEigenPair(y)
    var eigenPair2 = createEigenPair(z)
    
    // Sort the eigenpairs in descending order.
    if eigenPair1.w > eigenPair0.w {
      let temp = eigenPair1
      eigenPair1 = eigenPair0
      eigenPair0 = temp
    }
    if eigenPair2.w > eigenPair0.w {
      let temp = eigenPair2
      eigenPair2 = eigenPair0
      eigenPair0 = temp
    }
    if eigenPair2.w > eigenPair1.w {
      let temp = eigenPair2
      eigenPair2 = eigenPair1
      eigenPair1 = temp
    }
    
    // Check how close we are to the solution.
    let revisedValues = SIMD3(eigenPair0.w, eigenPair1.w, eigenPair2.w)
    var eigenValueError = revisedValues / eigenValues
    eigenValueError -= 1
    eigenValueError = SIMD3(
      eigenValueError[0].magnitude,
      eigenValueError[1].magnitude,
      eigenValueError[2].magnitude)
    
    // Orthonormalize with the Gram-Schmidt method.
    x = unsafeBitCast(eigenPair0, to: SIMD3<Double>.self)
    y = unsafeBitCast(eigenPair1, to: SIMD3<Double>.self)
    z = unsafeBitCast(eigenPair2, to: SIMD3<Double>.self)
    x /= revisedValues[0]
    y -= (y * x).sum() * x
    z -= (z * x).sum() * x
    
    // Handle the case of a 2-fold repeated root.
    if (y * y).sum() < 1e-3 * 1e-3 {
      guard (z * z).sum().magnitude > 1e-3,
            (z * x).sum().magnitude < 1e-3 else {
        fatalError("Could not use z as a reference to fix y.")
      }
      y = cross(leftVector: z, rightVector: x)
    }
    
    y /= (y * y).sum().squareRoot()
    z -= (z * y).sum() * y
    
    // Handle the case of a 2-fold repeated root.
    if (z * z).sum() < 1e-3 * 1e-3 {
      guard (y * y).sum().magnitude > 1e-3,
            (x * y).sum().magnitude < 1e-3 else {
        fatalError("Could not use y as a reference to fix z.")
      }
      z = cross(leftVector: x, rightVector: y)
    }
    
    z /= (z * z).sum().squareRoot()
    
    // Check how close we are to orthonormal.
    let orthogonalityError = SIMD3(
      (x * y).sum().magnitude,
      (x * z).sum().magnitude,
      (y * z).sum().magnitude)
    
    if trialID <= 10 {
      if orthogonalityError.max() < 1e-16,
         eigenValueError.max() < 1e-8 {
        break
      }
    } else if trialID < 40 {
      if orthogonalityError.max() < 1e-12,
         eigenValueError.max() < 1e-6 {
        break
      }
    } else if trialID < 160 {
      if orthogonalityError.max() < 1e-8,
         eigenValueError.max() < 1e-4 {
        break
      }
    } else {
      guard orthogonalityError.max() < 1e-4,
            eigenValueError.max() < 1e-2 else {
        return (nil, nil, """
          Failed to refine eigenpairs after 160 iterations:
          λ0 = \(eigenValues[0]) -> \(revisedValues[0]) v0 = \(x) error0 = \(eigenValueError[0])
          λ1 = \(eigenValues[1]) -> \(revisedValues[1]) v1 = \(y) error1 = \(eigenValueError[1])
          λ2 = \(eigenValues[2]) -> \(revisedValues[2]) v2 = \(z) error2 = \(eigenValueError[2])
          Orthogonality error: \(orthogonalityError)
          """)
      }
    }
  }
  
  // Flip the first two vectors so they're close to the original space's
  // cardinal axes. Flip the last vector so it forms a right-handed basis.
  if (x * SIMD3(1, 1, 1)).sum() < 0 {
    x = -x
  }
  if (y * SIMD3(1, 1, 1)).sum() < 0 {
    y = -y
  }
  z = cross(leftVector: x, rightVector: y)
  
  return (eigenValues, (x, y, z), nil)
}
