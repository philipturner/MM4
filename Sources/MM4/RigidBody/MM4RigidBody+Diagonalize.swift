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
  
  func x(k: Int) -> Double? {
    let phaseShift = Double(k) * (2 * .pi / 3)
    let cubeRoot = Complex(length: length, phase: phase + phaseShift)
    let numerator = Complex(b) + cubeRoot + Complex(Δ0) / cubeRoot
    
    let imaginaryProportion = numerator.imaginary / numerator.real
    guard imaginaryProportion.magnitude < 1e-8 else {
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

// The entered matrix must be diagonalizable. Otherwise, there will be a fatal
// error that describes why it could not be diagonalized.
func diagonalize(
  matrix: (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>)
) -> (
  eigenValues: SIMD3<Double>,
  eigenVectors: (SIMD3<Double>, SIMD3<Double>, SIMD3<Double>)
)? {
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
    // Could not factor the characteristic polynomial.
    return nil
  }
  
  // Sort the eigenvalues in descending order.
  var eigenValues = [root0, root1, root2]
  eigenValues.sort(by: { $0 > $1 })
  
  // Find the eigenvector corresponding to each eigenvalue.
  func createEigenVector(_ eigenValue: Double) -> SIMD3<Double>? {
    var B = (SIMD3<Double>(SIMD3<Float>(matrix.0)),
             SIMD3<Double>(SIMD3<Float>(matrix.1)),
             SIMD3<Double>(SIMD3<Float>(matrix.2)))
    B.0[0] -= eigenValue
    B.1[1] -= eigenValue
    B.2[2] -= eigenValue
    
    // Make 10 attempts to invert the B matrix. Each attempt slightly modifies
    // the diagonal, bringing it closer to something invertible.
    for _ in 0..<10 {
      guard let inverseB = invert(matrix: B) else {
        B.0[0] *= 1 + 1e-10
        B.1[1] *= 1 + 1e-10
        B.2[2] *= 1 + 1e-10
        continue
      }
      return normalize(vector: inverseB.2)
    }
    return nil
  }
  guard var x = createEigenVector(eigenValues[0]),
        var y = createEigenVector(eigenValues[1]),
        var z = createEigenVector(eigenValues[2]) else {
    // Failed to generate the eigenvectors from the eigenvalues.
    return nil
  }
  
  func createOrthogonalityError() -> Double {
    let xyError = (x * y).sum().magnitude
    let xzError = (x * z).sum().magnitude
    let yzError = (y * z).sum().magnitude
    return SIMD3(xyError, xzError, yzError).max()
  }
  
  // Ensure the eigenvectors are mutually perpendicular, then form a
  // right-handed basis from them.
  var orthogonalityError = createOrthogonalityError()
  
  let trialCount = 10
  for trialID in 1...trialCount {
    x = gemv(matrix: matrix, vector: x)
    y = gemv(matrix: matrix, vector: y)
    z = gemv(matrix: matrix, vector: z)
    
    var eigenValueError = SIMD3(
      (x * x).sum().squareRoot(),
      (y * y).sum().squareRoot(),
      (z * z).sum().squareRoot())
    eigenValueError /= SIMD3(eigenValues[0], eigenValues[1], eigenValues[2])
    eigenValueError -= 1
    eigenValueError *= eigenValueError
    
    x = normalize(vector: x)!
    y -= (y * x).sum() * x
    z -= (z * x).sum() * x
    y = normalize(vector: y)!
    z -= (z * y).sum() * y
    z = normalize(vector: z)!
    
    orthogonalityError = createOrthogonalityError()
    
    if orthogonalityError.magnitude < 1e-16,
       eigenValueError.max() < 1e-16 {
      break
    } else if trialID == trialCount {
      // Failed to refine eigenpairs.
      return nil
    }
  }
  
  // Flip the eigenvectors so they're as close as possible to the original
  // coordinate space's cardinal axes.
  if (x * SIMD3(1, 1, 1)).sum() < 0 {
    x = -x
  }
  if (y * SIMD3(1, 1, 1)).sum() < 0 {
    y = -y
  }
  let s1 = x[1] * y[2] - x[2] * y[1]
  let s2 = x[2] * y[0] - x[0] * y[2]
  let s3 = x[0] * y[1] - x[1] * y[0]
  let crossProduct = SIMD3(s1, s2, s3)
  if (crossProduct * z).sum() < 0 {
    z = -z
  }
  
  return (
    SIMD3(eigenValues[0], eigenValues[1], eigenValues[2]),
    (x, y, z))
}
