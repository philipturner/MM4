//
//  MM4Vector.swift
//  
//
//  Created by Philip Turner on 11/22/23.
//

let MM4VectorWidth: Int = 4
typealias MM4DoubleVector = SIMD4<Double>
typealias MM4FloatVector = SIMD4<Float>
typealias MM4Int8Vector = SIMD4<Int8>
typealias MM4Int32Vector = SIMD4<Int32>
typealias MM4UInt8Vector = SIMD4<UInt8>
typealias MM4UInt32Vector = SIMD4<UInt32>

extension MM4RigidBody {
  @inline(__always)
  func extractScalar(_ scalarID: Int, array: [MM4FloatVector]) -> SIMD3<Float> {
    let vID = scalarID / MM4VectorWidth
    let lane = scalarID &- vID &* MM4VectorWidth
    let x = array[vID &* 3 &+ 0][lane]
    let y = array[vID &* 3 &+ 1][lane]
    let z = array[vID &* 3 &+ 2][lane]
    return SIMD3(x, y, z)
  }
  
  @inline(__always)
  func swizzleFromVectorWidth<T: BinaryFloatingPoint>(
    _ vectorized: (MM4FloatVector, MM4FloatVector, MM4FloatVector),
    _ vID: Int, baseAddress: UnsafeMutablePointer<SIMD3<T>>
  ) {
    let (x, y, z) = vectorized
    if vID == atoms.vectorCount &- 1 {
      let remaining = atoms.count &- vID &* MM4VectorWidth
      for lane in 0..<remaining {
        var vector3: SIMD3<T> = .zero
        vector3.x = T(x[lane])
        vector3.y = T(y[lane])
        vector3.z = T(z[lane])
        baseAddress[vID &* MM4VectorWidth &+ lane] = vector3
      }
    } else {
      for lane in 0..<MM4VectorWidth {
        var vector3: SIMD3<T> = .zero
        vector3.x = T(x[lane])
        vector3.y = T(y[lane])
        vector3.z = T(z[lane])
        baseAddress[vID &* MM4VectorWidth &+ lane] = vector3
      }
    }
  }
  
  @inline(__always)
  func swizzleToVectorWidth<T: BinaryFloatingPoint>(
    _ vID: Int, baseAddress: UnsafePointer<SIMD3<T>>
  ) -> (MM4FloatVector, MM4FloatVector, MM4FloatVector) {
    var x: MM4FloatVector = .zero
    var y: MM4FloatVector = .zero
    var z: MM4FloatVector = .zero
    
    if vID == atoms.vectorCount &- 1 {
      let remaining = atoms.count &- vID &* MM4VectorWidth
      for lane in 0..<remaining {
        let vector3 = baseAddress[vID &* MM4VectorWidth &+ lane]
        x[lane] = Float(vector3.x)
        y[lane] = Float(vector3.y)
        z[lane] = Float(vector3.z)
      }
    } else {
      for lane in 0..<MM4VectorWidth {
        let vector3 = baseAddress[vID &* MM4VectorWidth &+ lane]
        x[lane] = Float(vector3.x)
        y[lane] = Float(vector3.y)
        z[lane] = Float(vector3.z)
      }
    }
    return (x, y, z)
  }
  
  @inline(__always)
  func withMasses<T>(_ closure: (UnsafePointer<MM4FloatVector>) -> T) -> T {
    self.masses.withUnsafeBufferPointer {
      let rawMasses = OpaquePointer($0.baseAddress)
      let vMasses = UnsafeRawPointer(rawMasses)!
        .assumingMemoryBound(to: MM4FloatVector.self)
      return closure(vMasses)
    }
  }
  
  @inline(__always)
  func withSegmentedLoop(chunk: Int, _ closure: (Range<Int>) -> Void) {
    var loopEnd = 0
    while loopEnd < atoms.vectorCount {
      let loopStart = loopEnd
      loopEnd = min(loopEnd &+ chunk, atoms.vectorCount)
      closure(loopStart..<loopEnd)
    }
  }
}
