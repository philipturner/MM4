//
//  MM4Vector.swift
//  MM4  
//
//  Created by Philip Turner on 11/22/23.
//

let MM4VectorWidth: Int = 4

typealias MM4FloatVector = SIMD4<Float>
typealias MM4DoubleVector = SIMD4<Double>
typealias MM4Int8Vector = SIMD4<Int8>
typealias MM4Int16Vector = SIMD4<Int16>
typealias MM4Int32Vector = SIMD4<Int32>
typealias MM4Int64Vector = SIMD4<Int64>
typealias MM4UInt8Vector = SIMD4<UInt8>
typealias MM4UInt16Vector = SIMD4<UInt16>
typealias MM4UInt32Vector = SIMD4<UInt32>
typealias MM4UInt64Vector = SIMD4<UInt64>

// Private APIs for creating thermal velocities.

let MM4VectorPairWidth: Int = 8
typealias MM4FloatVectorPair = SIMD8<Float>
typealias MM4UInt16VectorPair = SIMD8<UInt16>
