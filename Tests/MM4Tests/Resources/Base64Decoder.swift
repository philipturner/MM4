//
//  Base64Decoder.swift
//  MM4Tests
//
//  Created by Philip Turner on 1/16/24.
//

import Foundation

struct Base64Decoder {
  static func decodeAtoms(_ string: String) -> [SIMD4<Float>] {
    let options: Data.Base64DecodingOptions = [
      .ignoreUnknownCharacters
    ]
    guard let data = Data(base64Encoded: string, options: options) else {
      fatalError("Could not decode the data.")
    }
    guard data.count % 16 == 0 else {
      fatalError("Data did not have the right alignment.")
    }
    
    let rawMemory: UnsafeMutableBufferPointer<SIMD4<Float>> =
      .allocate(capacity: data.count / 16)
    let encodedBytes = data.copyBytes(to: rawMemory)
    guard encodedBytes == data.count else {
      fatalError("Did not encode the right number of bytes.")
    }
    
    let output = Array(rawMemory)
    rawMemory.deallocate()
    return output
  }
  
  static func decodeBonds(_ string: String) -> [SIMD2<UInt32>] {
    let options: Data.Base64DecodingOptions = [
      .ignoreUnknownCharacters
    ]
    guard let data = Data(base64Encoded: string, options: options) else {
      fatalError("Could not decode the data.")
    }
    guard data.count % 8 == 0 else {
      fatalError("Data did not have the right alignment.")
    }
    
    let rawMemory: UnsafeMutableBufferPointer<SIMD2<UInt32>> =
      .allocate(capacity: data.count / 8)
    let encodedBytes = data.copyBytes(to: rawMemory)
    guard encodedBytes == data.count else {
      fatalError("Did not encode the right number of bytes.")
    }
    
    let output = Array(rawMemory)
    rawMemory.deallocate()
    return output
  }
}
