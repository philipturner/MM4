//
//  Utilities.swift
//
//
//  Created by Philip Turner on 12/21/23.
//

import XCTest

func XCTAssertEqual<T: SIMD>(
  _ lhs: T,
  _ rhs: T,
  accuracy: T.Scalar,
  _ message: String = "",
  file: StaticString = #filePath,
  line: UInt = #line
)
where T.Scalar : FloatingPoint {
  for lane in 0..<T.scalarCount {
    XCTAssertEqual(
      lhs[lane], rhs[lane], accuracy: accuracy,
      "v[\(lane)] did not match for '\(message)'.",
      file: file, line: line)
  }
}
