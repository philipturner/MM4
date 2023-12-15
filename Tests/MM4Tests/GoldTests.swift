import XCTest
import MM4

final class GoldTests: XCTestCase {
  static func createGoldParameters() throws -> MM4Parameters {
    var paramsDesc = MM4ParametersDescriptor()
    paramsDesc.atomicNumbers = goldAtoms.map { UInt8($0.w) }
    paramsDesc.bonds = []
    return try MM4Parameters(descriptor: paramsDesc)
  }
  
  static func createGoldPositions() -> [SIMD3<Float>] {
    return goldAtoms.map { SIMD3<Float>($0.x, $0.y, $0.z) }
  }
  
  func testGold() throws {
    // Check that the number of atoms and their elements are the same.
    let params = try Self.createGoldParameters()
    XCTAssertEqual(goldAtoms.map { UInt8($0.w) }, params.atoms.atomicNumbers)
    XCTAssertEqual(goldAtoms.count, params.atoms.count)
    XCTAssertEqual(goldAtoms.indices, params.atoms.indices)
    
    // Check the atom code, atomic mass, and vdW parameters.
    for atomID in params.atoms.indices {
      let code = params.atoms.codes[atomID]
      XCTAssertEqual(code.rawValue, MM4AtomCode.gold.rawValue)
      
      let paramsMass = Double(params.atoms.masses[atomID])
      let mass: Double = 196.9665695 * MM4YgPerAmu
      XCTAssert(abs(paramsMass - mass) < 1e-5,
                "\(paramsMass) - \(mass)")
      
      let parameters = params.atoms.parameters[atomID]
      XCTAssertEqual(parameters.charge, 0)
      XCTAssertEqual(parameters.epsilon.default, Float(0.078))
      XCTAssertEqual(parameters.radius.default, Float(2.074))
    }
  }
}

private let goldAtoms: [SIMD4<Float>] = [
  SIMD4<Float>(0.408, 0.612, 0.204, 79),
  SIMD4<Float>(0.612, 0.408, 0.204, 79),
  SIMD4<Float>(0.612, 0.612, 0.000, 79),
  SIMD4<Float>(0.816, 0.612, 0.204, 79),
  SIMD4<Float>(0.612, 0.816, 0.204, 79),
  SIMD4<Float>(0.408, 0.204, 0.612, 79),
  SIMD4<Float>(0.612, 0.000, 0.612, 79),
  SIMD4<Float>(0.612, 0.204, 0.408, 79),
  SIMD4<Float>(0.816, 0.204, 0.612, 79),
  SIMD4<Float>(0.000, 0.612, 0.612, 79),
  SIMD4<Float>(0.204, 0.408, 0.612, 79),
  SIMD4<Float>(0.204, 0.612, 0.408, 79),
  SIMD4<Float>(0.408, 0.408, 0.408, 79),
  SIMD4<Float>(0.408, 0.612, 0.612, 79),
  SIMD4<Float>(0.612, 0.408, 0.612, 79),
  SIMD4<Float>(0.612, 0.612, 0.408, 79),
  SIMD4<Float>(0.816, 0.408, 0.408, 79),
  SIMD4<Float>(0.816, 0.612, 0.612, 79),
  SIMD4<Float>(1.020, 0.408, 0.612, 79),
  SIMD4<Float>(1.020, 0.612, 0.408, 79),
  SIMD4<Float>(1.223, 0.612, 0.612, 79),
  SIMD4<Float>(0.204, 0.816, 0.612, 79),
  SIMD4<Float>(0.408, 0.816, 0.408, 79),
  SIMD4<Float>(0.408, 1.020, 0.612, 79),
  SIMD4<Float>(0.612, 0.816, 0.612, 79),
  SIMD4<Float>(0.612, 1.020, 0.408, 79),
  SIMD4<Float>(0.816, 0.816, 0.408, 79),
  SIMD4<Float>(0.816, 1.020, 0.612, 79),
  SIMD4<Float>(1.020, 0.816, 0.612, 79),
  SIMD4<Float>(0.612, 1.223, 0.612, 79),
  SIMD4<Float>(0.612, 0.204, 0.816, 79),
  SIMD4<Float>(0.204, 0.612, 0.816, 79),
  SIMD4<Float>(0.408, 0.408, 0.816, 79),
  SIMD4<Float>(0.408, 0.612, 1.020, 79),
  SIMD4<Float>(0.612, 0.408, 1.020, 79),
  SIMD4<Float>(0.612, 0.612, 0.816, 79),
  SIMD4<Float>(0.816, 0.408, 0.816, 79),
  SIMD4<Float>(0.816, 0.612, 1.020, 79),
  SIMD4<Float>(1.020, 0.612, 0.816, 79),
  SIMD4<Float>(0.408, 0.816, 0.816, 79),
  SIMD4<Float>(0.612, 0.816, 1.020, 79),
  SIMD4<Float>(0.612, 1.020, 0.816, 79),
  SIMD4<Float>(0.816, 0.816, 0.816, 79),
  SIMD4<Float>(0.612, 0.612, 1.223, 79),
]
