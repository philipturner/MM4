import XCTest
import MM4

// Tests parameters for carbon-only molecules.
final class AdamantaneTests: XCTestCase {
  func testAdamantane() throws {
    try testAdamantaneVariant(atomCode: .alkaneCarbon)
  }
  
  func testSilaAdamantane() throws {
    try testAdamantaneVariant(atomCode: .silicon)
  }
}

private func testAdamantaneVariant(atomCode: MM4AtomCode) throws {
  let adamantane = Adamantane(atomCode: atomCode)
  var paramsDesc = MM4ParametersDescriptor()
  paramsDesc.atomicNumbers = adamantane.atomicNumbers
  paramsDesc.bonds = adamantane.bonds
  
  let params = try MM4Parameters(descriptor: paramsDesc)
  print()
  print("atomic numbers (Z):", params.atoms.atomicNumbers)
  print("atomic parameters (r, eps, Hred):", params.atoms.parameters.map {
    ($0.radius.default, $0.epsilon.default, $0.hydrogenReductionFactor)
  })
  print("atomic masses:", params.atoms.masses)
}
