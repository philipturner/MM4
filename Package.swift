// swift-tools-version: 5.9
// The swift-tools-version declares the minimum version of Swift required to build this package.

import PackageDescription

let package = Package(
  name: "MM4",
  products: [
    // Products define the executables and libraries a package produces, making them visible to other packages.
    .library(
      name: "MM4",
      targets: ["MM4"]),
  ],
  dependencies: [
    .package(url: "https://github.com/philipturner/swift-openmm", branch: "main"),
    .package(url: "https://github.com/apple/swift-docc-plugin", branch: "main"),
  ],
  targets: [
    // Targets are the basic building blocks of a package, defining a module or a test suite.
    // Targets can depend on other targets in this package and products from dependencies.
    .target(
      name: "MM4",
      dependencies: [
        .product(name: "OpenMM", package: "swift-openmm"),
      ]),
    .testTarget(
      name: "MM4Tests",
      dependencies: ["MM4"]),
  ]
)
