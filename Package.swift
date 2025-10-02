// swift-tools-version: 6.1

import PackageDescription

let package = Package(
  name: "MM4",
  products: [
    .library(
      name: "MM4",
      targets: ["MM4"]),
  ],
  dependencies: [
    .package(url: "https://github.com/apple/swift-atomics", .upToNextMajor(from: "1.3.0")),
    .package(url: "https://github.com/apple/swift-docc-plugin", branch: "main"),
    .package(url: "https://github.com/philipturner/swift-numerics", branch: "Quaternions"),
    .package(url: "https://github.com/philipturner/swift-openmm", branch: "main"),
  ],
  targets: [
    .target(
      name: "MM4",
      dependencies: [
        .product(name: "Atomics", package: "swift-atomics"),
        .product(name: "Numerics", package: "swift-numerics"),
        .product(name: "OpenMM", package: "swift-openmm"),
      ]),
    .testTarget(
      name: "MM4Tests",
      dependencies: ["MM4"]),
  ]
)
