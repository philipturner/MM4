//
//  Constants.swift
//
//
//  Created by Philip Turner on 10/14/23.
//

/// 0.001660539
///
/// millidyne-angstrom = attojoule
public let MM4AJPerKJPerMol: Double = 1000 / 6.02214076e23 / 1e-18

/// 1.660539
public let MM4ZJPerKJPerMol: Double = 1000 / 6.02214076e23 / 1e-21

/// 602.14
///
/// millidyne-angstrom = attojoule
public let MM4KJPerMolPerAJ: Double = 1e-18 * 6.02214076e23 / 1000

/// 0.602214
public let MM4KJPerMolPerZJ: Double = 1e-21 * 6.02214076e23 / 1000

/// 0.001
public let MM4AJPerZJ: Double = 0.001

/// 1000
public let MM4ZJPerAJ: Double = 1000

/// 0.2081943
///
/// Converts from debyes to elementary charge \* angstroms.
public let MM4EAngstromPerDebye: Double = 0.2081943

/// 4.803205
///
/// Converts from elementary charge \* angstroms to debyes.
public let MM4DebyePerEAngstrom: Double = 1 / 0.2081943

/// 0.008314
public let MM4BoltzInKJPerMolPerK: Double = 8.314462618 / 1000

/// 0.013806
public let MM4BoltzInZJPerK: Double = 8.314462618 / 1000 * MM4ZJPerKJPerMol

/// 1.660539
///
/// Converts from atomic mass units to yoctograms.
public let MM4YgPerAmu: Double = MM4ZJPerKJPerMol

/// 0.602214
///
/// Converts from yoctograms to atomic mass units.
public let MM4AmuPerYg: Double = MM4KJPerMolPerZJ
