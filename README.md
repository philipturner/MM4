# MM4

General-purpose simulator for molecular nanotechnology.

Documentation: [philipturner.github.io/MM4](https://philipturner.github.io/MM4)

### Simulation

Levels of theory:
- Quantum mechanics
- Molecular mechanics
- Rigid body mechanics

### Supported Atoms

Quantum mechanics supports any element or bonding topology. The following restrictions apply to molecular and rigid body mechanics.

| MM4 Atom Code | 6-ring | 5-ring | 4-ring | 3-ring |
| - | - | - | - | - |
| H             | 5   | n/a | n/a           | n/a           |
| C             | 1   | 123 | not supported | not supported |
| N (trivalent) | 8   | 8   | not supported | not supported |
| O             | 6   | 6   | not supported | not supported |
| F             | 11  | n/a | n/a           | n/a           |
| Si            | 19  | 19  | not supported | not supported |
| P (trivalent) | 25  | 25  | not supported | not supported |
| S             | 15  | 15  | not supported | not supported |
| Ge            | 31  | 31  | not supported | not supported |

The following are officially supported in the current release. Other atoms are experimental.

| MM4 Atom Code | 6-ring | 5-ring | 4-ring | 3-ring |
| - | - | - | - | - |
| H             | 5   | n/a | n/a           | n/a           |
| C             | 1   | 123 | not supported | not supported |

### Supported Bonds

Quantum mechanics supports any element or bonding topology. The following restrictions apply to molecular and rigid body mechanics.

| Element | H | C | N | O | F | Si | P | S | Ge |
| ------- | - | - | - | - | - | - | - | - | - |
| H       |   | X |   |   |   | X |   |   | X |
| C       | X | X | O | O | O | O | O | O | O |
| N       |   | O |   |   |   |   |   |   |   |
| O       |   | O |   |   |   |   |   |   |   |
| F       |   | O |   |   |   |   |   |   |   |
| Si      | X | O |   |   |   | X |   |   |   |
| P       |   | O |   |   |   |   |   |   |   |
| S       |   | O |   |   |   |   |   |   |   |
| Ge      | X | O |   |   |   |   |   |   | X |

The following are officially supported in the current release. Other bonds are experimental.

| Element | H | C |
| ------- | - | - |
| H       |   | X |
| C       | X | X |

Key:
- X = nonpolar sigma bond
- O = polar sigma bond

### Forced Motions

|         | Velocity             | Force           |
| ------- | -------------------- | --------------- |
| Linear  | anchor with velocity | external force\* |
| Angular | flywheel             | linear to rotary converter |

> \*There is currently an undocumented API restriction that prevents external forces from being set on a per-atom granularity. One can only use a per-rigid body granularity, with handles for selecting atoms that are affected. A future API change will fix this restriction.

### Units

The following internally consistent unit systems are used for MM4 and OpenMM. The `MM4` module defines several constants for converting between units.

|        | MM4   | OpenMM    |
| ------ | ----- | --------- |
| Angle  | rad   | rad       |
| Energy | zJ    | kJ/mol    |
| Force  | pN    | kJ/mol/nm |
| Mass   | yg\*  | amu       |
| Length | nm    | nm        |
| Speed  | nm/ps | nm/ps     |
| Time   | ps    | ps        |

> \*Mass is currently in amu, but should be changed to yoctograms for consistency with other units. Expect a breaking API change for all mass-related APIs like `hydrogenMassRepartitioning`.

### Releases

Current version: v1.0.0-beta0

v1.0.0
- Accurate simulation of 5-ring carbons
- Anchors
- External forces
- Experimental, untested support for non-carbon elements

Future versions:
- High-precision energy measurements
- Quantum mechanics
- Rigid body mechanics
- Tested support for non-carbon elements

## Tips

List:
- Compile this package in Swift release mode. Vectorized code is known to be extremely slow in debug mode. However, it may not be a bottleneck for small enough systems (under 1000 atoms).
