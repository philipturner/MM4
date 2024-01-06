# MM4

Molecular mechanics simulator for molecular nanotechnology.

Documentation: [philipturner.github.io/MM4](https://philipturner.github.io/MM4)

### Atoms

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

The following are officially supported in the current release. Other atoms are unsupported.

| MM4 Atom Code | 6-ring | 5-ring | 4-ring | 3-ring |
| - | - | - | - | - |
| H             | 5   | n/a | n/a           | n/a           |
| C             | 1   | 123 | not supported | not supported |
| Si            | 19  | 19  | n/a           | n/a           |

### Bonds

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

The following are officially supported in the current release. Other bonds are unsupported.

| Element | H | C | Si |
| ------- | - | - | - |
| H       |   | X | X |
| C       | X | X |   |
| Si      | x |   | X |

Key:
- X = nonpolar sigma bond
- O = polar sigma bond

### Forces

Available forces:
- bend ❌
- bend-bend ❌
- nonbonded ❌
  - van der Waals force ❌
  - overlap repulsion ❌
  - electrostatic force ❌
- stretch ❌
- stretch-bend ❌
- stretch-stretch ❌
- torsion ❌
- torsion-bend ❌
- torsion-stretch ❌

### Units

The following unit systems are used for MM4 and OpenMM. The `MM4` module defines several constants for converting between them.

| Unit   | MM4   | OpenMM    |
| ------ | ----- | --------- |
| Angle  | rad   | rad       |
| Energy | zJ    | kJ/mol    |
| Force  | pN    | kJ/mol/nm |
| Mass   | yg    | amu       |
| Length | nm    | nm        |
| Speed  | nm/ps | nm/ps     |
| Time   | ps    | ps        |

| Value in SI Units | SI Unit | MM4   | OpenMM    |
| ----------------- | ------- | ----- | --------- |
| Angle             | rad     | 1     | 1         |
| Energy            | J       | 1e-21 | 1.66e-21  |
| Force             | N       | 1e-12 | 1.66e-12  |
| Mass              | kg      | 1e-27 | 1.66e-27  |
| Length            | m       | 1e-9  | 1e-9      |
| Speed             | m/s     | 1000  | 1000      |
| Time              | s       | 1e-12 | 1e-12     |

MM4's unit system is internally consistent. Units for force and energy are derived from mass and velocity.

```
energy = 0.5 * m * v^2 = (10^-27) (10^3)^2 = 10^-21
force = dU / dx = (10^-21) / (10^-9) = 10^-12
```
