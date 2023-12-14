# MM4

General-purpose simulator for molecular nanotechnology.

Documentation: [philipturner.github.io/MM4](https://philipturner.github.io/MM4)

### Simulation

Levels of theory:
- Molecular mechanics
- Rigid body mechanics

### Supported Atoms

Rigid body mechanics supports any bonding topology. The following restrictions apply molecular mechanics.

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
| Au            | n/a | n/a | n/a           | n/a           |

The following are officially supported in the current release. Other atoms are unsupported.

| MM4 Atom Code | 6-ring | 5-ring | 4-ring | 3-ring |
| - | - | - | - | - |
| H             | 5   | n/a | n/a           | n/a           |
| C             | 1   | 123 | not supported | not supported |
| Si            | 19  | 19  | n/a           | n/a           |
| Au            | n/a | n/a | n/a           | n/a           |

Gold atoms should be anchors or simulated with rigid body mechanics. They only participate in nonbonded interactions. For example, MM4 can simulate the vdW interaction between graphene and gold. It can prevent a graphene sheet from physically colliding with a gold surface. The interaction force would be combined with a quantum mechanical simulator, such as GFN-FF, that handles pi-bonded carbons but not bulk metals.

### Supported Bonds

Rigid body mechanics supports any bonding topology. The following restrictions apply molecular mechanics.

| Element | H | C | N | O | F | Si | P | S | Ge | Au |
| ------- | - | - | - | - | - | - | - | - | - | - |
| H       |   | X |   |   |   | X |   |   | X |   |
| C       | X | X | O | O | O | O | O | O | O |   |
| N       |   | O |   |   |   |   |   |   |   |   |
| O       |   | O |   |   |   |   |   |   |   |   |
| F       |   | O |   |   |   |   |   |   |   |   |
| Si      | X | O |   |   |   | X |   |   |   |   |
| P       |   | O |   |   |   |   |   |   |   |   |
| S       |   | O |   |   |   |   |   |   |   |   |
| Ge      | X | O |   |   |   |   |   |   | X |   |
| Au      |   |   |   |   |   |   |   |   |   |   |

The following are officially supported in the current release. Other bonds are unsupported.

| Element | H | C | Si | Au |
| ------- | - | - | - | - |
| H       |   | X | X |   |
| C       | X | X |   |   |
| Si      | x |   | X |   |
| Au      |   |   |   |   |

Key:
- X = nonpolar sigma bond
- O = polar sigma bond

### Forced Motions

|         | Velocity             | Force           |
| ------- | -------------------- | --------------- |
| Linear  | anchor with velocity | external force  |
| Angular | flywheel             | linear to rotary motion converter |

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

### Releases

Current version: v1.0.0-beta0

v1.0.0
- Accurate simulation of 5-ring carbons
- Anchors
- External forces
- Support for hydrocarbons and hydrosilicons

Future versions:
- High-precision energy measurements
- Rigid body mechanics
- Support for non-carbon elements

## Tips

List:
- Compile this package in Swift release mode. Vectorized code is known to be extremely slow in debug mode. However, it may not be a bottleneck for small systems (under 1000 atoms).
