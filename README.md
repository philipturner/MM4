# MM4

General-purpose simulator for molecular nanotechnology.

Documentation: [philipturner.github.io/MM4](https://philipturner.github.io/MM4)

### Simulation

Levels of theory:
- Molecular dynamics
- Rigid body dynamics
- Hybrid approach for deformable bodies

### Supported Atoms

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

### Supported Bonds

Key:
- X = nonpolar sigma bond
- O = polar sigma bond

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

### Forced Motions

|         | Force           | Velocity             |
| ------- | --------------- | -------------------- |
| Linear  | external force  | anchor with velocity |
| Angular | external torque | flywheel             |

### Releases

Current version: not released yet

v1.0.0
- Accurate simulation of 6-ring carbons
- Anchors
- External forces

v1.1.0
- Accurate simulation of 5-ring carbons
- Rigid body dynamics

Future versions:
- External torques
- Hybrid approach for deformable bodies
- Recognizing angular velocity of objects with multiple collinear anchors
- Support for non-carbon elements
