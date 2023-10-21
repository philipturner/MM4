# MM4

General-purpose simulator for molecular nanotechnology.

Swift API: [philipturner.github.io/MM4](https://philipturner.github.io/MM4)

C API: [MM4.h](./Bindings/C/MM4.h)

Python API: in progress

## Parameters

### Supported Atoms

| MM4 Atom Code | 6-ring | 5-ring | 4-ring | 3-ring |
| - | - | - | - | - |
| H               | 5   | n/a | n/a           | n/a           |
| C (sp3)         | 1   | 123 | not supported | not supported |
| N (trivalent)   | 8   | 8   | not supported | not supported |
| O               | 6   | 6   | not supported | not supported |
| F               | 11  | n/a | n/a           | n/a           |
| Si              | 19  | 19  | not supported | not supported |
| P (trivalent)   | 25  | 25  | not supported | not supported |
| S               | 15  | 15  | not supported | not supported |
| Ge              | 31  | 31  | not supported | not supported |

### Supported Bonds

Key:
- X = nonpolar covalent bond
- O = polar covalent bond
- blank = not supported

| Element | H | C<sub>sp3</sub> | N | O | F | Si | P | S | Ge |
| --------------- | - | - | - | - | - | - | - | - | - |
| H               |   | X |   |   |   | X |   |   | X |
| C<sub>sp3</sub> | X | X | O | O | O | O | O | O | O |
| N               |   | O |   |   |   |   |   |   |   |
| O               |   | O |   |   |   |   |   |   |   |
| F               |   | O |   |   |   |   |   |   |   |
| Si              | X | O |   |   |   | X |   |   |   |
| P               |   | O |   |   |   |   |   |   |   |
| S               |   | O |   |   |   |   |   |   |   |
| Ge              | X | O |   |   |   |   |   |   | X |

## Roadmap

First development round (current):
- Internal workings figured out, Swift code written
- C API implemented

Second development round (in ~weeks):
- Inner workings tested
- Python API implemented
- Optimize the CPU-side cost of various operations

Third development round (in ~months):
- Optimize every force to the physical limits, estimate no more than 50% speedup
- Any necessary support for multi-GPU systems
