# MM4

Molecular Mechanics force field, version 4. The simulator used to create _Nanosystems (1992)_, but updated with modern ab initio parameters.

Documentation: [philipturner.github.io/MM4](https://philipturner.github.io/MM4)

### Supported Parameters

Atoms:

| Element | Ring Types |
| ------- | ---------- |
| H             | n/a  |
| C             | 5, 6 |
| Si            | 5, 6 |
| P (trivalent) | 5, 6 |
| S             | 5, 6 |
| Ge            | 5, 6 |

Bonds:

| Element | H | C | Si | Ge |
| ------- | - | - | -- | -- |
| C       | X | X | X  | X  |
| Si      | X | X | X  |    |
| P       |   | X |    |    |
| S       |   | X |    |    |
| Ge      | X | X |    | X  |

Forces:
- bend
- external
- nonbonded
  - van der Waals force
  - overlap repulsion
  - electrostatic force
- stretch
- stretch-bend

### Levels of Theory

MM4 offloads all molecular dynamics calculations to OpenMM. Rigid body dynamics must be integrated on the CPU by the library user. Any communication between CPU and GPU causes a latency bottleneck. This bottleneck manifests as a large $O(1)$ term in the polynomial for algorithmic complexity.

|  | Stable Time Step | Minimum Latency/Step | Maximum ns/day | Scaling | Force Computation | Integration |
| :-----------------: | :--------: | :--------: | :-----: | :-: | :-: | :-: |
| Molecular Dynamics (no cutoff) | 2.5 fs |  175 μs | 1250 ns/day | $O(n^2)$ | GPU | GPU |
| Molecular Dynamics             | 2.5 fs |  375 μs | 550 ns/day  | $O(n)$ |  GPU | GPU |
| Rigid Body Dynamics            | 80 fs  | 1500 μs | 4000 ns/day | $O(n)$ |  GPU | CPU |

For large atom counts and lower-end hardware, the $O(n)$ term will dominate. This term is around the compute cost of biomolecular force fields (e.g. AMBER). However, GPU hardware allows several thousand calculations to occur each clock cycle. This fact makes MM4 much faster than CPU-based simulators (GROMACS, LAMMPS) running the same type of force field.

> NOTE: There is currently a massive bottleneck in the $O(n)$ term for nonbonded forces. It makes MM4 roughly 3x slower than it should be. The current performance of MM4 w/ GPU could equate to GROMACS w/ CPU, until the bottleneck is fixed.

### Units

MM4 and OpenMM use slightly different unit systems. MM4 adheres to the SI system: nanometer, yoctogram, picosecond. Units for force and energy are derived from dimensional analysis.

```
energy = m * v^2 = yg * (nm/ps)^2 = zJ
force = dU / dx = zJ / nm = pN
```

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

## Linker

MM4 needs to link against OpenMM, which can be problematic for the Swift compiler. On Unix platforms, the easiest method is by setting the environment variable, `OPENMM_LIBRARY_PATH`. This activates a piece of code in the package manifest for `swift-openmm`. On Windows (especially through the VSCode terminal), this does not work. An alternative is to modify your VSCode project's package manifest. Copy `OpenMM.dll` into the same folder as `Package.swift`, then explicitly link it in the manifest:

```swift
targets: [
  .executableTarget(
    name: "CLI",
    dependencies: [
       // Also add other necessary dependencies, such as 'HDL'.
      .product(name: "MM4", package: "MM4"),
    ],
    linkerSettings: [
      // In this example, the computer's user is 'username'. The
      // package manifest is located in the folder 'workspace'.
      .unsafeFlags(["-LC:/Users/username/Documents/.../workspace"]),
      .linkedLibrary("OpenMM"),
    ]),
],
```
