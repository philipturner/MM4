# MM4

Molecular Mechanics force field, version 4. The simulator used to create _Nanosystems (1992)_, but updated with modern ab initio parameters.

Documentation: [philipturner.github.io/MM4](https://philipturner.github.io/MM4)

### Atoms

Officially supported:

| Element | Ring Types |
| ------- | ---------- |
| H             | n/a  |
| C             | 5, 6 |
| Si            | 5, 6 |
| P (trivalent) | 5, 6 |
| S             | 5, 6 |
| Ge            | 5, 6 |

Experimental:

| Element | Ring Types |
| ------- | ---------- |
| N (trivalent) | 5, 6 |
| O             | 5, 6 |
| F             | n/a  |

### Bonds

Officially supported:

| Element | H | C |
| ------- | - | - |
| C       | X | X |
| Si      | X | X |
| P       |   | X |
| S       |   | X |
| Ge      | X | X |

Experimental:

| Element | H | C |
| ------- | - | - |
| N       |   | X |
| O       |   | X |
| F       |   | X |

### Forces

Officially supported:
- bend
- external
- nonbonded
  - van der Waals force
  - overlap repulsion
  - electrostatic force
- stretch
- stretch-bend

Experimental:
- bend-bend
- stretch-stretch
- torsion
- torsion-bend
- torsion-stretch

### Levels of Theory

MM4 offloads all molecular dynamics calculations to OpenMM. Rigid body dynamics must be integrated on the CPU by the library user. Any communication between CPU and GPU causes a latency bottleneck. This bottleneck manifests as a large $O(1)$ term in the polynomial for algorithmic complexity.

|  | Stable Time Step | Minimum Latency/Step | Maximum ns/day | Scaling | Force Computation | Integration |
| :-----------------: | :--------: | :--------: | :-----: | :-: | :-: | :-: |
| Molecular Dynamics (w/o cutoff)        | 4.35 fs |  100 μs | 3200 ns/day | $O(n^2)$ | GPU | GPU |
| Molecular Dynamics (w/o neighbor list) | 4.35 fs |  200 μs | 1600 ns/day | $O(n^2)$ | GPU | GPU |
| Molecular Dynamics                     | 4.35 fs |  700 μs | 500 ns/day  | $O(n)$ |  GPU | GPU |
| Rigid Body Dynamics                    | 80 fs   | 1500 μs | 4000 ns/day | $O(n)$ |  GPU | CPU |

For large atom counts and lower-end hardware, the $O(n)$ term will dominate. It is about the compute cost of biomolecular force fields (e.g. AMBER). However, GPU hardware allows several thousand calculations to occur each clock cycle. This makes MM4 much faster than CPU-based simulators (GROMACS, LAMMPS) running the same type of force field.

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
