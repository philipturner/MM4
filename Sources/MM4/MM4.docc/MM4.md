# ``MM4``

Molecular Mechanics force field, version 4. The simulator used to create _Nanosystems (1992)_, but updated with modern ab initio parameters.

## Overview

The simulator supports the atoms enumerated by <doc:MM4AtomCode>, with some
restrictions on permitted bond types. Elements except C, Si, Ge must have
every covalent bond shared with a carbon atom. In addition, no dissimilar
non-H/C atoms may be separated by a bond once removed. Within these
restrictions, there exists a large variety of permissible structures, more
than with exclusively carbon and hydrogen.

All force terms from MM3, and most terms from MM4, are included, except:
- Electronegativity effect correction to equilibrium bond angle
- Electronegativity effect correction to bond stretching stiffness
- Torsion-torsion interaction
- Bend-torsion-bend interaction in bulk diamond

The simulator uses 32-bit single precision and an NVE ensemble. Single 
precision is often associated with poor energy conservation. However, the 
reason for using NVE is not energy conservation. The default time step is
so large that integration error dwarfs rounding error. Switching to 64-bit
double precision would not improve energy conservation.

> Note: Accuracy (whether an archer hits the right target) is orthogonal
  to precision (whether the arrows hit the same place every time).
  Conservation of energy only affects the number of significant figures of
  measured energy (precision). It does not make system dynamics more
  accurate, except in the case where energy systematically drifts upward.

Single precision is an implementation choice to make the simulator more
compatible with different GPU architectures. Most vendors have no or few
FP64 units (unlike CPU, where FP64 is widely supported). Mixed precision
can sometimes run the bulk of computations in FP32, with only the energy
computations in FP64. However, even mixed precision has implementation
issues, being difficult to run on GPUs with zero FP64 units.
