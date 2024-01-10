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

 All major platforms (Mac, Linux, Windows) and all major GPU architectures
 (Apple, AMD, Intel, Nvidia) are supported. This includes multi-GPU systems
 with up to around a dozen GPUs (e.g. a single supercomputer node). This is
 because it delegates all of the computation to OpenMM, an existing,
 well-tested framework for easily porting molecular simulation algorithms to
 multiple vendors. OpenMM avoids any O(n^3) algorithms and provides the
 flexibility to opt out of O(n^2) algorithms. For example, the
 [minimizer](<doc:MM4ForceField/minimize(tolerance:maxIterations:)>) 
 uses limited-memory BFGS with O(n) complexity.

 > Note: Throughout the entire software stack, from setup to minimization to
   simulation, this simulator is O(n). Linear scaling exists at every system
   size and material composition. This is a rarity among simulation
   algorithms, most of which are inherently O(n^2) or
   [employ sparse O(n^2) matrix factorizations](https://xtb-docs.readthedocs.io/en/latest/gfnff.html).
   Molecular mechanics is one of the few techniques that can simulate
   million-atom nanosystems on affordable hardware.

 This is an NVE simulator, using 32-bit single precision and a massive
 timestep. Energy fluctuates, but shouldn't systematically drift upward or
 downward (unless the timestep exceeds C-H stretching resonance frequency). 
 The default timestep (4 fs) is so massive that integration error dwarfs 
 numerical error. Even if the forcefield used 64-bit double or mixed precision,
 that would not improve its ability to conserve energy.

 > Note: Accuracy (whether an archer hits the right target) is orthogonal
   to precision (whether the arrows hit the same place every time).
   Conservation of energy only affects the number of significant figures of
   measured energy (precision). It does not make system dynamics more
   accurate, except in the case where energy systematically drifts upward.

 The integrator uses a multiple time-stepping (MTS) scheme. Cheaper bonded
 forces, such as bond-stretch and bond-bend, are only stable at ~2 fs
 without constraints. Expensive forces like torsions, nonbonded, and
 electrostatic can execute at double the timestep. The value you enter for
 [`timeStep`](<doc:MM4ForceField/timeStep>)
 specifies the execution rate of expensive forces. Always
 assume the C-H stretching forces execute at half the specified timestep.
 For example, in the note below, bond stretching forces don't actually
 execute at the quoted '2 fs'.

 > Note: To maximize the simulation speed, hydrogen mass repartitioning (HMR)
   is enabled by default. This method makes hydrogens heaver and makes
   non-hydrogen atoms lighter, to decrease the C-H/Si-H/Ge-H stretching
   frequency. You can, and should, disable HMR for energy-conserving
   simulations where timestep falls below 2 fs. This can be done by setting
 <doc:MM4ParametersDescriptor/hydrogenMassScale> to `1`.
   Note that changing the repartitioning will shift the center of mass,
   which is where bulk angular momentum is applied.

 Single precision is an implementation choice to make the simulator more
 compatible with different GPU architectures. Most vendors have no or few
 FP64 units (unlike CPU, where FP64 is widely supported). Mixed precision
 can sometimes run the bulk of computations in FP32, with only the energy
 computations in FP64. However, even mixed precision has implementation
 issues, being difficult to run on GPUs with zero FP64 units.

 The best method to maximize energy conservation is a timestep of 1 fs.
 This requires 4x as much computation for the same trajectory length,
 compared to the default (4 fs). It executes expensive forces at 1 fs and
 cheap forces at 0.5 fs. Integration error scales O(h^2), so this should be
 around 16x more precise. Any smaller timestep, and rounding error from
 FP32 will overtake integration error.
