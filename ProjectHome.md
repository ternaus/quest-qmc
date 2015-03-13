QUantum Electron Simulation Toolbox (QUEST) is a Fortran 90/95 package that implements the Determinant Quantum Monte Carlo (DQMC) method for quantum electron simulations. The original versions of DQMC programs, developed by the condensed matter theory group at UCSB including R. L Sugar, D. J. Scalapino, S. R. White, and E. Y. Loh, and maintained by R. Scalettar, have been extensively used to study magnetism, metal-insulator transitions, and superconductivity by the Hubbard model.


[Download and Install](https://code.google.com/p/quest-qmc/wiki/Download)

[How to start the simulation](https://code.google.com/p/quest-qmc/wiki/HowToUse)

[List of publications](https://code.google.com/p/quest-qmc/wiki/publications)

---

Stable version 1.4.9 released.

[quest-1.4.9.tgz](https://drive.google.com/file/d/0ByATqxQWgjs9cW9DMkxLTTE4MlE/edit?usp=sharing)

  * BLAS/LAPACK was replaced by OpenBLAS
  * python templates for several lattices added.

---

People who worked on this release:
  * Chia-Chen Chang
  * Vladimir Iglovikov


---

Know issues:

  * Gives wrong values for the specific heat.
  * Portland compiler is not supported.
  * GPU is not tested.
  * MPI/OpenMP is not tested.


---


---

Stable version 1.4.4 released.

[quest-1.4.4.tgz](https://drive.google.com/file/d/0ByATqxQWgjs9VnVqOW1uZ0I2Y1k/edit?usp=sharing)

  * Definition of the kinetic and total energy changed.
    1. Kinetic energy now includes all quadratic terms of the Hamiltonian.
    1. Total energy now matches most intuitive definition.
  * New measurements added:
    1. Potential energy
    1. Hopping energy
    1. Double occupancy
    1. Square of the magnetization.
  * Symmetry factor added to the output.
  * BLAS/LAPACK routines updated to the version 3.5.0
  * Added documentation about specific heat derivation and implementation in the code.
  * Bug fix: U/2 term added in the specific heat
  * Bug fix: added missing terms in density-density correlation function.
  * Bug fix: wrong allocation of the array size in the dqmc\_gtau.F90 fixed.

---

People who worked on this release:
  * Simone Chiesa
  * Chia-Chen Chang
  * Vladimir Iglovikov

---

Know issues:

  * Gives wrong values for the specific heat.
  * Portland compiler is not supported.
  * GPU is not tested.
  * MPI/OpenMP is not tested.
  * Checkerboard decomposition for the matrix multiplication is not tested.
