

# Downloading and Installation. #

The last version of QUEST is available as [quest-1.4.9.tgz](https://drive.google.com/file/d/0ByATqxQWgjs9cW9DMkxLTTE4MlE/edit?usp=sharing). This tar file can be extracted by
```
   tar -xzf quest-1.4.9.tgz
```
which will create a directory
```
   quest-qmc
```
To compile the library in most Linux systems with the GNU compilers, just run the command:
```
   make lapack
   make
```
A test program is available in the EXAMPLES directory, including several input files (.in extension). To do a quick check of your QUEST installation, first switch to the EXAMPLES/test directory and then issue the command:
```
   ./test small.in
```
This will generate the output file small.out.
# Multicore processors #

This tar file includes a basic version of the required BLAS and LAPACK libraries. To obtain better performance, QUEST must be linked with an optimized version of these libraries. In order to benefit from parallelism in multicore processors, the BLAS library (and optionally LAPACK) must be multithreaded (such as ATLAS or Intel's MKL). Usually, these libraries are provided by the vendor in high performance computer systems. Please consult the system documentation for details and edit accordingly the BLASLAPACK variable in the make.inc file.
# GPU acceleration #

QUEST also contains experimental support for GPU accelerators using the Nvidia CUDA development tools and the [MAGMA](http://icl.cs.utk.edu/magma/) library. This code is not active by default. To activate the GPU support, uncomment all the CUDA and MAGMA related variables in make.inc file and activate the flag DQMC\_CUDA in PRG\_FLAGS.

# Previous stable versions of the code. #
[Downloads](https://drive.google.com/folderview?id=0ByATqxQWgjs9SWpPMmFfSzk0Q0k&usp=sharing#list)

# Developers version of the code. #
The QUEST _mercurial_ repository is open for read-only for our users to be able to get the latest bug fixed.
```
hg clone https://code.google.com/p/quest-qmc/
```
# Changelog: #

## 1.4.9 ##
  * BLAS/LAPACK was replaced by OpenBLAS
  * python templates for several lattices added.
## 1.4.4 ##
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

## 1.4.0 ##
  * Number of bugfixes.
  * Last Intel compiler's support that was broken in previous versions was repaired.
  * QUEST does not print HSF variables into the file anymore. 5-15% of the speed gain.
  * Changed name of the file describing symmetries of the lattice.
  * Makefile examples for most commonly used compiler/libraries configurations with optimal compiler flags added.
  * Examples of the lattice geometries added.


## 1.3.0 ##

  * Support for time dependent measurement in the geom example.
  * Input file name specified in command line. The output file name is a input parameter.
  * Old non-working examples are removed. The generic geom example should be able to compute all these simulations.

## 1.2.0 ##

  * The ASQRD method is now the default stratification algorithm. An error will stop the simulation if parameter L is not a multiple of north.
  * More robust computation recycling. Active when nwrap = north and fixwrap = 1.
  * Full GPU implementation of the ASQRD method (requires MAGMA).

## 1.1.0 ##

  * New alternative algorithm for stratification (computation of G): ASQRD method with prepivoting.
  * This method is much faster because does most of the computations with regular QR instead of QR with pivoting.
Also, this implementation uses OpenMP for parallelism in multicore machines (a multithreaded LAPACK library such as Intel MKL is required).
However this method could be less stable than the old stratification procedure in extreme cases.
  * To activate the ASQRD method define the macro DQMC\_ASQRD at the end of the make.inc file.
The source code of this implementation is dqmc\_gfun.F90 and computeg.c++
  * Computation recycling: when nwrap = north a lot of computations can be saved from one Metropolis step to be used in the following steps.
This can be activated by the parameter fixwrap = 1 and L should be a multiple of nwrap. This feature is still experimental.
  * GPU implementation of matrix products using the CUBLAS library and custom CUDA code (multb.cu and multb.h). Depending on the GPU hardware this could be faster than using only the CPU for large lattices. Define the macro DQMC\_CUDA to activate this feature.
  * Built-in profiling support for the most relevant computational kernels, activated by the macro DQMC\_PROFILE.
  * BLAS/LAPACK and C++ routines: profile.h, profile.c++, blaslapack.h and blaslapack.c++
  * Fortran kernels: dqmc\_gfun.F90 and dqmc\_matb.F90
  * Added BLAS/LAPACK source files so QUEST is now self-contained.
However for improved performance an external optimized BLAS/LAPACK library is still required.