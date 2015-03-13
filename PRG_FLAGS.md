# Introduction. #
There are different compiler flags that are used in the compilation of the latest version of the QUEST.

# Description. #
  1. Matrix multiplication flags:
    * -DDQMC\_ASQRD: Method of the stratification (computation of G), which does most of the computations with regular QR, instead of QR with pivoting.
    * -DDQMC\_CKB: Use checkboard method for generating matrix B
> The above flags should not be used simultaneously.
  1. -D\_PREFIX - uses prefix matrix multiplication. (**Not tested**)
  1. -D\_OpenBC - open boundary conditions (**Not tested**)
  1. -DDQMC\_PROFILE - Flag used for profiling. ( **Should be removed in the next release.**) Using [Gprof](http://en.wikipedia.org/wiki/Gprof) or [Valgrind](http://en.wikipedia.org/wiki/Valgrind) looks like a much better idea.)
  1. -DDQMC\_CUDA - This flag is used if we want to use GPU for come parts of the computations. This flag is mandatory if [MAGMA](http://icl.utk.edu/magma/) is used instead of LAPACK/BLAS. (**Not tested**)
  1. -D\_SXX - **No clue**
  1. -D\_QMC\_MPI - Multicore support. (**Not tested**).