## Why do we need to speedup TDM ##

In profile of "ggeom", the subroutine _DQMC\_TDM1\_Meas_ takes the major computational time.

E.g. Try the input arguments
| ns | L | north | nwarm | npass |
|:---|:--|:------|:------|:------|
| 20 | 100 | 10 | 100 | 200 |

with processors @2.97GHz.

The CPU\_time in seconds in a single core sequential execution in MKL:

| Total time | TDM | Others |
|:-----------|:----|:-------|
| 21153 | **14883** | 6270 |

With optimized multi-threaded BLAS/LAPACK
the performance can be improved in a multi-core system by just naively setup the multiple threads.
```
export OMP_NUM_THREADS=n
```
where n is the number of threads.
| n | Total time | TDM | Others |
|:--|:-----------|:----|:-------|
| 6 | 11358 | **9521** | 1837 |
| 12 | 10140 | **8670** | 1470 |

**Problem:** The speedup of TDM part is limited with the growth of number of threads.

## The performance after use Parallel DQMC\_TDM1\_Meas ##
There is an optional way to compute TDMs by a OpenMP parallel subroutine.

With the same input arguments, we can obtain a better performance on TDM part and further reduce the total time for the whole simulation.
| n | Total time | TDM | Others |
|:--|:-----------|:----|:-------|
| 6 | 4588 | **2772** | 1816 |
| 12 | 3126 | **1671** | 1455 |

## How to use the parallel TDM ##
The parallel TDM is implemented by OpenMP. You can specify the number of threads through the command line. For example,

  * 1. max threads:
```
./ggeom in -p
```

  * 2. 1 <= n <= max, n is the number of threads you specify:
```
./ggeom in -p n
```

  * 3. original single thread
```
./ggeom in
```

You can check the number of cores in your system by
```
lscpu
```
Note that in most of the cases, the max number of threads equals the number of cores in your system.

## Problem linking with OpenBLAS ##
The current code works well with MKL.
While linking with OpenBLAS, it occurs an error as:
'BLAS : Program is Terminated. Because you tried to allocate too many memory regions.'
which needs to be fixed.