Different libraries, different compilers, Which combination give fastest code?

# Introduction #

This package can be used with different libraries, different compilers and compiler settings. Below, I will try to give some numbers that will tell us what is better, what is worse.

# Main result #

| **CPU** | **Compiler** | **Library** | **Matrix multiplication** | **small** | **medium** | **large** |
|:--------|:-------------|:------------|:--------------------------|:----------|:-----------|:----------|
| i7-3520M @ 2.90GHz |gcc 4.8.2 | BLAS/LAPACK 3.5.0 | -DDQMC\_ASQRD | 103 | 293 | 4334|
| i7-3520M @ 2.90GHz |gcc 4.8.2 | OpenBLAS 2.10 RC2 | -DDQMC\_ASQRD |  38 | 96 | 1190 |
| i7-3520M @ 2.90GHz|icc 14.0.2 | BLAS/LAPACK 3.5.0 | -DDQMC\_ASQRD | 53  | 159 | 2609 |
| i7-3520M @ 2.90GHz |icc 14.0.2 | OpenBLAS 2.10 RC2 | -DDQMC\_ASQRD |  36 |  93 |  |
| i7-3520M @ 2.90GHz|icc 14.0.2 | Intel MKL 11.1.0.080 | -DDQMC\_ASQRD |35 | 84 | 1189|
| i7-3520M @ 2.90GHz|pg 14.6-0 | BLAS/LAPACK 3.5.0 | -DDQMC\_ASQRD | 83 | 248 |  |
| Intel Pentium 4 @ 3.20GHz| icc 14.0.1 | BLAS/LAPACK 3.5.0 | -DDQMC\_ASQRD | 153|456 |7120 |
| Intel Pentium 4 @ 3.20GHz| icc 14.0.1 | Intel MKL 11.1.0.080 | -DDQMC\_ASQRD | 108| 301| 4286|
| Intel Pentium 4 @ 3.20GHz| gcc 4.4.6 | BLAS/LAPACK 3.5.0 | -DDQMC\_ASQRD | 287| 858| 13348 |
| Intel Pentium 4 @ 3.20GHz| gcc 4.4.6 | OpenBLAS 2.10 RC2 | -DDQMC\_ASQRD | 112 |303 | 4153 |
| i5-3470 @ 3.20GHz| gcc  4.8.2 | BLAS/LAPACK 3.5.0 | -DDQMC\_ASQRD |80 | 233|3535 |
| i5-3470 @ 3.20GHz| gcc  4.8.2 + optimizations | BLAS/LAPACK 3.5.0 | -DDQMC\_ASQRD |65 |194 | 3214 |
| i5-3470 @ 3.20GHz| gcc  4.8.2 | Intel MKL 11.1.0.080 | -DDQMC\_ASQRD |29 |77 |923 |
| Intel Xeon E5-1620 @ 3.60GHz| gcc  4.6.3 | BLAS/LAPACK 3.4.2 | -DDQMC\_ASQRD |94 | 265 |  |
| Intel Xeon E5-1620 @ 3.60GHz| gcc  4.6.3 + optimizations | BLAS/LAPACK 3.4.2 | -DDQMC\_ASQRD |68 | 204 |  |
| Intel Xeon E5-1620 @ 3.60GHz| gcc  4.6.3 + optimizations | OpenBLAS 2.10 RC | -DDQMC\_ASQRD | 37 | 88 |  |
| Intel Xeon E5-1620 @ 3.60GHz| icc 14.0.1 | Intel MKL 11.1.0.080 | -DDQMC\_ASQRD | 29 | 77 | 905 |

# Details. #

## CPU ##
cpu name and frequency was obtained using

less /proc/cpuinfo

## Compilers. ##
Fortran compilers:
  * gfortran
  * Intel fortran
  * Portland Group compiler

C++ compilers:
  * g++
  * Intel C++
  * pgc++


To measure the speed I will use programs from the EXAMPLE/test folder.
i.e.
  1. ./test large.in
  1. ./test median.in
  1. ./test small.in

---

gcc compiler flags:

FC\_FLAGS = -m64 -O3 -funroll-loops

CXX\_FLAGS = -m64 -O3 -funroll-loops

---

gcc + optimizations compiler flags:

FC\_FLAGS = -m64 -O3 -funroll-loops -msse2 -msse3 -msse4 -msse4.1 -mavx

CXX\_FLAGS = -m64 -O3 -funroll-loops -msse2 -msse3 -msse4 -msse4.1 -mavx

---

Intel compiler flags:

FC\_FLAGS =  -m64 -O3 -unroll

CXX\_FLAGS = -m64 -O3 -unroll

---

Portland Group compiler flags:

FC\_FLAGS   = -fast –Minform=inform -Mipa=fast,inline

CXX\_FLAGS   = -fast –Minform=inform -Mipa=fast,inline

---
