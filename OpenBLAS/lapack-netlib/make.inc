FORTRAN     = ifort
OPTS        = -O2 -openmp -fPIC
POPTS       = -O2 -openmp -fPIC -pg -fno-inline
NOOPT       = -O2 -openmp -fPIC -O0
PNOOPT      = -O2 -openmp -fPIC -pg -fno-inline -O0
LOADOPTS    = -O2 -openmp -fPIC 
CC          = clang
override CFLAGS      = -O2 -Wall -m64 -DF_INTERFACE_INTEL -fPIC -DNO_WARMUP -DMAX_CPU_NUMBER=1 -DASMNAME=_ -DASMFNAME=__ -DNAME=_ -DCNAME= -DCHAR_NAME="_" -DCHAR_CNAME="" -DNO_AFFINITY -I. -DHAVE_LAPACK_CONFIG_H
ARCH        = ar
ARCHFLAGS   = -ru
RANLIB      = ranlib
LAPACKLIB   = ../libopenblas_haswell-r0.2.10.rc2.a
TMGLIB      = ../libopenblas_haswell-r0.2.10.rc2.a
BLASLIB     = ../../../libopenblas_haswell-r0.2.10.rc2.a
LAPACKELIB  = ../libopenblas_haswell-r0.2.10.rc2.a
LAPACKLIB_P = ../libopenblas_haswell-r0.2.10.rc2_p.a
SUFFIX      = o
PSUFFIX     = po
CEXTRALIB   = 
TIMER       = NONE
LOADER      = ifort
SHELL = /bin/sh
PLAT = _LINUX
DRVOPTS  = $(OPTS)
ARCHFLAGS= -ru
#RANLIB   = ranlib

