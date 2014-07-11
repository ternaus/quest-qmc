QUEST_DIR = .

include make.inc

#all: example_mkl

#all : example_

all : example_OpenBLAS

example_OpenBLAS : libopenblas libdqmc
	(cd EXAMPLE; $(MAKE))

example_: liblapack libblas libdqmc
	(cd EXAMPLE; $(MAKE))

example_mkl: libdqmc
	$(MAKE) -C EXAMPLE	

libblas:
	(cd BLAS; $(MAKE))

liblapack:
	(cd LAPACK; $(MAKE))

libdqmc:
	$(MAKE) -C SRC

libopenblas:
	(cd OpenBLAS; $(MAKE))

clean:
	(cd BLAS; $(MAKE) clean)
	(cd LAPACK; $(MAKE) clean)
	(cd SRC; $(MAKE) clean)
	(cd EXAMPLE; $(MAKE) clean)
	(rm -f $(DQMCLIB))
	(cd OpenBLAS; $(MAKE) clean)

