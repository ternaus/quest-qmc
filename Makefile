QUEST_DIR = .

include make.inc

all: example_mkl

#all : example_

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

clean:
	(cd BLAS; $(MAKE) clean)
	(cd LAPACK; $(MAKE) clean)
	(cd SRC; $(MAKE) clean)
	(cd EXAMPLE; $(MAKE) clean)
	(rm -f $(DQMCLIB))

