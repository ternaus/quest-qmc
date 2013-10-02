############################################################################
#
#  Program:         QUEST V1.4
#  Module:          Makefile
#  Purpose:         Top-level Makefile
#  Modified:        10/24/2012
#
############################################################################

QUEST_DIR = .

include make.inc

#all: libblas liblapack lib example_

all : lib example_
libblas:
	(cd BLAS; $(MAKE))

liblapack:
	(cd LAPACK; $(MAKE))

lib:
	(cd SRC; $(MAKE))

example_:
	(cd EXAMPLE; $(MAKE))

clean:
	(cd BLAS; $(MAKE) clean)
	(cd LAPACK; $(MAKE) clean)
	(cd SRC; $(MAKE) clean)
	(cd EXAMPLE; $(MAKE) clean)
	(rm -f $(DQMCLIB))

