############################################################################
#  QUEST Makefile - Simplified for MacOS
############################################################################
QUEST_DIR = $(shell pwd)

# Compiler settings
FC        = gfortran
CXX       = g++
FC_FLAGS  = -fopenmp -Wall -O3 -funroll-loops
CXX_FLAGS = -Wall -O3 -funroll-loops

# Libraries
CXXLIB    = -lstdc++
libOpenBLAS = $(QUEST_DIR)/OpenBLAS/libopenblas.a
LAPACKLIB = $(libOpenBLAS)

# Archiver
ARCH    = ar
ARFLAG  = cr
RANLIB  = ranlib

# Main library
DQMCLIB = libdqmc.a

# Required libraries for driver routines
LIB = $(CXXLIB) $(LAPACKLIB)

.PHONY: all libdqmc example lapack clean

# Default target
all: libdqmc example

# Check if OpenBLAS exists before building main targets
libdqmc: check_openblas
	$(MAKE) -C SRC

example: libdqmc
	$(MAKE) -C EXAMPLE

# OpenBLAS handling
lapack:
	$(MAKE) -C $(QUEST_DIR)/OpenBLAS

check_openblas:
	@if [ ! -f $(libOpenBLAS) ]; then \
		echo "OpenBLAS not found. Building it first..." ; \
		$(MAKE) lapack ; \
	fi

# Cleanup
clean:
	$(MAKE) -C $(QUEST_DIR)/OpenBLAS clean
	$(MAKE) -C $(QUEST_DIR)/SRC clean
	$(MAKE) -C $(QUEST_DIR)/EXAMPLE clean
	$(RM) $(QUEST_DIR)/$(DQMCLIB)

# Help target
help:
	@echo "Available targets:"
	@echo "  all     - Build everything (default)"
	@echo "  libdqmc - Build DQMC library"
	@echo "  example - Build examples"
	@echo "  lapack  - Build OpenBLAS"
	@echo "  clean   - Clean all built files"
	@echo "  help    - Show this help message"

# Export variables for sub-makefiles
export FC FC_FLAGS CXX CXX_FLAGS ARCH ARFLAG RANLIB DQMCLIB LIB
