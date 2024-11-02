############################################################################
#  QUEST Makefile - Using Apple Accelerate Framework
############################################################################
QUEST_DIR = $(shell pwd)
export QUEST_DIR

# Compiler settings
FC        = gfortran
CXX       = g++
FC_FLAGS  = -fopenmp -Wall -O3 -funroll-loops -fexternal-blas
CXX_FLAGS = -Wall -O3 -funroll-loops

# Libraries
CXXLIB    = -lstdc++
LAPACKLIB = -framework Accelerate

# Archiver
ARCH    = ar
ARFLAG  = cr
RANLIB  = ranlib

# Main library
DQMCLIB = libdqmc.a

# Required libraries for driver routines
LIB = $(CXXLIB) $(LAPACKLIB)

.PHONY: all libdqmc example clean

# Default target
all: libdqmc example

# Build main library
libdqmc:
	$(MAKE) -C SRC

example: libdqmc
	$(MAKE) -C EXAMPLE

# Cleanup
clean:
	$(MAKE) -C $(QUEST_DIR)/SRC clean
	$(MAKE) -C $(QUEST_DIR)/EXAMPLE clean
	$(RM) $(QUEST_DIR)/$(DQMCLIB)

# Help target
help:
	@echo "Available targets:"
	@echo "  all     - Build everything (default)"
	@echo "  libdqmc - Build DQMC library"
	@echo "  example - Build examples"
	@echo "  clean   - Clean all built files"
	@echo "  help    - Show this help message"

# Export variables for sub-makefiles
export FC FC_FLAGS CXX CXX_FLAGS ARCH ARFLAG RANLIB DQMCLIB LIB
