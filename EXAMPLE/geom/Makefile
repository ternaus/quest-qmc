FLAGS=$(FC_FLAGS) -I../../SRC

.SUFFIXES : .mod .o .F90

# Default target
all: ggeom

# Build the geometry example
ggeom: ggeom.F90 ../../$(DQMCLIB)
	$(FC) $(FLAGS) -o ggeom ggeom.F90 ../../$(DQMCLIB) $(LIB)

clean:
	rm -f ggeom *.o *.mod
