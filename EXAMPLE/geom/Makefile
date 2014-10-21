FLAGS=$(FC_FLAGS) -I$(QUEST_DIR)/SRC

.SUFFIXES : .mod .o .F90

test:   $(QUEST_DIR)/$(DQMCLIB) ggeom.F90
	$(FC) $(FLAGS) -o ggeom ggeom.F90 $(QUEST_DIR)/$(DQMCLIB) $(LIB)

clean:
	rm -f ggeom
