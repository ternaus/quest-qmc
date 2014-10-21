FLAGS=$(FC_FLAGS) $(PRG_FLAGS) -I$(QUEST_DIR)/SRC

.SUFFIXES : .mod .o .F90

mpi:    $(QUEST_DIR)/$(DQMCLIB) verify.F90
	$(FC) $(FLAGS) -o verify verify.F90 $(QUEST_DIR)/$(DQMCLIB) $(LIB) 

clean:
	rm -f verify

