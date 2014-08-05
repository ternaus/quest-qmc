.SUFFIXES : .mod .o .F90

test:    $(QUEST_DIR)/$(DQMCLIB) test.F90
	$(FC) $(FLAGS) -I$(QUEST_DIR)/SRC -o test  test.F90 $(QUEST_DIR)/$(DQMCLIB) $(LIB) 

clean:
	rm -f test

