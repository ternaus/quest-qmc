QUEST_DIR = ..

include $(QUEST_DIR)/make.inc

.SUFFIXES : .mod .o .F90

FLAGS = $(FC_FLAGS) $(PRG_FLAGS)

%.o: 	%.F90
	$(FC) $(FLAGS) $(SPRNG_INC) -c $<

%.o: 	%.c++ blaslapack.h profile.h kernel_cpu.h kernel_gpu.h
	$(CXX) $(PRG_FLAGS) $(CXX_FLAGS) -c $<

%.o:	%.cu blaslapack.h profile.h kernel_gpu.h
	$(NVCC) $(PRG_FLAGS) $(CU_FLAGS) -c $<

SRC =	dqmc_mpi.F90 dqmc_util.F90 dqmc_cfg.F90 \
	dqmc_geom_param.F90 dqmc_latt.F90 dqmc_bonds.F90 dqmc_reclatt.F90 \
	dqmc_hamilt.F90 dqmc_symm.F90 \
	dqmc_wspace.F90 dqmc_struct.F90 dqmc_geom_wrap.F90 dqmc_matb.F90 dqmc_ckb.F90 \
	dqmc_seqb.F90 dqmc_gfun.F90 dqmc_phy0.F90 \
	dqmc_phy2.F90 dqmc_hubbard.F90 dqmc_gtau.F90 \
	dqmc_tdm1.F90 dqmc_tdm2.F90 dqmc_kbonds.F90 dqmc_2dperl.F90 #dqmc_simqmc.F90

ifdef NVCC
OBJ = 	$(SRC:.F90=.o) computeg.o blaslapack.o profile.o kernel_cpu.o kernel_gpu.o
else
OBJ = 	$(SRC:.F90=.o) computeg.o blaslapack.o profile.o kernel_cpu.o
endif

INC = 	$(SRC:.F90=.mod)

lib: 	$(OBJ)
	$(ARCH) $(ARFLAG) ../$(DQMCLIB) $(OBJ)
	$(RANLIB) ../$(DQMCLIB)

clean:
	rm -f *.o *.mod 
