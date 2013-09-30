QUEST_DIR = ..
include $(QUEST_DIR)/make.inc

OBJ = dasum.o daxpy.o dcabs1.o dcopy.o ddot.o dgbmv.o dgemm.o dgemv.o dger.o \
      dnrm2.o drot.o drotg.o drotm.o drotmg.o dsbmv.o dscal.o dsdot.o dspmv.o \
      dspr2.o dspr.o dswap.o dsymm.o dsymv.o dsyr2.o dsyr2k.o dsyr.o dsyrk.o \
      dtbmv.o dtbsv.o dtpmv.o dtpsv.o dtrmm.o dtrmv.o dtrsm.o dtrsv.o dzasum.o \
      dznrm2.o xerbla.o zaxpy.o zcopy.o zdotc.o zdotu.o zdrot.o zdscal.o zgbmv.o \
      zgemm.o zgemv.o zgerc.o zgeru.o zhbmv.o zhemm.o zhemv.o zher2.o zher2k.o \
      zher.o zherk.o zhpmv.o zhpr2.o zhpr.o zrotg.o zscal.o zswap.o zsymm.o \
      zsyr2k.o zsyrk.o ztbmv.o ztbsv.o ztpmv.o ztpsv.o ztrmm.o ztrmv.o ztrsm.o \
      ztrsv.o lsame.o idamax.o izamax.o

lib: 	$(OBJ)
	$(ARCH) $(ARFLAG) ../libblas.a $(OBJ)
	$(RANLIB) ../libblas.a

%.o: 	%.f
	$(FC) $(FC_FLAGS) -c $<

clean:
	rm -f ../libblas.a *.o *genmod*
	
