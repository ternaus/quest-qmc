QUEST_DIR = ..
include $(QUEST_DIR)/make.inc

OBJ =   dlamch.o ilaenv.o iparmq.o ieeeck.o dlarnv.o dlaruv.o \
	dsyev.o dlansy.o dsytrd.o dorgtr.o dsteqr.o dlatrd.o dsytd2.o dorgql.o \
	dlasr.o dorg2l.o dgesv.o \
	dgeqp3.o dlaqps.o dlaqp2.o \
	dorgqr.o dorg2r.o \
	dormqr.o dorm2r.o \
	dgetrf.o dgetf2.o \
	dgetri.o dtrtri.o dtrti2.o \
	dgetrs.o \
	dgejsv.o dgelqf.o dgesvj.o dlacpy.o dlaset.o dlaswp.o dormlq.o dpocon.o \
	dgsvj0.o dgsvj1.o dorml2.o dlatrs.o dlacn2.o drscl.o dlabad.o \
	dgelq2.o \
	dgeqrf.o dgeqr2.o dlarfb.o dgeqr2.o dlarft.o dlarfg.o dlarf.o iladlr.o \
	iladlc.o \
	zheev.o zlanhe.o zlascl.o zhetrd.o dsterf.o zungtr.o zsteqr.o \
	dlae2.o dlaev2.o dlanst.o dlapy2.o dlartg.o dlascl.o dlasrt.o zlaset.o \
	zlasr.o disnan.o dlassq.o zungql.o zungqr.o zhetd2.o zlatrd.o zlassq.o \
	zlacgv.o zlarfg.o zlarfb.o zung2r.o zung2l.o dlaisnan.o zlarf.o ilazlr.o \
	ztrtri.o zgetf2.o zlaswp.o ztrti2.o \
	ilazlc.o zladiv.o dlapy3.o zlarft.o dladiv.o zgetri.o zgetrf.o
      
lib: 	$(OBJ)
	$(ARCH) $(ARFLAG) ../liblapack.a $(OBJ)
	$(RANLIB) ../liblapack.a

# Automatic copy from LAPACK source
#.PRECIOUS: %.f
#%.f:
#	cp /home/aetdom/tmp/lapack-3.3.1/SRC/$@ .

# DLAMCH must be compiled without optimizations
dlamch.o: dlamch.f
	$(FC) $(NOOPT_FLAGS) -c $<

%.o: 	%.f
	$(FC) $(FC_FLAGS) -c $<

clean:
	rm -f ../liblapack.a *.o

