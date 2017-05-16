FC=ifort
FC=gfortran
#FFLAG1= -mcmodel=medium 
#-ffixed-line-length-none
FFLAG=-fbounds-check
executable1=stack
executable2=Stack_all
executable3=stack_all_test
objects1=stack.o hilbert.o zfour.o S_trans.o diftc_fast.o clogc.o sacio.o
objects2=stack_all.o hilbert.o zfour.o S_trans.o diftc_fast.o clogc.o sacio.o do_pws.o do_tfpws.o filter.o taperf.o
objects3=stack_all_test.o hilbert.o zfour.o S_trans.o diftc_fast.o clogc.o sacio.o do_pws.o do_tfpws.o
objects4=lstack.o sacio.o
objects5=do_S_trans.o hilbert.o zfour.o S_transform.o diftc_fast.o clogc.o sacio.o
objects6=do_S_transform.f90 hilbert.o zfour.o S_transform.o diftc_fast.o clogc.o sacio.o
objects7=do_S_transform_2016_11_17.f90 hilbert.o zfour.o S_transform.o diftc_fast.o clogc.o sacio.o
all: sacio.mod $(executable1) $(executable2) $(executable3) lstack do_S_trans do_S_transform do_S_transform_2016_11_17
sacio.mod:sacio.f90
	$(FC) $^ -c
%.o:%.f90
	$(FC) $(FFLAG) $< -c
$(executable1):$(objects1)
	$(FC) $(FFLAG1) $^ -o $@
$(executable2):$(objects2)
	$(FC) $(FFLAG1) $^ -o $@
$(executable3):$(objects3)
	$(FC) $(FFLAG1) $^ -o $@
lstack:$(objects4)
	$(FC) $(FFLAG1) $^ -o $@
do_S_trans:$(objects5)
	$(FC) $(FFLAG) $^ -o $@
do_S_transform:$(objects6)
	$(FC) $(FFLAG) $^ -o $@
do_S_transform_2016_11_17:$(objects7)
	$(FC) $(FFLAG) $^ -o $@
install:
	cp $(executable1) $(executable2) $(executable3) lstack do_S_trans do_S_transform do_S_transform_2016_11_17 ../bin
uninstall:
	-rm ~/bin/$(excecutable)
clean:
	-rm *.o *.mod
