OBJ= coul90.o input.o numerov.o compute.o\
   main.o



BASE := $(shell expr $(CURDIR) : "\(.*\)/.*")
VERDATE := $(shell git log -1 --format=%cd  )
VERREV := $(shell git log -1 --pretty=format:"%h")
COMPDATE :=$(shell date)


LIB = -llapack -lblas
F90 = gfortran
FFLAGS = 

COMPILE_OPT1=    -fcheck=all -O3   
COMPILE_OPT2= -ffixed-line-length-0
COMPILE_OPT3= 
COMPILE_OPT4= -Wunused-variable
COMPILE_OPT5=
COMPILE_OPT6=
COMPILE_OPT7=
COMPILE_OPT8=

FOPT = $(COMPILE_OPT1) $(COMPILE_OPT2) $(COMPILE_OPT3) $(COMPILE_OPT4) \
       $(COMPILE_OPT5) $(COMPILE_OPT6) $(COMPILE_OPT7) $(COMPILE_OPT8) \
      -DVERDATE="'$(VERDATE)'" -DVERREV="'$(VERREV)'" \
      -DMAKEF90="'$(F90)'" -DMAKEF77="'$(FC)'" \
      -DMAKEFOPT1="'$(COMPILE_OPT1)'" \
      -DMAKEFOPT2="'$(COMPILE_OPT2)'" \
      -DMAKEFOPT3="'$(COMPILE_OPT3)'" \
      -DMAKEFOPT4="'$(COMPILE_OPT4)'" \
      -DMAKEFOPT5="'$(COMPILE_OPT5)'" \
      -DMAKEFOPT6="'$(COMPILE_OPT6)'" \
      -DMAKEFOPT7="'$(COMPILE_OPT7)'" \
      -DMAKEFOPT8="'$(COMPILE_OPT8)'" \
      -DMAKELIBSTD1="'$(LIBSTD1)'" \
      -DMAKELIBSTD2="'$(LIBSTD2)'" \
      -DMAKELIBSTD3="'$(LIBSTD3)'" \
      -DMAKELIBSTD4="'$(LIBSTD4)'" \
      -DMAKELIBSTD5="'$(LIBSTD5)'" \
      -DMAKELIBSTD6="'$(LIBSTD6)'" \
      -DMAKELIBSTD7="'$(LIBSTD7)'" \
      -DMAKELIBSTD8="'$(LIBSTD8)'" \
      -DCOMPDATE="'$(COMPDATE)'"


.SUFFIXES: .F90 .f90 .f95 .F95

all: bound


bound:  $(OBJ)
	$(F90)   -o npscattering $(FOPT) $(FFLAGS) $(OBJ) $(LIB)




.F90.o          :
	$(F90) $(FFLAGS) $(FOPT)  -c $<

.f95.o          :
	$(F90) $(FFLAGS) $(FOPT)  -c $<

.F.o          :
	$(F90) $(FFLAGS) $(FOPT)  -c $<

.f.o          :
	$(F90) $(FFLAGS) $(FOPT)  -c $<


clean: 
	rm -f npscattering $(objectsbound)  *.mod core *.o 



