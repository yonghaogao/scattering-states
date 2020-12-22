OBJ= coul90.o \
   main.o


LIB = -mkl
F90 = ifort
FFLAGS = -g

.SUFFIXES: .F90 .f90 .f95 .F95

all: bound


bound:  $(OBJ)
	$(F90) -o npscattering $(FFLAGS) $(OBJ) $(LIB)




.F90.o          :
	$(F90) $(FFLAGS) -c $<

.f95.o          :
	$(F90) $(FFLAGS) -c $<

.F.o          :
	$(F90) $(FFLAGS) -c $<

.f.o          :
	$(F90) $(FFLAGS) -c $<


clean: 
	rm -f npscattering $(objectsbound)  *.mod core *.o 



