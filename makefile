OBJ= coul90.o \
    main.o


LIB = -L ../lapack-3.9.0  -llapack  -lrefblas
F90 = gfortran
FFLAGS = -O2 -Wtabs   -ffixed-line-length-0  

.SUFFIXES: .F90 .f90 .f95 .F95

all: bound


bound:  $(OBJ)
	$(F90) -o npbound $(FFLAGS) $(OBJ) $(LIB)




.F90.o          :
	$(F90) $(FFLAGS) -c $<

.f95.o          :
	$(F90) $(FFLAGS) -c $<

.F.o          :
	$(F90) $(FFLAGS) -c $<

.f.o          :
	$(F90) $(FFLAGS) -c $<


clean: 
	rm -f npbound $(objectsbound)  *.mod core *.o 


