##----------------------------------------------------------
include ./Mkinclude

OBJECT1=bound_exp

MODS=prep_bdy_exp.o

##----------------------------------------------------------
## make rules

all: main

main: prep_bdy_exp.o
	$(FC) $(FCOPTIM) $(FFLAGS) -o $(OBJECT1) $(MODS) $(LIBS)

#f90.o:
prep_bdy_exp.o: prep_boundary_experiment.f90
	$(FC) $(FCOPTIM) $(FFLAGS) $(INCLUDES) -o $@ -c $<

clean:
	rm -f *.o *.mod boundary_exp_00* $(OBJECT1) $(OBJECT2)
##----------------------------------------------------------
