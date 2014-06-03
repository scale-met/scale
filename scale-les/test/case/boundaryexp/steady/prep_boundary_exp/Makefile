##----------------------------------------------------------
include ./Mkinclude

OBJECT1=bound_exp

MODS=init_bdy_exp.o

##----------------------------------------------------------
## make rules

all: main

main: init_bdy_exp.o
	$(FC) $(FCOPTIM) $(FFLAGS) -o $(OBJECT1) $(MODS) $(LIBS)

#f90.o:
init_bdy_exp.o: init_boundary_experiment.f90
	$(FC) $(FCOPTIM) $(FFLAGS) $(INCLUDES) -o $@ -c $<

clean:
	rm -f *.o *.mod $(OBJECT1) $(OBJECT2)
##----------------------------------------------------------
