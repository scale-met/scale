################################################################################
#
# Makefile for each test project
# %> make clean       (cleanup)
# %> make depends     (generate depends.info and include)
# %> make
#
# if you have the problem about compilation order,
# please try the following commands;
#
# %> make -n       (display all procedure for check)
# or
# %> make -B       (ignore the timestamp)
#
################################################################################
include ../../sysdep/Makedef.$(SCALE_SYS)

VPATH = ../../common:../../io:../../grid:../../atmos:../../ocean

PROJ_EXE = init_boundary

all: makedir makebin

makedir:
	mkdir -p ../../../bin/$(SCALE_SYS)

makebin: $(PROJ_EXE)
	install $(PROJ_EXE) ../../../bin/$(SCALE_SYS)

$(PROJ_EXE): $(PROJ_EXE).o
	$(LD) $(LFLAGS) -o $(PROJ_EXE) $(PROJ_EXE).o $(MODS)


.SUFFIXES: .o .f90 .f .c
.f90.o:
	$(FC)  $(FFLAGS) -o $@ -c $<

.f.o:
	$(F77) $(FFLAGS) -o $@ -c $<

.c.o:
	$(CC)  $(CFLAGS) -o $@ -c $<

%.o: %.mod


.PHONY: depends freeze clean copylocal mktar

clean:
	rm -f $(PROJ_EXE) depends.info *.mod *.o *~

freeze: clean depends copylocal mktar

copylocal:
	@for i in $(MODSRC) ; do \
		cp $$i . ; \
	done
	cp ../../sysdep/Makedef.$(SCALE_SYS) .

mktar:
	tar -czf ../$(PROJ_EXE).tgz ./

depends:
	@../mkdepends.pl $(PROJ_EXE) $(VPATH)

-include depends.info
