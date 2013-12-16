# Makefile for compiling Cloud Resolving Model
# created by K.Suzuki (2005/09/26)
SYSDEP_DIR  = ../../../sysdep

include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)

.SUFFIXES:
.SUFFIXES: .o .f90 .f

# Definition of Macro

FC1 = $(FC)
FC2 = $(FC)

TARGET = mkpara
OBJS =  mkpara.o        \
        mod_const.o     \
	mod_spline.o    \
	mod_spline_2d.o \
        TINVSS.o

# Suffix Rule

.f90.o:
	$(FC1) $(FFLAGS) -o $@ -c $<

.f.o:
	$(FC2) $(FFLAGS) -o $@ -c $<

# Definition of Target

$(TARGET): $(OBJS)
	$(FC1) $(FFLAGS) $(OBJS) -o $@

# Dependencies

mkpara.o	: mkpara.f90  \
                  mod_const.o mod_spline.o mod_spline_2d.o TINVSS.o
mod_const.o	: mod_const.f90
mod_spline.o	: mod_spline.f90 TINVSS.o
mod_spline_2d.o	: mod_spline_2d.f90 TINVSS.o
TINVSS.o	: TINVSS.f

# Cleaning

clean:
	rm -rf *.o *.mod *.f90~ $(TARGET)
