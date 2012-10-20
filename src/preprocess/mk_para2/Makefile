# Makefile for compiling Cloud Resolving Model
# created by K.Suzuki (2005/09/26)

.SUFFIXES:
.SUFFIXES: .o .f90 .f

# Definition of Macro

#FC1 = f90
#FC2 = f77
FC1 = gfortran-mp-4.7
FC2 = gfortran-mp-4.7

# No Option
FFLAGS =

# IBM fortran on Mac Quad
# FFLAGS = -O5

# Absoft fortran on Mac
# FFLAGS = -O

# Dec fortran on alpha
# FFLAGS = -fast -tune host -O4 -assume byterecl -convert big_endian

# gfortran on Mac
FFLAGS = -O3 -fconvert=big-endian

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
