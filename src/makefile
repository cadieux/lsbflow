#---
# lsb.x  makefile
#
#--- 

#--- Macro definitions

# Fortran compiler
# FC = gfortran
FC = mpif90.openmpi
# FC = mpif90.mpich2
# FC = tau_f90.sh -tau_options=-optVerbose
# FC = ifort
# Compiler flags debug/run
# gfortran:
# CFLAGS = -fdefault-real-8 -ffree-line-length-none -Wall -fcheck=all -O0 -g
# CFLAGS = -fdefault-real-8 -ffree-line-length-none -fcheck=all -O0 -g
CFLAGS = -fdefault-real-8 -ffree-line-length-none
# CFLAGS = -ffree-line-length-none
# ifort:
# CFLAGS = -r8 -llapack -lfftw3
# CFLAGS = -r8
# The LDFLAGS variable sets flags for library linker
# LDFLAGS = -llapack -lfftw3
LDFLAGS = -I/usr/include/ -L/usr/lib -llapack -lblas -lfftw3 -lm
# LDFLAGS = -I/usr/usc/fftw/3.3.3/include -L/usr/usc/fftw/3.3.3/lib -llapack -lfftw3 -lm

SOURCES = main.f90 init.f90 common.f90 mpi_setup.f90 iofiles.f90 modfftw.f90 horfft.f90 chebyshev.f90 derivatives.f90 orrsomm.f90 helmholtz.f90 navierstokes.f90 io.f90  filters.f90 mpi_io.f90 dealiasing.f90 sgsterm.f90
FSOURCES =
OBJECTS = $(SOURCES:.f90=.o)
OBJECTS += $(FSOURCES:.f=.o)
TARGET = ../run/lsb.x

#--- Target dependencies

$(TARGET) : $(OBJECTS)
	$(FC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

#--- Specific dependencies

mpi_setup.o : common.o

modfftw.o : common.o

horfft.o : common.o mpi_setup.o modfftw.o

io.o : common.o iofiles.o horfft.o

mpi_io.o : common.o iofiles.o horfft.o

chebyshev.o : common.o modfftw.o

derivatives.o : common.o chebyshev.o

orrsomm.o : common.o iofiles.o

helmholtz.o : common.o

filters.o : common.o horfft.o

dealiasing.o : common.o horfft.o

sgsterm.o : common.o mpi_setup.o

init.o : common.o mpi_setup.o iofiles.o horfft.o chebyshev.o orrsomm.o helmholtz.o derivatives.o io.o mpi_io.o filters.o dealiasing.o

navierstokes.o : common.o derivatives.o horfft.o helmholtz.o io.o iofiles.o init.o  filters.o dealiasing.o sgsterm.o

main.o : common.o io.o mpi_io.o init.o navierstokes.o

#--- General rules

%.o: %.f90
	$(FC) $(CFLAGS) -c $< $(LDFLAGS)

%.o: %.f
	$(FC) $(CFLAGS) -c $< $(LDFLAGS)

#--- Utility targets
default:: opt

all:: $(TARGET)

opt::
	@make all "CFLAGS = -O2  $(CFLAGS)"

# debug ::
# 	@make all "CFLAGS = -O0 -debug all -check all -traceback -warn all $(CFLAGS)"

prof ::
	@make all "CFLAGS = -O2 -pg $(CFLAGS)"

debug ::
	@make all "CFLAGS = -O0 -fbounds-check -fbacktrace -Wall -Wextra $(CFLAGS)"

clean::
	@\rm -f *.o *.mod *genmod.f90 $(TARGET) 
