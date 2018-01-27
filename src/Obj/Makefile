# -*- Makefile -*-
SHELL=/bin/sh

EXE=../EXE.x

F90=gfortran-mp-7

FFLAGS= -Wall -O

all: $(EXE)

VPATH=../
S90 = tools.f90 \
      KMC_routines.f90 \
      Init_system.f90 \
      Kernel.f90

O90=$(S90:.f90=.o)

$(EXE): $(O90)
	$(F90) $(FFLAGS) $(O90) -o $(EXE)

%.o: %.f90
	$(F90) -c $(FFLAGS) $<

clean:
	rm -f *.o *.mod $(EXE)
