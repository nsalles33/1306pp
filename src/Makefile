
F90=gfortran-mp-7

FFLAGS= -Wall

all: KMC.x

S90 = KMC_routines.f90 \
      Kernel.f90

O90=$(S90:.f90=.o)

KMC.x: $(O90)
	$(F90) $(FFLAGS) $(O90) -o KMC.x

%.o: %.f90
	$(F90) -c $(FFLAGS) $<

clean:
	rm -f *.o *.mod KMC.x
