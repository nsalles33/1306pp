SHELL=/bin/sh

CC=gcc-mp-7
F90=gfortran-mp-7

FFLAGS= -Wall -cpp -shared -DSHARED_LIB -fPIC 

all: ex01 ex02

ex01:
	$(F90) $(FFLAGS) -o lib_ex01.so ../KMC_routines.f90 ../tools.f90 ex01_diffvac.f90

ex02:
	$(F90) $(FFLAGS) -o lib_ex02.so ../KMC_routines.f90 ../tools.f90 ex02_adsrbdesrb.f90

clean:
	rm *.o *.mod *.so


