SHELL=/bin/sh

CC=gcc-mp-7
F90=gfortran-mp-7

FFLAGS= -Wall -cpp -shared -DSHARED_LIB -fPIC 
CFLAGS= -fPIC -shared

all: sys01 sys02 sys01c

sys01:
	$(F90) $(FFLAGS) -o lib_sys01.so ../lev00_new_types.f90 ../lev00_module.f90 ../lev00_tools.f90 sys01_diffvac.f90

sys02:
	$(F90) $(FFLAGS) -o lib_sys02.so ../lev00_new_types.f90 ../lev00_module.f90 ../lev00_tools.f90 sys02_adsrbdesrb.f90

sys01c:
	$(F90) -c -DSHARED_LIB -fPIC ../lev00_new_types.f90 
	$(CC) $(CFLAGS) -o lib_sys01c.so lev00_new_types.o sys01_diffvac.c

clean:
	rm *.o *.mod *.so


