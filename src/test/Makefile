# -*- Makefile -*-
SHELL=/bin/sh

EXE=../EXE.x

F90=gfortran-mp-7
CC=gcc-mp-7

FFLAGS= -Wall -O -fdefault-real-8 $(FDEBUG) $(LSHARED)
FDEBUG= -g -fbounds-check
LSHARED= -cpp -ldl -DSHARED_LIB

CFLAGS= -fPIC

all: flib clib

flib: main lib-f

clib: main lib-c

VPATH=./
S90 = module.f90 buffer.f90 lib_dlopen.f90 test.f90

O90=$(S90:.f90=.o)

main: $(O90)
	$(F90) $(FFLAGS) $(O90)

%.o: %.f90
	$(F90) -c $(FFLAGS) $<

lib-f:
	$(F90) -o lib.so -shared -fPIC module.f90 lib_shared.f90

lib-c:
	$(CC) -c -fPIC lib_shared.c
	$(CC) -o lib.so -shared -fPIC module.o lib_shared.o 

clean:
	rm -f *.o *.mod *.so a.out
