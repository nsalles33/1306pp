compiler = gfortran
objects = routines.o maintest.o
flags = -O
double = -fdefault-real-8

all: maintest.x

maintest.x : $(objects)
	$(compiler) $(double) -o maintest.x $(objects)

%.o: %.f90
	$(compiler) $(double) ${flags} -c $<

clean: 
	rm -f *.o
	rm -f *.mod
	rm -f *.x
