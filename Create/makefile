#
# Makefile for the Fortran 90 code, Create_File.f90, used to develop
# forcing functions for the stream temperature model, RBM.  Output
# is from DHSVM Version 3.1.1
#
# Start of the makefile
#
f90comp=gfortran
# Defining variables
#
objects = Create_File.o Julian.o
# Makefile
Create_File: $(objects)
	$(f90comp) -o Create_File $(objects)
Julian.o: Julian.f90
	$(f90comp) -c Julian.f90
Create_File.o: Create_File.f90
	$(f90comp) -c Create_File.f90

# Cleaning everything
clean:
	rm Create_File
	rm $(objects)
