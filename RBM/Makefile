########################################################################
### rbm_DHSVM.f makefile ############################################### 
########################################################################
#
# Semi-Lagrangian stream temperature model
#
# The compiler
FC = gfortran

#If compiling on SUN and LINUX use  (remember -O)
FFLAGS = -O3 -C 

FCFLAGS = -O3

HFILES =	RBM.fi

OBJECTS =	RBM.o

exe:	$(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o RBM

RBM.o: RBM.f

clean:
	/bin/rm *.o

%.o: %.f
	$(FC) $(FCFLAGS) -c  $<
