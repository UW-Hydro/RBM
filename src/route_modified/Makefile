########################################################################
### rout.f makefile #################################################### 
########################################################################
#
# Routing algorithm written originally by D. Lohmann
# and modified by J. Yearsley for purposes of preparing files necessary
# to run the semi-Lagrangian stream temperature model, RBM.
#
# 
#This program uses the non-standard Fortran argument GETARG
#
FFLAGS = -O -C 
#
# Uses the GFORTRAN compiler
# 
FC=gfortran

HFILES=    	parameter.h

OBJECTS=	rout.o 			\
		make_convolution.o	\
                init_routines.o		\
                read_routines.o		\
                write_routines.o	\
                unit_hyd_routines.o

exe:		$(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o rout

rout.o:				rout.f
make_convolution.o:		make_convolution.f
init_routines.o:		init_routines.f
read_routines.o:		read_routines.f
write_routines.o:		write_routines.f
unit_hyd_routines.o:		unit_hyd_routines.f

clean:
	/bin/rm *.o
