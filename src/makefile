#
# Makefile for the grid-based semi-Lagrangian water temperature model, RBM10_VIC
#
# Start of the makefile
#
# Defining variables
#
objects = rbm10_VIC.o Begin.o Systmm.o Particle_Track.o \
          Energy.o Julian.o tntrp.o Read_Forcing.o \
          Block_Energy.o Block_Hydro.o Block_Network.o \
          Block_Reservoir.o \
          density_subroutine.o flow_subroutine.o \
          reservoir_subroutine_implicit.o \
          reservoir_single_subroutine.o \
          Water_Balance.o Write.o
f90comp = gfortran #-O2#-O0 -g
# Makefile
rbm10_VIC: $(objects)
	$(f90comp) -o rbm10_VIC $(objects)
Block_Energy.o: Block_Energy.f90
	$(f90comp) -c Block_Energy.f90
block_energy.mod: Block_Energy.f90
	$(f90comp) -c Block_Energy.f90
Block_Hydro.o: Block_Hydro.f90
	$(f90comp) -c Block_Hydro.f90
block_hydro.mod: Block_Hydro.o Block_Hydro.f90
	$(f90comp) -c Block_Hydro.f90
Block_Network.o: Block_Network.f90
	$(f90comp) -c Block_Network.f90
block_network.mod: Block_Network.o Block_Network.f90
	$(f90comp) -c Block_Network.f90
Block_Reservoir.o: Block_Reservoir.f90
	$(f90comp) -c Block_Reservoir.f90
block_reservoir.mod: Block_Reservoir.o Block_Reservoir.f90
	$(f90comp) -c Block_Reservoir.f90
Begin.o: block_energy.mod block_hydro.mod block_network.mod block_reservoir.mod Begin.f90
	$(f90comp) -c Begin.f90
Read_Forcing.o: block_energy.mod block_hydro.mod block_network.mod block_reservoir.mod Read_Forcing.f90
	$(f90comp) -c Read_Forcing.f90 
Systmm.o: block_network.mod block_energy.mod block_hydro.mod block_reservoir.mod Systmm.f90
	$(f90comp) -c Systmm.f90
Energy.o: block_energy.mod Energy.f90
	$(f90comp) -c Energy.f90
Particle_Track.o: block_hydro.mod block_network.mod block_reservoir.mod Particle_Track.f90
	$(f90comp) -c Particle_Track.f90
Write.o: Write.f90
	$(f90comp) -c Write.f90
density_subroutine.o: density_subroutine.f90
	$(f90comp) -c density_subroutine.f90
flow_subroutine.o: block_hydro.mod block_network.mod block_reservoir.mod flow_subroutine.f90
	$(f90comp) -c flow_subroutine.f90
reservoir_subroutine_implicit.o: block_network.mod block_reservoir.mod reservoir_subroutine_implicit.f90
	$(f90comp) -c reservoir_subroutine_implicit.f90
reservoir_single_subroutine.o: block_network.mod block_reservoir.mod reservoir_single_subroutine.f90
	$(f90comp) -c reservoir_single_subroutine.f90
Julian.o: Julian.f90
	$(f90comp) -c Julian.f90
tntrp.o: tntrp.f90
	$(f90comp) -c tntrp.f90
Water_Balance.o: block_hydro.mod block_network.mod Water_Balance.f90
	$(f90comp) -c Water_Balance.f90
rbm10_VIC.o: rbm10_VIC.f90
	$(f90comp) -c rbm10_VIC.f90

# Cleaning everything
clean:
	rm block_energy.mod block_hydro.mod block_network.mod block_reservoir.mod \
           rbm10_VIC
	rm $(objects)
