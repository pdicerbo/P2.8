Assignment 3 - Parallel tempering
==================================

One problem frequently encountered in molecular simulations is the so-called sampling problem.
When high energy barriers between states exist, simulations tend to get trapped in local minima.
To overcome this problem, a large number of enhanced sampling techniques have been proposed.
One of the most common one is parallel tempering (or replica exchange), in which different replicas
of the system are run at different temperatures, thereby accelerating the crossing of the barriers.
The coordinates of neighbor replicas are swapped using a Metropolis criterion.

1. Implement a framework for multi-replica simulations. You should split MPI_COMM_WORLD in
   subgroups so that each simulation runs on multiple processors. Modify the program in such a way
   that it accepts multiple input files (e.g. ``mpirun -np 4 ../src/simplemd.x in0 in1`` should
   run two simulations with 4 processors each). In this manner you can generate input files
   ``in0`` and ``in1`` with a script so that they write on different output files and do not interfere
   between each other. This will allow different independent simulations to be run.
   Use a different ``idum`` for each of them so that trajectories are different.

2. Implement in simplemd a parallel tempering alogrithm, in which each processor runs a replica and
   exchanges between them are attempted every m steps. The exchange probability has to follow the Metropolis-Hastings criterion.
   The right place to add it is just after (or before) the energy file is written, that is at the end of the MD loop.
   You should do the following things:
   - Add an input option "exchangestride" to control how frequently exchanges should be tried
   - Everytime an exchange is tried, choose a partner replica. E.g., on even exchanges, replica 0 should try to exchange with replica 1, replica 2 with replica 3, etc. On even exchanges, replica 1 should try to exchange with replica 2, etc. The following should make the job:
   ::

     int partner=irep+(((istep+1)/exchangestride+irep)%2)*2-1;
     if(partner<0) partner=0;
     if(partner>=nrep) partner=nrep-1;

  - Send to the partner replica the value of configurational energy and temperature.
    range.
  - Compute the acceptance using the Metropolis criterion (see lecture).
  - Compare the acceptance with the random number to decide if the exchange should be accepted or not.
  - Notice that the two replicas should use the same random number to take a consistent decision.
    To this aim, you could either send the random number together with configurational energy and temperature
    or transfer a boolean at this stage.
  - In case the exchange is accepted, swap coordinates and velocities among the two replicas. Also scale the velocities
    according to the temperature of the two replicas.

3. Verify that the parallel tempering implementation is working correctly. To this aim:

   - Run 16 simulations at temperatures from 1 to 1.5 following a geometric schedule.
   - Compute the average value of the potential energy (column 4 from energies.dat file)
   - Verify that if you do not make any exchange (e.g. set exchangestride to a huge number)
     or if you do frequent exchanges (e.g. set exchangestride to 10) you get compatible numbers.
   - Compute the average acceptance rate and check how it depends on the size of the system.









