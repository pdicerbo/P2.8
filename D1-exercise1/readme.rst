Assignment 1 - Parallelize force calculation
============================================

The calculation of forces in Molecular Dynamics simulations of a system of Lennard-Jones particles
involves the calculation of all pairwise distances at each time step of the simulation.
This in principle requires a double loop over all particles.
Later we will use neighbor lists. For the moment, disable them by setting ``listcutoff``
to a huge number (e.g. 1e6).

TASK A

1. Modify the code and the makefile so as to use MPI library. Look for all ``fopen()``
   statements and make sure that processes which are not root (rank!=0) are opening 
   ``/dev/null`` so that there are no duplicated lines in the output.

2. Parallelize the computation of the forces within the double loop using MPI.
   A possible approch is to modify the loop on particle pairs and makes it such that
   every processor only includes some of the pairs. At the end of the loop you can use
   ``MPI_Allreduce()`` to sum the forces over the processors.

3. Measure the scalability of the code, that is how the speeds depend on the number
   of processors at fixed number of atoms.

Verlet lists (neighbor lists) are used to speed up the calculation of pairwise distances and forces.
They are already implemented in the code, so you can see their effect by just setting
``listcutoff`` to a reasonable number (e.g. 3). Simulation will be much faster.

Although neighbor lists are not calculated every step, their calculation also require a loop over
all pairs. As a consequence, if your code has neighbor lists implemented, it is much more
convenient to directly parallelize their calculation with MPI.
Notice that by doing that each process will contain only a part of the list of neighbors.
This will decrease the memory usage per process.
Moreover, since each process will only contain a portion of the pairs of interacting particles,
the force calculation will be implicitly parallelized.
Also notice that steps when the neighbor list is updated are significantly slower than other steps.
Thus, for doing benchmarks you should use a large number of steps (e.g. ``nsteps 1000``)

TASK B

1. Parallelize the computation of the neighbor lists. Each process should only
   store a subset of the pairs of interacting atoms.
2. Measure the scalability of the code.

