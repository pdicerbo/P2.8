Introduction
=============

Take some time to familiarize with the simplemd code in ``src/simplemd.cpp``
In particular, make sure you understand the functions *compute_forces*, *compute_list* and *check_list*.
You can compile it typing ``make all``.

Also try to use it to run short MD simulations. Look in the input/ directory, where you will find
a sample input file. You can run simplemd by typing
::

  ../src/simplemd.x in

Keywords are explained in the input file. The ones that are relevant for the assignments are:

- ``listcutoff``

  An appropriate value for a calculation with neighbor lists is 3.
  Set it to a huge number (e.g. 1e6) to disable neighbor list.

- ``nsteps``

  Use it to tune the length of the simulation.
  Notice that benchmarks will typically need at least a few hundreds of steps to
  be meaningful.

You can generate starting configurations using the command:
::

  ./lattice 10 > crystal.xyz

Here 10 means that your xyz file will contain 10*10*10*4 atoms.
You will use this to check how the performances of the code scale with the number of particles.

Moreover, you can have a look to the ``energies.dat`` file. It has several columns:

*step-number* *simulation-time* *instantaneous-temperature* *configurational-energy* *total-energy* *conserved-quantity*

You can use it to quickly check if two simulations are identical. E.g., if you optimize the code and run it
again, the new ``energies.dat`` file should be identical to the old one.

