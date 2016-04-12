Assignment 2 - Code Linked cells
================================

As the system size increases, the size of the neighbor lists becomes very large, and their construction becomes computationally expensive, involving a double loop over all pairs and a logical testing of their distance.

An important trick used in molecular simulations is the use of linked cells.
Notice that often linked cells are used to compute forces directly. Since the program you are using
already implements neighbor lists, it is convenient to use linked cells to parallelize the update of
the neighbor list. This will make the update scale linearly in the number of particles.

1. Write a parallel code implementing the linked cells. To do it you should modify
   the routine that updates the neighbor list. In particular you should
  
   - Divide the simulation box in a number of domains such that only neighboring
     domains are interacting.
   - Make a loop over the atoms assigning each of them to one of the domains
   - Make a double loop over domains considering only pairs of domains that are close to each
     other.
   - Make a loop over the particles of neighboring domains, perhaps adding the pair to the neighbor
     list.

2. Measure the scalability of the code. 

3. Keeping the number of processors fixed, measure the scalability of the code
   when you change the number of atoms in the simulation box.
