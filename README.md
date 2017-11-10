# FEM-Explicit2D
This is a simple finite element 2D explicit time integration code. The problem is a square domain with one end. fixed and one end under a dynamic stretch over time. 

This is a very draft version. I only made sure the program can run. The immediate next step is:
 1. To add output subroutines so I can visualize the result in other software like Paraview.
 2. To add comments to the code (IMPORTANT)
 3. To organize the subroutines so that the pass of data is more clearer. 

At this stage, to run the code, I use gfortran compiler on University of Cincinnati's HPC system Homer. The only library the program uses right now is BLAS. In the short future it will also need to use LAPACK. I am also planning to make the program parallezible in future versions using OpenMP and Open-MPI. It is necessary when the problem becomes big and computationally expensive. 
