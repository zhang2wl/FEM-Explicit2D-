# FEM-Explicit2D
This Fortran 90 code provides a framework of solving dynamic solid mechanics problems. It is not a comprehensive nonlinear FEM solver but only provides information about how I did my FEM parallel programming. 

Features:
 1. Providing the element, node, and BC data in LS-DYNA input format (which can be obtained from the free software LS-Prepost), one can obain the displacement, velocity, acceleration and stress over time. 
 2. Only 2D quad element is supported at this moment. 
 3. MPI is used for parallel computing. Parallization is applied to element internal force (R_int) calculation. Each processor only calculate the R_int vector on a sub-domain, and all the other processors send the R_int to one processor for it to calculate acceleration for the next time step. For the problem tested (10,000 elements), the speed up in the parallalization part is 8 times using 10 processors.
 4. An open-source library METIS is used to do domain composition
 5. A highly optimized linear algibra library BLAS is used to do some basic linear algibra calculation.
 6. vtk format is used to output data to be viewed in Paraview.

Future updates:
 1. To optimize the coding style to separate input and computing part.
 2. Add nonlinearity to the problem
 3. Update the comments
 4. Optimize the parallization algorithm if I have time. 
