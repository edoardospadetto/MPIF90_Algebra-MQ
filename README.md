# Hamiltonian parallel diagonalization using ScaLAPACK and SLEPc 
Journey across MPI (Message Protocol Interface), ScaLAPACK and SLEPc using Fortran 90. 

This program computes in parallel the eigenvalues and eigenvectors of an Hamiltonian. In our project we used the Heisenber Model Hamiltonian and the Transverse Field Ising Model.
The [first version](scalapack_version) is based on ScaLAPACK, which mainly works with dense matrices. The [second version](slepc_version) is based on SLEPc, which is optimized for sparse matrices. 

 
- Added a matrix multiplication code to perform Matmul in HPC environment


