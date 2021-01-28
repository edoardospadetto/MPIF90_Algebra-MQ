# MPIF90_examples
Journey across MPI (Message Protocol Interface) using Fortran 90
---
## 1) Hello World example

Just the MPI Hello World with 16 processes.

hello_world.sh
---
## 2) MPI SEND MPI RECEIVE 
Some considerations: 

All the variables are declared in different processes in this code, but stores different values .
The values are setted with SEND and RECEIVE routines. 
