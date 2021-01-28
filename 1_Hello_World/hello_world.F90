! Hello World MPI
! Author: Edoardo Spadetto 2021

program hello_mpi 

	use mpi
	
	implicit none 
	
	! #### MPI Variables ####
	integer :: ierr ! MPI error handler 
	integer :: rank ! the process ID
	integer :: nprocs ! number of processes
	
	! Initialize MPI 
	call MPI_INIT(ierr)
	
	call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr) 
	
	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
	! Finalize MPI 
	
	print*, "Hello Wolrd from proc : " , rank , "of " , nprocs
	call MPI_FINALIZE(ierr)

end program
