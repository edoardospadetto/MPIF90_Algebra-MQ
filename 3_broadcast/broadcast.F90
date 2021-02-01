!Broadcast MPI
! Author: Edoardo Spadetto 2021

program send_recv

	use mpi
	
	implicit none 
	
	! #### MPI Variables ####
	integer :: ierr ! MPI error handler 
	integer :: rank ! the process ID
	integer :: nprocs ! number of processes
	
	!####SEND RECV variables
	integer, dimension (MPI_STATUS_SIZE) :: status1 !status of the send receive calls
	
	character*(MPI_MAX_PROCESSOR_NAME) :: hostname
	integer :: namesize 
	
	
	!## TEST DATA to send and receive ###
	integer :: data0 =0, data1=0
	
	! Initialize MPI 
	call MPI_INIT(ierr)
	

	call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr) 
	
	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)


	call MPI_GET_PROCESSOR_NAME(hostname, namesize, ierr)
	
	
	
	!MPI_BCAST shares a variable from one process to others.
	!If we want to see the variable we must call the function inside
	!all the processes that want the variable. 
	!If the BCAST is not called in the "sender" the other that asked for the 
	!variable keep waiting.
	
	!MPI_BCAST( the data , size of the data, MPI_type, fromwho, MPI_COMM_WORLD,ierr )
	

	if(rank .eq. 0) then 
		data0 = 50 
	call MPI_BCAST(data0,1, MPI_INT, 0 ,MPI_COMM_WORLD,ierr )
		 
	end if
	
	call MPI_BCAST(data0,1, MPI_INT, 0, MPI_COMM_WORLD,ierr )
	
	
	print*, "data0", rank , data0 
		!call MPI_RECV(data1, 1, MPI_INT, 0, 1, MPI_COMM_WORLD,status1, ierr )
		
	
        ! Finalize MPI 
	call MPI_FINALIZE(ierr)

end program
