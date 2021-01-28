! Hello World MPI
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
	
	! communicate size i.e number of process, 
	! MPI_COMM_WORLD is not the only communicator in MPI
	call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr) 
	
	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

	!get hostname for each process
	call MPI_GET_PROCESSOR_NAME(hostname, namesize, ierr)
	!hostname(1:namesize) the host name is smaller than the max size, to keep it small 
	!that sintax is preferred 
	print*, "Hello I am  : " , hostname(1:namesize), "rank " , rank, "nprocs" , nprocs
	
	
	!invoke single processes
	! in C++ send is like this 
	! SYNOPSIS
       !int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)

	!INPUT PARAMETERS
	 !      buf    - initial address of send buffer (choice), in recv where to store the data
	 !     count  - number of elements in send buffer (nonnegative integer)
	 !      datatype - datatype of each send buffer element (handle)	      
	 !      dest   - rank of destination (integer) / in recv represent the rank of the sender
	 !      tag    - message tag (integer), used for to be sure it is the right data
	 !      comm   - communicator (handle)
	 !+ ierr in fortran 
	 
	 !IN RECV we have also a status variable

	if(rank .eq. 0) then 
		data0 = 50 
		print*, "data0 in proc w rank " , rank, "is " ,data0
		
		!Then i send to the other process the data 
		call MPI_SEND(data0 , 1, MPI_INT, 1, 1, MPI_COMM_WORLD, ierr)
		! Receive back 
		print*, "sended from 0 "
		!wrong params keep waiting 4 ever
		call MPI_RECV(data0, 1, MPI_INT, 1, 2, MPI_COMM_WORLD,status1, ierr )
		print*, "received in 0"
		print*, "datao is " , data0 
		
		print*, "Finished process " , rank
		print*, " "
	end if 
	
	if (rank .eq. 1) then 
		!I noticed that we are using variable that exists for every process, also 
		! rank = 1 owns data 0, let's see what prints when i ask for that value.
		
		call MPI_RECV(data1, 1, MPI_INT, 0, 1, MPI_COMM_WORLD,status1, ierr )
		print*, "RECEIVED in 1"
		print*, "Rnk = " , rank, "see the received , data1" , & 
			data1, "and also the other variable as i mentioned data0 = " , data0
		data1=100 
		print*, "Rnk = " , rank, " modified data1 " , data1
		call MPI_SEND(data1 , 1, MPI_INT, 0, 2, MPI_COMM_WORLD, ierr)
		print*, "sended from 1"
		print*, "Finished process " , rank
		print*, " "
		
	end if
        ! Finalize MPI 
	call MPI_FINALIZE(ierr)

end program
