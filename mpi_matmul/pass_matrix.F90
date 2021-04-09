module cbp_DistributeHandler
	use mpi
	implicit none 
	contains
	
	!Just print a local Matrix on terminal
	    SUBROUTINE printmatrix(mat)
	      IMPLICIT NONE
	      INTEGER :: ii,jj
	      INTEGER, DIMENSION(:,:), INTENT(IN) :: mat
	      INTEGER, DIMENSION(SIZE(SHAPE(mat))) :: sizem
	      sizem=SHAPE(mat)

	      PRINT*, ""
	      DO ii=1,sizem(1)
		DO jj=1,sizem(2)
		  WRITE(*, fmt="(I2, tr2)", advance="no") mat(ii,jj)
		END DO
		PRINT*, ""
	      END DO


	    END SUBROUTINE printmatrix
	
	
	!Computes the rank of the process block with coordinate merow+down, mecol+side
	!If the process grid pass the limits it restart from the other side
	!Pacman Effect :3
	function next2me(merow, mecol, down,side, totrows, totcols) result (idx)
		implicit none
		integer :: merow, mecol, side, down, totrows, totcols, itcol, itrow
		integer :: idx
	
		itcol = (mecol+side) 
		itrow = (merow+down)
		
		if (itrow .gt. totrows)  then	
			itrow = mod(itrow-1,totrows)+1
		else if( itrow .le. 0  ) then 
			itrow = totrows- mod(itrow, totrows)
		end if 
		
		if (itcol .gt. totcols)then
			itcol = mod(itcol-1,totcols)+1
		else if( itcol .le. 0  ) then 
			itcol = totcols -mod(itcol, totcols)
		end if 
		
		!get to the rank
		!print*, merow,mecol, "||", down, side,"||" ,itrow, itcol
		idx = (itcol-1)*totrows + itrow-1
		
		return
		
	end function
	
	
	!Inititalizes the process grid
	!Block Cij is located alwais in the ij process
	!Aik is in the i row
	!Bkj in the j column.
	!It selects the k to assing at each ij grid block
	
	subroutine init_process_grid(process_grid, rowproc, colproc,rank)
		implicit none
		integer , dimension(:,:) :: process_grid
		integer , dimension(:), allocatable :: setofpieces
		integer :: ii, jj, kk, rowproc, colproc, ierr, rank
		
		allocate(setofpieces(min(rowproc, colproc)))
		
	
		
		process_grid = 0
		do ii = 1, max(rowproc,colproc)
			do jj = 1, max(rowproc,colproc) !columns
				process_grid(ii,jj) = mod( (ii+jj-2) , max(rowproc,colproc) )+1
			
				if( process_grid(ii,jj) .gt. min(colproc,rowproc) ) then
					process_grid(ii,jj) = 0 !Invalid Pieces
				end if
			end do 
		
		end do
		
		do kk = 1, min(rowproc, colproc)
			setofpieces(kk) = kk 
		end do 
		
		do ii = 1, min(rowproc, colproc)
		do jj = 1, max(rowproc, colproc)
		
		if(rowproc .gt. colproc) then
			if(process_grid(jj,ii) .eq. 0) then
			do kk = 1, min(rowproc, colproc)
				if(.not. any(abs(process_grid(jj,1:colproc)) == setofpieces(kk)) ) then 
					process_grid(jj,ii) = - setofpieces(kk)
				end if
			end do
			end if
		
		else if(rowproc .lt. colproc) then
			if(process_grid(ii,jj) .eq. 0) then
			do kk = 1, min(rowproc, colproc)
				if(.not. any(abs(process_grid(1:rowproc, jj)) == setofpieces(kk)) ) then 
					process_grid(ii,jj) = - setofpieces(kk)
				end if
			end do
			end if
		end if
		
		
		end do 
		end do
		
		
		deallocate(setofpieces)
		
		!Dispatch Process Grid to all processes
		if(rank .eq. 0 ) then
		call printmatrix(process_grid(:rowproc, :colproc))
		call MPI_BCAST(process_grid,size(process_grid,dim=1)**2, MPI_INTEGER, 0 ,MPI_COMM_WORLD,ierr ) 
		else
		call MPI_BCAST(process_grid,size(process_grid,dim=1)**2, MPI_INTEGER, 0 ,MPI_COMM_WORLD,ierr )
		end if
		
	end subroutine 
	
	
	!Next step of the process grid to perform matmul
	subroutine update_process_grid(process_grid, rowproc, colproc, rank)
	        use mpi
		implicit none
		integer , dimension(:,:) :: process_grid
		integer :: ld, fromc, jj, colproc, rowproc, rank
		integer, dimension(:), allocatable :: temp_ 
		integer ::ierr
		
	
		ld = size(process_grid, 1)
		
		if(rank .eq. 0 ) then 
		allocate(temp_(ld));
		
		
		if (colproc .lt. rowproc) then
			
			temp_ = process_grid(:,1)
		
			do jj = 1, ld-1
				fromc = mod(jj,ld) +1
				process_grid(:,jj) = process_grid(:,fromc)  	
			end do 
				process_grid(:,ld) = temp_
			
			
			
			do jj = 1, ld
				if (process_grid(jj,colproc) .ne. 0 ) then
					if( findval(process_grid(jj,1:colproc), -abs(process_grid(jj,colproc)) ) .gt. 0)	then 
					process_grid(jj ,findval(process_grid(jj,1:colproc), -abs(process_grid(jj,colproc)) ) ) = 0
					end if
				end if
				if (process_grid(jj,ld) .ne. 0 ) then 
					if( findval(process_grid(jj,1:colproc), 0)  .gt. 0)	then
					process_grid(jj ,findval(process_grid(jj,1:colproc), 0) ) = - abs(process_grid(jj,ld))
					end if
					!process_grid(jj,ld) =  0
				end if 
				
			end do
		
		else if (colproc .gt. rowproc) then
		
			temp_ = process_grid(1,:)
		
			do jj = 1, ld-1
				fromc = mod(jj,ld) +1
				process_grid(jj,:) = process_grid(fromc,:)  	
			end do 
				process_grid(ld,:) = temp_
			
		
			do jj = 1, ld
				if (process_grid(rowproc,jj) .ne. 0 ) then
					if( findval(process_grid(1:rowproc, jj), -abs(process_grid(rowproc, jj))) .gt. 0) then 
					process_grid(findval(process_grid(1:rowproc, jj), -abs(process_grid(rowproc, jj))), jj ) = 0
					end if
				end if
				if (process_grid(ld,jj) .ne. 0 ) then 
				       if( findval(process_grid(1:rowproc,jj), 0)  .gt. 0)	then
					process_grid(findval(process_grid(1:rowproc,jj), 0), jj ) = - abs(process_grid(ld,jj))
					end if
				end if 
				
			end do
		else 
			temp_ = process_grid(1,:)
			
			do jj = 1, ld-1
				fromc = mod(jj,ld) +1
				process_grid(jj,:) = process_grid(fromc,:)  	
			end do 
				process_grid(ld,:) = temp_
		end if
			deallocate(temp_);
			
			call MPI_BCAST(process_grid, size(process_grid,dim=1)**2, MPI_INTEGER, 0 ,MPI_COMM_WORLD,ierr )
	
		else 
			call MPI_BCAST(process_grid,size(process_grid,dim=1)**2, MPI_INTEGER, 0 ,MPI_COMM_WORLD,ierr )
			
		end if
		
		 
	end subroutine
	
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	
	!given row and column in the process grid computes the rank
	function grid2rank(row, col , rowproc, colproc) result(idx)
		implicit none
		integer :: row, col, idx, rowproc, colproc
	
		idx = (col-1)*rowproc + row-1
		
		return
		
	end function
	
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	!Check error of mpi
	subroutine mpic_checkerror(ierr)
	implicit none
	character(512) ::what
	integer :: ierr, ierr0,ll
	if(ierr .ne. 0 ) then
	call MPI_ERROR_STRING(ierr, what, ll, ierr0)
	print*, trim(what)
	end if
	end subroutine
	
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	!Send and Receive a data.
	function MPI_comm_wrap(thedata,merow,mecol,otherrow,othercol,rowproc,colproc,mod,tag) result (request)
		implicit none
		integer, dimension(:) :: thedata
		integer :: dims , otherrow,othercol, request, ierr, tag, merow, mecol, rowproc, colproc
		character*1 :: mod
		
		dims = size(thedata, dim = 1)
		!print *,"dimension : ",  dims
		
		
		if (mod .eq. 's') then
			call MPI_ISEND( thedata, dims, MPI_INTEGER,&
				 	grid2rank(otherrow,othercol,rowproc, colproc) , & 
				 	4*grid2rank(otherrow,othercol,rowproc,colproc)+tag ,&  !tag is the  dest
				 	MPI_COMM_WORLD, request, ierr)
			call mpic_checkerror(ierr)
		else if (mod .eq. 'r') then
			call MPI_IRECV( thedata, dims, MPI_INTEGER,&
				 	grid2rank(otherrow,othercol,rowproc,colproc) , & 
				 	4*grid2rank(merow,mecol,rowproc,colproc)+tag ,&  !tag is the receiver, so me
				 	MPI_COMM_WORLD, request, ierr)
			call mpic_checkerror(ierr)
		else 
			print*, "Mod error"
		end if
		
	
		return
	end function
	
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@òòòòòò@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		
		
!Local Matrix -> 2D array of geenric dimensions
!MatrixDescriptor -> 
!Merow-mecol-> Process grid coordinates of the process who perform the action
!Otherrow-othercol->Process grid coordinates of eho receives (if i am sender), of send (if i am receiver)
!Mod -> 's' send 'r' receive
!Tag -> A number from 0 to 10 Don't worry its made unique for this step. 

		
	subroutine MPI_comm_wrap_matmul(LocalMatrix, MatrixDescriptor ,& 
				      merow,mecol,otherrow,othercol,rowproc,colproc,& 
				      mod,tag ,requestM, requestD)
		
		implicit none
		complex, dimension(:,:) :: LocalMatrix
		integer, dimension(4) :: MatrixDescriptor
		
		integer :: dims(2) , otherrow,othercol, ierr, tag, merow, mecol, rowproc, colproc
		integer :: requestM , requestD
		character*1 :: mod
		
		dims(1) = size(LocalMatrix, dim = 1)
		dims(2) = size(LocalMatrix, dim = 2)
		!print *,"dimension : ",  dims
		
		
		if (mod .eq. 's') then
		
			call MPI_ISEND( LocalMatrix, dims(1)*dims(2), MPI_COMPLEX,&
				 	grid2rank(otherrow,othercol,rowproc, colproc) , & 
				 	10*grid2rank(otherrow,othercol,rowproc,colproc)+2*tag ,&  !tag is the  dest
				 	MPI_COMM_WORLD, requestM, ierr)
				 	
			call MPI_ISEND( MatrixDescriptor, 4, MPI_INTEGER,&
				 	grid2rank(otherrow,othercol,rowproc, colproc) , & 
				 	10*grid2rank(otherrow,othercol,rowproc,colproc)+2*tag+1 ,&  !tag is the  dest
				 	MPI_COMM_WORLD, requestD, ierr)
				 	
			call mpic_checkerror(ierr)
		else if (mod .eq. 'r') then
			
			call MPI_IRECV( LocalMatrix, dims(1)*dims(2), MPI_COMPLEX,&
				 	grid2rank(otherrow,othercol,rowproc,colproc) , & 
				 	10*grid2rank(merow,mecol,rowproc,colproc)+2*tag ,&  !tag is the receiver, so me
				 	MPI_COMM_WORLD, requestM, ierr)
				 	
			call MPI_IRECV( MatrixDescriptor, 4, MPI_INTEGER,&
				 	grid2rank(otherrow,othercol,rowproc,colproc) , & 
				 	10*grid2rank(merow,mecol,rowproc,colproc)+2*tag+1 ,&  !tag is the receiver, so me
				 	MPI_COMM_WORLD, requestD, ierr)
			
			
			call mpic_checkerror(ierr)
		else 
			print*, "Mod error"
		end if
		
	
		
	end subroutine
	
	
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@òòòòòò@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!Various cases, sa send A, ra recv A etc..
!If rowproc > colproc you do not need copies of B. it is guaranteed that at least one process has one. hence negative id 
!are not dispatched, the same for A in the case rowproc < colproc. 
!The process ij looking at the new process grid understand who to send to  or receive from the piece of this step.
! ra is receive a , sa send a etc etc.. 
	subroutine Find_coop(old_proc_grid, new_proc_grid,idxrow,idxcol,otherrow, othercol, rowproc, colproc, mod )
	implicit none
	integer , dimension (:,:) :: old_proc_grid, new_proc_grid
	integer :: idxrow, idxcol, otherrow, othercol, rowproc, colproc
	logical :: saferow, safecol
	character*2:: mod
	
	!A PART 
	!A has to stay in the same column.
	if(mod .eq. "sa") then
		otherrow = idxrow
		if(rowproc .gt. colproc) then
			othercol = findval(abs(new_proc_grid(idxrow, 1:colproc)), abs(old_proc_grid(idxrow,idxcol))) 
			
		else 
			if (old_proc_grid(idxrow,idxcol) .gt. 0) then 
				othercol = findval(new_proc_grid(idxrow, 1:colproc), old_proc_grid(idxrow,idxcol))
			else 
				otherrow =0
				othercol =0 
			end if
			!if(new_proc_grid(othercol, othercol) .lt. 0 ) then 
			!	otherrow =0
			!	othercol =0 
			!end if 
			if (othercol .lt. 0) then
			print*, "error coop not found for sending A "
			end if 
		end if
		
		
		
			
		
	else if(mod .eq. "ra") then
		otherrow = idxrow
		if(rowproc .gt. colproc) then
			othercol = findval(abs(old_proc_grid(idxrow, 1:colproc)), abs(new_proc_grid(idxrow,idxcol))) 
			
		else 
			if (new_proc_grid(idxrow,idxcol) .gt. 0) then
				othercol = findval(old_proc_grid(idxrow, 1:colproc), new_proc_grid(idxrow,idxcol)) 
			else 
				otherrow = 0
				othercol = 0
			end if
			!if(new_proc_grid(othercol, othercol) .lt. 0 ) then 
			!	otherrow =0
			!	othercol =0 
			!end if 
			if (othercol .lt. 0) then
			print*, "error coop not found for receiving A "
			end if 
		end if
		
		
		
		
	!B part
	!B stays in the same row
	else if (mod .eq. "sb") then
		othercol = idxcol
		if(rowproc .lt. colproc) then
			otherrow = findval(abs(new_proc_grid(1:rowproc, idxcol)), abs(old_proc_grid(idxrow,idxcol))) 
			if (otherrow .lt. 0) then
			print*, "error coop not found for receiving A "
			end if 
		else 
			if (old_proc_grid(idxrow,idxcol) .gt. 0) then 
				otherrow = findval(new_proc_grid(1:rowproc, idxcol), old_proc_grid(idxrow,idxcol)) 
			else 
				otherrow = 0
				othercol = 0
			end if
			
			if (otherrow .lt. 0) then
			print*, "error coop not found for receiving A "
			end if 
			
		end if
		
		
	
	else if(mod .eq. "rb") then
		othercol = idxcol
		if(rowproc .lt. colproc) then
			otherrow = findval(abs(old_proc_grid(1:rowproc, idxcol)), abs(new_proc_grid(idxrow,idxcol)))
			 
		else 
			if (new_proc_grid(idxrow,idxcol) .gt. 0) then
				otherrow = findval(old_proc_grid(1:rowproc, idxcol), new_proc_grid(idxrow,idxcol))
			else 
				otherrow =0
				othercol =0 
			end if
			
			if (otherrow .lt. 0) then
			print*, "error coop not found for receiving B "
			end if 
			
		end if
		
		
	else 
		print*, "not valid mod to coop with a process"
	end if
	

	
	end subroutine




!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@òòòòòò@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!Find first value equal to the requested in an array, if not found resturns -1	
	function findval(array, val) result(ii)
	integer , dimension(:) :: array
	integer :: val,ii
	logical :: good 
	
	good = .false.
	do ii = 1, size(array)
		if(array(ii) .eq. val) then
		good = .true. 
		return
		exit
		end if 
	end do

	ii = -1

	return
		
	
	end function
	


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@òòòòòò@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine ComputeMatmul(localA, localB, localC, descriptorA, descriptorB , descriptorC, rowproc, colproc, rank, process_grid)
	complex, dimension(:,:,:):: localA , localB 
	complex, dimension(:,:):: localC
	integer , dimension(4,2) :: descriptorA, descriptorB
	integer :: descriptorC(4)
	
	integer , dimension(:,:)  :: process_grid
	integer :: rowproc, colproc, rank 
	integer ::nprocs
	
	
	
	integer , dimension(:,:) , allocatable :: old_process_grid
	
	integer :: idxrow, idxcol
	integer :: switch, withwhorow, withwhocol
	integer :: requestsendA, requestrecvA ,requestsendB, requestrecvB 
	integer :: requestsendDA, requestrecvDA ,requestsendDB, requestrecvDB 
	logical :: sa, ra, sb , rb
	
	integer, dimension(MPI_STATUS_SIZE) :: statussendA,statussendB, statusrecvA,statusrecvB
	integer, dimension(MPI_STATUS_SIZE) :: statussendDA,statussendDB, statusrecvDA,statusrecvDB
	
	integer :: ii,jj, pace
	integer :: ierr ! MPI error handler 
	!Here make the checks
	logical :: testrb
	
	integer , dimension(:,:), allocatable :: stat 
		
		
		
	
	ierr = 0
	allocate (old_process_grid(max(rowproc,colproc), max(rowproc, colproc)))
	
	
	!end here
	nprocs = colproc*rowproc
	pace = max(colproc,rowproc)
	switch = 1
	
	allocate(stat(pace,8))
	
	!Each process understand its identity
	idxcol= rank/rowproc+1			!grid assigned by cols
	idxrow= mod(rank,rowproc)+1
	
	localC = cmplx(0.0,0.0)
	
	if(process_grid(idxrow,idxcol) .gt. 0  ) then 

	localC(:descriptorC(1),&
	       :descriptorC(2))&  
				= & 
	localC(:descriptorA(1,switch),&
	       :descriptorB(2,switch))&  
			+ matmul(& 
	localA(:descriptorA(1,switch),:descriptorA(2,switch),switch),& 
	localB(:descriptorB(1,switch),:descriptorB(2,switch),switch))
	
	!if(descriptorA(2,switch) .ne. descriptorB(1,switch) ) then
	!	write(*,*)idxrow, idxcol,process_grid(idxrow,idxcol),"|" ,descriptorA(3:,switch) ,"|", descriptorB(3:,switch)
	!	end if
	
	!print*,idxrow, idxcol, "|",  descriptorA(4,switch),descriptorB(3,switch) & 
	!	     		,  descriptorA(4,switch)+descriptorA(2,switch) -1 & 
	!	     		,  descriptorB(3,switch) + descriptorB(1,switch)-1, "|"& 
	!	     		,  descriptorA(3,switch) , descriptorB(4,switch), "|" 
		      
	!print*, "r",idxrow, idxcol,process_grid(idxrow,idxcol),& 
				!sum(matmul(& 
			!localA(:descriptorA(1,switch),:descriptorA(2,switch),switch),& 
			!localB(:descriptorB(1,switch),:descriptorB(2,switch),switch))), & 
			!sum(localA(:descriptorA(1,switch),:descriptorA(2,switch),switch)),& 
			!sum(localB(:descriptorB(1,switch),:descriptorB(2,switch),switch))	
			!descriptorA(1,switch)-1,& 
			!	descriptorB(2,switch)-1,& 
			!	descriptorA(2,switch)-1 
			
	end if
	
	
	
	
	
	do ii = 2,pace
		
		!Step UPDATE-------------------------------
		old_process_grid = process_grid
		
		sa = .false.
		sb = .false. 
		ra = .false. 
		rb = .false.
		
		
		call update_process_grid(process_grid, rowproc, colproc,rank);
		
		!print*, idxrow, idxcol, 1,  'a |', descriptorA(:,switch)
		!print*, idxrow, idxcol, 1,  'b |', descriptorB(:,switch)
		
		if(rank.eq.0) then
		call printmatrix(process_grid(:rowproc, :colproc))
		end if
		
		call MPI_BARRIER(MPI_COMM_WORLD, ierr)
		!FIND COOP and SEND
		!############################################################################################################!
		!send B -------------------------------------------------------------------------------------------------------
		call find_coop(old_process_grid, process_grid, idxrow, idxcol, withwhorow, withwhocol, rowproc, colproc, "sb")
		if ((withwhorow .gt. 0 ).and. (withwhocol .gt. 0 )) then
			if((withwhorow .eq. idxrow) .and. (withwhocol .eq. idxcol)) then 
				LocalB(:,:,mod(switch,2)+1) = LocalB(:,:,switch)
				descriptorB(:,mod(switch,2)+1) = descriptorB(:,switch)
				!print*, "Already here sendB"
			else 
			sb = .true.
			call  MPI_comm_wrap_matmul(LocalB(:,:,switch),descriptorB(:,switch),idxrow,idxcol,& 
				withwhorow, withwhocol,rowproc,colproc,'s',1, requestsendB, requestsendDB)
			!print*,'sb', idxrow, idxcol, withwhorow, withwhocol, ii,  '| b |', descriptorB(:,switch)
			end if 
		!	print*, "sending"
		end if
		!############################################################################################################!
		
		
		!############################################################################################################!
		!send A--------------------------------------------------------------------------------------------------------
		call find_coop(old_process_grid, process_grid, idxrow, idxcol, withwhorow, withwhocol, rowproc, colproc, "sa")
		!print*, idxrow, idxcol, "->", withwhorow, withwhocol ,"sa"
		if ((withwhorow .gt. 0 ).and. (withwhocol .gt. 0 )) then
		
			if((withwhorow .eq. idxrow) .and. (withwhocol .eq. idxcol)) then 
				LocalA(:,:,mod(switch,2)+1) = LocalA(:,:,switch)
				descriptorA(:,mod(switch,2)+1) = descriptorA(:,switch)
				!print*, "Already Here send A"
			else 
			sa = .true.
			!print*,'sa', idxrow, idxcol, withwhorow, withwhocol, ii , '| a |', descriptorA(:,switch) 
			call  MPI_comm_wrap_matmul(LocalA(:,:,switch),descriptorA(:,switch),idxrow,idxcol,&
				withwhorow, withwhocol,rowproc,colproc,'s',0,requestsendA, requestsendDA) 
			end if
		end if
		!############################################################################################################!
		
		
		
		!@@@@@@@@@@@@@@@@@@@@@@
		switch = mod(switch,2)+1
		!@@@@@@@@@@@@@@@@@@@@@@
		
		!FIND COOP and RECV
		!############################################################################################################!
		!recv B--------------------------------------------------------------------------------------------------------
		call find_coop(old_process_grid, process_grid, idxrow, idxcol, withwhorow, withwhocol, rowproc, colproc, "rb")
		if ((withwhorow .gt. 0 ).and. (withwhocol .gt. 0 )) then
			if(.not. ( (withwhorow .eq. idxrow) .and. (withwhocol .eq. idxcol) ) ) then 
				rb = .true.
				
				call MPI_comm_wrap_matmul(LocalB(:,:,switch),descriptorB(:,switch),idxrow,idxcol,& 
					withwhorow, withwhocol,rowproc,colproc,'r',1, requestrecvB, requestrecvDB)
				
				!print*,'rb',withwhorow, withwhocol, idxrow, idxcol,  ii ,  '| b |'!, descriptorB(:,switch)
		
			end if
		end if
		!############################################################################################################!
		
		
		!############################################################################################################!
		!recv A--------------------------------------------------------------------------------------------------------
		call find_coop(old_process_grid, process_grid, idxrow, idxcol, withwhorow, withwhocol, rowproc, colproc, "ra")
		
		if ((withwhorow .gt. 0 ).and. (withwhocol .gt. 0 )) then
			
			if(.not. ( (withwhorow .eq. idxrow) .and. (withwhocol .eq. idxcol) ) ) then 
			ra = .true.
			
			call MPI_comm_wrap_matmul(LocalA(:,:,switch),descriptorA(:,switch),idxrow,idxcol,& 
				withwhorow, withwhocol,rowproc,colproc,'r',0, requestrecvA, requestrecvDA)
			
			!print*,'ra',withwhorow, withwhocol, idxrow, idxcol,  ii, '| a |'!, descriptorA(:,switch) 
	
			
			end if
		end if
		!############################################################################################################!
		
		
		
		if(sb) then
			call MPI_WAIT(requestsendB, statussendB, ierr )
			call MPI_WAIT(requestsendDB, statussendDB, ierr )
		end if 
		if(sa) then
			call MPI_WAIT(requestsendA, statussendA, ierr )
			call MPI_WAIT(requestsendDA, statussendDA, ierr )
		end if 
		if (ra) then 
			call MPI_WAIT(requestrecvA, statusrecvA, ierr )
			call MPI_WAIT(requestrecvDA, statusrecvDA, ierr )
		end if
		if (rb) then 
			call MPI_WAIT(requestrecvB, statusrecvB, ierr )
			call MPI_WAIT(requestrecvDB, statusrecvDB, ierr )
		end if
		!-------------------------End Step Update ----------------------------------
		!MATMUL
		if(process_grid(idxrow,idxcol) .gt. 0  ) then 
	
		!print*,idxrow, idxcol, "|",  descriptorA(4,switch),descriptorB(3,switch)& 
		!     		,  descriptorA(4,switch)+descriptorA(2,switch) -1 & 
		!     		,  descriptorB(3,switch) + descriptorB(1,switch) -1 ,"|" , & 
		!     		  descriptorA(3,switch) , descriptorB(4,switch) 
		if(descriptorA(2,switch) .ne. descriptorB(1,switch) ) then
		!write(*,*)idxrow, idxcol,process_grid(idxrow,idxcol),ii,"|" ,descriptorA(3:,switch) ,"|", descriptorB(3:,switch)
		end if
			localC(:descriptorC(1),&
	       	       :descriptorC(2))&  
				= & 
			localC(:descriptorC(1),&
	       	       :descriptorC(2))&  
			+ matmul(& 
			localA(:descriptorA(1,switch),:descriptorA(2,switch),switch),& 
			localB(:descriptorB(1,switch),:descriptorB(2,switch),switch))
			!print*, "r",idxrow, idxcol,process_grid(idxrow,idxcol),& 
			!	sum(matmul(& 
			!localA(:descriptorA(1,switch),:descriptorA(2,switch),switch),& 
			!localB(:descriptorB(1,switch),:descriptorB(2,switch),switch))), & 
			!sum(localA(:descriptorA(1,switch),:descriptorA(2,switch),switch)),& 
			!sum(localB(:descriptorB(1,switch),:descriptorB(2,switch),switch))	
			!descriptorA(1,switch)-1,& 
			!	descriptorB(2,switch)-1,& 
			!	descriptorA(2,switch)-1 !
		end if
		
		!--------------------------------------------------------------------------
	end do 
	
	
	
	deallocate(stat)
	
	deallocate(old_process_grid)
	
	


end subroutine	

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!Process #0 contains the global Matrix, request the descriptor from all the others, and according to it 
!send back the correct chunk of the global matrix

subroutine dispatch_global_matrix(A, localA, descriptor, rowproc, colproc,& 
		 idxrow, idxcol,rank,process_grid)
		use mpi
		implicit none
		!Use the descriptor and the process grid slepy dumb idiot.
		complex, dimension(:,:) :: A
		complex, dimension(:,:) :: localA
		integer, dimension(:,:) :: process_grid
		integer, dimension(4) :: descriptor, extdescriptor
		integer :: rowproc, colproc, idxrow, idxcol, rank, dimrow,dimcol
		integer, dimension(MPI_STATUS_SIZE) :: status_
		integer :: limr0, limr1, limc0, limc1, dimr, dimc
		integer :: ierr, ii,jj
		
	
		if(rank .eq. 0) then 
		
		do ii = 1, rowproc
			do jj = 1, colproc
				if(ii*jj .ne. 1) then 
					
					!if(process_grid(ii,jj) .gt. 0) then
					
					
					call MPI_RECV(extdescriptor, 4, MPI_INTEGER, & 
						      grid2rank(ii,jj,rowproc,colproc) & 
						     ,3*(ii*(colproc+1)+jj)+1 , MPI_COMM_WORLD, status_, ierr)
					call mpic_checkerror(ierr)
					!print* , extdescriptor, ii, jj,  process_grid(ii, jj)
					!if(process_grid(ii,jj) .gt. 0) then
					if(product(extdescriptor) .gt. 0 ) then
					!print*, "test g", ii , jj,  A(extdescriptor(3), extdescriptor(4)),extdescriptor(1)*extdescriptor(2) !here the problem
					call MPI_SSEND(A(extdescriptor(3):extdescriptor(3)+extdescriptor(1)-1 ,& 
							extdescriptor(4):extdescriptor(4)+extdescriptor(2)-1) , extdescriptor(1)*extdescriptor(2),& 
						 MPI_COMPLEX, grid2rank(ii,jj,rowproc,colproc) ,& 
						 ii*(colproc+1)+jj , MPI_COMM_WORLD, ierr)
					
					call mpic_checkerror(ierr)
						
					end if
				else 
					localA = A(descriptor(3): descriptor(3)+descriptor(1)-1,& 
						   descriptor(4):descriptor(4)+descriptor(2)-1)
				end if
				
			end do 
		end do 
		else 
		
		!if(process_grid(idxrow,idxcol) .gt. 0) then
		
		
		call MPI_SSEND( descriptor , 4,& 
					 MPI_INTEGER, 0 ,& 
					 3*(idxrow*(colproc+1)+idxcol)+1 , MPI_COMM_WORLD, ierr)
		
		call mpic_checkerror(ierr)
		!print* , descriptor, idxrow, idxcol, process_grid(idxrow, idxcol)
		!if(process_grid(ii,jj) .gt. 0) then
		if(product(descriptor) .gt. 0 ) then
		!if(process_grid(idxrow,idxcol) .gt. 0) then
		call MPI_RECV(localA(:descriptor(1),:descriptor(2)), descriptor(1)*descriptor(2), MPI_COMPLEX, 0 & 
			,idxrow*(colproc+1)+idxcol , MPI_COMM_WORLD, status_, ierr)
		
		call mpic_checkerror(ierr)
		
		
	
		end if

		
		end if
		
		
				
		
	
	
end subroutine
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine GetBackMatrix(localC, descriptor, C, rowproc, colproc,& 
		 idxrow, idxcol,rank,process_grid)
	 
		use mpi
		implicit none
		
		
		complex, dimension(:,:) :: C
		complex, dimension(:,:) :: localC
		complex, dimension(:,:),allocatable :: tempc 
		integer, dimension(:,:) :: process_grid
		complex :: res
		integer, dimension(4) :: descriptor, extdescriptor
		integer :: rowproc, colproc, idxrow, idxcol, rank, dimrow,dimcol, dimArow,dimAcol
		integer, dimension(MPI_STATUS_SIZE) :: status_
		integer :: limr0, limr1, limc0, limc1, dimr, dimc
		integer :: ierr, ii,jj
		
		res = cmplx(0.0,0.0)
		if(rank .eq. 0) then 
		
		do ii = 1, rowproc
			do jj = 1, colproc
				if(ii*jj .ne. 1) then 
					
					
						call MPI_RECV(extdescriptor, 4, MPI_INTEGER, & 
							      grid2rank(ii,jj,rowproc,colproc), & 
							      3*(ii*(colproc+1)+jj) , MPI_COMM_WORLD, status_, ierr)
						
						call mpic_checkerror(ierr)
							   
						allocate (tempc(extdescriptor(1),extdescriptor(2))) 
						tempc = cmplx(0.0,0.0)
						call MPI_RECV(tempc,& 
								extdescriptor(1)*extdescriptor(2), MPI_COMPLEX, & 
								grid2rank(ii,jj,rowproc,colproc), & 
								3*(ii*(colproc+1)+jj)+1 , MPI_COMM_WORLD,status_, ierr)
						
						call mpic_checkerror(ierr)		
						C(extdescriptor(3):extdescriptor(3)+extdescriptor(1)-1 ,& 
								extdescriptor(4):extdescriptor(4)+extdescriptor(2)-1) =& 
								tempc
						
						deallocate (tempc)
						
					
				else 
					
					C(descriptor(3): descriptor(3)+descriptor(1)-1,& 
						   descriptor(4):descriptor(4)+descriptor(2)-1) = localC
					
				end if
				
			end do 
		end do 
	
		else 
		
		
		
			call MPI_SSEND( descriptor , 4,& 
					MPI_INTEGER, 0 ,& 
					3*(idxrow*(colproc+1)+idxcol) , MPI_COMM_WORLD, ierr)
			
			call mpic_checkerror(ierr)
			
			
			call MPI_SSEND(localC(:descriptor(1), :descriptor(2)) , descriptor(1)*descriptor(2),& 
				 MPI_COMPLEX, 0 ,& 
				 3*(idxrow*(colproc+1)+idxcol)+1 , MPI_COMM_WORLD, ierr)
			
			
			
			call mpic_checkerror(ierr)	
		

		
		end if
		
		

end subroutine

	

	
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

end module



program passm
	use cbp_DistributeHandler
	use mpi
	implicit none
	! #### MPI Variables ####
	
	integer :: rank ! the process ID
	integer :: nprocs ! number of processes
	integer , dimension(:,:) , allocatable :: process_grid 
	integer :: rowproc, colproc
	integer :: ierr , idxcol, idxrow, ii, jj,kk
	!Global
	!real , dimension(:,:),  allocatable :: tmpa , tmpb, tmpc
	complex , dimension(:,:) , allocatable :: A, B, C,D
	integer :: dimA(2), dimB(2), dimC(2)
	
	integer , dimension(4,2):: descriptorA(4,2), descriptorB(4,2)
	integer :: descriptorC(4)
	integer :: localdimA(2), localdimB(2), localdimC(2)
	integer :: testr(2), testc(2) , testf(2)
	

	
	
	
	
	!debug
	
	complex , allocatable, dimension(:,:,:) :: localA, localB
	complex , allocatable, dimension(:,:) :: localC
	
	rowproc = 8
	colproc = 8
	
	
	call MPI_INIT(ierr)
	

	
	
	! communicate size i.e number of process, 
	! MPI_COMM_WORLD is not the only communicator in MPI
	
	call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr) 
	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
	
	

	dimA=(/523,200/)
	dimB=(/200,100/)
	dimC=(/dimA(1), dimB(2)/)
	
	if(rank .eq. 0 ) then 
	
	allocate(A(dimA(1), dimA(2)))
	allocate(B(dimB(1), dimB(2)))
	allocate(C(dimC(1), dimC(2)))
	allocate(D(dimC(1), dimC(2)))
	!allocate(tmpa(dimA(1), dimA(2)))
	!allocate(tmpb(dimB(1), dimB(2)))
	!descriptorA = 0
	!descriptorB = 0 
	!descriptorC = 0 
	
	!Fill A, B with random numbers
	call random_number(A%re)
	call random_number(A%im)
	call random_number(B%re)
	call random_number(B%im)
	D = matmul(A,B)
	print*, " "
	!print*, ceiling(real(dimA(1))/real(rowproc)), ceiling(real(dimB(2))
	do ii = 1, rowproc
		do jj = 1,colproc
			do kk = 1, min(rowproc,colproc)
				testr = (/(ii-1)*ceiling(real(dimA(1))/real(rowproc))+1,& 
					  ii*ceiling(real(dimA(1))/real(rowproc))/)
				testc = (/(jj-1)*ceiling(real(dimB(2))/real(colproc))+1,& 
					  jj*ceiling(real(dimB(2))/real(colproc))/)
				testf = (/(kk-1)*ceiling(real(dimB(1))/real(min(colproc,rowproc)))+1,& 
					  kk*ceiling(real(dimB(1))/real(min(colproc,rowproc)))/)
				if(testc(2) .gt. dimB(2)) then
					testc(2) = dimB(2)
				end if
				if(testr(2) .gt. dimA(1)) then 
					testr(2) = dimA(1)
				end if
				if(testf(2) .gt. dimB(1)) then 
					testf(2) = dimB(1)
				end if
				!print*, testr, testc
				print*, ii, jj,kk, sum(matmul(A(testr(1):testr(2),testf(1):testf(2)),& 
							   B(testf(1):testf(2),testc(1):testc(2)) )), & 
						    sum(A(testr(1):testr(2),testf(1):testf(2))), & 
						    sum(B(testf(1):testf(2),testc(1):testc(2)))
						    !testr(2)-testr(1), testc(2)-testc(1), testf(2)-testf(1)!
			end do
		end do 
	end do
	
	end if
	
	idxcol= rank/rowproc+1			!grid assigned by cols
	idxrow= mod(rank,rowproc)+1
	
	!-----------------------------------------------------------------
	!DIVIDE MATRIX IN THE PROCESS GRID
	allocate (process_grid(max(rowproc,colproc), max(rowproc,colproc)))
	
	call init_process_grid(process_grid, rowproc, colproc,rank)
	
	
	if(nprocs .ne. rowproc*colproc) then
		print *, "invalid process grid, fatal error ", nprocs 
		STOP 
	else 
		localdimA(1) = ceiling(real(dimA(1))/real(rowproc))
		localdimB(2) = ceiling(real(dimB(2))/real(colproc))
		localdimC(1) = localdimA(1)
		localdimC(2) = localdimB(2)
		localdimA(2) = ceiling(real(dimA(2))/real(min(colproc,rowproc)))
		localdimB(1) = localdimA(2)
		
		!print*,"        " ,localdimA, "|", localdimB , "|", localdimC 	
		descriptorA(:,1) =  (/localdimA(1), localdimA(2), (idxrow-1)*localdimA(1)+1, & 
				(abs(process_grid(idxrow,idxcol))-1)*localdimA(2)+1 /) 
	 
		descriptorB(:,1) =  (/localdimB(1), localdimB(2),& 
				 (abs(process_grid(idxrow,idxcol))-1)*localdimB(1)+1, (idxcol-1)*localdimB(2)+1 /) 
		
		descriptorC =  (/localdimC(1), localdimC(2), (idxrow-1)*localdimC(1)+1,& 
				 (idxcol-1)*localdimC(2)+1/) 
		
		!print*,idxrow, "|", idxcol,  dimA(2) , descriptorA(4,1)
		!---LIMIT CASES---
		
		if(idxrow .eq. rowproc) then
			!if(descriptorA(1,1) .lt. dimA(1) - descriptorA(3,1) +1) then 
			!	print*, "error chunk too big"
			!	STOP
			!else 
				descriptorA(1,1) = dimA(1) - descriptorA(3,1) +1
				descriptorC(1) = dimC(1) - descriptorC(3) +1
			!end if
		end if
		if(idxcol .eq. colproc) then
			!if(dimB(2) - descriptorB(4,1) +1 .gt. descriptorB(2,1)) then 
			!	print*, "error chunk too big"
			!	STOP
			!else 
				descriptorB(2,1) = dimB(2) - descriptorB(4,1) +1
				descriptorC(2) = dimC(2) -descriptorC(4)+1
			!end if
		
		end if 
		if(abs(process_grid(idxrow,idxcol)) .eq. min(colproc,rowproc)) then 
			!if(descriptorA(2,1) .lt. dimA(2) - descriptorA(4,1) +1) then 
			!	print*, "error chunk too big"
			!	STOP
			!else 
				descriptorA(2,1) = dimA(2) - descriptorA(4,1) +1
				descriptorB(1,1) = dimB(1) - descriptorB(3,1) +1
			!end if
			if(descriptorA(2,1) - descriptorB(1,1) .ne. 0) then
				print*, "error"
			end if 
		end if
		
		if( process_grid(idxrow, idxcol) .le. 0 ) then
			if (rowproc .lt. colproc) then
				descriptorA  = 0 
			else if (colproc .lt. rowproc) then
				descriptorB  = 0
			end if
			
		end if
		!-------------------!
	
	end if
	

	!print*, "dd" , idxrow, idxcol,   "| ", descriptorA(:,1)
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	!###############DEBUG##################
	!if(rank .eq. 0 ) then
	!
	!	do ii = 1, rowproc
	!		do jj = 1, colproc
	!		 B( (abs(process_grid(ii,jj))-1)*localdimB(1)+1,(jj-1)*localdimB(2)+1 )  = cmplx( descriptorB(1,1), descriptorB(2,1) ) 
			! if (process_grid(ii,jj) .gt. 0 ) then 
			
			
	!		 A( (ii-1)*localdimA(1)+1,(abs(process_grid(ii,jj))-1)*localdimA(2)+1 )  = cmplx( descriptorA(1,1), descriptorA(2,1) ) 
			 !!print *, "init",ii, jj,  (ii-1)*localdimA(1)+1, process_grid(ii,jj)*localdimA(2)+1
			 !end if
	!		end do
	!	end do 
	
	!end if
	
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
				
	!#####################################
	!And allocate space for the local matrices
	!Doubled the space, in the 3dim, one for computation and recv, the other queued to send.
	
	allocate(localA(localdimA(1),localdimA(2),2))
	allocate(localB(localdimB(1),localdimB(2),2))
	allocate(localC(localdimC(1),localdimC(2)))
	
	localC = cmplx(0.0,0.0)
	localA = cmplx(0.0,0.0)
	localB = cmplx(0.0,0.0)
	
		!if(process_grid(idxrow,idxcol) .gt. 0) then
		!print*, idxrow, idxcol, '|', descriptorA(1,1), descriptorA(2,1),"-", & 
		!	descriptorB(1,1), descriptorB(2,1), "-", descriptorC(1), descriptorC(2), '|', sum(localC)
		!end if
	
	
	!Dispatch global matrix to processes ---------------------------
	!There is a bug in this function, 
	!I notice that the local matrices are not what i expect from my idea
	! In fact setting their first element (of the local) to a number(setting the global) 
	!
	call dispatch_global_matrix(A, localA(:,:,1), descriptorA(:,1), & 
		rowproc, colproc, idxrow, idxcol, rank, process_grid)
	call dispatch_global_matrix(B, localB(:,:,1), descriptorB(:,1), & 
		 rowproc, colproc, idxrow, idxcol, rank, process_grid)
		 
	
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	
	!if ((process_grid(idxrow,idxcol) .gt. 0) .and. (idxcol .lt. colproc) & 
	!		.and. (idxrow .lt. rowproc) ) then 
	!print*, sum(localA(:,:,1))/(descriptorA(1,1)*descriptorA(2,1)) , idxrow, idxcol
	!end if
	
	call ComputeMatmul(localA, localB, localC, descriptorA, descriptorB , descriptorC, rowproc, colproc, rank, process_grid)
	
	call GetBackMatrix(localC, descriptorC, C, rowproc, colproc,& 
		 			idxrow, idxcol,rank,process_grid)
	
	
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if(rank .eq. 0 ) then
	print*, sum(C-D), sum(C), sum(D)
	
	end if
	
	
				
	 

	
	!Once initialized i can start iterating over pace
	
	!HERE CALL MATMUL 
	
	!if(rank .eq.0) then 
	!call printmatrix(process_grid(1:rowproc, 1:colproc))
	!end if
	
	!if(rank .eq. 0 ) then 
	!call printmatrix(process_grid)
	!end if
	
	!seem that the local arrays are corrupted..
	deallocate (localA) 
	deallocate (localB) 
	deallocate (localC)
	
	
	if(rank .eq. 0 ) then 
	deallocate (A) 
	deallocate (B) 
	deallocate (C)
	deallocate (D)
	end if
	deallocate (process_grid)
	call MPI_FINALIZE(ierr)
end program 
