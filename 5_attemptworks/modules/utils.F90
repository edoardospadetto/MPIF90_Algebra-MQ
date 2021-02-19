module utils

contains

	!average value of single sping given op
	function getavgvalue(state,direction,N) result (avgs)

		implicit none

		character*1:: direction 
		double complex, dimension(:):: state
		double complex, dimension(2,2) :: px, py, pz
		integer:: ii, jj ,kk, testa,testb, N
		double complex, dimension(N) :: avgs 
		
		!print*, norm2(state%re)**2+norm2(state%im)**2
		avgs=dcmplx(0.d0,0.d0)
		
		px = dcmplx(0.d0,0.d0)
		py = dcmplx(0.d0,0.d0)
		pz = dcmplx(0.d0,0.d0) 
		
		px(1,2) = dcmplx(1.d0,0.d0)
		py(1,2) = dcmplx(0.d0,-1.d0)
		pz(1,1) = dcmplx(1.d0,0.d0) 
		
		px(2,1) = dcmplx(1.d0,0.d0)
		py(2,1) = dcmplx(0.d0,1.d0)
		pz(2,2) = dcmplx(-1.d0,0.d0) 
		
		if (mod(size(state, dim=1),2 ) .ne. 0) then
			print*, "error, not multi spin sistem" 
		end if

		!print *, state

		if(direction .eq. 'x') then 
			do ii = 0, size(state, dim=1) -1
				do jj = 0, size(state, dim=1) -1
					do kk = 0, N-1	
						if( (xor(ii,jj) .eq. 0) .xor.  (xor(ii,jj) .eq. 2**kk)  ) then 
							testa=mod(ii/(2**kk),2)+1
							testb=mod(jj/(2**kk),2)+1
							!print*, testa,testb
							avgs(kk+1) = conjg(state(ii+1))*state(jj+1)*px(testa,testb)	+ avgs(kk+1)  	
						end if 
					end do 
				end do 
			end do
		else if (direction .eq. 'y') then
			do ii = 0, size(state, dim=1) -1
				do jj = 0, size(state, dim=1) -1
					do kk = 0, N-1	
						if( (xor(ii,jj) .eq. 0) .xor.  (xor(ii,jj) .eq. 2**kk)  ) then 
							testa=mod(ii/(2**kk),2)+1
							testb=mod(jj/(2**kk),2)+1	
							avgs(kk+1) = conjg(state(ii+1))*state(jj+1)*py(testa,testb)	+ avgs(kk+1)  	
						end if 	
					end do 
				end do 
			end do
		else if (direction .eq. 'z') then
			do ii = 0, size(state, dim=1) -1
				do jj = 0, size(state, dim=1) -1
					do kk = 0, N-1	
						if( (xor(ii,jj) .eq. 0) .xor.  (xor(ii,jj) .eq. 2**kk)  ) then 
							testa=mod(ii/(2**kk),2)+1
							testb=mod(jj/(2**kk),2)+1
							avgs(kk+1) = conjg(state(ii+1))*state(jj+1)*pz(testa,testb)	+ avgs(kk+1)  	
						end if 
					end do 
				end do 
			end do
		else 
			print*, "Error: invalid dir."
	end if
	
	return
	
end function
	
	

end module
