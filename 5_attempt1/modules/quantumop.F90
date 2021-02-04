module quantumop
    use debugqic
    use printer
        !variables of the simulation
        double precision :: xmin,xmax,e
        integer:: dim



    contains
        ! set variables of the simulation
        subroutine setenv(x1,x2,L)
            implicit none
            double precision :: x1,x2
            integer :: L
            xmin = x1
            xmax = x2
            dim  = L
            e = (x2-x1)/real(2*L+1)
        end subroutine setenv
        !###############################################
        ! compute Forward and backward fourier transform
        function zfft(input, N, direction) result(output)
        !###############################################
            implicit none

            integer :: N,direction,ii
            integer*8 :: plan
            double precision :: scale = 1.0
            double complex, dimension(N), intent(INOUT):: input
            double complex, dimension(N)::output,input1
            !print*, "yyy"


            !print*, size(input)
            !input1 = input/sum(input)
            !do ii = 1, n
            !    print*, input(ii)
            !end do



            call breakifn("FFT:  Input must be 1D double complex array", size(shape(input))  .eq. 1, .true. )
            call breakifn("FFT: Output must be 1D double complex array", size(shape(output)) .eq. 1, .true. )
            n = size(input, dim=1)
            call breakifn("FFT: Input and output must have the same dimension", size(input) .eq.size(output))


            if (direction .eq. 1) then
                call dfftw_plan_dft_1d(plan,N,input1,output,1,0)
                scale = dble(1/dble(N))
            else if (direction .eq. -1) then
                call dfftw_plan_dft_1d(plan,N,input1,output,-1,0)
            else
                call breakifn("Invalid Direction of the FFT direct +1, inverse -1", .false., .true.)
            endif
            input1 = input


           call dfftw_execute_dft(plan, input1, output)
          ! call dfftw_execute_dft(plan, input1, output)
           !print*, plan
           call dfftw_destroy_plan(plan)
           !print*, "oo"
           output = output*scale
           !do ii = 1, n
            !   print*, output(ii)
           !end do


        end function

        !###############################################
        ! compute kinetic energy therm
        function K_1d(cnst) result(K)
        !###############################################
            implicit none
            integer :: N,ii

            double  complex :: cnst
            double  complex :: K(2*dim+1,2*dim+1)
            K= (.0,.0)

            do ii = 2, 2*dim
                    K(ii,ii) =(-2.0,.0)
                    K(ii+1,ii)=(1.0,.0)
                    K(ii-1,ii)=(1.0,.0)
            end do
            K(1,1) =(-2.0,.0)
            K(2*dim+1,2*dim+1) =(-2.0,.0)
            K(2,1)=(1.0,.0)
            K(2*dim,2*dim+1)=(1.0,.0)

            K=cnst*K*(1/e**2)
        end function K_1d
        !###############################################
        !compute the potential energy therm
        function x2_1d(cnst) result(x2)
        !###############################################
            implicit none
            integer :: ii
            double complex :: x2(2*dim+1,2*dim+1), cnst
            double complex :: tmp
            x2 = (0.0,0.0)
            do ii = 1, 2*dim+1
                tmp%im=0.0
                tmp%re = e*(ii)+xmin
                x2(ii,ii) = tmp**2
            end do
            x2=cnst*x2
        end function x2_1d

        !###############################################
        !compute and output eigenvectors real part and probabilyty density function
        subroutine eigvecout(x1,x2,hamiltonian,which)
        !###############################################
            use matrixqic
            implicit none
            double complex, dimension(:,:) :: hamiltonian
            character*20 :: numb
            integer :: which, L,ii
            double precision :: x1, x2,eps

            L = size(hamiltonian, 1)
            eps= (x2-x1)/real(L)

            write(numb,"(I3,'_',I4,'_',E8.2)") which,L,x2
            numb = trim(adjustl(numb))
            !compute prob
            hamiltonian(:, which) = hamiltonian(:, which)/sqrt(norm2(hamiltonian(:, which)%re)**2&
                                                        +norm2(hamiltonian(:, which)%im)**2)
            !output data
            do ii = 1,L
                call wddo(trim("./results/real")//trim(numb)//".txt" , x1 + (x2-x1)/real(L)*ii, (hamiltonian(ii,which)%re/sqrt(e)))
                call wddo(trim("./results/prob")//trim(numb)//".txt" , x1 + (x2-x1)/real(L)*ii, &
                        (realpart(hamiltonian(ii,which))**2+imagpart(hamiltonian(ii,which))**2)/e )
            end do
        end subroutine eigvecout

        !###############################################
        ! Compute and output stats about eigvals
        subroutine deltaeigvalhs(x1,x2,eigv,cnst)
        !###############################################
            use matrixqic
            implicit none
            double precision, dimension(:) :: eigv
            integer :: L,ii
            double precision :: x1, x2,eps, dsum, cnst
            character*20 :: numb, nb

            L = size(eigv, 1)
            eps= (x2-x1)/real(L)

            write(numb,"(I4,'_',E8.2)") L,abs(x1)
            numb = trim(adjustl(numb))
            write(nb,"('_',E8.2)") abs(x1)
            nb = trim(adjustl(nb))
            dsum =0.
            !compute s300
            do ii = 1, 300
                dsum= dsum + (eigv(ii)-cnst*(ii-0.5))**2
            end do
                !output s300
                call widdo("./results/score"//trim(nb)//".txt", L, eps, dsum/(300.0)**2)
                dsum =0
                !compute stot (not used in the report)
                do ii = 1, size(eigv)

                    dsum= dsum + (eigv(ii)-cnst*(ii-0.5))**2
                    call widdo("./results/deltaeig"//trim(numb)//".txt", ii, eigv(ii), (eigv(ii)-cnst*(ii-0.5))**2)

                end do
                !compute and output s20
                dsum =0
                do ii = 1, 20

                    dsum= dsum + (eigv(ii)-cnst*(ii-0.5))**2

                end do
                call widdo("./results/score20_"//trim(nb)//".txt", L, eps, dsum/(10.0**2))
        end subroutine deltaeigvalhs


        !##################################################
        function mpediv2(now,cnst,dt, bigT) result(V)
        !###################################################
        implicit none
        ! i want exp(integral(x - t/T)**2)
        !approximate middle point the integral
        ! exp(dt(x - t/T)**2) luckily the matrix is diagonal
        !i got already dim defined from the enviromnet
        !used a vector instead of a mutrix full of zeros
        double complex, dimension(2*dim+1):: V
        double precision:: now,cnst,dt,bigT,x, velo
        integer:: ii
        !if (bigt .le. 0.00001) then
            velo =dble(1)/bigt
        !else
            !velo=0
        !end if

        V=0.0
        do ii = 1, 2*dim+1
            x = xmin+  (xmax-xmin)*ii/(2*dim+1)
            V(ii)= exp(cmplx(0.D0,-0.5*cnst*dt*(x-now*velo)**2))


            !V(ii)%re = cos(x*x*dt+dt*now*((velo)**2)*(now+dt)-(velo)*x*dt*(2*now+dt))
            !V(ii)%im = sin(x*x*dt+dt*now*((velo)**2)*(now+dt)-(velo)*x*dt*(2*now+dt))
        end do


        end function

        function mke(cnst,Ezero,dt) result(Kf)
        !###################################################
        implicit none
        ! i want temporal evolution for the momentum
        !but jnowing that the eigenfucntion of the operator
        !are the harmonics
        ! i can use ftt , and then a diagonal matrix.
        !and the invert the fft
        !remember E = hkk/2m
        !then E0 = h/2m *1*1
        !maybe you have to correct the constants later, for the energy
        !used a vector instead of a matrix full of zeros
        double complex, dimension(2*dim+1):: Kf
        double precision :: Ezero, cnst,dt
        integer:: ii,N
        double precision :: PI=4.D0*DATAN(1.D0), dp
        Kf=0.0
        dp = 2*PI/(xmax-xmin) !maybee
        N=2*dim+1
        do ii = 1, dim

            E=Ezero*(dp*(ii-1))**2
            !write(*,*) E,ii
            Kf(ii) = exp(cmplx(0.D0,-cnst*E*dt))

        end do
        do ii = dim+1, N

            E=(dp*(2*dim+1-ii))**2
            Kf(ii) = exp(cmplx(0.D0,-cnst*E*dt))

        end do


        Kf = Kf

        end function

        !##################################################
        function exseventemporalevolution(vec,bigt,dt,time,w) result(vec1)
        !###################################################
        implicit none
        !the form should be T/2 fft-1 K fft T/2 psi
        double complex, dimension(2*dim+1) :: vec,vec1
        double complex, dimension(2*dim+1):: Vte, Kte
        double precision :: bigt,dt,time,w
        integer:: ii
        !compute K
        !mke(cnst,Ezero,dt)
        Kte = mke(dble(0.5),dble(0.5) ,dt)
        !compute T/2
        ! mpediv2(now,cnst,dt, bigT)
        Vte =mpediv2(time,dble(w*w*.5), dt, bigt)
        vec1=0.0
        !print*, time
        vec1=vec

        do ii=1,2*dim+1
            vec1(ii) = vec(ii)*Vte(ii)
        end do
        !print*, size(vec1)

        vec1=vec
        vec1= zfft(vec1, 2*dim+1, 1)
        do ii=1,2*dim+1
            vec1(ii) = Kte(ii)*vec1(ii)            !print*, vec1(ii)
        end do
        vec1= zfft(vec1, 2*dim+1, -1)




        do ii=1,2*dim+1
            vec1(ii) = Vte(ii)*vec1(ii)!print*, vec1(ii)
        end do



        !zfft(input, N, direction) result(output)







        end function



end module quantumop
