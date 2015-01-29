module init

implicit none

contains

subroutine initialize()
    ! initialize mpi, fft's, vertical discretization
    ! and run preprocessor to find w+, w_ for w = w0 + b_w_ + b+w+
    use mpicom
    use dim
    use flags
    use runparam
    use parameters
    use time
    use iofiles
    use formatstrings
    use derivatives
    use modhorfft
    use modfftw
    use core
    use tavg
    use paral
    use io
    use parallel_io
    use grid
    use filters
    use dealias
    ! local vars
    integer :: AllocateStatus
    real :: cput1, cput2
    ! measure cpu time
    cput1 = mpi_wtime()

    ! calculate all useful dimensions
    call init_dim(nproch,nprocv)
    ! debug: make sure all procs got inputs
    if (debug>=2) write(*,*) 'myid=',myid,'nxpp=',nxpp

    ! print to screen useful information about run
    if (debug>=2 .and. myid==0) then
        write(*,*) 'resolution: nx, ny, nz =', nx, ny, nz
        write(*,*) 'num MPI procs =',numprocs
        write(*,*) 'domain decomposition in phys space:'
        write(*,*) 'nx, nypl, nzp = ', nx, nypl, nzp
        write(*,*) 'domain decomposition in fourier space:'
        write(*,*) 'nxpl, ny, nzp = ', nxpl, ny, nzp
    end if

    ! set formatstrings for debug output
    if (debug>=1) call set_formatstrings()

    ! setup mpi groups, comms, and data type
    call setup_mpi()
    ! setup all necessary qties for FFTs
    call setup_fftw(1)
    ! setup all necessary qties for dealiasing
    if (dealiasflag==1) call setup_dealiasing(1)


    ! allocate tseries
    if (.not. allocated(tseries)) allocate(tseries(nsteps))

    ! allocate arrays for solver
    ! a) physical space arrays
    if (.not. allocated(ubar)) allocate(ubar(nzp), ubarp(nzp), ubarpp(nzp), wbar(nzp), wbarp(nzp), wbarpp(nzp), tbar(nzp), tbarp(nzp), tbarpp(nzp), stat=AllocateStatus)
    if (AllocateStatus /= 0) write(*,*) "**Not Enough Memory - Phys. Space 3D arrays allocation **"

    ! b) fourier space arrays
    if (.not. allocated(uf)) allocate(uf(nxpl,ny,nzp), vf(nxpl,ny,nzp), wf(nxpl,ny,nzp),  unif(nxpl,ny,nzp,3), unmif(nxpl,ny,nzp,3), nlnif(nxpl,ny,nzp,3), nlnmif(nxpl,ny,nzp,3), pf(nxpl,ny,nzp), pnf(nxpl,ny,nzp), tempf(nxpl,ny,nzp), tnf(nxpl,ny,nzp), tnmf(nxpl,ny,nzp), nltnf(nxpl,ny,nzp), nltnmf(nxpl,ny,nzp), stat=AllocateStatus)
    if (AllocateStatus /= 0) write(*,*) "**Not Enough Memory - Fourier Space 3D arrays allocation **"

    its = 0

    if (nstart==1) then ! start from restart file
        ! read in all values for important quantities
!         call restart()
        call read_grid()
        call read_paral()

        if ( io_mpi == 1 ) then
            call mpi_read_restart_file()
        else 
            call read_restart_file()
        end if
        its = isum
        if (myid==0 .and. debug>=1) write(*,'(A,G12.5)') 'restart dt =',dt
        ! setup vertical derivative matrix
        call setup_zder()
    
    else ! starting from scratch
        ! setup grid based on input parameters
        call setup_grid()
        call write_grid() ! write grid to file

        ! setup vertical derivative matrix
        call setup_zder()

        ! set initial conditions for ubar, tbar, u, v, w, temp
        call setup_ic()
        call write_paral() ! write mean parallel flow qties to file

        ! set timestep dt based on input cfl, grid, and initial conditions
        call setdt(uf,vf,wf,cfl)
        if (myid==0 .and. debug>=1) write(*,*) 'based input cfl, grid, and initial dt =', dt
    end if

    if (linstab>=1) call setup_linstab()


    ! allocate avg arrays
    if ( otfavg==1 .and. istartavg<isum+nsteps) then
        if ( .not. allocated(uf_tavg) ) allocate(uf_tavg(nxpl,ny,nzp), vf_tavg(nxpl,ny,nzp), wf_tavg(nxpl,ny,nzp), tempf_tavg(nxpl,ny,nzp), pf_tavg(nxpl,ny,nzp), stat=AllocateStatus)
        if (AllocateStatus /= 0) then
            write(*,*) "**Not Enough Memory - time average arrays allocation **"
        end if
        ! zero avg arrays
        uf_tavg = 0.; vf_tavg = 0.; wf_tavg = 0.; tempf_tavg = 0.; pf_tavg = 0.;
    end if

    ! measure cpu time
    cput2 = mpi_wtime()

    if (myid==0) write(*,'(//,A,G14.6)') 'time reqd to initialize =',cput2-cput1


end subroutine initialize


subroutine setup_grid()
    use dim
    use flags
    use formatstrings
    use mpicom, only: myid
    use runparam
    use parameters
    use grid
    use io
    ! vars
    integer :: i,j,k
    real :: pi,twopi
    if (.not. allocated(xpts)) allocate(xpts(nx),ypts(ny),y(nypl),zpts(nzp),wavx(nxpl),wavy(ny),wavxx(nxpl),wavyy(ny),gf(nzp),gf2(nzp),xi(nzp))
    ! define some constants
    pi = 4.*atan(1.)
    twopi = 8.*atan(1.)
    ! change parameters for linear stability test
    if (linstab>=1) then
        xlen = twopi*6.02*sqrt(xnu*x0/u0)/astar
        ylen = xlen
    end if

    ! make horizontal grids
    dx = xlen/(nx - 1.)
    do i=1,nx
        xpts(i) = x0 + (i-1)*dx
    end do

    dy = ylen/(ny - 1.)
    do j=1,ny
        ypts(j) = (j-1)*dy
    end do
    do j=1,nypl
        y(j) = (myid*nypl + j-1)*dy
    end do

    ! set up wave numbers
    alpha = twopi/xlen
    do i=1,nxhpl
        wavx(2*i-1) = (myid*nxhpl + i - 1.)*alpha
    end do
    do i=2,nxpl,2
        wavx(i) = wavx(i-1)
    end do
    wavxx(:) = wavx(:)*wavx(:)

    beta = twopi/ylen
    do j=1,nyh
        wavy(ny-j+1) = -j*beta
        wavy(j) = (j-1.)*beta
    end do
    do j=1,ny
        wavyy(j) = wavy(j)*wavy(j)
    end do

    ! set up vertical grid
    if (debug>=2 .and. myid==0) write(*,*) 'zlen=', zlen, 'mapping=',mapping

    select case(mapping)
        case(0) ! channel flow mapping where zlen is channel width (usually 1)
            xl = zlen
            do k=1,nzp
                xi(k) = xl/2.0*cos(pi*(k-1.)/nz) ! xi = [-xl/2,xl/2]
                zpts(k) = xi(k) ! gf = 1/(dz/dxi)
                gf(k) = 1.
            end do

        case(1) ! semi-infinite domain for boundary layer flow
            xl = 1.
            aa = 1.
            bb = 1.875*sqrt(xnu*x0/u0)
!             bb = 1.875*sqrt( xnu*( x0 + 50.0*xlen )/u0 )
            do k=2,nzp
                xi(k) = xl/2.0*( 1. + cos(pi*(k-1.)/nz) ) ! xi = [0,xl]
                zpts(k) = bb*xi(k)/(aa - xi(k))
                ! gf = (dz/dxi)^-1
                gf(k) = (aa - xi(k))**2.0/(bb*aa) !bb*aa/(zpts(k) + bb)**2.0 
                gf2(k) = (aa - xi(k))**3./(2.0*bb*aa)
            end do
            xi(1) = 1.
            zpts(1) = 1.E+20 ! approximating infinity (can be increased)
            zlen = zpts(1)
!             gf(1) = bb*aa/(zpts(1) + bb)**2.0
            gf(1) = 0.0! 1/infty = 0

        case(2) ! finite algebraic mapping for boundary layer flow
            xl = 1.
            ! larger bb -> less pts in BL at inflow
            ! it's important to have it large enough for outflow BL to be well-resolved
            bb = 1.875*sqrt( xnu*( x0 + blfac*xlen )/u0 )
!             bb = 1.875*sqrt( xnu*( x0 + xlen )/u0 ) 
            aa = 1. + bb/zlen
            do k=1,nzp
                xi(k) = xl/2.0*( 1. + cos(pi*(k-1.)/real(nz)) ) ! xi = [0,xl]
                zpts(k) = bb*xi(k)/( aa - xi(k) )
                ! gf = (dz/dxi)^-1
                gf(k) = (aa - xi(k))**2.0/(aa*bb) !bb*aa/(zpts(k) + bb)**2.0
                gf2(k) = (aa - xi(k))**3./(2.0*aa*bb)
    !                 if (cheb==0) gf(k) = (aa - xi(k))**2.0/( bb*0.5*pi/xl*sin(pi*real(k-1.)/real(nz)) )
            end do

        case(3) ! channel flow from z=[0,zlen]
            xl = zlen
            do k=1,nzp
                xi(k) =  xl/2.0*cos(pi*(k-1.)/real(nz))  ! xi = [0,xl]
                zpts(k) = xl/2.0 + xi(k)
                gf(k) = 1.
            end do

!         case(4) ! finite algebraic mapping for boundary layer flow
!             xl = 1.
!             ! larger bb -> less pts in BL at inflow
!             ! it's important to have it large enough for outflow BL to be well-resolved
!             bb = 1.875*sqrt( xnu*( x0 + 8.0*xlen )/u0 )
! !             bb = 1.875*sqrt( xnu*( x0 + xlen )/u0 ) 
!             aa = 1. + bb/zlen
!             do k=1,nzp
!                 xi(k) = xl/2.0*( 1. + cos(pi*(k-1.)/real(nz)) ) ! xi = [0,xl]
!                 zpts(k) = bb*xi(k)/( aa - xi(k) )
!                 ! gf = (dz/dxi)^-1
!                 gf(k) = (aa - xi(k))**2.0/(aa*bb) !bb*aa/(zpts(k) + bb)**2.0
!                 gf2(k) = (aa - xi(k))**3./(2.0*aa*bb)
!             end do

    end select

    if (debug>=2) then
        if (myid==0) write(*,*) 'twopi=',twopi,', xlen=',xlen,', alpha=',alpha
        write(*,*) 'myid =',myid,', wavx ='
        write(*,*) wavx(1:nxpl)
    end if

    if (debug>=1 .and. myid==0) then
        write(*,*) 'xpts ='
        write(*,fm1) xpts(1:nx)
        write(*,*) 'ypts ='
        write(*,fm2) ypts(1:ny)
        write(*,*) 'aa, bb'
        write(*,*) aa, bb
        write(*,*) 'zpts, gf'
        write(*,'(2(G12.4,'',''))') (zpts(k),gf(k),k=1,nzp)
    end if

end subroutine setup_grid


subroutine setup_zder()
    use flags, only: debug
    use dim
    use chebyshev
    use mats
    use helmholtz

    ! allocate zder matrix
    if (.not. allocated(d)) allocate(dor(nzp,nzp),d(nzp,nzp),d2(nzp,nzp))
    ! use Chebyshev discretization in z 
    call build_cheb_zder_matrix(dor,d,d2)
    
    ! diagonalize derivative matrix for helmholtz & poisson solvers
    ! only interior points for Dirichlet boundary conditions
    if (.not. allocated(diag)) allocate( diag(nzm), smat(nzm,nzm), simat(nzm,nzm) )
    call diagonalize_matrix(d2(2:nz,2:nz),nzm,diag,smat,simat)

end subroutine setup_zder


subroutine setup_ic()
    ! set initial conditions for u,v,w, ubar, and tbar
    use dim
    use mpicom, only: myid
    use flags, only: debug
    use runparam
    use parameters
    use grid
    use formatstrings
    use paral
    use core
    use modhorfft
    use mats
    use time, only: dt
    use derivatives
    use prng
    ! local vars
    integer :: AllocateStatus
    integer :: i, k, n, iter , j
    real :: pi, rot
    real :: eps, okr, dokr
    real :: amp, alph
    real, allocatable, dimension(:) :: blasfl_zpts
    real, allocatable, dimension(:,:) :: blasfl_dor, blasfl_d, blasfl_d2

    real, allocatable, dimension(:,:,:) :: u, v, w, temp

    allocate( u(nxpp,nypl,nzp), v(nxpp,nypl,nzp), w(nxpp,nypl,nzp), temp(nxpp,nypl,nzp), stat=AllocateStatus )
    if (AllocateStatus /= 0) then
        write(*,*) "**Not Enough Memory - setup_ic arrays allocation **"
    end if

    ! set some useful parameters using inputs
    xkap=xnu/prtl
    pi=4.*atan(1.)


    select case(ic)
        case(0) 
        ! 0: ubar=blasius, tbar=thigh - gradt*zpts(k), u=v=w=T=0 
            rot=sqrt(tayl)*xnu/(2.0*zlen**2) ! leave Taylor number = 0 in input, rotation not implemented
            rotx=0.0
            roty=0.0
            rotz=rot
            gradt = (thigh-tlow)/zlen 
            re = x0*u0/xnu
!             do k = 1, nzp
!                 tbar(k) = thigh - gradt*zpts(k)
! !                 tbarp(k) = gradt
!                 ! test temp initial condition
! !                 temp(:,:,k) = 0.5*(thigh-tlow)*sin(2.0*pi*zpts(k)/zlen)*exp(-zpts(k))
!             end do
!             call dpdz(1,1,nzp,tbar,tbarp,1)
!             call dpdz(1,1,nzp,tbar,tbarpp,2)
!             do i = 1, nx
!                 temp(i,:,:) = temp(i,:,:)*sin(10.0*pi*xpts(i)/xlen)
!             end do
            tbar = 0.0
            tbarp = 0.0
            tbarpp = 0.0

            ! use  Von Karman-Pohlhausen approximate Blasius soln
            ! for ubar - only use if blasfl won't converge
            dokr = 4.64 ! constant set for Blasius
            okr = dokr*sqrt(x0*xnu/u0)
            ubar = 1.
            ubarp = 0.0
            ubarpp = 0.0
            wbar = 3./16.*dokr
            do k=1,nzp
                eps = zpts(k)/okr
                if (eps<=1.) then
                    ubar(k) = 1.5*eps - 0.5*eps**3.
                    ubarp(k) = 1.5 - 1.5*eps**2.0
                    ubarpp(k) = -3.*eps
    !                     eps2 = eps*4.64
                    wbar(k) = 0.5*dokr*(0.75*eps**2.0 - 3./8.*eps**4.)
                end if
            end do
                
            ! parameters to solve for Blasius flow
            eps = 1.E-11
            iter = 100
            alph = 0.33!*0.001/xnu
            okr = sqrt(x0*xnu/u0) ! non-dimensionalization length
            allocate( blasfl_zpts(nzp), blasfl_dor(nzp,nzp), blasfl_d(nzp,nzp), blasfl_d2(nzp,nzp) )
            ! non-dimensionalize vertical grid and its derivative matrices
            blasfl_zpts = zpts/okr
            blasfl_dor = dor
            blasfl_d = d*okr
            blasfl_d2 = matmul(blasfl_d, blasfl_d)
!                 if (mapping==2) then
!                     do k=1,nzp
!                         blasfl_zpts(k) = bb/okr*xi(k)/(1.+bb/zlen/okr - xi(k))
!                         blasfl_d(k,:) = dor(k,:)*( 1. + bb/zlen/okr - xi(k) )**2.0/( bb/okr*(1. + bb/zlen/okr) )
!                     end do
!                     blasfl_dor = dor
!                     blasfl_d2 = matmul(blasfl_d, blasfl_d)
!                 end if

            if (debug>=2 .and. myid==0) then
                write(*,'(//,A)') 'blasfl_zpts='
                write(*,'(G12.4,'','')') blasfl_zpts(1:nzp)
            end if

            ! account for small ceilings
    !             if (mapping>=3 .and. zpts(1) < 1.E3) alph = alph*10.0
            ! solve for Blasius flow
            call blasfl(eps,iter,alph,nzp,blasfl_zpts,blasfl_dor,blasfl_d,blasfl_d2,ubar,wbar)
            deallocate(blasfl_zpts, blasfl_d, blasfl_d2)

            ! dimensionalize results
            ubar = u0*ubar
!             ubar(nzp) = 0.0! velocity at wall must be zero
            wbar = wbar*sqrt(xnu*u0/x0)
            wbar(1) = wbar(2)
            wbar(nzp) = 0.0
            !             if (mapping>=3) wbar = wbar/1.875/okr
    !             if (cheb==1) then
!                 ubarp = matmul(d,ubar)
    !                 ubarpp = matmul(d,ubarp)
!                 ubarpp = matmul(d2,ubar)
!                 wbarp = matmul(d,wbar)
!                 wbarpp = matmul(d2,wbar)
    !             end if
            call dpdz(1,1,nzp,ubar,ubarp,1)
            call dpdz(1,1,nzp,ubar,ubarpp,2)
            call dpdz(1,1,nzp,wbar,wbarp,1)
            call dpdz(1,1,nzp,wbar,wbarpp,2)
!             
            u = 0.0
            v = 0.0
            w = 0.0
            temp = 0.0
            ! random initial perturbations for T
!             amp = 0.01
!             n = nx*nypl*nzp
!             call init_random_seed(n)
!             call random_number(u(:,:,:))
!             call init_random_seed(n)
!             call random_number(v(:,:,:))
!             do i=1,nx-1,2
!                 temp(i,:,:) = amp*sqrt( -2.0*log(u(i,:,:)) )*cos( 2.0*pi*u(i,:,:) )
!                 temp(i+1,:,:) = amp*sqrt( -2.0*log(u(i,:,:)) )*sin( 2.0*pi*v(i,:,:) )
!             end do
        

        case(1) 
        ! 1: make your own initial conditions for ubar, tbar, u, v, w, temp in physical space
            ! e.g. you could copy the above and uncomment the random temp perturbations

    end select

    ! set dt based on input cfl
    ! call maxcfl / setdt
    dt = cfl*dx/max(maxval(abs(ubar+maxval(abs(u)))),1.)

    ! print initial conditions for debugging
    if (debug>=1 .and. myid==0) then
        write(*,'(//,A)') '-----------------initial conditions----------------'
        write(*,*) 'dt=',dt
        write(*,'(//,A)') 'zpts, ubar, ubarp, ubarpp, wbar, tbar ='
        write(*,'(6(G12.4,'',''))') (zpts(k),ubar(k),ubarp(k),ubarpp(k), wbar(k), tbar(k),k=1,nzp)
        if (debug>=2) then
            write(*,'(//,A)') 'u(x,z) initial ='
            write(*,fm1) (u(1:nx,1,k),k=1,nzp,4)
            write(*,'(//,A)') 'v(x,z) initial ='
            write(*,fm1) (v(1:nx,1,k),k=1,nzp,4)
            write(*,'(//,A)') 'w(x,z) initial ='
            write(*,fm1) (w(1:nx,1,k),k=1,nzp,4)
            write(*,'(//,A)') 'temp(x,z) initial ='
            write(*,fm1) (temp(1:nx,1,k),k=1,nzp,4)
        end if
    end if

    ! transform initial conditions to spectral space
    ! so that it is ready for solver
    call horfft(u,uf,-1)
    call horfft(v,vf,-1)
    call horfft(w,wf,-1)
    call horfft(temp,tempf,-1)
    unif = 0.0
    unmif = 0.0
    nlnif = 0.0 
    nlnmif = 0.0
    tnf = 0.0
    tnmf = 0.0
    nltnf = 0.0
    nltnmf = 0.0
    pf = 0.0
    pnf = 0.0

    deallocate( u, v, w, temp, stat=AllocateStatus )
    if (AllocateStatus /= 0) then
        write(*,*) "**Error setup_ic arrays deallocation **"
    end if


end subroutine setup_ic


subroutine blasfl(eps,iter,alph,nzp,zpts,dor,d,d2,ubar,wbar)
    ! computes solution to Blasius flow ubar, wbar given:
    ! eps: tolerance, 
    ! iter: max num of iterations,
    ! alph: [0.33, 10.] for zlen = [infty, O(10)]
    ! nzp: num of pts in vertical grid
    ! zpts: vertical grid points, with clustering of points near wall
    ! dor: unmapped vertical derivative matrix d/dz
    ! d: mapped vertical derivative matrix d/dz
    ! d2: mapped second derivative matrix d^2/dz^2
    use mpicom, only: myid
    ! i/o
    real, intent(in) :: eps, alph
    integer, intent(in) :: iter, nzp
    real, dimension(nzp), intent(in) :: zpts
    real, dimension(nzp,nzp), intent(in) :: dor, d, d2
    real, dimension(nzp), intent(out) :: ubar, wbar
    ! local vars
    integer :: ii, i, j, AllocateStatus, info, nzpp
    integer, allocatable, dimension(:) :: ipiv
    real, allocatable, dimension(:,:) :: d3, A
    real, allocatable, dimension(:) :: f, df, rhs, temp1, temp2, rat
    real :: err
    nzpp = nzp+1
    allocate( A(nzpp,nzpp), d3(nzp,nzp), f(nzp), df(nzpp), rhs(nzp), temp1(nzp), temp2(nzp), rat(nzp), ipiv(nzpp), stat=AllocateStatus )
    if (AllocateStatus /= 0) then
        write(*,*) "**Not Enough Memory - blasfl array allocation **"
    end if

    ! initial guess
    f = ( -1. + exp(-zpts(1:nzp)*alph) )/alph
    ! calculate third derivative matrix
    d3 = matmul(d2,d)

    if (myid==0) write(*,'(//,A,I5)') 'max iter in blasfl:',iter

    ! iterative loop to get better solution 
    ! to parallel boundary layer equation
    rat = 1./( 1. + zpts(1:nzp) )
    do ii=1,iter
        ! calc RHS = d^2f/dz^2
        rhs = matmul(d2,f)
        ! calc LHS matrix = 0.5 f*d^2/dz^2 + d^3/dz^3 + 0.5 eta d^2/dz^2
        do i=1,nzp
            do j=1,nzp
                A(i,j) = 0.5*rat(i)*f(i)*d2(i,j) + rat(i)*d3(i,j) + 0.5*(zpts(i)*rat(i))*d2(i,j)
            end do
            A(i,i) = A(i,i) + 0.5*rat(i)*rhs(i)
        end do
        ! take care of all boundaries of matrix A
        A(1:nzp,nzpp) = d(1:nzp,nzp)
        A(1,1:nzp) = dor(1,1:nzp)
        A(nzp,1:nzp) = 0.0
        A(nzpp,1:nzp) = d(nzp,1:nzp)
        A(1,nzpp) = 0.0
        A(nzp,nzpp) = 0.0
        A(nzpp,nzpp) = 0.0
        A(nzp,nzp) = 1.0
        temp1 = matmul(d3,f)
        temp2 = matmul(d2,f)
        df(1:nzp) = -0.5*rat(:)*f(:)*temp2(:) - 0.5*zpts(:)*rat(:)*temp2(:) - rat(:)*temp1(:)
        df(1) = 0.0
        df(nzp) = 0.0
        df(nzpp) = 0.0
        ! computes the solution to a real system of linear equations A * X = B
        call dgesv( nzpp, 1, A, nzpp, ipiv, df, nzpp, info )
        ! add solution df to f
        f = f + df(1:nzp)
        ! calculate 1-norm of iteration df
        err = sum(abs(df(1:nzp)))
        ! check convergence
        if (err<=eps) exit
    end do
    ! report convergence or write(*,*) program if not convergent
    if (ii < iter) then
        if (myid==0) write(*,'(A,I5)') 'blasfl: converged at iter=',ii
    else
        write(*,*) '**error: blasfl convergence no reached**'
    end if

    ! create ubar and wbar based on solution f
    !     if (mapping==1 .or. mapping==2) then
    ubar = 1.0 + matmul(d,f)
    wbar = 0.5*( zpts(1:nzp)*ubar(1:nzp) - (f(1:nzp) + zpts(1:nzp)) )
    !     else
    !         ubar =1. + 0.5*matmul(d,f)
    !         wbar = 0.5*( zpts(1:nzp)*ubar(1:nzp) - (0.5*f(1:nzp) + zpts(1:nzp)) )
    !     end if

end subroutine blasfl


subroutine setup_linstab()
    ! solve for most unstable Orr-Sommerfeld mode and inject it into initial fields
    use dim
    use core
    use runparam, only: nsteps
    use parameters, only: u0, cguess, xnu, x0
    use grid, only: alpha, beta, wavx
    use orrsomm
    use mpicom, only: myid
    use flags, only: debug
    use io
    use paral, only: ubarp, ubarpp ! to remove
    ! vars for stability test
    integer :: iglb
    real :: alphax, betax, eps, reylin
    ! parameters for orr-sommerfeld solver
    eps = 1.E-10 ! normalized mode energy
    iglb = -1 ! -1: all eigenvalues, 0: only the most unstable
    alphax = alpha ! wavelength of domain
    betax = 0.0! zero for 2D problems
    omega = alphax*cguess*u0 ! guess for the most unstable mode (freq)
    reylin = (6.02*sqrt(xnu*x0/u0)*u0)/xnu ! Reynolds number based on BL thickness
    if (debug>=1 .and. myid==0) then
        write(*,*)'alphax=',alphax
        write(*,*) 'omega_guess=',omega
        write(*,*) 'd*=',1.72*sqrt(xnu*x0/u0)
        write(*,*) 'Re_d_99.9=',reylin
    end if
!     allocate( uvec(nzpp2), vvec(nzpp2), wvec(nzpp2) )
    if (.not. allocated(uvec)) allocate( uvec(nzp), vvec(nzp), wvec(nzp) )
    call solve_orrsomm(uvec,vvec,wvec,eps,alphax,betax,omega,iglb)
    if (debug>=1 .and. myid==0) write(*,*) 'most unstable omega:',AIMAG(omega)
    
    ! add most unstable mode to initial velocities
!     if (wavx(3)==2.0*alpha) then
    if (myid==0) then
!         uf(3,1,:) = 3.E-6*ubarpp(:)
!         uf(4,1,:) = 1.E-5*ubarpp(:)
!         wf(3,1,:) = -1.E-6*alpha*ubarp(:)
!         wf(4,1,:) = 1.E-6*alpha*ubarp(:)
        uf(3,1,:) = uf(3,1,:) + real(uvec(:))
        uf(4,1,:) = uf(4,1,:) + aimag(uvec(:))
        vf(3,1,:) = vf(3,1,:) + real(vvec(:))
        vf(4,1,:) = vf(4,1,:) + aimag(vvec(:))
        wf(3,1,:) = wf(3,1,:) + real(wvec(:))
        wf(4,1,:) = wf(4,1,:) + aimag(wvec(:))
    endif
    ! print perturbation for debugging
    if (debug>=1) then
        if (myid==0) write(*,'(//,A)') '-------- orr-sommerfeld perturbation --------'
        call printuvw(uf,vf,wf)
    end if

    deallocate( uvec, vvec, wvec )

end subroutine setup_linstab


subroutine setdt(uf,vf,wf,cfl)
    ! set timestep based on desired cfl
    use dim
    use modhorfft
    use paral
    use grid
    use time
    ! i/o
    real, dimension(nxpl,ny,nzp), intent(inout) :: uf,vf,wf
    real, intent(in) :: cfl
    ! local vars
    integer :: k, AllocateStatus
    real, allocatable, dimension(:,:,:) :: u, v, w
    real :: dtx, dty, dtz
    ! allocate vars
    allocate( u(nxpp,nypl,nzp), v(nxpp,nypl,nzp), w(nxpp,nypl,nzp), stat=AllocateStatus )
    if (AllocateStatus /= 0) then
        write(*,*) "**Not Enough Memory - setdt array allocation **"
    end if
    ! transfer vel to phys space
    call horfft(u,uf,1)
    call horfft(v,vf,1)
    call horfft(w,wf,1)
    ! add mean flow
    do k=1,nzp
        u(:,:,k) = u(:,:,k) + ubar(k)
    end do
    ! transform w into dz/abs(w)
    do k=1,nz
        w(:,:,k) = (zpts(k) - zpts(k+1))/abs(w(:,:,k))
    end do
    w(:,:,nzp) = 0.0
    dtx = cfl*dx/maxval(abs(u))
    dty = cfl*dy/maxval(abs(v))
    dtz = cfl*maxval(abs(w))

    dt = min( dtx, dty, dtz, 1. )
    dt2 = dt
    dt3 = dt

    deallocate( u, v, w )

end subroutine setdt


subroutine finalize()
    use core
    use tavg
    use flags, only: nstart,debug
    use runparam, only: dealiasflag
    use dealias
    use modfftw
    use mats
    use paral
    use helmholtz
    use grid
    use time
    use mpicom
    use iofiles
    integer :: AllocateStatus

    if (allocated(tseries)) deallocate(tseries)
    call setup_fftw(0)
    if (dealiasflag==1) call setup_dealiasing(0)
    ! deallocate all core arrays
!     if (allocated(uf)) deallocate(uf,vf,wf,tempf,pf,fxf,fyf,fzf,fpf,unm,vnm,wnm,tnf,stat=AllocateStatus)
    if (allocated(uf)) deallocate(uf, vf, wf, tempf, pf, pnf, unmif, unif, nlnmif, nlnif, tnf, tnmf, nltnf, nltnmf, stat=AllocateStatus)
    if (AllocateStatus /= 0) then
        write(*,*) "**Error Deallocating - core arrays in finalize **"
    endif
    ! deallocate all paral arrays
    if (allocated(ubar)) deallocate(ubar, ubarp, ubarpp, wbar, wbarp, wbarpp, tbar, tbarp, tbarpp, stat=AllocateStatus)
    if (AllocateStatus /= 0) then
        write(*,*) "**Error Deallocating - paral arrays in finalize **"
    endif
    ! deallocate all grid arrays
    if (nstart==0 .and. allocated(xi)) deallocate(xi, gf2, stat=AllocateStatus)
    if (allocated(wavx)) deallocate(wavx, wavy, wavxx, wavyy, xpts, ypts, zpts, y, gf, stat=AllocateStatus)
    if (AllocateStatus /= 0) then
        write(*,*) "**Error Deallocating - grid arrays in finalize **"
    endif
    ! deallocate all mats arrays
    if (allocated(diag)) deallocate(diag, dor, d, d2, smat, simat, stat=AllocateStatus)
    if (AllocateStatus /= 0) then
        write(*,*) "**Error Deallocating - mats arrays in finalize **"
    endif
    ! deallocate all tavg arrays
    if (allocated(uf_tavg)) deallocate( uf_tavg, vf_tavg, wf_tavg, tempf_tavg, pf_tavg, stat=AllocateStatus)
    if (AllocateStatus /= 0) write(*,*) "**Error Deallocating - tavg arrays in finalize **"

    call mpi_barrier(comm,ierr) ! to be safe

!     if (debug/=0 .and. myid==0) then
!         close(iunit_res) ! close file
!     end if

    call mpi_finalize(ierr) ! to be safe

end subroutine finalize
    

end module init
