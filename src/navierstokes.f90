module navierstokes
! contains all main N-S solver subroutines

implicit none
public
save

contains

subroutine timestep()
    ! performs integration in time of full 3D Navier-Stokes equations
    ! as described in J Andrezj Domaradzki, "An Analytic Green's Functions
    ! Method in Pseudo-Spectral Navier-Stokes Solvers for Boundary Layer 
    ! and Channel Flows", JCP, Vol. 88, No.1, 1990
    use dim
    use flags, only: debug
    use core
    use paral
    use time
    use runparam!, only: ifrz, iomod, dtadapt, sgsflag, ifilt, modfilt
    use grid, only: wavx!, dx, alpha
    use filters
    use derivatives
    use modhorfft
    use sgsinput, only: sgsmodel
    ! to remove after debugging
    use io
    use mpicom!, only: myid, ierr
    real :: diverr!, newdt!, cpu_t1, cpu_t2


    ! physical filtering for Truncated Navier-Stokes approach
    if ( ifilt>=1 .and. mod(istep,modfilt)==0 ) then
        if (sgsflag==1 .and. sgsmodel==4) then
            call tns_automatic_filtering(uf,vf,wf,tempf,pf)
        else
            call filterall3pt(uf,vf,wf,tempf,pf)
!         call filterallspec(uf,vf,wf,0)
        end if
    end if

    ! complete time step update    
    call navierstokes_update(uf, vf, wf, unif, unmif, nlnif, nlnmif, pf, pnf, tempf, tnf, tnmf, nltnf, nltnmf)
!     call nl_visc_press_update_zh(uf, vf, wf, unif, unmif, nlnif, nlnmif, pf, pnf,  sgsflag)
    
    if(myid==0 .and. debug>=2) write(*,*) "N-S update successful" 

    ! normalize FFTs (may be unnecessary)
    call norm(uf,vf,wf,tempf)

    if(myid==0 .and. debug>=2) write(*,*) "fft norm successful" 

    if (ifrz==1 .and. wavx(1)==0) then
        ! keep constant mean vel: remove avg velocity to counter viscous decay
        uf(1:2,1,:) = 0.0
        vf(1:2,1,:) = 0.0
        wf(1:2,1,:) = 0.0
!         tempf(1:2,1,:) = 0.0
    end if

    ! update time
    t = t + dt
    its = its + 1
    isum = its
    tseries(istep) = t

    dt3 = dt2
    dt2 = dt

    if (mod(istep,iomod)==0 .and. debug>=1) then
        call divcheck(uf,vf,wf,diverr)
        if (debug>=2) call printvars(uf,vf,wf,tempf,pf) ! add press calc to this
    end if

    if(myid==0 .and. debug>=2) write(*,*) "divcheck successful" 

    ! adaptive time-stepping to ensure CFL < 0.45
    if ( dtadapt==1 .and. (mod(istep,iomod)==0 .or. mod(istep,5)==0) ) then
        ! check for NaN and abort if so
        if (any(isnan(uf)) .or. any(isnan(vf)) .or. any(isnan(wf)) .or. any(isnan(tempf)) ) then
            if (myid==0) write(*,*) 'NaN detected at t=',t,'iteration #:',its,'. Abort'
            call mpi_finalize(ierr)
        end if
        call maxcalc(uf,vf,wf,tempf,1.0,dt,dt2,dt3)

        if(myid==0 .and. debug>=2) write(*,*) "maxcalc successful" 
    end if
    

end subroutine timestep

subroutine force(u,v,w,fx,fy,fz)
    ! compute body force
    ! built to return flow to inlet conditions
    ! using fringe/sponge layer method
    ! aka Rayleigh damping
    use dim
    use time, only: dt, istep
    use parameters, only: u0,xnu,xlen
    use paral, only: ubar
    use grid
    use runparam, only: fsize, lfsize, sfsize
    ! to remove
    use flags, only: debug
    use io
    use formatstrings
    use mpicom, only: myid
    ! i\o
    real, dimension(nxpp,nypl,nzp), intent(in) :: u,v,w
    real, dimension(nxpp,nypl,nzp), intent(out) :: fx,fy,fz
    ! local vars
    integer :: i, j, k, xshift, ftype
    real :: epsx,dxs,dxl,dxr,kappa,eta0,eta,pih,s0,s,h
    real :: rfsize, dxf, sfs, lfs
!     real :: wav2
    integer :: AllocateStatus
    real, allocatable, dimension(:) :: r, lambda

    !
    ! COMPUTE FORCING COEFFICIENT LAMBDA(X)    
    !
    ! fringe layer parameters
    pih=2.0*atan(1.0)
    epsx = 1.E-5 !dt/2.0 ! free parameter which should be < dt
    h = zpts(1)
    ! fsize = size of fringe layer in percent of domain; 15% is minimum
    sfs = sfsize*fsize
    lfs = lfsize*fsize

    rfsize = fsize-sfs-lfs        ! size of right portion of fringe
    dxs = sfs*xlen        ! length to ramp up forcing over
    dxl = lfs*xlen        ! length of fringe on the left
    dxr = rfsize*xlen        ! length of fringe on the right
    dxf = dxs+dxl+dxr        ! total length of fringe
    kappa = 2.0*u0/(dxf)*(-log(epsx))    ! init strength parameter kappa
    !     if (myid==0 .and. debug>=1) write(*,*) 'fringe kappa =', kappa
    ! Herbst & Henningson exponential damping function
    allocate(r(nxpp), lambda(nxpp), stat=AllocateStatus)
    if (AllocateStatus/=0) write(*,*) "** forcing: allocation error **"
    r = 0.0
    lambda = 0.0

!     ftype = 0
    ftype = 1

    if (ftype==0) then
        do i=1,nx
            eta0 = ( xpts(i) - (xpts(nx) - dxf) )/dxs
            eta = ( xpts(i) - (xpts(nx) - dxl) )/dxl
            s0 = 1./( 1. + exp(1./(eta0 - 1.) + 1./eta0) )
            s = 1./( 1. + exp(1./(eta - 1.) + 1./eta) )
            if (eta0 < 0.) s0 = 0.0
            if (eta0 > 1.) s0 = 1.0
            if (eta < 0.) s = 0.0
            if (eta > 1.) s = 1.0
            r(i) = kappa*(s0 - s)
        end do
        r(nx) = 0.0
    else if (ftype==1) then
        do i=1,nx
            eta0 = ( xpts(i) - (xpts(nx) - dxf + dxs) )/(dxs/2.0)
            eta = ( xpts(i) - (xpts(nx) - dxl) )/(dxl/2.5)
            s0 = exp(-eta0**2)
            s = exp(-eta**2)
            if (eta0>0.0.and. eta<0.) s0 = 1.
            if (eta<0.) s = 0.0
            if (eta>=0.) s0 = 0.0
            r(i) = kappa*(s0 + s)
        end do
    else
        write(*,*) "**ftype error**"
    end if
    ! now shift fringe so that it ends around x/zpts(1)=0.5
    xshift = 0
    do i=1,nx/4
        if (xpts(i)<xpts(1)+dxl/1.5 .and. ftype==1) xshift = i
        if (xpts(i)<xpts(1)+dxl/2.25 .and. ftype==0) xshift = i
!         if (xpts(i)*u0/xnu<=1.5E5) xshift = i
    end do
    ! Write final R(K) to LAMBDA
    if (xshift>1) then
        lambda(xshift+1:nx) = r(1:nx-xshift)
        lambda(1:xshift) = r(nx-xshift+1:nx)
    else
        lambda = r
    end if
    if (debug>=1 .and. myid==0 .and. istep==1) write(*,'(//,A,/,'//fm1//')') 'lambda(x)=',lambda(1:nx)
    !
    ! CALCULATE FORCING TERMS
    !
    fx = u
    fy = v
    fz = w

    ! forcing back to zero perturbation
    do k=1,nzp
        do j=1,nypl
            fx(:,j,k) = lambda(:)*( ubar(k) - fx(:,j,k) )
!             fx(:,j,k) = lambda(:)*( 0.0 - fx(:,j,k) )
            fy(:,j,k) = lambda(:)*( 0.0 - fy(:,j,k) )
            fz(:,j,k) = lambda(:)*( 0.0 - fz(:,j,k) )
!             fz(:,j,k) = 0.0
        end do
    end do


    ! deallocate arrays
    deallocate( r, lambda, stat=AllocateStatus)
    if (AllocateStatus/=0) write(*,*) "** forcing: deallocation error **"

end subroutine force


subroutine bound(bcf,amp)
    ! computes boundary condition
    use dim
    use time, only: t,dt
    use runparam, only: istartbound, samp, fsize, sfsize, lfsize, bctype
    use grid, only: xpts, y, zpts, wavx
    use parameters, only: xlen, u0, xnu, x0
    use modhorfft
    ! to remove
    use flags, only: debug
    use mpicom!, only: myid
    use formatstrings
    ! i/o
    real, intent(out) :: amp
    real, dimension(nxpl,ny,nzp), intent(out) :: bcf
    ! local vars
    integer :: i, j
    real :: xloc, xwidth, h, eta, eta0, bamp, bwidth, bfac ! Spalart
    real :: pi, A, A_half, phi, f, g, x1, x2, xm, omega, lambda_y ! Taraneh
    real :: ref_time, theta
    real, allocatable, dimension(:,:,:) :: bc

    pi = 4.0*atan(1.0)
    ! allocate phys space array for boundary conditions
    allocate( bc(nxpp,nypl,nzp) )
    bc = 0.0

    if (bctype==1) then ! Spalart suction from top
        h = zpts(1) ! Spalart ceiling
        xloc = (1.E5*xnu/u0) ! Spalart suction location based on Re_x
        if (xpts(nx/2)<xloc) then
            ! ensure suction is well within domain and outside forcing
            xloc = x0 + 3./7.*xlen
        end if
        xwidth = 0.24*h
        amp = samp*u0/(0.24*sqrt(pi))!*1.25
!         bamp = amp/1.5 ! use for e^(x^2) bell curve inflow
        bamp = samp/2.25 !0.3/(0.5*fsize*xlen/h) ! use for square wave input using Henningson fcn

        ref_time = xlen/u0!*2. !xlen/u0*4.
        if ( t - istartbound <= ref_time ) then
            amp = amp*sin( 0.5*pi*(t - istartbound)/ref_time )**2.0
            bamp = bamp*sin( 0.5*pi*(t - istartbound)/ref_time )**2.0
        else if (t - istartbound > ref_time) then
            amp = amp*1.
            bamp = bamp*1.
        else
            amp = 0.0
            bamp = 0.0
        end if
  
!         bwidth = 0.35*fsize*xlen ! use 20 points to smooth to expected value
        bwidth = xwidth*2.25
        x1 = xpts(nx) - 2.*bwidth
!         x1 = xpts(nx)  - (0.75*fsize*xlen) ! - 0.2*fsize*xlen
!         x2 = xpts(nx) - xlen*1.1/10.
!         x2 = xpts(nx)  - fsize*xlen - 2.*xwidth !xpts(nx)  - bwidth ! - 0.2*fsize*xlen ! - lfsize*xlen*(1.-1./2.25)

        do i=1,nx
            ! suction
            eta = ( xpts(i) - xloc ) / xwidth
            bc(i,:,1) = amp*exp( -eta**2.0 )
            ! blowing in fringe
            ! just use forcing profile instead of this
            eta0 = ( xpts(i) - x1 ) / bwidth
!             eta = ( xpts(i) - x2 ) / bwidth
!             f = 1./( 1. + exp( 1./(eta0 - 1.) + 1./eta0 ) )
!             g = 1./( 1. + exp( 1./(eta - 1.) + 1./eta ) )
!             if (eta0 < 0.) f = 0.0
!             if (eta0 > 1.) f = 1.0
!             if (eta < 0.) g = 0.0
!             if (eta > 1.) g = 1.0
!             bc(i,:,nzp) = bamp*(f - g)
!             eta0 = ( xpts(i) - (xpts(nx) - 3.25*bwidth) ) / bwidth
!             eta0 = ( xpts(i) - (x0+(1. - 0.6*fsize)*xlen) ) / bwidth
!             eta0 = ( xpts(i) - x2 ) / xwidth
!             bc(i,:,1) = bc(i,:,1) - amp*exp( -eta0**2.0 )
            bc(i,:,nzp) = bamp*exp( -eta0**4.0 )
            ! test case
    !         bc(i,:,1) = -samp*( exp(-((xpts(i)-xloc-xwidth/1.5)/(xwidth/4.))**2.0) )- exp(-((xpts(i)-xloc+xwidth/1.5)/(xwidth/4.))**2.0) )
        end do

    else if (bctype==2) then ! Taraneh blowing/suction from wall: H-Type
        ! h-type transition parameters
        A = 0.002*u0 ! primary amplitude
        A_half = 0.0001*u0 ! subharmonic amplitude
        phi = 0.0! subharmonic phase shift
        x1 = 1.655E5*xnu/u0 ! start of disturbance strip
        x2 = 1.81E5*xnu/u0 ! end of disturbance strip
        xm = ( x1 + x2 ) / 2.0 ! x at middle of strip
        ref_time = x0/u0
        omega = 1.24E-4*ref_time ! Tollmien-Schlichting wave frequency in time
        lambda_y = 0.151367188*x0 ! spanwise subharmonic wavelength
        theta = t - istartbound
        
        ! suction/blowing ansatz s described in Sayadi et al 2011 CTR brief
        do j=1,nypl
            do i=1,nx
                eta = 0.0
                if ( xpts(i) >= x1 .and. xpts(i) <= xm ) then
                    eta = ( xpts(i) - x1 ) / ( xm - x1 )
                end if

                if ( xpts(i) >= xm .and. xpts(i) <= x2 ) then
                    eta = ( x2 - xpts(i) ) / ( x2 - xm )
                end if

                if (xpts(i)>=x1 .and. xpts(i)<=x2) then
                    f = 15.1875*eta**5 - 35.4375*eta**4 + 20.25*eta**3
                    if ( xpts(i)>xm ) f = -f
                    g = cos( 2.0*pi*y(j) / lambda_y )
                    bc(i,j,nzp) = A*f*sin(omega*theta) + A_half*f*g*cos(omega/2.0*theta + phi)
                end if
            end do
        end do
    end if
    ! debug
!     if (t==istartbound .and. t<istartbound+dt .and. debug>=1 .and. myid==0) then
!         write(*,'(//,A,/,'//fm1//')') 'w(x,1,H)=',bc(1:nx,1,1)!,j=1,nypl)
!         write(*,'(//,A,/,'//fm1//')') 'w(x,1,0)=',bc(1:nx,1,nzp)!,j=1,nypl)
!         write(*,*) 'w(kx=0,ky=0,1)=',bcf(1:2,1,1)
!         write(*,*) 'w(kx=0,ky=0,nzp)=',bcf(1:2,1,nzp)
!     end if
    call horfft(bc,bcf,-1) ! transfer to spectral space

    if (bctype==1) then
        ! ensure the divergence is respected by figuring out the exact outflow
        if (wavx(1)==0.) then
            amp = bcf(1,1,1)
            bamp = bcf(1,1,nzp)
            bfac = 0.
            if (bamp/=0) bfac = amp/bamp
        end if
    !     call mpi_barrier(comm,ierr)
        call mpi_bcast(bfac,1,mpi_double_precision,0,comm,ierr)
        bcf(:,:,nzp) = bfac*bcf(:,:,nzp) ! make inflow equal outflow
    end if


    deallocate( bc )

end subroutine bound


subroutine maxcalc(uf,vf,wf,tempf,dtmax,dt,dt2,dt3)
    ! Calculate minimum and maximum values of (u,v,w,temp) 
    ! and print out their locations
    ! Set the timestep based on max calculated CFL if outside bounds
    use dim
    use mpicom
    use modhorfft
    use grid
    use time, only: t, isum
    use paral
    use iofiles
    use runparam, only: ifilt, modfilt, fsize, istartbound, iomod
    ! i/o
    real, dimension(nxpl,ny,nzp), intent(inout) :: uf,vf,wf,tempf
    real, intent(in) :: dtmax
    real, intent(inout) :: dt, dt2, dt3
    ! local vars
    integer :: k, AllocateStatus, findex
    real, allocatable, dimension(:,:,:) :: u,v,w,temp
    real :: cflzmin_lim, cflzmax_lim, cflxmin_lim, cflxmax_lim, cflymin_lim, cflymax_lim
    real :: cflxmax, cflymax, cflzmax, cflxmax_glob, cflymax_glob, cflzmax_glob
    real :: umin, umax0, vmin, vmax, wmin, wmax, tmin, tmax
    real :: umin_glob, umax_glob, vmin_glob, vmax_glob, wmin_glob, wmax_glob, tmin_glob, tmax_glob
    real, parameter :: CDT_CONST=1.25
    real :: timestart, timeend, cdt
    logical change_dt
    logical exceed_limit, below_limit

!     timestart = MPI_WTIME()
    call cpu_time(timestart)

    ! Upper & Lower bounds of allowable vertical CFL#s
    ! and lower limit for x-direction CFL# (actually this is
    ! half the upper bound, because when doubled it still remains
    ! less than 0.2).
!     cflzmin_lim = 0.45; cflzmax_lim = 2.00
!     cflxmin_lim = 0.08; cflxmax_lim = 0.45
!     cflymin_lim = 0.08; cflymax_lim = 0.45
    cflzmin_lim = 0.25; cflzmax_lim = 0.45
    cflxmin_lim = 0.10; cflxmax_lim = 0.25
    cflymin_lim = 0.12; cflymax_lim = 0.25

    allocate( u(nxpp,nypl,nzp), v(nxpp,nypl,nzp), w(nxpp,nypl,nzp), temp(nxpp,nypl,nzp) , stat=AllocateStatus )
    if (AllocateStatus /= 0) then
        write(*,*) "**Not Enough Memory - maxcalc array allocation **"
    end if

    ! Transform data to physical space
    call horfft(u,uf,1)
    call horfft(v,vf,1)
    call horfft(w,wf,1)
    call horfft(temp,tempf,1)

    ! Compute the min/max velocities and maximum CFL
    ! in each processor. Perform REDUCE operation when done
    vmax = maxval(v)
    vmin = minval(v)
    wmax = maxval(w)
    wmin = minval(w)
    do k=1,nzp
        u(:,:,k) = u(:,:,k) + ubar(k)
        temp(:,:,k) = temp(:,:,k) + tbar(k)
    end do
    umax0 = maxval(u)
    umin = minval(u)
    tmin = minval(temp)
    tmax = maxval(temp)

!     cflxmax = maxval(abs(u))*dt/dx
!     cflymax = maxval(abs(v))*dt/dy
    ! since grid is nonuniform in z-direction, CFL-z is dependent on z-loc
    ! use average of cell size above and below since w can be positive or negative 
    do k=1,nz
        w(:,:,k) = w(:,:,k)*dt/(zpts(k)-zpts(k+1))
!         w(:,:,k) = w(:,:,k)*dt/( 0.5*(zpts(k-1)-zpts(k+1)) )
    end do
    ! at top we use cell below because there is no cell above
    w(:,:,1) = w(:,:,1)*dt/( zpts(1) - zpts(2) )
    ! at bottom we use cell above because there is no cell below
    w(:,:,nzp) = w(:,:,nzp)*dt/( zpts(nz) - zpts(nzp) )

    cflxmax = maxval(abs(u))*dt/dx
    cflymax = maxval(abs(v))*dt/dy
!     cflzmax = maxval(abs(w))

    ! calculate maxcfl only in physical region to avoid slowing down simulation for non-physical things that happen in fringe layer
    if (t>istartbound) then
        findex = int((1.-0.8*fsize)*nx) !int((1.-1.*fsize)*nx) 
!     cflxmax = maxval( abs( u(1:findex,:,:) ) )*dt/dx
!     cflymax = maxval( abs( v(1:findex,:,:) ) )*dt/dy
        cflzmax = maxval( abs( w(1:findex,:,:) ) )
    else
        cflzmax = maxval(abs(w))
    end if


    ! perform reduce operation to compute max CFL# and min/max velocities/scalar
    ! could speed this up by keeping only local copy and if any proc has d
    call MPI_ALLREDUCE(umax0,umax_glob,1,mpi_double_precision,MPI_MAX,comm,ierr)
    call MPI_ALLREDUCE(umin,umin_glob,1,mpi_double_precision,MPI_MIN,comm,ierr)

    call MPI_ALLREDUCE(vmax,vmax_glob,1,mpi_double_precision,MPI_MAX,comm,ierr)
    call MPI_ALLREDUCE(vmin,vmin_glob,1,mpi_double_precision,MPI_MIN,comm,ierr)

    call MPI_ALLREDUCE(wmax,wmax_glob,1,mpi_double_precision,MPI_MAX,comm,ierr)
    call MPI_ALLREDUCE(wmin,wmin_glob,1,mpi_double_precision,MPI_MIN,comm,ierr)

    call MPI_ALLREDUCE(tmax,tmax_glob,1,mpi_double_precision,MPI_MAX,comm,ierr)
    call MPI_ALLREDUCE(tmin,tmin_glob,1,mpi_double_precision,MPI_MIN,comm,ierr)

    call MPI_ALLREDUCE(cflxmax,cflxmax_glob,1,mpi_double_precision,MPI_MAX,comm,ierr)
    call MPI_ALLREDUCE(cflymax,cflymax_glob,1,mpi_double_precision,MPI_MAX,comm,ierr)
    call MPI_ALLREDUCE(cflzmax,cflzmax_glob,1,mpi_double_precision,MPI_MAX,comm,ierr)


     exceed_limit = ( (cflzmax_glob > cflzmax_lim) .or.    (cflxmax_glob > cflxmax_lim) .or.    (cflymax_glob > cflymax_lim) )
     below_limit    = ( (cflzmax_glob < cflzmin_lim) .and.    (cflxmax_glob < cflxmin_lim) .and.    (cflymax_glob < cflymin_lim) .and. (dt < dtmax) )

     change_dt = ( exceed_limit .or. below_limit )

    if (change_dt) then

        if ( exceed_limit ) cdt = CDT_CONST        ! Reduce dt by 1/cdt
        if ( below_limit ) cdt = 1.0/CDT_CONST ! Increase dt by 1/cdt

        dt2 = dt3
        dt2 = dt
        dt = dt/cdt

    endif 

    deallocate( u, v, w, temp )

!     timeend = MPI_WTIME()
    call cpu_time(timeend)

    ! display abs max values
    if ( mod(isum,iomod) == 0 .or. change_dt ) then
        if (myid == 0) then
            write(*,'(A,I6,A,F12.4)')'Max CFL calc at iter #', isum,' t=',t
            write(*,*)'Max|u+ubar|=',max(abs(umin_glob),abs(umax_glob))
            write(*,*)'Max|v|    =',max(abs(vmin_glob),abs(vmax_glob))
            write(*,*)'Max|w|    =',max(abs(wmin_glob),abs(wmax_glob))
            write(*,*)'Max|rho|=',max(abs(tmin_glob),abs(tmax_glob))
            ! display all cfl #s
            write(*,*) 'Max. x-CFL#: ',cflxmax_glob
            write(*,*) 'Max. y-CFL#: ',cflymax_glob
            write(*,*) 'Max. z-CFL#: ',cflzmax_glob
            if (change_dt) write(*,*) 'TIMESTEP *REDUCED* by a factor of 4/5 to dt=',dt
            write(*,*) 'Time required for min/max/cfl computation', timeend-timestart
        endif
    end if

end subroutine maxcalc


subroutine nlvel_skewsym(uif, nltermif, sgsflag)
    use dim
    use sgs_models
    use time
    use flags, only: nstart
    use paral, only: ubar, ubarp
    use runparam, only: iomod, forcing, dealiasflag
    use derivatives
    use modhorfft
    use dealias
    ! to remove
!     use filters
    use mpicom, only: myid
    use flags, only: debug
    ! i/o
    real, dimension(nxpl,ny,nzp,3), intent(inout) :: uif
    real, dimension(nxpl,ny,nzp,3), intent(out) :: nltermif
    integer, intent(in) :: sgsflag
    ! local vars
    integer :: k, l, m, nq, AllocateStatus
    real :: cpu_t1, cpu_t2
    real, allocatable, dimension(:,:,:,:) :: ui, nltermi, fi
    real, allocatable, dimension(:,:,:,:,:) :: duidxjf, duidxj, uiuj, tauij, hijf, dhijdxjf, ujduidxj

    nq = 3
    ! allocate arrays
    if ( .not. allocated(ui) ) allocate( ui(nxpp,nypl,nzp,nq), nltermi(nxpp,nypl,nzp,nq), duidxj(nxpp,nypl,nzp,nq,nq), duidxjf(nxpl,ny,nzp,nq,nq), ujduidxj(nxpp,nypl,nzp,nq,nq), uiuj(nxpp,nypl,nzp,nq,nq), tauij(nxpp,nypl,nzp,nq,nq), hijf(nxpl,ny,nzp,nq,nq), dhijdxjf(nxpl,ny,nzp,nq,nq), fi(nxpp,nypl,nzp,nq), stat=AllocateStatus )
    if ( AllocateStatus/=0 ) write(*,*) "** opt_nlvel: allocation error **"

    
    ! zero other arrays
    nltermi = 0.0
    duidxjf = 0.0
    uiuj = 0.0
    tauij = 0.0 
    hijf = 0.0
    dhijdxjf = 0.0

    ! compute derivatives duidxjf
    do l = 1, nq
        ! compute all du_i/dx_j
        call dfdx(uif(:,:,:,l),duidxjf(:,:,:,l,1))
        call dfdy(uif(:,:,:,l),duidxjf(:,:,:,l,2))
        call dfdz(uif(:,:,:,l),duidxjf(:,:,:,l,3))
    end do

    ! bring req'd velocities and derivatives to physical space
    do m = 1, nq
        do l = 1, nq
            call horfft(duidxj(:,:,:,l,m),duidxjf(:,:,:,l,m),1)
        end do
        call horfft(ui(:,:,:,m),uif(:,:,:,m),1)
    end do

    ! add base inflow velocity and derivative
    do k = 1, nzp
        ui(:,:,k,1) = ui(:,:,k,1) + ubar(k)
        duidxj(:,:,k,1,3) = duidxj(:,:,k,1,3) + ubarp(k)
    end do

    do l = 1, nq
        do m = 1, nq          
            if ( dealiasflag==1 ) then
                ! nltermi = -0.5*u_j d(u_i)/dx_j
                call dealiasing_2d_fftw( ui(:,:,:,m), duidxj(:,:,:,l,m), ujduidxj(:,:,:,l,m) )
                nltermi(:,:,:,l) = nltermi(:,:,:,l) - 0.5*ujduidxj(:,:,:,l,m)
                ! uiuj
                call dealiasing_2d_fftw( ui(:,:,:,m), ui(:,:,:,l), uiuj(:,:,:,m,l) )
            else
                ! nltermi = -0.5*u_j d(u_i)/dx_j
                nltermi(:,:,:,l) = nltermi(:,:,:,l) - 0.5*ui(:,:,:,m)*duidxj(:,:,:,l,m)
                ! uiuj
                uiuj(:,:,:,m,l) = ui(:,:,:,m)*ui(:,:,:,l)
            end if
        end do
    end do


    ! compute force terms
    if ( forcing==1 ) then
        call force(ui(:,:,:,1),ui(:,:,:,2),ui(:,:,:,3),fi(:,:,:,1),fi(:,:,:,2),fi(:,:,:,3))
        ! add force terms to nltermi
        do l = 1, nq
            nltermi(:,:,:,l) = nltermi(:,:,:,l) + fi(:,:,:,l)
        end do
    end if


    ! compute sgs terms in physical space
    if ( sgsflag==1 ) then
        call cpu_time(cpu_t1)
        call tau_sgs(ui,uiuj,duidxj,tauij)
        if (myid==0 .and. debug>=2) write(*,*) "tau_sgs step complete"
        call cpu_time(cpu_t2)
        if (myid==0 .and. debug>=1 .and. mod(istep,iomod)==0) write(*,*) 'time spent computing SGS terms = ',cpu_t2 - cpu_t1
    end if


    ! h_ij = 0.5*(u_i u_j) + tau_ij  (may make it non-symmetric)
    uiuj = 0.5*uiuj + tauij

    if (myid==0 .and. debug>=2) write(*,*) "uiuj + tauij done"
    
    do m = 1, nq
        ! take nltermi to fourier space
        call horfft(nltermi(:,:,:,m),nltermif(:,:,:,m),-1)
        do l = 1, nq
            ! take uiuj to spectral space to take divergence
            call horfft(uiuj(:,:,:,l,m),hijf(:,:,:,l,m),-1)
        end do
    end do

    if (myid==0 .and. debug>=2) write(*,*) "transfer to spectral space complete"

    do l = 1, nq
        ! compute d( 0.5*uiuj + tauij )/dx_j which is not necessarily symmetric
        call dfdx(hijf(:,:,:,l,1),dhijdxjf(:,:,:,l,1))
        call dfdy(hijf(:,:,:,l,2),dhijdxjf(:,:,:,l,2))
        call dfdz(hijf(:,:,:,l,3),dhijdxjf(:,:,:,l,3))
    end do

    if (myid==0 .and. debug>=2) write(*,*) "hijf derivatives complete"


    ! FINAL STEP: nltermif = [-0.5*u_j d(u_i)/dx_j + F_i] - d( 0.5*uiuj + tauij )/dx_j
    do l = 1, nq
        do m = 1, nq
            nltermif(:,:,:,l) = nltermif(:,:,:,l) - dhijdxjf(:,:,:,l,m)
!             nltermif(:,:,:,l) = nltermif(:,:,:,l) - dhijdxjf(:,:,:,m,l)
        end do
    end do

    if (myid==0 .and. debug>=2) write(*,*) "nlvel complete"

    ! non-linear term has been computed in skew-symmetric form with full velocity
    ! now update each velocity in time according to time integration scheme

    ! deallocate arrays
    if ( allocated(ui) ) deallocate( ui, nltermi, fi, duidxjf, duidxj, uiuj, tauij, hijf, dhijdxjf, ujduidxj, stat=AllocateStatus )
    if ( AllocateStatus/=0 ) write(*,*) "** opt_nlvel: deallocation error **"

end subroutine nlvel_skewsym


subroutine nlterm_skewsym(uif, tempf, nltermif, nltempf, sgsflag)
    use dim
    use sgs_models
    use time
    use flags, only: nstart
    use paral, only: ubar, ubarp, tbar, tbarp
    use runparam, only: iomod, forcing, dealiasflag
    use derivatives
    use modhorfft
    use dealias
    ! to remove
    use filters
    use mpicom, only: myid
    use flags, only: debug
    ! i/o
    real, dimension(nxpl,ny,nzp,3), intent(inout) :: uif
    real, dimension(nxpl,ny,nzp,3), intent(out) :: nltermif
    real, dimension(nxpl,ny,nzp), intent(inout) :: tempf
    real, dimension(nxpl,ny,nzp), intent(out) :: nltempf
    integer, intent(in) :: sgsflag
    ! local vars
    integer :: k, l, m, nq, AllocateStatus
    real :: cpu_t1, cpu_t2
    real, allocatable, dimension(:,:,:,:) :: ui, nltermi, fi
    real, allocatable, dimension(:,:,:,:,:) :: duidxjf, duidxj, uiuj, tauij, hijf, dhijdxjf, ujduidxj
    real, allocatable, dimension(:,:,:) :: temp, nltemp
    real, allocatable, dimension(:,:,:,:) :: dtdxif, dtdxi, uitemp

    nq = 3
    ! allocate arrays
    if ( .not. allocated(ui) ) allocate( ui(nxpp,nypl,nzp,nq), nltermi(nxpp,nypl,nzp,nq), duidxj(nxpp,nypl,nzp,nq,nq), duidxjf(nxpl,ny,nzp,nq,nq), ujduidxj(nxpp,nypl,nzp,nq,nq), uiuj(nxpp,nypl,nzp,nq,nq), tauij(nxpp,nypl,nzp,nq,nq), hijf(nxpl,ny,nzp,nq,nq), dhijdxjf(nxpl,ny,nzp,nq,nq), fi(nxpp,nypl,nzp,nq), temp(nxpp,nypl,nzp), dtdxif(nxpl,ny,nzp,nq), dtdxi(nxpp,nypl,nzp,nq), uitemp(nxpp,nypl,nzp,nq), nltemp(nxpp,nypl,nzp), stat=AllocateStatus )
    if ( AllocateStatus/=0 ) write(*,*) "** opt_nlvel: allocation error **"

    
    ! zero other arrays
    nltermi = 0.0
    duidxjf = 0.0
    uiuj = 0.0
    tauij = 0.0 
    hijf = 0.0
    dhijdxjf = 0.0

    nltemp = 0.0
    nltempf = 0.0
    temp = 0.0
    dtdxif = 0.0
    dtdxi = 0.0



    ! compute derivatives duidxjf
    do l = 1, nq
        ! compute all du_i/dx_j
        call dfdx(uif(:,:,:,l),duidxjf(:,:,:,l,1))
        call dfdy(uif(:,:,:,l),duidxjf(:,:,:,l,2))
        call dfdz(uif(:,:,:,l),duidxjf(:,:,:,l,3))
    end do

    ! bring req'd velocities and derivatives to physical space
    do m = 1, nq
        do l = 1, nq
            call horfft(duidxj(:,:,:,l,m),duidxjf(:,:,:,l,m),1)
        end do
        call horfft(ui(:,:,:,m),uif(:,:,:,m),1)
    end do

    
    do k = 1, nzp
        ! add base inflow velocity and derivative
        ui(:,:,k,1) = ui(:,:,k,1) + ubar(k)
        duidxj(:,:,k,1,3) = duidxj(:,:,k,1,3) + ubarp(k)
    end do

    ! velocity convection form u_j du_i/dx_j and u_i u_j
    do l = 1, nq
        do m = 1, nq          
            if ( dealiasflag==1 ) then
                ! nltermi = -0.5*u_j d(u_i)/dx_j
                call dealiasing_2d_fftw( ui(:,:,:,m), duidxj(:,:,:,l,m), ujduidxj(:,:,:,l,m) )
                nltermi(:,:,:,l) = nltermi(:,:,:,l) - 0.5*ujduidxj(:,:,:,l,m)
                ! uiuj
                call dealiasing_2d_fftw( ui(:,:,:,m), ui(:,:,:,l), uiuj(:,:,:,m,l) )
            else
                ! nltermi = -0.5*u_j d(u_i)/dx_j
                nltermi(:,:,:,l) = nltermi(:,:,:,l) - 0.5*ui(:,:,:,m)*duidxj(:,:,:,l,m)
                ! uiuj
                uiuj(:,:,:,m,l) = ui(:,:,:,m)*ui(:,:,:,l)
            end if
        end do
    end do

    ! compute sponge terms
    if ( forcing==1 ) then
        call force(ui(:,:,:,1),ui(:,:,:,2),ui(:,:,:,3),fi(:,:,:,1),fi(:,:,:,2),fi(:,:,:,3))
        ! add force terms to nltermi
        do l = 1, nq
            nltermi(:,:,:,l) = nltermi(:,:,:,l) + fi(:,:,:,l)
        end do
    end if


    ! compute sgs terms in physical space
    if ( sgsflag==1 ) then
        call cpu_time(cpu_t1)
        call tau_sgs(ui,uiuj,duidxj,tauij)
        call cpu_time(cpu_t2)
        if (myid==0 .and. debug>=1 .and. mod(istep,iomod)==0) write(*,*) 'time spent computing SGS terms = ',cpu_t2 - cpu_t1
    end if


    ! h_ij = 0.5*(u_i u_j) + tau_ij  (may make it non-symmetric)
    uiuj = 0.5*uiuj + tauij

    
    do m = 1, nq
        ! take nltermi to fourier space
        call horfft(nltermi(:,:,:,m),nltermif(:,:,:,m),-1)
        do l = 1, nq
            ! take uiuj to spectral space to take divergence
            call horfft(uiuj(:,:,:,l,m),hijf(:,:,:,l,m),-1)
        end do
    end do

    do l = 1, nq
        ! compute d( 0.5*uiuj + tauij )/dx_j which is not necessarily symmetric
        call dfdx(hijf(:,:,:,l,1),dhijdxjf(:,:,:,l,1))
        call dfdy(hijf(:,:,:,l,2),dhijdxjf(:,:,:,l,2))
        call dfdz(hijf(:,:,:,l,3),dhijdxjf(:,:,:,l,3))
    end do


    ! FINAL STEP: nltermif = [-0.5*u_j d(u_i)/dx_j + F_i] - d( 0.5*uiuj + tauij )/dx_j
    do l = 1, nq
        do m = 1, nq
            nltermif(:,:,:,l) = nltermif(:,:,:,l) - dhijdxjf(:,:,:,l,m)
!             nltermif(:,:,:,l) = nltermif(:,:,:,l) - dhijdxjf(:,:,:,m,l)
        end do
    end do

    ! non-linear vel term has been computed in skew-symmetric form with full velocity
    ! now update each velocity in time according to time integration scheme

    ! if non-zero temperature, compute derivatives
    if ( any( abs(tempf)>0.0 ) ) then
        ! compute all dT/dx_i
        call dfdx(tempf(:,:,:),dtdxif(:,:,:,1))
        call dfdy(tempf(:,:,:),dtdxif(:,:,:,2))
        call dfdz(tempf(:,:,:),dtdxif(:,:,:,3))
        ! bring all derivatives to physical space
        do l = 1, nq
            call horfft(dtdxi(:,:,:,l),dtdxif(:,:,:,l),1)
        end do
        ! bring temp to physical space
        call horfft(temp,tempf,1)
    end if

    ! add base temperature and derivative
    do k = 1, nzp
        temp(:,:,k) = temp(:,:,k) + tbar(k)
        dtdxi(:,:,:,3) = dtdxi(:,:,:,3) + tbarp(k)
    end do
    ! temperature convection form and u_iT
    do l = 1, nq
        if ( dealiasflag==1 ) then
            ! reuse ujduidxj = u_i dt/dx_i
            call dealiasing_2d_fftw( ui(:,:,:,l), dtdxi(:,:,:,l), ujduidxj(:,:,:,l,1) )
            nltemp = nltemp - 0.5*ujduidxj(:,:,:,l,1)
            ! uiuj
            call dealiasing_2d_fftw( ui(:,:,:,l), temp, uitemp(:,:,:,l) )
        else
            nltemp = nltemp - 0.5*ui(:,:,:,l)*dtdxi(:,:,:,l)
            uitemp(:,:,:,l) = ui(:,:,:,l)*temp
        end if
    end do

    ! take temperature terms to spectral space
    call horfft(nltemp,nltempf,-1)
    do l = 1, nq
        call horfft(uitemp(:,:,:,l),hijf(:,:,:,l,1),-1)        
    end do
    ! compute d(u_i T)/dx_i
    call dfdx(hijf(:,:,:,1,1),dtdxif(:,:,:,1))
    call dfdy(hijf(:,:,:,2,1),dtdxif(:,:,:,2))
    call dfdz(hijf(:,:,:,3,1),dtdxif(:,:,:,3))
    ! add each to the nonlinear temp term
    do l = 1, nq
        nltempf = nltempf - 0.5*dtdxif(:,:,:,l)
    end do

    ! non-linear temp term has been computed in skew-symmetric form with full temp
    ! now update temp in time according to time integration scheme


    ! deallocate arrays
    if ( allocated(ui) ) deallocate( ui, nltermi, fi, duidxjf, duidxj, uiuj, tauij, hijf, dhijdxjf, ujduidxj, temp, dtdxi, dtdxif, nltemp, stat=AllocateStatus )
    if ( AllocateStatus/=0 ) write(*,*) "** nlterm_skewsym: deallocation error **"

end subroutine nlterm_skewsym



subroutine nl_visc_press_update_zh(uf, vf, wf, unif, unmif, nlnif, nlnmif, pf, pnf,  sgsflag)
    ! computes non-linear term in skew symmetric form using a BDF3 scheme
    ! performs implicit viscous step via Helmholtz eqn (1 for each direction)
    ! uses Zang-Hussaini boundary conditions to improve tangential velocity at wall
    ! computes pressure correction using usual dirichlet BC's 
    use dim
    use time
    use flags, only: nstart
    use parameters, only: xnu
    use paral, only: ubarpp
    use modhorfft
    use mats
    use helmholtz
    use runparam, only: bctype, forcing, istartbound, iomod, dealiasflag
    use grid, only: wavx, wavxx, wavyy
    use derivatives
    use filters
    ! to remove
    use mpicom, only: myid
    use flags, only: debug
    use io
    ! i/o
    real, dimension(nxpl,ny,nzp), intent(inout) :: uf, vf, wf, pf, pnf
    real, dimension(nxpl,ny,nzp,3), intent(inout) :: unif, unmif, nlnif, nlnmif
    integer, intent(in) :: sgsflag
    ! local vars
    integer :: i, j, k, AllocateStatus, nq
    real, dimension(:,:,:,:), allocatable :: uif, nltermif, tmp, rhsf
    real, dimension(:,:,:), allocatable :: du0dz2, du0dz2f
    real :: dtnui, amp, wav2, pfac(2)
    ! pressure vars
    real, allocatable, dimension(:,:,:) :: divf, divdzf, dwdzf, dpf, psf
!     real, dimension(:,:,:,:), allocatable :: gradpf
    
    nq = 3

    ! allocate nl_visc term arrays
    if (.not. allocated(uif)) allocate(uif(nxpl,ny,nzp,nq), nltermif(nxpl,ny,nzp,nq), tmp(nxpl,ny,nzp,nq), rhsf(nxpl,ny,nzp,nq), du0dz2(nxpp,nypl,nzp), du0dz2f(nxpl,ny,nzp), stat=AllocateStatus)
    if ( AllocateStatus/=0 ) write(*,*) "** nl_visc_update_bdf3: allocation error **"

    ! allocate pressure arrays
    allocate( divf(nxpl,ny,nzp), divdzf(nxpl,ny,nzp), dwdzf(nxpl,ny,nzp), dpf(nxpl,ny,nzp), psf(nxpl,ny,nzp), stat=AllocateStatus)
    if (AllocateStatus /= 0) then
        write(*,*) "**Not Enough Memory - pressure arrays allocation **"
    end if


    ! assign velocities
    uif(:,:,:,1) = uf
    uif(:,:,:,2) = vf
    uif(:,:,:,3) = wf

    
    !!!!!!!!!!!!!!!!!!!! Non-linear term calculation !!!!!!!!!!!!!!!!!!!!!!!

    ! sharp spectral filter instead of dealiasing before nonlinear evaluation
    if ( dealiasflag>=2 ) then
        call filterqspec(uif,nq)
    end if

    ! compute nonlinear term in skew-symmetric form
    call nlvel_skewsym(uif,nltermif,sgsflag)

    if (myid==0 .and. debug>=2) write(*,*) "nlvel step complete"

    ! re-assign non-filtered velocities
    if ( dealiasflag==2 ) then
        uif(:,:,:,1) = uf
        uif(:,:,:,2) = vf
        uif(:,:,:,3) = wf
    end if

    ! temporarily store nltermif
    tmp = nltermif

    ! set AB3/SS3 coefficients
    if (istep<=3 .or. dt/=dt2 .or. dt2/=dt3) call abcoef(nstart,istep,dt,dt2,dt3,abfac)
    
    ! AB3/SS3 update of non-linear term
    nltermif = abfac(1)*nltermif + abfac(2)*nlnif + abfac(3)*nlnmif

    ! get BDF3 coefficients
    if (istep<=3 .or. dt/=dt2 .or. dt2/=dt3) call bdfcoef(nstart,istep,dt,dt2,dt3,bdf,bdfv)

    if ( istep>3 .and. myid==0 .and. (dt/=dt2 .or. dt2/=dt3) ) write(*,*) 'bdfv =', bdfv, 'bdf =', bdf

    
    
    !!!!!!!!!!!!!!!!!!!! RHS of viscous-nonlinear step !!!!!!!!!!!!!!!!!!!!!!!

    ! ( D^2 - k^2 - 1./(nu*dt) )u* = -1./(nu*dt)*(bdf(1)*u^n + bdf(2)*u^n-1 + bdf(3)*u^n-2) - 1/nu*(abfac(1)*N^n + abfac(2)*N^n-1 + abfac(3)*N^n-2)
    dtnui = 1./(xnu*dt)
    rhsf = -dtnui*(dt*nltermif + bdf(1)*uif + bdf(2)*unif + bdf(3)*unmif)

    ! add mean D^2 u (ubarpp) to first equation
    if (sum(abs(ubarpp))>0.) then
        du0dz2 = 0.0
        do k=1,nzp
            du0dz2(1:nx,:,k) = ubarpp(k)
        end do
        call horfft(du0dz2,du0dz2f,-1)
    else
        du0dz2f = 0.0
    end if
    rhsf(:,:,:,1) = rhsf(:,:,:,1) - du0dz2f

    ! now that RHS is complete
    ! relegate n-1 to n-2 level
    nlnmif = nlnif
    unmif = unif
    ! relegate n to n-1 level
    nlnif = tmp
    unif = uif

    

    !!!!!!!!!!!!!!!!!!!! BCs for viscous-nonlinear step !!!!!!!!!!!!!!!!!!!!!!!

    ! set Dirichlet boundary conditions for combined viscous non-linear step
    ! w - vertical velocity
    if ( (bctype>0 .and. t>=istartbound) ) then
        call bound(uif(:,:,:,3), amp)
!     else if ( forcing==1 ) then
!         uif(:,:,nzp,3) = 0.0 ! zero at wall
!         !D^2(w) -> 0 at top free flow
!         do j = 1, ny
!             do i = 1, nxpl
!                 wav2 = wavxx(i)+wavyy(j)
!                 uif(i,j,1,3) = -rhsf(i,j,1,3) / ( wav2 + dtnui/bdfv )
!             end do
!         end do
    else
        uif(:,:,1,3) = 0.0 ! homogeneous BC's
        uif(:,:,nzp,3) = 0.0
    end if
    
    ! set boundary conditions for horizontal velocities u, v
    uif(:,:,nzp,2) = 0.0 ! u, v (z=0) = 0 at wall
    uif(:,:,1,2) = 0.0 ! v = 0 at top
!     uif(:,:,1,1) = 0.0 ! u = 0 at top
    
    ! Zang-Hussaini BC for u* = dt*bdfv*dp^n/dx
    if (istep<2) then
        pfac(1:2) = 0.0
    else if (istep==2) then
        pfac(1) = 1.0
        pfac(2) = 0.0
    else
        pfac(1) =  2.0
        pfac(2) = -1.0
    end if
    ! add pressure correction term
    psf = pfac(1)*pf + pfac(2)*pnf !+ abfac(3)*pnmf
    ! take d(p*)/dx
    call dfdx(psf,dpf)
    pnf = pf ! relegate past pressure to current level
!     uif(:,:,1,1) = dt*bdfv*dpf(:,:,1)
    uif(:,:,nzp,1) = dt*bdfv*dpf(:,:,nzp)
    ! D^2(u,v) -> 0 at top free flow
    do j = 1, ny
        do i = 1, nxpl
            wav2 = wavxx(i)+wavyy(j)
!             uif(i,j,1,:) = -rhsf(i,j,1,:) / ( wav2 + dtnui/bdfv )
!             uif(i,j,1,1:2) = -rhsf(i,j,1,1:2) / ( wav2 + dtnui/bdfv )
            uif(i,j,1,1) = -rhsf(i,j,1,1) / ( wav2 + dtnui/bdfv )
        end do
    end do

    ! let w(0,H) = dt*nu*bdfv*d^2w/dz^2|(0,H)
    call d2fdz2(uif(:,:,:,3),dwdzf)
    uif(:,:,1,3) = uif(:,:,1,3) + xnu*dt*bdfv*dwdzf(:,:,1)
!     uif(:,:,nzp,3) = uif(:,:,nzp,3) + xnu*dt*bdfv*dwdzf(:,:,nzp)
    ! try w(H) = dp/dz(H)
!     call dfdz(psf,dpf)
!     uif(:,:,1,3) = uif(:,:,1,3) + dt*bdfv*dpf(:,:,1)
!     uif(:,:,nzp,3) = uif(:,:,nzp,3) + dt*bdfv*dpf(:,:,nzp)

    


    !!!!!!!!!!!!!!!!!!!! Solve viscous-nonlinear Helmholtz eqns !!!!!!!!!!!!!!!!!!!!!!!
    
    ! solve each Helmholtz equation with Dirichlet BCs
    do i = 1, nq
        if ( (bctype>0 .and. t>=istartbound) .or. forcing==1 ) then
            call solve_nhom_helmholtz(rhsf(:,:,:,i),uif(:,:,:,i),dtnui/bdfv,diag,smat,simat)
        else
            call solve_helmholtz(rhsf(:,:,:,i),dtnui/bdfv,nzm,diag,smat,simat)
        end if
    end do

    
    ! update velocities
!     uif = rhsf
    uf = rhsf(:,:,:,1)
    vf = rhsf(:,:,:,2)
    wf = rhsf(:,:,:,3)

    if (myid==0 .and. debug>=2) write(*,*) "viscous step complete"
    
    !!!!!!!!!!!!!!!!!!!! RHS of pressure step !!!!!!!!!!!!!!!!!!!!!!!

    ! form RHS of: (D^2-k^2)w^n+1 = -D(ikx u* + iky v*)^n - k^2 w*^n
    ! calc d/dz(ikx u +    ikyv)
!     call dfdz(uf,dpf)
!     call dfdz(vf,dwdzf)
!     call divuv(dpf,dwdzf,divdzf)
    call divuv(uf,vf,divf) ! divf = ikx u + iky v
    call dfdz(divf,divdzf) ! divdzf = D[ikx u + iky v]
    ! ozf = -D[ikx u + iky v] - k^2 w
    do j=1,ny
        do i=1,nxpl
            wav2 = wavxx(i) + wavyy(j)
            rhsf(i,j,:,3)= -divdzf(i,j,:) - wav2*wf(i,j,:)
        end do
    end do
    if (debug>=2) then
        if (myid==0) write(*,'(//,A)') 'w*, D[ikx u* + iky v*], D[ikx u* + iky v*] + k^2 w*} + dt*fzf'
        call printuvw(wf,divdzf,rhsf(:,:,:,3))
    end if


    !!!!!!!!!!!!!!!!!!!! BC's for pressure step !!!!!!!!!!!!!!!!!!!!!!!

    ! set boundary conditions for w^(n+1)
    if ( bctype>0 .and. t>=istartbound ) then
        call bound(wf,amp)
!     else if (forcing==1) then
!         wf(:,:,nzp) = 0.0
    else
        wf(:,:,1) = 0.0
        wf(:,:,nzp) = 0.0
    end if

!     ! extra boundary conditions from FLOWGUN to improve error at wall
!     call dfdz(uf - unif(:,:,:,1),dpf)
!     call dfdz(vf - unif(:,:,:,2),dwdzf)
!     call divuv(dpf,dwdzf,divdzf)
! !     call divuv(rhsf(:,:,:,1),rhsf(:,:,:,2),dpf) ! divf = ikx u + iky v
! !     call dfdz(dpf,divdzf) ! divdzf = D[ikx u + iky v]
!     ! ozf = -D[ikx u + iky v] - k^2 w
!     do j=1,ny
!         do i=1,nxpl
!             wav2 = wavxx(i)+wavyy(j)
!             wf(i,j,1) = wf(i,j,1) - xnu*dt*bdfv*( divdzf(i,j,1) + wav2*(wf(j,i,1) - unif(i,j,1,3)) )
!             wf(i,j,nzp) = wf(i,j,nzp) - xnu*dt*bdfv*( divdzf(i,j,nzp) + wav2*(wf(j,i,nzp) - unif(i,j,nzp,3)) )
!         end do
!     end do



    !!!!!!!!!!!!!!!!!!!! Solve Poisson eqn in pressure step !!!!!!!!!!!!!!!!!!!!!!!
    
    ! Solve poisson equation for w^n+1
    ! (D^2-k^2)w^n+1 = -D(ikx u + iky v)^n - k^2 w^n
    call solve_nhom_helmholtz(rhsf(:,:,:,3),wf,0.0,diag,smat,simat)

    wf = rhsf(:,:,:,3)
    ! phi*dt*bdfv = -[Dw^n+1 + ikx u^n + iky v^n]/k^2
    call dfdz(wf,dwdzf)
    do j=1,ny
        do i=1,nxpl
            wav2 = wavxx(i)+wavyy(j)
            if (wav2==0.) then
                pf(i,j,:) = 0.0 ! can't divide by zero, mean pressure is arbitrary
            else
                pf(i,j,:)= -( dwdzf(i,j,:) + divf(i,j,:) )/(wav2*dt*bdfv)
            end if
        end do
    end do
    ! compute u^n+1 and v^n+1
    call dfdx(pf,dpf)
    uf = uf - dt*bdfv*dpf
    
    if (debug>=2 .and. mod(istep,iomod)==0) then
        if (myid==0) write(*,'(//,A)') 'dw**/dz, dynp**, dp/dx'
        call printuvw(dwdzf,pf,dpf)
    !     call printuvw(wf,pf,uf)
    end if

    call dfdy(pf,dpf)
    vf = vf - dt*bdfv*dpf

!     ! get pressure alone: phi/(dt*bdfv) = (P^n+1 - P* + nu*div(u*))
!     ! P^n+1 = phi/(dt*bdfv) + P* - nu*div(u*)
!     call div(uif(:,:,:,1),uif(:,:,:,2),uif(:,:,:,3),divf)
!     pf = pf + psf - xnu*divf

    if (myid==0 .and. debug>=2) write(*,*) "pressure step complete"

    ! deallocate arrays
    if ( allocated(uif) ) deallocate( uif, nltermif, tmp, rhsf, psf, du0dz2, du0dz2f, divf, divdzf, dwdzf, dpf, stat=AllocateStatus )
    if ( AllocateStatus/=0 ) write(*,*) "** nl_visc_update_bdf3: deallocation error **"

end subroutine nl_visc_press_update_zh



subroutine navierstokes_update(uf, vf, wf, unif, unmif, nlnif, nlnmif, pf, pnf, tempf, tnf, tnmf, nltnf, nltnmf)
    ! computes non-linear term in skew symmetric form using a BDF3 scheme
    ! performs implicit viscous step via Helmholtz eqn (1 for each direction)
    ! uses Zang-Hussaini boundary conditions to improve tangential velocity at wall
    ! computes pressure correction using usual dirichlet BC's 
    use dim
    use time
    use flags, only: nstart
    use parameters
    use paral, only: ubarpp, tbarpp
    use modhorfft
    use mats
    use helmholtz
    use runparam
    use grid
    use derivatives
    use filters
    use sgsinput, only: a,b
    ! to remove
    use paral, only: ubar
    use mpicom, only: myid
    use flags, only: debug
    use io
    ! i/o
    real, dimension(nxpl,ny,nzp), intent(inout) :: uf, vf, wf, pf, pnf
    real, dimension(nxpl,ny,nzp), intent(inout) :: tempf, tnf, tnmf, nltnf, nltnmf
    real, dimension(nxpl,ny,nzp,3), intent(inout) :: unif, unmif, nlnif, nlnmif
!     integer, intent(in) :: sgsflag
    ! local vars
    integer :: i, j, k, AllocateStatus, nq
    real :: dtnui, amp, wav2, pfac(2), dtkapi

    real, dimension(:,:,:,:), allocatable :: uif, nltermif, tmp, rhsf
    real, dimension(:,:,:), allocatable :: du0dz2, du0dz2f, nltempf, tempsf, divf, divdzf, dwdzf, dpf, psf
    
    nq = 3

    ! allocate nl_visc term arrays
    if ( .not. allocated(uif) ) allocate( uif(nxpl,ny,nzp,nq), nltermif(nxpl,ny,nzp,nq), tmp(nxpl,ny,nzp,nq), rhsf(nxpl,ny,nzp,nq), du0dz2(nxpp,nypl,nzp), du0dz2f(nxpl,ny,nzp), nltempf(nxpl,ny,nzp), tempsf(nxpl,ny,nzp), divf(nxpl,ny,nzp), divdzf(nxpl,ny,nzp), dwdzf(nxpl,ny,nzp), dpf(nxpl,ny,nzp), psf(nxpl,ny,nzp), stat=AllocateStatus )
    if ( AllocateStatus/=0 ) write(*,*) "** navierstokes_update: allocation error **"


    if(myid==0 .and. debug>=2) write(*,*) "array allocation successful"

    ! assign velocities
    uif(:,:,:,1) = uf
    uif(:,:,:,2) = vf
    uif(:,:,:,3) = wf

    ! assign temperature
    tempsf = tempf

    
    !!!!!!!!!!!!!!!!!!!! Non-linear term calculation !!!!!!!!!!!!!!!!!!!!!!!

    ! sharp spectral filter instead of dealiasing before nonlinear evaluation
    if ( dealiasflag>=2 ) then
        call filterqspec(uif,nq)
        if ( any(abs(tempsf)>1.0E-10) ) call filterqspec(tempsf,1)
        if(myid==0 .and. debug>=2) write(*,*) "dealias filter successful"
    end if

    ! compute nonlinear term in skew-symmetric form
    call nlterm_skewsym(uif, tempsf, nltermif, nltempf, sgsflag)

    if(myid==0 .and. debug>=2) write(*,*) "nlterm successful"

    ! re-assign non-filtered velocities
    if ( dealiasflag==2 ) then
        uif(:,:,:,1) = uf
        uif(:,:,:,2) = vf
        uif(:,:,:,3) = wf
    end if

    ! temporarily store nltermif
    tmp = nltermif
    tempsf = nltempf

    ! set AB3/SS3 coefficients
    if (istep<=3 .or. dt/=dt2 .or. dt2/=dt3) call abcoef(nstart,istep,dt,dt2,dt3,abfac)
    
    ! AB3/SS3 update of non-linear terms
    nltermif = abfac(1)*nltermif + abfac(2)*nlnif + abfac(3)*nlnmif
    nltempf = abfac(1)*nltempf + abfac(2)*nltnf + abfac(3)*nltnmf

    ! get BDF3 coefficients
    if (istep<=3 .or. dt/=dt2 .or. dt2/=dt3) call bdfcoef(nstart,istep,dt,dt2,dt3,bdf,bdfv)
    if ( debug>=1 .and. myid==0 .and. istep>3 .and. (dt/=dt2 .or. dt2/=dt3) ) write(*,*) 'bdfv =', bdfv, 'bdf =', bdf



    !!!!!!!!!!!!!!!!!!!! RHS of viscous-nonlinear TEMPERATURE step !!!!!!!!!!!!!!!!!!!!!

    ! ( D^2 - k^2 - 1./(kap*dt) )T^n+1 = -1./(kap*dt)*(bdf(1)*T^n + bdf(2)*T^n-1 + bdf(3)*T^n-2) - 1/nu*(abfac(1)*N(T)^n + abfac(2)*N(T)^n-1 + abfac(3)*N(T)^n-2)
    xkap = xnu/prtl
    dtkapi = 1./(xkap*dt)
    rhsf(:,:,:,1) = -dtkapi*( dt*nltempf + bdf(1)*tempf + bdf(2)*tnf + bdf(3)*tnmf )

    ! add mean D^2 T (tbarpp) to first equation
    if (sum(abs(tbarpp))>0.) then
        du0dz2 = 0.0
        do k=1,nzp
            du0dz2(1:nx,:,k) = tbarpp(k)
        end do
        call horfft(du0dz2,du0dz2f,-1)
    else
        du0dz2f = 0.0
    end if
    rhsf(:,:,:,1) = rhsf(:,:,:,1) - du0dz2f

    ! solve homogeneous Helmholtz eqn for temperature
    call solve_helmholtz(rhsf(:,:,:,1),dtkapi/bdfv,nzm,diag,smat,simat)
    
    ! update previous time steps
    nltnmf = nltnf
    tnmf = tnf
    nltnf = tempsf
    tnf = tempf
    tempf = rhsf(:,:,:,1)

    if(myid==0 .and. debug>=2) write(*,*) "helmholtz temp successful"   

    !!!!!!!!!!!!!!!!!!!! RHS of viscous-nonlinear VELOCITY step !!!!!!!!!!!!!!!!!!!!!!!

    ! ( D^2 - k^2 - 1./(nu*dt) )u* = -1./(nu*dt)*(bdf(1)*u^n + bdf(2)*u^n-1 + bdf(3)*u^n-2) - 1/nu*(abfac(1)*N^n + abfac(2)*N^n-1 + abfac(3)*N^n-2)
    dtnui = 1./(xnu*dt)
    rhsf = 0.0
    rhsf = -dtnui*( dt*nltermif + bdf(1)*uif + bdf(2)*unif + bdf(3)*unmif )

    ! add mean D^2 u (ubarpp) to first equation
    if (sum(abs(ubarpp))>0.) then
        du0dz2 = 0.0
        do k=1,nzp
            du0dz2(1:nx,:,k) = ubarpp(k)
        end do
        call horfft(du0dz2,du0dz2f,-1)
    else
        du0dz2f = 0.0
    end if
    rhsf(:,:,:,1) = rhsf(:,:,:,1) - du0dz2f

    if ( any(abs(tempf)>1.0E-10) ) then
        ! add buoyancy force
        rhsf(:,:,:,1) = rhsf(:,:,:,1) - coex*gx*tempf/xnu
        rhsf(:,:,:,2) = rhsf(:,:,:,2) - coex*gy*tempf/xnu
        rhsf(:,:,:,3) = rhsf(:,:,:,3) - coex*gz*tempf/xnu
    end if

    if(myid==0 .and. debug>=2) write(*,*) "RHS vel successful" 

    ! relegate n-1 to n-2 level
    nlnmif = nlnif
    unmif = unif
    ! relegate n to n-1 level
    nlnif = tmp
    unif = uif

    !!!!!!!!!!!!!!!!!!!! BCs for viscous-nonlinear step !!!!!!!!!!!!!!!!!!!!!!!

    ! set Dirichlet boundary conditions for combined viscous non-linear step
    
    ! set boundary conditions for horizontal velocities u, v
    
    ! Zang-Hussaini BC for u* = dt*bdfv*dp^n/dx
    if (istep==1) then
        pfac(1:2) = 0.0
    else if (istep==2) then
        pfac(1) = 1.0
        pfac(2) = 0.0
    else
        pfac(1) =  2.0
        pfac(2) = -1.0
    end if
    ! add pressure correction term
    psf = pfac(1)*pf + pfac(2)*pnf !+ abfac(3)*pnmf
    pnf = pf ! relegate past pressure to current level
    ! take d(p*)/dx
    call dfdx(psf,dpf)
!     uif(:,:,1,1) = uif(:,:,2,1)
!     uif(:,:,1,1) = dt*bdfv*dpf(:,:,1)
!     uif(:,:,nzp,1) = dt*bdfv*dpf(:,:,nzp)
    uif(:,:,nzp,1) = dt*bdfv*dpf(:,:,nzp)
    ! D^2(u,v) -> 0 at top free flow
    do j = 1, ny
        do i = 1, nxpl
            wav2 = wavxx(i)+wavyy(j)
!             uif(i,j,1,:) = -rhsf(i,j,1,:) / ( wav2 + dtnui/bdfv )
            uif(i,j,1,1:2) = -rhsf(i,j,1,1:2) / ( wav2 + dtnui/bdfv )
!             uif(i,j,1,1) = -rhsf(i,j,1,1) / ( wav2 + dtnui/bdfv )
        end do
    end do

    ! Zang-Hussaini BC for w* = dt*nu*bdfv*d^2w^n/dz^2
        call d2fdz2(uif(:,:,:,3),dwdzf)    
    ! try w(H) = dt*bdfv*dp/dz(H) which is supposedly unstable
!     call dfdz(psf,dpf)
    if ( (bctype>0 .and. t>=istartbound) ) then
        call bound(uif(:,:,:,3), amp) ! imposed suction
        ! add Zang-Hussaini BC correction
        uif(:,:,1,3) = uif(:,:,1,3) + xnu*dt*bdfv*dwdzf(:,:,1)
!         uif(:,:,nzp,3) = uif(:,:,nzp,3) + xnu*dt*bdfv*dwdzf(:,:,nzp)
        ! add supposedly unstable BC w(H) = dt*bdfv*dp/dz(H)
!         uif(:,:,1,3) = uif(:,:,1,3) + dt*bdfv*dpf(:,:,1)
!         uif(:,:,nzp,3) = uif(:,:,nzp,3) + dt*bdfv*dpf(:,:,nzp)
    else
        uif(:,:,1,3) = xnu*dt*bdfv*dwdzf(:,:,1)
        uif(:,:,nzp,3) = xnu*dt*bdfv*dwdzf(:,:,nzp)
!         uif(:,:,nzp,3) = 0.0
!         uif(:,:,1,3) = dt*bdfv*dpf(:,:,1)
!         uif(:,:,nzp,3) = dt*bdfv*dpf(:,:,nzp)
    end if

    ! Zang-Hussaini BC for v* = dt*bdfv*dp^n/dy
    ! take d(p*)/dy
    call dfdy(psf,dpf)
!     uif(:,:,1,2) = dt*bdfv*dpf(:,:,1)
    uif(:,:,nzp,2) = dt*bdfv*dpf(:,:,nzp)
!     uif(:,:,1,2) = 0.0
!     uif(:,:,nzp,2) = 0.0
    ! or compute what BC should be from continuity
!     call dfdx(uif(:,:,:,1),dpf)
!     call dfdz(uif(:,:,:,3),dwdzf)
!     uif(:,:,1,2) = - dpf(:,:,1) - dwdzf(:,:,1)
!     uif(:,:,nzp,2) = - dpf(:,:,nzp) - dwdzf(:,:,nzp) 


    if(myid==0 .and. debug>=2) write(*,*) "vel BC's successful" 

    !!!!!!!!!!!!!!!!!!!! Solve viscous-nonlinear Helmholtz eqns !!!!!!!!!!!!!!!!!!!!!!!
    
    ! solve each Helmholtz equation with Dirichlet BCs
    do i = 1, nq
        if ( (bctype>0 .and. t>=istartbound) .or. forcing==1 ) then
            call solve_nhom_helmholtz(rhsf(:,:,:,i),uif(:,:,:,i),dtnui/bdfv,diag,smat,simat)
        else
            call solve_helmholtz(rhsf(:,:,:,i),dtnui/bdfv,nzm,diag,smat,simat)
        end if
    end do
    
    
!     uif = rhsf
    ! update velocities
    uf = rhsf(:,:,:,1)
    vf = rhsf(:,:,:,2)
    wf = rhsf(:,:,:,3)

    rhsf = 0.0

    if(myid==0 .and. debug>=2) write(*,*) "vel Helmholtz successful" 
    
    !!!!!!!!!!!!!!!!!!!! RHS of pressure step !!!!!!!!!!!!!!!!!!!!!!!

    ! form RHS of: (D^2-k^2)w^n+1 = -D(ikx u* + iky v*)^n - k^2 w*^n
    ! calc d/dz(ikx u +    ikyv)
    call dfdz(uf,dpf)
    call dfdz(vf,dwdzf)
    call divuv(dpf,dwdzf,divdzf)
    call divuv(uf,vf,divf) ! divf = ikx u + iky v
!     call dfdz(divf,divdzf) ! divdzf = D[ikx u + iky v]
    ! ozf = -D[ikx u + iky v] - k^2 w
    do j=1,ny
        do i=1,nxpl
            wav2 = wavxx(i) + wavyy(j)
            rhsf(i,j,:,3)= -divdzf(i,j,:) - wav2*wf(i,j,:)
        end do
    end do
    if (debug>=2) then
        if (myid==0) write(*,'(//,A)') 'w*, D[ikx u* + iky v*], D[ikx u* + iky v*] + k^2 w*} + dt*fzf'
        call printuvw(wf,divdzf,rhsf(:,:,:,3))
    end if


    !!!!!!!!!!!!!!!!!!!! BC's for pressure step !!!!!!!!!!!!!!!!!!!!!!!

    ! set boundary conditions for w^(n+1)
    if ( bctype>0 .and. t>=istartbound ) then
        call bound(wf,amp)
!     else if (forcing==1) then
!         wf(:,:,nzp) = 0.0
    else
        wf(:,:,1) = 0.0
        wf(:,:,nzp) = 0.0
    end if

!     ! extra boundary conditions from FLOWGUN to improve error at wall
!     call dfdz(uf - unif(:,:,:,1),dpf)
!     call dfdz(vf - unif(:,:,:,2),dwdzf)
!     call divuv(dpf,dwdzf,divdzf)
! !     call divuv(rhsf(:,:,:,1),rhsf(:,:,:,2),dpf) ! divf = ikx u + iky v
! !     call dfdz(dpf,divdzf) ! divdzf = D[ikx u + iky v]
!     ! ozf = -D[ikx u + iky v] - k^2 w
!     do j=1,ny
!         do i=1,nxpl
!             wav2 = wavxx(i)+wavyy(j)
!             wf(i,j,1) = wf(i,j,1) - xnu*dt*bdfv*( divdzf(i,j,1) + wav2*(wf(j,i,1) - unif(i,j,1,3)) )
!             wf(i,j,nzp) = wf(i,j,nzp) - xnu*dt*bdfv*( divdzf(i,j,nzp) + wav2*(wf(j,i,nzp) - unif(i,j,nzp,3)) )
!         end do
!     end do



    !!!!!!!!!!!!!!!!!!!! Solve Poisson eqn in pressure step !!!!!!!!!!!!!!!!!!!!!!!
    
    ! Solve poisson equation for w^n+1
    ! (D^2-k^2)w^n+1 = -D(ikx u + iky v)^n - k^2 w^n
    call solve_nhom_helmholtz(rhsf(:,:,:,3),wf,0.0,diag,smat,simat)

    wf = rhsf(:,:,:,3)
    ! phi*dt*bdfv = -[Dw^n+1 + ikx u^n + iky v^n]/k^2
    call dfdz(wf,dwdzf)
    do j=1,ny
        do i=1,nxpl
            wav2 = wavxx(i)+wavyy(j)
            if (wav2==0.) then
                pf(i,j,:) = 0.0 ! can't divide by zero, mean pressure is arbitrary
            else
                pf(i,j,:)= -( dwdzf(i,j,:) + divf(i,j,:) )/(wav2*dt*bdfv)
            end if
        end do
    end do
    ! compute u^n+1 and v^n+1
    call dfdx(pf,dpf)
    uf = uf - dt*bdfv*dpf
    
    if (debug>=2 .and. mod(istep,iomod)==0) then
        if (myid==0) write(*,'(//,A)') 'dw**/dz, dynp**, dp/dx'
        call printuvw(dwdzf,pf,dpf)
    !     call printuvw(wf,pf,uf)
    end if

    call dfdy(pf,dpf)
    vf = vf - dt*bdfv*dpf

    if(myid==0 .and. debug>=2) write(*,*) "pressure step successful" 

    ! deallocate arrays
    if ( allocated(uif) )deallocate( uif, nltermif, tmp, rhsf, du0dz2, du0dz2f, nltempf, tempsf, divf, divdzf, dwdzf, dpf, psf, stat=AllocateStatus )
    if ( AllocateStatus/=0 ) write(*,*) "** nl_visc_update_bdf3: deallocation error **"

end subroutine navierstokes_update



subroutine abcoef(nstart, istep, dt1, dt2, dt3, abfac)
    ! computes appropriate coefficients for non-linear step
    integer, intent(in) :: nstart, istep
    real, intent(in) :: dt1, dt2, dt3
    real, dimension(3), intent(out) :: abfac

    ! first step starting from scratch
    if ( nstart==0 ) then
        if ( istep==1 ) then
            ! Euler time step
            abfac(1) =   1.0
            abfac(2:3) = 0.0
        else if ( istep==2 ) then
            !SS2
            abfac(1) =  2.0
            abfac(2) = -1.0
            abfac(3) =  0.0
            ! AB2
!             abfac(1) =  1.5
!             abfac(2) = -0.5
!             abfac(3) =  0.0
        ! time step not constant
!         else if (dt1/=dt2 .or. dt2/=dt3) then
!             abfac(1) = abvar1(dt1,dt2,dt3)
!             abfac(2) = abvar2(dt1,dt2,dt3)
!             abfac(3) = abvar3(dt1,dt2,dt3)
        else
            abfac(1) =  3.0
            abfac(2) = -3.0
            abfac(3) =  1.0
        end if
    ! time step not constant
!     else if (dt1/=dt2 .or. dt2/=dt3) then
!         abfac(1) = abvar1(dt1,dt2,dt3)
!         abfac(2) = abvar2(dt1,dt2,dt3)
!         abfac(3) = abvar3(dt1,dt2,dt3)
    else
        abfac(1) =  3.0
        abfac(2) = -3.0
        abfac(3) =  1.0
    end if

end subroutine abcoef


subroutine bdfcoef(nstart, istep, dt1, dt2, dt3, bdf, bdfv)
    ! computes appropriate coefficients for non-linear step
    integer, intent(in) :: nstart, istep
    real, intent(in) :: dt1, dt2, dt3
    real, dimension(3), intent(out) :: bdf
    real, intent(out) :: bdfv

    ! first step starting from scratch
    if ( nstart==0 ) then
        if ( istep==1 ) then
            ! Euler time step
            bdfv = 1.0
            bdf(1) = 1.0
            bdf(2:3) = 0.0
        else if ( istep==2 ) then
            ! BDF2
            bdfv = 2./3.
            bdf(1) =  2.0
            bdf(2) = -0.5
            bdf(3) =  0.0
        ! time step not constant
        else if (dt1/=dt2 .or. dt2/=dt3) then
            bdfv = 1./bdfvarv(dt1,dt2,dt3)
            bdf(1) = bdfvar1(dt1,dt2,dt3)
            bdf(2) = bdfvar2(dt1,dt2,dt3)
            bdf(3) = bdfvar3(dt1,dt2,dt3)
        else
            bdfv = 6./11.
            bdf(1) =  3.0
            bdf(2) = -1.5
            bdf(3) =  1./3.
        end if
    ! time step not constant
    else if (dt1/=dt2 .or. dt2/=dt3) then
        bdfv = 1./bdfvarv(dt1,dt2,dt3)
        bdf(1) = bdfvar1(dt1,dt2,dt3)
        bdf(2) = bdfvar2(dt1,dt2,dt3)
        bdf(3) = bdfvar3(dt1,dt2,dt3)
    else
        bdfv = 6./11.
        bdf(1) =  3.0
        bdf(2) = -1.5
        bdf(3) =  1./3.
    end if

end subroutine bdfcoef


! Adams-Bashforth 3rd order variable time step coefficients

real function abvar1(dt1,dt2,dt3)
    real, intent(in) :: dt1,dt2,dt3
    real :: dt123, dt23
    dt123 = 2.*dt1 + 6.*dt2 + 3.*dt3
    dt23 = dt2 + dt3
    abvar1 = 1. + dt1*dt123/(6.*dt2*dt23)
end function abvar1


real function abvar2(dt1,dt2,dt3)
    real, intent(in) :: dt1,dt2,dt3
    real :: dt123
    dt123 = 2.*dt1 + 3.*dt2 + 3.*dt3
    abvar2 = -dt1*dt123/(6.*dt2*dt3)
end function abvar2


real function abvar3(dt1,dt2,dt3)
    real, intent(in) :: dt1,dt2,dt3
    real :: dt12, dt23
    dt12 = 2.*dt1 + 3.*dt2
    dt23 = dt2 + dt3
    abvar3 = dt1*dt12/(6.*dt3*dt23)
end function abvar3


! BDF3 variable time step coefficients

real function bdfvarv(dt1,dt2,dt3)
    real, intent(in) :: dt1,dt2,dt3
    real :: dt12, dt123
    dt12 = dt1+dt2
    dt123 = dt1+dt2+dt3 
    bdfvarv = 1. + dt1/dt12 + dt1/dt123
end function bdfvarv


real function bdfvar1(dt1,dt2,dt3)
    real, intent(in) :: dt1,dt2,dt3
    real :: dt12, dt23, dt123
    dt12 = dt1+dt2
    dt23 = dt2+dt3
    dt123 = dt1+dt2+dt3
    bdfvar1 = -1. - dt1/dt2 - dt1/dt123 - dt1*dt1/(dt123*dt2) - dt1*dt1*dt12/(dt123*dt23*dt2) 
    bdfvar1 = -bdfvar1
end function bdfvar1


real function bdfvar2(dt1,dt2,dt3)
    real, intent(in) :: dt1,dt2,dt3
    real :: dt12, dt23, dt123
    dt12 = dt1+dt2
    dt23 = dt2+dt3
    dt123 = dt1+dt2+dt3
    bdfvar2 = dt1*dt1/(dt12*dt2) + dt1*dt1/(dt2*dt123) + dt1*dt1*dt12/(dt2*dt23*dt123) + dt1*dt1*dt12/(dt3*dt23*dt123)
    bdfvar2 = -bdfvar2 
end function bdfvar2


real function bdfvar3(dt1,dt2,dt3)
    real, intent(in) :: dt1,dt2,dt3
    real :: dt12, dt23, dt123
    dt12 = dt1+dt2
    dt23 = dt2+dt3
    dt123 = dt1+dt2+dt3
    bdfvar3 = -dt1*dt1*dt12/(dt3*dt23*dt123)            
    bdfvar3 = -bdfvar3
end function bdfvar3


end module navierstokes
