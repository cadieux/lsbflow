! module double
! ! use, intrinsic :: iso_fortran_env, dp=>real64, idp=>int64
! implicit none
! public
! save
! integer, parameter :: dp = selected_real_kind(15,307)
! integer, parameter :: idp = 8
! end module

module dim

implicit none
public
save

integer :: nx, ny, nz, nzp, nxpl, nxhpl, nypl, nzpl ! main dimensions
! extra stuff used by certain subroutines
integer :: nzp2, nzm, nxh, nyh, ndimx, ndimy, ndimz, ndimz2, ndimx2, nxp, nyp, nxpp, nx2h, nzpp, nzpp2, nyplp, nxin, nyin, nzpin, nyplin
! for reading inputs
namelist/dimen/ nx, ny, nz

contains
    subroutine init_dim(nproch,nprocv)
        ! i/o
        integer, intent(in) :: nproch,nprocv

        nzp = nz+1
        nxp=nx+1
        nyp=ny+1

        nzp2=2*nzp
        nzm=nz-1
        nxh=nx/2
        nyh=ny/2

        ndimx=nx
        ndimy=ny
        ndimz=nzp

        ndimz2=2*ndimz
        ndimx2=2*ndimx

        nxpp=nx+2
        nx2h=nxpp/2
        nzpp=nzp+1
        nzpp2=nzp+2
        
        nzpl = nzp/nprocv
        nypl = ny/nproch
        nyplp = nypl + 1
        nxpl = nx/nproch
        nxhpl = nxpl/2

    end subroutine init_dim

end module dim


module flags

implicit none
public

integer :: nstart, iwrite, debug, io_mpi
namelist/flaglist/ nstart, iwrite, debug, io_mpi

end module flags


module runparam

implicit none
public
save
integer :: nsteps, iomod, otfavg, istartavg, mapping, ic, bctype, forcing, linstab, ifrz, dtadapt, sgsflag, ifilt, modfilt, interpflag, dealiasflag
real :: cfl, blfac, fsize, lfsize, sfsize, samp, astar, istartbound

namelist/run/ cfl, dtadapt, nsteps, iomod, otfavg, istartavg, mapping, blfac, ic, bctype, istartbound, samp, interpflag, forcing, fsize, lfsize, sfsize, linstab, astar, ifrz, sgsflag, ifilt, modfilt, dealiasflag

end module runparam


module parameters

implicit none
public
save

real :: xlen, ylen, zlen, x0, u0, xnu, prtl, rayl, ri, thigh, tlow, coex, gx, gy, gz
complex :: cguess
real :: xkap, deld, re, tayl, uup, guup, gradt, rotx, roty, rotz
namelist/param/ xlen, ylen, zlen, x0, u0, xnu, prtl, rayl, tayl, ri, thigh, tlow, coex, gx, gy, gz, cguess

end module parameters


module sgsinput

implicit none
public
save

integer :: sgsmodel
real :: a,b, CS
real :: kc, filt_order, filt_amp

namelist/sgsin/ sgsmodel, a, b, CS, kc, filt_order, filt_amp

end module sgsinput


module grid

implicit none
public
save    

real, allocatable, dimension(:) :: wavx, wavy, wavxx, wavyy
real, allocatable, dimension(:) ::  xpts, ypts, zpts, y, xin, yin, yin_loc, zin
real, allocatable, dimension(:) :: gf, gf2, xi, gfin
real :: alpha, beta, dx, dy, dz, xl, bfac, aa, bb

end module grid


module wallunit

implicit none
public
save
real :: Refr
real, allocatable, dimension(:) :: zplus
real, allocatable, dimension(:,:) :: zp
real, allocatable, dimension(:,:,:) :: uplus, upt, vpt, wpt

end module wallunit


module time

implicit none
public
save
integer :: its, istep, isum
real :: t, dt, dt2, dt3, abfac(3), bdf(3), bdfv
real, allocatable, dimension(:) :: tseries

end module time


module paral

implicit none
public
save
! for parallel boundary layer or mean flow qties    
    real, allocatable, dimension(:) :: ubar, ubarp, ubarpp, wbar, wbarp, wbarpp, tbar, tbarp, tbarpp

end module paral


module core

implicit none
public
save
    ! velocities and nonlinear term at time n, (n-1)
    real, allocatable, dimension(:,:,:) :: uf, vf, wf
    real, allocatable, dimension(:,:,:,:) :: unif, unmif, nlnif, nlnmif
    ! pressure
    real, allocatable, dimension(:,:,:) :: pf, pnf
    ! temp at previous time & nonlinear term at time n, (n-1)
    real, allocatable, dimension(:,:,:) :: tempf, tnf, tnmf, nltnf, nltnmf
    
end module core


module tavg

implicit none
public
save
    ! time-averaged velocities and temperature
    real, allocatable, dimension(:,:,:) :: uf_tavg,vf_tavg,wf_tavg,tempf_tavg,pf_tavg

end module tavg


module mats

implicit none
public
save
    ! stores arrays used for Helmholtz solver
    real, allocatable, dimension(:) :: diag, fdiag
    real, allocatable, dimension(:,:) :: dor, d, d2, smat, simat, fsmat, fsimat
end module mats


module formatstrings

implicit none
public
save
    ! to automatically write nicely to screen
    character(len=5) :: n1p,n2p,n3p,n1f,n2f,n3f
    character(len=100) :: fm1,fm2,fm3,fm1f,fm2f,fm3f

contains
    subroutine set_formatstrings()
        use dim
        use flags, only: debug
!         use mpicom, only: myid
        write(n1p,'(I5)') nx
        write(n2p,'(I5)') nypl
        write(n3p,'(I5)') nzpl
        write(n1f,'(I5)') nxpl
        write(n2f,'(I5)') ny
        write(n3f,'(I5)') nzpl
        fm1 = '('//trim(n1p)//'(G12.5,'',''))'
        fm2 = '('//trim(n2p)//'(G12.5,'',''))'
        fm3 = '('//trim(n3p)//'(G12.5,'',''))'
        fm1f = '('//trim(n1f)//'(G12.5,'',''))'
        fm2f = '('//trim(n2f)//'(G12.5,'',''))'
        fm3f = '('//trim(n3f)//'(G12.5,'',''))'
        ! debug
!         if (debug>=2 .and. myid==0) then
!             write(*,*) 'formatstrings test:'
!             write(*,*) 'nx, nypl, nzpl = ',n1p,n2p,n3p
!             write(*,*) 'nxpl, ny, nzpl = ',n1f,n2f,n3f
!         end if
    end subroutine set_formatstrings

end module formatstrings


module prng

implicit none
public
save

contains

subroutine init_random_seed(n)
    ! creates a pseudo-random seed for random_number() of size n
!     use iso_fortran_env!, only: int64
!     use ifport
    integer, allocatable :: rseed(:)
    integer :: i, n, un, istat, dt(8), pid, t(2), s
    integer :: count, tms
             
    call random_seed(size = n)
    allocate(rseed(n))
    ! First try if the OS provides a random number generator
    un = 12
    open(unit=un, file="/dev/urandom", access="stream",form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
        read(un) rseed
        close(un)
    else
        ! Fallback to XOR:ing the current time and pid. The PID is
        ! useful in case one launches multiple instances of the same
        ! program in parallel.
        call system_clock(count)
        if (count /= 0) then
            t = transfer(count, t)
        else
            call date_and_time(values=dt)
            tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
                  + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
                  + dt(3) * 24 * 60 * 60 * 60 * 1000 &
                  + dt(5) * 60 * 60 * 1000 &
                  + dt(6) * 60 * 1000 + dt(7) * 1000 &
                  + dt(8)
            t = transfer(tms, t)
        end if
        s = ieor(t(1), t(2))
        pid = getpid() + 1099279 ! Add a prime
        s = ieor(s, pid)
        if (n >= 3) then
            rseed(1) = t(1) + 36269
            rseed(2) = t(2) + 72551
            rseed(3) = pid
            if (n > 3) then
                rseed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
            end if
        else
            rseed = s + 37 * (/ (i, i = 0, n - 1 ) /)
        end if
    end if
    call random_seed(put=rseed)
end subroutine init_random_seed

end module prng