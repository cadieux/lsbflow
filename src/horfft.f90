module modhorfft

implicit none

contains

subroutine test_fft(u,avgerror)
! test accuracy of fft for a given function
use flags, only: debug
use dim, only: nx, nxpp, nxpl, ny, nypl, nzpl
use mpicom, only: myid, numprocs, comm, ierr
! i/o
real, dimension(nxpp,nypl,nzpl), intent(in) :: u
real, intent(out) :: avgerror
! vars
integer :: i,j,k,ii
real, allocatable, dimension(:,:,:) :: qf, q

allocate(qf(nxpl,ny,nzpl),q(nxpp,nypl,nzpl))

q = u
call horfft(q,qf,-1)
q = 0.0
call horfft(q,qf,1)

avgerror = 0.0
do k=1,nzpl
    do j=1,nypl
        do i=1,nx
            avgerror = avgerror + abs( u(i,j,k)-q(i,j,k) )/(nx*nypl*nzpl)
        end do
    end do
end do

if (debug>=1 .and. myid==0) then
    write(*,*) 'avg error introduced by FFT =', avgerror
end if

if (debug>=3) then
    if (myid==0) then
        write(*,*)    '----------------------------------------------------'
        write(*,*)    'u(:,:,1) in test_fft:'
    end if

    do ii=1,numprocs
        if (myid+1 == ii) then
            write(*,'(256(G12.4,'',''))') ( (u(1,j,k),i=1,nx),j=1,nypl )
        else
            call sleep(1)
        end if
    end do
    call MPI_BARRIER(comm,ierr)
endif


deallocate(qf,q)

end subroutine test_fft


subroutine horfft(ain,ainf,is)
! SUBROUTINE THAT PERFORMS 2-D FFT OVER ALL nzpl XY PLANES ON A
! PROCESSOR.
!
! is = -1: PHYSICAL TO SPECTRAL SPACE
! is =    1: SPECTRAL TO PHYSICAL SPACE
!
! THIS SUBROUTINE IS A FFTW VERSION OF "horfft".
!
! PROGRAMMED BY: T.SAKAI 10-16-2012 (MODIFIED FROM ORIGINAL horfft)

use dim, only: nx,nxpp,nypl,nzpl,nxpl,ny
use mpicom
use MODFFTW


! i/o
integer, intent(IN) :: is
real, dimension(nxpp,nypl,nzpl), intent(inout) :: ain
real, dimension(nxpl,ny,nzpl), intent(inout) :: ainf

! Scratch arrays
real, allocatable, dimension(:,:,:) :: aintmp, aintmpf

! Local variables
integer :: k
real :: scalefactor, eps

allocate( aintmp(nxpp,nypl,nzpl), aintmpf(nxpl,ny,nzpl) )
aintmp=0.0
aintmpf=0.0
eps = 1.E-20 ! constant to remove noise coming from FFT

if ( is == -1 ) then

    aintmp = ain

    !-------------------
    !-1-D FORWARD FFTs-|
    !-------------------

    ! Perform 1-D FFT in x-direction (real-to-Complex) for each
    ! k-plane in processor
    do k = 1,nzpl
        call dfftw_execute_dft_r2c(PLAN_X_FWD,aintmp(1,1,k),aintmp(1,1,k))
    end do

    ! Global data transposition
    call dataflip_xtoy(aintmp,ainf)

    ! Perform 1-D FFT in y-direction (Complex-to-Complex) for each
    ! k-plane in processor
    do k = 1,nzpl
        call dfftw_execute_dft(PLAN_Y_FWD,ainf(1,1,k),ainf(1,1,k))
    end do

    ! Scale the result
    scalefactor = 1.0/real(nx*ny)
    ainf = ainf*scalefactor

else

    aintmpf = ainf

    !-------------------
    !-1-D INVERSE FFTs-|
    !-------------------

    ! Perform 1-D FFT in y-direction (Complex-to-Complex)
    ! for each k-plane in processor
    do k = 1,nzpl
        call dfftw_execute_dft(PLAN_Y_BACKWD,aintmpf(1,1,k),aintmpf(1,1,k))
    end do

    ! Global data transposition
    call dataflip_ytox(aintmp,aintmpf)

    ! Perform 1-D INVERSE FFT in x-direction (Complex-to-real)
    ! for each k-plane in processor
    do k = 1,nzpl
        call dfftw_execute_dft_c2r(PLAN_X_BACKWD,aintmp(1,1,k),aintmp(1,1,k))
    end do

    ain = aintmp

end if

deallocate( aintmp, aintmpf )

end subroutine horfft


subroutine norm(uf,vf,wf,tempf)
use dim
use grid, only: wavx
! inputs/outputs
real, dimension(nxpl,ny,nzpl), intent(inout) :: uf,vf,wf,tempf
! local vars
integer :: nyhp, i,j,k, jp


!    nyh = ny/2
nyhp = nyh + 1
! intended only to affect k_x=0
if (wavx(1)==0) then
    ! remove any imaginary part to the integral (avg) over entire domain
    do k = 1,nzp        
        uf(2,1,k) = 0.0
        vf(2,1,k) = 0.0
        wf(2,1,k) = 0.0
        tempf(2,1,k) = 0.0
    enddo
endif
! if this is a 2-D problem, don't normalize
if (ny > 2) then
! intended only to affect k_x=0
    if (wavx(1)==0) then
        ! Normalize by taking average of FFT value for ky and conjugate of -ky
        ! e.g. 0.5*( q(kx=0,ky,z) + conjg(q(kx=0,-ky,z)) )
        do j = 2,nyh
            jp = ny + 2 - j
            do k = 1,nzp
                uf(1,j,k) = .5*( uf(1,j,k) + uf(1,jp,k) )
                uf(2,j,k) = .5*( uf(2,j,k) - uf(2,jp,k) )
                vf(1,j,k) = .5*( vf(1,j,k) + vf(1,jp,k) )
                vf(2,j,k) = .5*( vf(2,j,k) - vf(2,jp,k) )
                wf(1,j,k) = .5*( wf(1,j,k) + wf(1,jp,k) )
                wf(2,j,k) = .5*( wf(2,j,k) - wf(2,jp,k) )
                tempf(1,j,k) = .5*( tempf(1,j,k) + tempf(1,jp,k) )
                tempf(2,j,k) = .5*( tempf(2,j,k) - tempf(2,jp,k) )
            enddo
        enddo
        ! Copy over in correct order the normalized values into -ky half
        ! s.t. ky = 0, b, 2b, ...,(nyh-1)b,0,-(nyh-1)b,-(nyh-2)b,...,2b,b
        ! and q(kx,-ky,z) = conjg(q(kx=0,ky,z))
        do j = 2,nyh
            jp = ny + 2 - j
            do k = 1,nzp
                uf(1,jp,k) =  uf(1,j,k)
                uf(2,jp,k) = -uf(2,j,k)
                vf(1,jp,k) =  vf(1,j,k)
                vf(2,jp,k) = -vf(2,j,k)
                wf(1,jp,k) =  wf(1,j,k)
                wf(2,jp,k) = -wf(2,j,k)
                tempf(1,jp,k) =  tempf(1,j,k)
                tempf(2,jp,k) = -tempf(2,j,k)
            enddo
        enddo

    endif ! (wavx(1)==0)
    ! remove any value at ky=-nyh*beta which should really be ky=0
    do k = 1,nzp
        do i = 1,nxpl
            uf(i,nyhp,k) = 0.0
            vf(i,nyhp,k) = 0.0
            wf(i,nyhp,k) = 0.0
            tempf(i,nyhp,k) = 0.0
        enddo
    enddo

endif ! (ny > 2)


end subroutine norm


!**********************************************************
!
!-This file contains subroutines used to flip data across
!-processors on the horizontal plane in such a way that
!-one can transition from domain decomposition in the x-direction
!-to domain decomposition in the y-direction and vice versa.
!
!-Developed by PD: January 2005.
!
!-These subroutines are patterned after the similar ones
!-in Kraig Winters' spectral code (See JAOT article).
!
!-PD-6/16/08--C A R E F U L: This is a new version of the routines
!-developed at Cornell by PD.
!-The data transposition is no longer done with immediate and
!-ready communication modes for the sends & receives. Such an approach
!-did not factor in switch topology and contention issues in the MPI
!-library routines. As a result, on the ARL-MJM cluster with
!-either Open-MPI or Infiniserve-MPI (and not MPICH which was
!-problem free in the past) time-out problems would occur during
!-code execution.
!-
!-Following Kraig Winters' advice we use MPI_ALL_TO_ALL global
!-communication routines to perform the data transposition.
!-
!-NOTE-PD-6/16/08: The local array buffer is practically useless right now.
!-However, I have kept it in the code because dataflip is called by
!-a couple of other routines, namely RILEY_DBK, REGRID & EXTRAPOLATE
!-(in both postprocessor and main solver).
!-This array should eventually be removed !
!***********************************************************



subroutine dataflip_xtoy(x,xf)
!************************************************************
! Subroutine that reorders data from domain decomposition in
! in x-direction to d.d. in y-direction
!************************************************************
use dim, only: nx,nxpp,nypl,nzpl,nxpl,ny
use mpicom
!-Input variable x( , , ): nzpl 2-D slices of data partitioned normal to y-
!-direction and
!-Output variable xf( , ,): nzpl 2-D slices of data partitioned normal to x-direction

real, dimension(nxpp,nypl,nzpl), intent(in) :: x
real, dimension(nxpl,ny,nzpl), intent(out) :: xf
real, dimension(:,:,:,:), allocatable :: xtempin, xtempout
!     real, dimension(nxpl,nzpl,nypl,nproch) :: xtempin, xtempout
!     integer :: numbytes
!     integer :: status_array(MPI_STATUS_SIZE)

!     integer :: iproc,vslab
!     integer :: istart,jstart,source,dest

integer :: i,j,k
integer :: AllocateStatus

integer :: size_of_block, iglob, jglob, iproc

!-Allocate local arrays
!- SEND BUFFER for MPI_ALLTOALL
!- RECV BUFFER for MPI_ALLTOALL
allocate( xtempin(nxpl,nzpl,nypl,nproch),xtempout(nxpl,nzpl,nypl,nproch),stat = AllocateStatus)
!
 if (AllocateStatus /= 0) then
    write(*,*) "**Not Enough Memory - DATAFLIP_XTOY**"
 end if

!     status_array = 0


!-Pack input array into send buffer.
!-Order of fastest incrasing indices: i,k,j,iproc (in send buffer)
!-So that upon reception only a simple index change will only be needed.
!-(i.e. received data will be contiguous in j-direction when stored in memory)
do iproc=1,nproch
    do k=1,nzpl
        do j=1,nypl
            do i=1,nxpl

                 iglob = (iproc-1)*nxpl + i

                 xtempin(i,k,j,iproc) = x(iglob,j,k)

            end do
        end do
    end do
end do

!-Size of blocks to be transmitted to other processors
size_of_block = nxpl*nypl*nzpl

!***TAK 5-23-2012: ADD SYNCRONIZATION BEFORE GLOBAL COLLECTIVE
!    OPERATION. THIS IS TO ENSURE THE PORTABILITY.
!call MPI_BARRIER(comm,ierr)

!-Perform data transposition
call MPI_ALLTOALL(xtempin,size_of_block,mpi_double_precision,xtempout,size_of_block,mpi_double_precision,comm,ierr)

if (ierr /= 0) then
    if (myid == 0) write(*,*) 'MPI_ALLTOALL error in DATAFLIP_XTOY)'
    call MPI_FINALIZE(ierr)
end if

!-Now copy the receive buffer from xtempout to the output array
do iproc=1,nproch
    do k=1,nzpl
        do j=1,nypl
            do i=1,nxpl

                 jglob = (iproc-1)*nypl + j

                 xf(i,jglob,k) = xtempout(i,k,j,iproc)

            end do
        end do
    end do
end do

!-Add ta barrier to be on the safe side
!call MPI_BARRIER(comm,ierr)

!-De-Allocate local array
deallocate(xtempin,xtempout,stat=AllocateStatus)
if (AllocateStatus /= 0) then
    write(*,*) "**Error Deallocating - DATAFLIP_XTOY**"
end if

end subroutine dataflip_xtoy


subroutine dataflip_ytox(x,xf)
!************************************************************
! Subroutine that reorders data from domain decomposition in
! in y-direction to d.d. in x-direction
!************************************************************
use dim, only: nx,nxpp,nypl,nzpl,nxpl,ny
use mpicom
!-Input variable xf( , ,): nzpl 2-D slices of data partitioned normal to x-direction
!-Output variable x( , , ): nzpl 2-D slices of data partitioned normal to y-direction.
real, dimension(nxpp,nypl,nzpl), intent(out) :: x
real, dimension(nxpl,ny,nzpl), intent(in) :: xf
real, dimension(:,:,:,:), allocatable :: xtempin,xtempout
!     real, dimension(nypl,nzpl,nxpl,nproch) :: xtempin,xtempout
!     integer status_array(MPI_STATUS_SIZE)

!     integer :: iproc,vslab
!     integer :: istart,jstart,source,dest,numwords,tag

integer :: i,j,k
integer :: AllocateStatus

integer :: size_of_block, iglob, jglob, iproc

!-Allocate local arrays
!- SEND BUFFER for MPI_ALLTOALL
!- RECV BUFFER for MPI_ALLTOALL
allocate( xtempin(nypl,nzpl,nxpl,nproch), xtempout(nypl,nzpl,nxpl,nproch),stat = AllocateStatus)

if (AllocateStatus /= 0) then
    write(*,*) "**Not Enough Memory - DATAFLIP_YTOX**"
end if

!     status_array = 0

!-Pack input array into send buffer.
!-Order of fastest incrasing indices: j,k,i,iproc (in send buffer)
!-So that upon reception only a simple index change will only be needed.
!-(i.e. received data will be contiguous in i-direction when stored in memory)
!-All processors post nonblocking buffered sends & receives
!-(i.e. each processor goes
!-about its business) and set-up    buffers for the sends

do iproc=1,nproch
    do k=1,nzpl
        do j=1,nypl
            do i=1,nxpl

                jglob = (iproc-1)*nypl + j

                 xtempin(j,k,i,iproc) = xf(i,jglob,k)

            end do
        end do
    end do
end do

!-Size of blocks to be transmitted to other processors
size_of_block = nxpl*nypl*nzpl

!***TAK 5-23-2012: ADD SYNCRONIZATION BEFORE GLOBAL COLLECTIVE
!    OPERATION. THIS IS TO ENSURE THE PORTABILITY.
!call MPI_BARRIER(comm,ierr)

!-Perform data transposition
call MPI_ALLTOALL(xtempin,size_of_block,mpi_double_precision,xtempout,size_of_block,mpi_double_precision,comm,ierr)

if (ierr /= 0) then
    if (myid == 0) write(*,*) 'MPI_ALLTOALL error in DATAFLIP_YTOX'
    call MPI_FINALIZE(ierr)
end if

!-Now copy the receive buffer from xtempout to the output array
do iproc=1,nproch
    do k=1,nzpl
        do j=1,nypl
            do i=1,nxpl

                iglob = (iproc-1)*nxpl + i

                x(iglob,j,k) = xtempout(j,k,i,iproc)

            end do
        end do
    end do
end do

!-Add a barrier to be on the safe side
!call MPI_BARRIER(comm,ierr)

!-De-Allocate local array
deallocate(xtempin,xtempout,stat=AllocateStatus)
if (AllocateStatus /= 0) then
    write(*,*) "**Error Deallocating - DATAFLIP_YTOX**"
end if

end subroutine dataflip_ytox



end module modhorfft