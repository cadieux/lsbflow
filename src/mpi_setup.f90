module mpicom
! use mpi

implicit none
include 'mpif.h'
public
save

integer :: myid, ierr, numprocs, comm
integer :: nprocs,nproch,nprocv,hproc,vproc

contains

subroutine init_mpi()
    ! MPI init
    comm = mpi_comm_world
    call mpi_init(ierr)
    call mpi_comm_rank(comm,myid,ierr)
    call mpi_comm_size(comm,numprocs,ierr)
    ! use the number of available mpi processes as final nprocs
    nprocs = numprocs
    ! only do domain decomposition in horizontal
    nproch = nprocs
    ! domain decomposition in vertical is disabled because it was too slow
    nprocv = 1
end subroutine init_mpi


subroutine setup_mpi()

    ! Check for errors in domain decomposition
    call decomp_check()

end subroutine setup_mpi


subroutine decomp_check()
    ! check for domain decomposition issues
    use dim, only: nx, ny, nzp, nxh, nyh, nxpl, nxhpl, nypl, nzpl
    character (len=50) :: f1, f2, f3

    f1 = ' '
    f2 = ' '
    f3 = ' '
    if ( nxhpl*2 /= nxpl .or. nxpl*nproch /= nx .or. nxh*2 /= nx ) then
        ierr = 1
        f1 = ' in the x direction '
    end if
    if ( nypl*nproch /= ny .or. nyh*2 /= ny ) then
        ierr = 1
        f2 = ', in the y direction '
    end if
    if ( nzpl*nprocv /= nzp ) then
        ierr = 1
        f3 = ', in the z direction '
    end if

    if (ierr == 1) then
        if (myid == 0) write(*,*) 'ERROR: incorrect domain decompositon'//trim(f1)//trim(f2)//trim(f3)
        call MPI_BARRIER(comm,ierr)
        call MPI_FINALIZE(ierr)
    end if
end subroutine decomp_check


end module mpicom