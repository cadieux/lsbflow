module parallel_io

implicit none
public
save

contains

subroutine mpi_write_restart_file()
    ! proper parallel io - all threads/procs write data to single file
    ! without any communication / non-blocking
    use mpicom
    use dim
    use core
    use time
    use runparam, only: iomod, nsteps
    use parameters
    use iofiles
    ! local vars
    integer :: i,j,k,ii, l, ivar, nvars, err, fhandle, res, bufsize, wordsize, prebufsize, prev_itn!, offset
    integer(kind=mpi_offset_kind) :: offset
!     logical :: i_exist
    character (len=32) :: filename, prev_filename
    real(KIND=4), dimension(:,:), allocatable :: buf
    real(KIND=4), dimension(:), allocatable :: prebuf
    real :: timeend, timestart
    ! mpi-io error checking
    character*(MPI_MAX_ERROR_STRING) :: error_string
    integer ::  ierror, resultlen

!    if (myid == 0) write (*,*) 'WRITING PARALLEL I/O RESTART FILE'
!     timestart = MPI_WTIME()
    call cpu_time(timestart)
    
!     wordsize = 8 ! memory size of integer and real
    prebufsize = 9 ! size of prefix file information
    nvars = 20 !8 ! number of variables to store
    res = nx*ny*nzp ! total size of variables
    bufsize = nxpl*ny*nzpl ! size of each variable

    if (.not. allocated(buf)) allocate(buf(bufsize,nvars), prebuf(prebufsize), stat=err)
    if (err /= 0) print *, "buf: Allocation request denied"
    
    ! prebuf contains useful information about reading the file
    prebuf = [ real(isum), real(t), real(dt), real(dt2), real(dt3), real(nx), real(ny), real(nzp), real(nxpl) ]

    ! copy all velocities into buffer
    ii = 1
    do i = 1, nxpl
        do j = 1, ny
            do k = 1, nzpl 
                ! store in z, y, x fashion to have contiguous blocks 
                ! to allow different numbers of processors to retrieve same data without communication
                buf(ii,1) = uf(i,j,k)
                buf(ii,2) = vf(i,j,k)
                buf(ii,3) = wf(i,j,k)
                buf(ii,4) = tempf(i,j,k)
                do l = 1, 3
                    buf(ii,4+l) = unif(i,j,k,l)
                    buf(ii,7+l) = unmif(i,j,k,l)
                    buf(ii,10+l) = nlnif(i,j,k,l)
                    buf(ii,13+l) = nlnif(i,j,k,l)
                end do
!                 buf(ii,5) = unm(i,j,k)
!                 buf(ii,6) = vnm(i,j,k)
!                 buf(ii,7) = wnm(i,j,k)
                buf(ii,17) = tnf(i,j,k)
                buf(ii,18) = tnmf(i,j,k)
                buf(ii,19) = nltnf(i,j,k)
                buf(ii,20) = nltnmf(i,j,k)
                ii = ii + 1
            end do
        end do
    end do

    ! set word size when writing
    wordsize = sizeof(buf(1,1))
    
    ! rename previous restart file if present - too slow!
!     filename = restart_file ! restart file name
!     if (myid==0) then
!         inquire(file=filename, exist=i_exist) ! check if a restart file exists 
!         if (i_exist .eqv. .true.) then
!             prev_itn = int((isum-iomod)/iomod)
!             ! if it exists and was likely created by this sim
!             ! then change its name appropriately
!             if (prev_itn > 0) then 
!                 write(prev_filename,'(I5.5)') prev_itn
!                 prev_filename = 'start'//trim(prev_filename)//'.dat'
!             else ! it was probably from a different sim with different parameters
!                 ! just call it prev_start.dat
!                 prev_filename = 'prev_'//filename
!             end if
!             call rename(filename,prev_filename) ! rename the existing start file
!         end if
!     end if
    
!     call mpi_barrier(comm,ierr) ! to avoid over-writing previous restart file

    ! just name the start file based on iomod and isum
    if (istep==nsteps) then
        filename = 'start.dat' !restart_file
    else
        prev_itn = int(isum/iomod)
        write(prev_filename,'(I5.5)') prev_itn
        filename = 'start'//trim(prev_filename)//'.dat'
    end if
    
    ! open mpi file to write on
    call mpi_file_open(comm, trim(filename), mpi_mode_wronly + mpi_mode_create, mpi_info_null, fhandle, ierr)
    if (ierr/=0) then
        call mpi_error_string(ierr, error_string, resultlen, ierror)
        if (myid==0) write(*,*) "MPI_FILE_OPEN Error: ",error_string," on ",filename
    end if
    
    ! all procs write the prebuf to file header
    offset = 0
!     call mpi_file_set_view(fhandle, offset, mpi_real, mpi_integer, 'external32', mpi_info_null, ierr)
    call mpi_file_set_view(fhandle, offset, mpi_real, mpi_integer, 'native', mpi_info_null, ierr)
    if (ierr/=0) then
        call mpi_error_string(ierr, error_string, resultlen, ierror)
        if (myid==0) write(*,*) "MPI_FILE_SET_VIEW prebuf Error: ",error_string," on ",filename
    end if

    call mpi_file_write(fhandle, prebuf, prebufsize, mpi_real, mpi_status_ignore, ierr)
    if (ierr/=0) then
        call mpi_error_string(ierr, error_string, resultlen, ierror)
        if (myid==0) write(*,*) "mpi_file_write prebuf Error: ",error_string," on ",filename
    end if

    ! each proc now writes its portion of each variable in order
    do ivar=1,nvars
        offset = (prebufsize + (ivar-1)*res + myid * bufsize ) * wordsize
!         call mpi_file_set_view(fhandle, offset, mpi_real, mpi_integer, 'external32', mpi_info_null, ierr)
        call mpi_file_set_view(fhandle, offset, mpi_real, mpi_integer, 'native', mpi_info_null, ierr)
        if (ierr/=0) then
            call mpi_error_string(ierr, error_string, resultlen, ierror)
            if (myid==0) write(*,*) "MPI_FILE_SET_VIEW buf Error: ",error_string," on ",filename
        end if

        call mpi_file_write(fhandle, buf(:,ivar), bufsize, mpi_real, mpi_status_ignore, ierr)
        if (ierr/=0) then
            call mpi_error_string(ierr, error_string, resultlen, ierror)
            if (myid==0) write(*,*) "mpi_file_write buf Error: ",error_string," on ",filename
        end if
    end do
    
    call mpi_file_close(fhandle, ierr)
    if (ierr/=0) then
        call mpi_error_string(ierr, error_string, resultlen, ierror)
        if (myid==0) write(*,*) "MPI_FILE_CLOSE Error: ",error_string," on ",filename
    end if

    if (allocated(buf)) deallocate(buf, prebuf, stat=err)
    if (err /= 0) print *, "buf: Deallocation request denied"

!     timeend = MPI_WTIME()
    call cpu_time(timeend)

    if (myid==0) write(*,'(A,G12.5,A)') 'parallel restart file written in ',timeend - timestart, ' seconds.'

end subroutine mpi_write_restart_file


subroutine mpi_read_restart_file()
    ! proper parallel io - all threads/procs read data from single file
    ! without any communication
    use mpicom
    use flags, only: debug
    use dim
    use core
    use time
    use runparam, only: iomod
    use parameters
    use iofiles
    use io
    use formatstrings
    ! local vars
    integer :: i,j,k,ii, l, ivar, nvars, err, fhandle, res, bufsize, wordsize, prebufsize
    integer :: nx_in, ny_in, nzp_in, nxpl_in
    integer(kind=mpi_offset_kind) :: offset
    logical :: i_exist
    character (len=9) :: filename
    real(KIND=4), dimension(:,:), allocatable :: buf
    real(KIND=4), dimension(:), allocatable :: prebuf
    real :: timeend, timestart
    ! mpi-io error checking
    character*(MPI_MAX_ERROR_STRING) :: error_string
    integer ::  ierror, resultlen

!     timestart = MPI_WTIME()
    call cpu_time(timestart)

!     wordsize = 4 ! memory size of integer and real
    prebufsize = 9 ! size of prefix file information
    nvars = 20 !8 ! number of variables to store
    res = nx*ny*nzp ! total size of variables
    bufsize = nxpl*ny*nzpl ! size of each variable
    
    if (myid == 0) write(*,'(A)') 'READING PARALLEL I/O RESTART FILE'
    
    if (.not. allocated(buf)) allocate(buf(bufsize,nvars), prebuf(prebufsize), stat=err)
    if (err /= 0) print *, "buf: Allocation request denied"

    wordsize = sizeof(buf(1,1))
    
    filename = 'start.dat' !restart_file
    ! check that file exists and if not, exit gracefully
    inquire(file=filename, exist=i_exist) ! check if a restart file exists
    if (i_exist .neqv. .true.) write(*,*) "** Error: no restart file 'start.dat' available **"
    ! open file using mpi
    call mpi_file_open(comm, filename, mpi_mode_rdonly, mpi_info_null, fhandle, ierr)
    if (ierr/=0) then
        call mpi_error_string(ierr, error_string, resultlen, ierror)
        if (myid==0) write(*,*) "MPI_FILE_OPEN Error: ",error_string," on ",filename
    end if
    
    ! all procs read the prebuf file header
    offset = 0
!     call mpi_file_set_view(fhandle, offset, mpi_real, mpi_integer, 'external32', mpi_info_null, ierr)
    call mpi_file_set_view(fhandle, offset, mpi_real, mpi_integer, 'native', mpi_info_null, ierr)
    if (ierr/=0) then
        call mpi_error_string(ierr, error_string, resultlen, ierror)
        if (myid==0) write(*,*) "MPI_FILE_SET_VIEW prebuf Error: ",error_string," on ",filename
    end if

    call mpi_file_read(fhandle, prebuf, prebufsize, mpi_real, mpi_status_ignore, ierr)
    if (ierr/=0) then
        call mpi_error_string(ierr, error_string, resultlen, ierror)
        if (myid==0) write(*,*) "mpi_file_read prebuf Error: ",error_string," on ",filename
    end if

    ! assign prebuf to associated variables
    isum = nint( prebuf(1) )
    t = prebuf(2)
    dt = prebuf(3)
    dt2 = prebuf(4)
    dt3 = prebuf(5)
    nx_in = nint( prebuf(6) )
    ny_in = nint( prebuf(7) )
    nzp_in = nint( prebuf(8) )
    nxpl_in = nint( prebuf(9) )
    
    ! error checking
    if (nx_in/=nx .or. ny_in/=ny .or. nzp_in/=nzp) then
        if ( myid==0 ) then
            write(*,*) "**error reading restart file: dimension mismatch**"
            write(*,*) '(nx_in ny_in nzp_in) = (',nx_in, ny_in, nzp_in,')'
            write(*,*) '(nx ny nzp) = (',nx,ny,nzp,')'
        end if
        call mpi_barrier(comm,ierr)
        call mpi_finalize(ierr)
        write(*,*) 
    end if
    
    ! each proc now reads its portion of each variable in order
    do ivar=1,nvars
        offset = (prebufsize + (ivar-1)*res + myid * bufsize ) * wordsize
!         call mpi_file_set_view(fhandle, offset, mpi_real, mpi_integer, 'external32', mpi_info_null, ierr)
        call mpi_file_set_view(fhandle, offset, mpi_real, mpi_integer, 'native', mpi_info_null, ierr)
        if (ierr/=0) then
            call mpi_error_string(ierr, error_string, resultlen, ierror)
            if (myid==0) write(*,*) "MPI_FILE_SET_VIEW buf Error: ",error_string," on ",filename
        end if

        call mpi_file_read(fhandle, buf(:,ivar), bufsize, mpi_real, mpi_status_ignore, ierr)
        if (ierr/=0) then
            call mpi_error_string(ierr, error_string, resultlen, ierror)
            if (myid==0) write(*,*) "mpi_file_read buf Error: ",error_string," on ",filename
        end if
    end do
    
!    call mpi_barrier(comm,ierr)
    
    ! organize read data into variables in correct order
    ii = 1
    do i = 1, nxpl
        do j = 1, ny
            do k = 1, nzpl 
                ! store in z, y, x fashion to have contiguous blocks 
                ! to allow different numbers of processors to retrieve same data without communication
                ! pay a cost in memory access being inefficient but no communication
                uf(i,j,k) = buf(ii,1)
                vf(i,j,k) = buf(ii,2)
                wf(i,j,k) = buf(ii,3)
                tempf(i,j,k) = buf(ii,4)
!                 unm(i,j,k) = buf(ii,5)
!                 vnm(i,j,k) = buf(ii,6)
!                 wnm(i,j,k) = buf(ii,7)
!                 tnf(i,j,k) = buf(ii,8)
                do l = 1, 3
                    unif(i,j,k,l) = buf(ii,4+l)
                    unmif(i,j,k,l) = buf(ii,7+l)
                    nlnif(i,j,k,l) = buf(ii,10+l)
                    nlnmif(i,j,k,l) = buf(ii,13+l)
                end do
                tnf(i,j,k) = buf(ii,17)
                tnmf(i,j,k) = buf(ii,18)
                nltnf(i,j,k) = buf(ii,19)
                nltnmf(i,j,k) = buf(ii,20)
                ii = ii + 1
            end do
        end do
    end do
    

    if (allocated(buf)) deallocate(buf, prebuf, stat=err)
    if (err /= 0) print *, "buf: Deallocation request denied"
    
    call mpi_file_close(fhandle, ierr) ! close restart file
    if (ierr/=0) then
        call mpi_error_string(ierr, error_string, resultlen, ierror)
        if (myid==0) write(*,*) "MPI_FILE_CLOSE Error: ",error_string," on ",filename
    end if

!     timeend = MPI_WTIME()
    call cpu_time(timeend)

    if (debug>=2) then
        if (myid==0) write(*,'(/,A)') 'uf, vf, wf from restart file'
!        if (myid==0) then
!            write(*,'[/,A,G12.4,A)') '-----------uf(kx,0,z) at t=',t,'-------------'
!            write(*,fm1f) ((uf(i,1,k),i=1,nxpl),k=1,nzp)
!            write(*,'[/,A,G12.4,A)') '-----------vf(kx,0,z) at t=',t,'-------------'
!            write(*,fm1f) ((vf(i,1,k),i=1,nxpl),k=1,nzp)
!            write(*,'[/,A,G12.4,A)') '-----------wf(kx,0,z) at t=',t,'-------------'
!            write(*,fm1f) ((wf(i,1,k),i=1,nxpl),k=1,nzp)
!        end if
        call printuvw(uf,vf,wf)
    end if

    if (myid==0) write(*,'(A,A,A,G12.5,A)') 'restart file ',trim(filename),' read in ',timeend - timestart, ' seconds.'
    
!    call mpi_barrier(comm,ierr)

end subroutine mpi_read_restart_file


subroutine mpi_snapshot(uf,vf,wf,tempf,pf,name)
    use mpicom
    use dim
!    use core
    use time
    use runparam, only: iomod
    use parameters
    use paral, only: ubar, tbar
    use modhorfft
    use iofiles
    ! i/o
    integer, intent(in) :: name
    real, dimension(nxpl,ny,nzpl), intent(inout) :: uf,vf,wf,tempf,pf
    ! local vars
    character (len=40) :: fout, filenum
    integer(kind=mpi_offset_kind) :: offset
!    integer*4 :: offset
    integer :: i,j,k,ii, ivar, err, fhandle, prebufsize, bufsize, wordsize, res, nvars
    logical :: i_exist
    real(KIND=4), dimension(:,:), allocatable :: buf
    real(KIND=4), dimension(:), allocatable :: prebuf
    real :: timeend, timestart 
    integer :: AllocateStatus
    real, allocatable, dimension(:,:,:) :: u, v, w, temp, dynp
    ! mpi-io error checking
    character*(MPI_MAX_ERROR_STRING) :: error_string
    integer ::  ierror, resultlen

!     timestart = MPI_WTIME()
    call cpu_time(timestart)
    
    prebufsize = 7
    nvars = 5
    res = nx*ny*nzp
    bufsize = nx*nypl*nzp ! size of array each processor writes
    
!    write(*,'(A,I5)') 'myid =',myid

    if (.not. allocated(buf)) allocate(buf(bufsize,nvars), prebuf(prebufsize), stat=err)
    if (err /= 0) print *, "buf: Allocation request denied"
    
    ! allocate physical arrays
    allocate( u(nxpp,nypl,nzpl), v(nxpp,nypl,nzpl), w(nxpp,nypl,nzpl), temp(nxpp,nypl,nzpl), dynp(nxpp,nypl,nzpl), stat=AllocateStatus )
    if (AllocateStatus /= 0) then
        write(*,*) "**Not Enough Memory - mpi_snapshot core arrays allocation **"
    end if

    ! transfer all qties to phys space
    call horfft(u,uf,1)
    call horfft(v,vf,1)
    call horfft(w,wf,1)
    call horfft(temp,tempf,1)
    call horfft(dynp,pf,1)
    
    ! add base flow velocity to u
    if (name/=2) then
        do k=1,nzp
            u(:,:,k) = u(:,:,k) + ubar(k)
            temp(:,:,k) = temp(:,:,k) + tbar(k)
        end do
    end if

    ! write filename using current iteration and iomod
    write(filenum,'(I5.5)') int(isum/iomod)
    if (name==1) then
        fout = 'qavg'//trim(filenum)//'.dat'
    else if (name==2) then
        fout = 'q_ref.dat'
    else
        fout = 'q'//trim(filenum)//'.dat'
    end if
!     if (myid == 0) write (*,*) 'WRITING PARALLEL I/O PHYSICAL FIELDS in '//fout
     
     ! if file already exist, rename instead of overwriting
     if (myid==0) then
        inquire(file=fout, exist=i_exist) ! check if file already exists 
        if (i_exist .eqv. .true.) call rename(fout, 'prev_'//fout) ! rename instead of overwriting
    end if
        
    ! prebuf contains useful information about reading the file
    prebuf = [ real(isum), real(t), real(dt), real(nx), real(ny), real(nzp), real(nypl) ]
    
    ! copy all velocities into buffer
    ii = 1
    do j = 1, nypl
        do k = 1, nzp
            do i = 1, nx
                ! store in x, z, y fashion to have contiguous blocks 
                ! to allow different numbers of processors to retrieve same data without communication
                buf(ii,1) = u(i,j,k)
                buf(ii,2) = v(i,j,k)
                buf(ii,3) = w(i,j,k)
                buf(ii,4) = temp(i,j,k)
                buf(ii,5) = dynp(i,j,k)
                ii = ii + 1
            end do
        end do
    end do

    call mpi_barrier(comm,ierr) ! to ensure file is not overwritten

    ! open mpi file to write on
    call mpi_file_open(comm, fout, mpi_mode_wronly + mpi_mode_create, mpi_info_null, fhandle, ierr)
    if (ierr/=0) then
        call mpi_error_string(ierr, error_string, resultlen, ierror)
        if (myid==0) write(*,*) "MPI_FILE_OPEN Error: ",error_string," on ", fout
    end if
    ! all procs write the prebuf to file header
    offset = 0
    call mpi_file_set_view(fhandle, offset, mpi_real, mpi_integer, 'native', mpi_info_null, ierr)
    if (ierr/=0) then
        call mpi_error_string(ierr, error_string, resultlen, ierror) 
        if (myid==0) write(*,*) "MPI_FILE_SET_VIEW prebuf Error: ",error_string," on ", fout
    end if

    call mpi_file_write(fhandle, prebuf, prebufsize, mpi_real, mpi_status_ignore, ierr)
    if (ierr/=0) then
        call mpi_error_string(ierr, error_string, resultlen, ierror) 
        if (myid==0) write(*,*) "mpi_file_write prebuf Error: ",error_string," on ", fout
    end if
    ! each proc now writes its portion of each variable in order
    wordsize = sizeof(buf(1,1))
!    offset = ( prebufsize + myid*bufsize )*wordsize
!    write(*,'(A,I5,A,I8)') 'myid =',myid,'offset =',offset
    do ivar=1,nvars
        offset = (prebufsize + (ivar-1)*res + myid * bufsize ) * wordsize
        
!        write(*,'(A,I5,A,I5,A,I8)') 'myid =',myid,', ivar =',ivar,', offset =',offset
        
        call mpi_file_set_view(fhandle, offset, mpi_real, mpi_integer, 'native', mpi_info_null, ierr)
        if (ierr/=0)  then
            call mpi_error_string(ierr, error_string, resultlen, ierror) 
            if (myid==0) write(*,*) "MPI_FILE_SET_VIEW buf Error: ",error_string," on ", fout
        end if

        call mpi_file_write(fhandle, buf(:,ivar), bufsize, mpi_real, mpi_status_ignore, ierr)
        if (ierr/=0) then
            call mpi_error_string(ierr, error_string, resultlen, ierror) 
            if (myid==0) write(*,*) "mpi_file_write buf Error: ",error_string," on ", fout
        end if
!        offset = offset + res*wordsize
    end do

    call mpi_file_close(fhandle, ierr)
    if (ierr/=0) then
        call mpi_error_string(ierr, error_string, resultlen, ierror) 
        if (myid==0) write(*,*) "MPI_FILE_CLOSE Error: ",error_string," on ", fout
    end if

    if (allocated(buf)) deallocate(buf, prebuf, u, v, w, dynp, temp, stat=err)
    if (err /= 0) print *, "buf: Deallocation request denied"

!     timeend = MPI_WTIME()
    call cpu_time(timeend)

    if (myid==0) write(*,'(A,G12.5,A)') 'mpi_snapshot file '//trim(fout)//' written in ',timeend - timestart, ' seconds.'


end subroutine mpi_snapshot


subroutine mpi_read_snapshot(q,n1,n2pl,n3,n2,name)
    use mpicom
    use dim
!    use core
    use time, only: dt
    use runparam, only: iomod
    use parameters
    use paral, only: ubar
    use modhorfft
    use iofiles
    ! i/o
    integer, intent(in) :: n1,n2,n3,n2pl,name
    real, dimension(n1,n2pl,n3,5), intent(out) :: q
    ! local vars
    character (len=40) :: fout, filenum
    integer(kind=mpi_offset_kind) :: offset
!    integer*4 :: offset
    integer :: i,j,k,l,ii, ivar, err, fhandle, prebufsize, bufsize, wordsize, res, nvars
    integer :: nx_in, ny_in, nzp_in, nypl_in, isum_in
    logical :: i_exist
    real(KIND=4), dimension(:,:), allocatable :: buf
    real(KIND=4), dimension(:), allocatable :: prebuf
    real :: timeend, timestart , t_in, dt_in
    integer :: AllocateStatus
    ! mpi-io error checking
    character*(MPI_MAX_ERROR_STRING) :: error_string
    integer ::  ierror, resultlen

!     timestart = MPI_WTIME()
    call cpu_time(timestart)
    
    prebufsize = 7
    nvars = 5
    res = nx*ny*nzp
    bufsize = nx*nypl*nzp ! size of array each processor writes
    
!    write(*,'(A,I5)') 'myid =',myid

    if (.not. allocated(buf)) allocate(buf(bufsize,nvars), prebuf(prebufsize), stat=err)
    if (err /= 0) print *, "buf: Allocation request denied"

    ! write filename using current iteration and iomod
    if (name==1) then
        fout = 'qavg_ref.dat'
    else
        fout = 'q_ref.dat'
    end if
!     if (myid == 0) write (*,*) 'WRITING PARALLEL I/O PHYSICAL FIELDS in '//fout
     
       
    ! open mpi file to read in
    call mpi_file_open(comm, fout, mpi_mode_rdonly, mpi_info_null, fhandle, ierr)
    if (ierr/=0) then
        call mpi_error_string(ierr, error_string, resultlen, ierror)
        if (myid==0) write(*,*) "MPI_FILE_OPEN Error: ",error_string," on ", fout
    end if
    ! all procs write the prebuf to file header
    offset = 0
    call mpi_file_set_view(fhandle, offset, mpi_real, mpi_integer, 'native', mpi_info_null, ierr)
    if (ierr/=0) then
        call mpi_error_string(ierr, error_string, resultlen, ierror) 
        if (myid==0) write(*,*) "MPI_FILE_SET_VIEW prebuf Error: ",error_string," on ", fout
    end if

    call mpi_file_read(fhandle, prebuf, prebufsize, mpi_real, mpi_status_ignore, ierr)
    if (ierr/=0) then
        call mpi_error_string(ierr, error_string, resultlen, ierror) 
        if (myid==0) write(*,*) "mpi_file_write prebuf Error: ",error_string," on ", fout
    end if
    
    ! assign prebuf to associated variables
    isum_in = nint( prebuf(1) )
    t_in = prebuf(2)
    dt_in = prebuf(3)
    nx_in = nint( prebuf(4) )
    ny_in = nint( prebuf(5) )
    nzp_in = nint( prebuf(6) )
    nypl_in = nint( prebuf(7) )

    ! error checking
    if (nx_in/=n1 .or. ny_in/=n2 .or. nzp_in/=nzp) then
        if ( myid==0 ) then
            write(*,*) "**error reading q_ref file: dimension mismatch**"
            write(*,*) '(nx_in ny_in nzp_in) = (',nx_in, ny_in, nzp_in,')'
            write(*,*) '(n1 n2 n3) = (',n1,n2,n3,')'
        end if
        call mpi_barrier(comm,ierr)
        call mpi_finalize(ierr)
        write(*,*) 
    end if

    ! if everything is great, set dt here
    dt = dt_in

    ! each proc now writes its portion of each variable in order
    wordsize = sizeof(buf(1,1))
!    offset = ( prebufsize + myid*bufsize )*wordsize
!    write(*,'(A,I5,A,I8)') 'myid =',myid,'offset =',offset
    do ivar=1,nvars
        offset = (prebufsize + (ivar-1)*res + myid * bufsize ) * wordsize
        
!        write(*,'(A,I5,A,I5,A,I8)') 'myid =',myid,', ivar =',ivar,', offset =',offset
        
        call mpi_file_set_view(fhandle, offset, mpi_real, mpi_integer, 'native', mpi_info_null, ierr)
        if (ierr/=0)  then
            call mpi_error_string(ierr, error_string, resultlen, ierror) 
            if (myid==0) write(*,*) "MPI_FILE_SET_VIEW buf Error: ",error_string," on ", fout
        end if

        call mpi_file_read(fhandle, buf(:,ivar), bufsize, mpi_real, mpi_status_ignore, ierr)
        if (ierr/=0) then
            call mpi_error_string(ierr, error_string, resultlen, ierror) 
            if (myid==0) write(*,*) "mpi_file_write buf Error: ",error_string," on ", fout
        end if
!        offset = offset + res*wordsize
    end do

    call mpi_file_close(fhandle, ierr)
    if (ierr/=0) then
        call mpi_error_string(ierr, error_string, resultlen, ierror) 
        if (myid==0) write(*,*) "MPI_FILE_CLOSE Error: ",error_string," on ", fout
    end if


    ii = 1
    do j = 1, nyplin
        do k = 1, nzp_in
            do i = 1, nx_in
                do l=1,nvars
                    ! stored in x,z,y fashion to have contiguous blocks 
                    ! to allow different numbers of processors to retrieve same data without communication
                    ! pay a cost in memory access being inefficient but no communication
                    q(i,j,k,l) = buf(ii,l)
                end do
                ii = ii + 1
            end do
        end do
    end do


    if (allocated(buf)) deallocate(buf, prebuf, stat=err)
    if (err /= 0) print *, "buf: Deallocation request denied"

!     timeend = MPI_WTIME()
    call cpu_time(timeend)

    if (myid==0) write(*,'(A,G12.5,A)') 'mpi_snapshot file '//trim(fout)//' read in ',timeend - timestart, ' seconds.'


end subroutine mpi_read_snapshot

end module parallel_io
