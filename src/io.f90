module io

implicit none
public
save

contains

subroutine write_endofsim_stats(cpu_t1,cpu_t2)
    use time
    use runparam
    use mpicom,only: myid
    ! i/o
    real, intent(in) :: cpu_t1, cpu_t2
    ! local vars
    real :: cpu_time
    integer :: cpu_days, cpu_hours, cpu_minutes, cpu_seconds

    cpu_time = cpu_t2 - cpu_t1
    cpu_days = int(cpu_time/60./60./24.)
    cpu_time = cpu_time - cpu_days*60.*60.*24.
    cpu_hours = int( cpu_time/60./60.0)
    cpu_time = cpu_time - cpu_hours*60.*60.0
    cpu_minutes = int(cpu_time/60.)
    cpu_time = cpu_time - cpu_minutes*60.0
    cpu_seconds = int( cpu_time )
    
    if (myid==0) then
        write(*,'(//,A)') '------------- End of Simulation Stats --------------'
        write(*,'(4(A,I2))') 'cpu time =',cpu_days,' days ',cpu_hours,':',cpu_minutes,':',cpu_seconds
        write(*,'(A,I10)') '# of its performed =', nsteps
        write(*,'(A,G12.4)') 'avg its / sec =',nsteps/(cpu_t2-cpu_t1)
        write(*,'(A,G14.6)') 'time simulated =', tseries(nsteps) - tseries(1)
        write(*,'(A,I10)') 'total # of its performed =', isum
        write(*,'(A,G14.6)') 'final time step dt =', dt
        write(*,'(A,G14.6)') 'total simulation time =', t
    end if
end subroutine write_endofsim_stats


subroutine printuvw(q1f,q2f,q3f)
    ! prints 3 given variable at time t in physical space to screen
    use dim
    use mpicom, only: myid
    use formatstrings
    use time, only: t
    use modhorfft
    ! i/o
    real, dimension(nxpl,ny,nzpl), intent(inout) :: q1f,q2f,q3f
    ! local vars
    integer :: i,k
    real, allocatable, dimension(:,:,:) :: q1, q2, q3
    allocate( q1(nxpp,nypl,nzpl), q2(nxpp,nypl,nzpl), q3(nxpp,nypl,nzpl) )
    ! take all vars to phys space
    call horfft(q1,q1f,1)
    call horfft(q2,q2f,1)
    call horfft(q3,q3f,1)
    ! print them
    if (myid==0) then
        write(*,'(//,A,G12.4,A)') '-----------u(x,1,z) at t=',t,'-------------'
        write(*,fm1) ((q1(i,1,k),i=1,nx),k=1,nzp)
        write(*,'(//,A,G12.4,A)') '-----------v(x,1,z) at t=',t,'-------------'
        write(*,fm1) ((q2(i,1,k),i=1,nx),k=1,nzp)
        write(*,'(//,A,G12.4,A)') '-----------w(x,1,z) at t=',t,'-------------'
        write(*,fm1) ((q3(i,1,k),i=1,nx),k=1,nzp)
    end if
    deallocate(q1,q2,q3)
end subroutine printuvw


subroutine printvars(q1f,q2f,q3f,q4f,q5f)
    ! prints given variables at time t in physical space to screen
    use dim
    use mpicom, only: myid
    use formatstrings
    use time, only: t
    use modhorfft
    ! i/o
    real, dimension(nxpl,ny,nzpl), intent(inout) :: q1f,q2f,q3f,q4f,q5f
    !     real, intent(in) :: t
    ! local vars
    integer :: i,k
    real, allocatable, dimension(:,:,:) :: q1, q2, q3, q4, q5
    allocate( q1(nxpp,nypl,nzpl), q2(nxpp,nypl,nzpl), q3(nxpp,nypl,nzpl), q4(nxpp,nypl,nzpl), q5(nxpp,nypl,nzpl) )
    ! take all vars to phys space
    call horfft(q1,q1f,1)
    call horfft(q2,q2f,1)
    call horfft(q3,q3f,1)
    call horfft(q4,q4f,1)
    call horfft(q5,q5f,1)

!     q5 = q5 - 0.5*( q1*2.0 + q2**2.0 + q3**2.0 ) ! du = 0.5*|v|^2
    ! print them
    if (myid==0) then
        write(*,'(//,A,G12.4,A)') '-----------u(x,1,z) at t=',t,'-------------'
        write(*,fm1) ((q1(i,1,k),i=1,nx),k=1,nzp)
        write(*,'(//,A,G12.4,A)') '-----------v(x,1,z) at t=',t,'-------------'
        write(*,fm1) ((q2(i,1,k),i=1,nx),k=1,nzp)
        write(*,'(//,A,G12.4,A)') '-----------w(x,1,z) at t=',t,'-------------'
        write(*,fm1) ((q3(i,1,k),i=1,nx),k=1,nzp)
        write(*,'(//,A,G12.4,A)') '-----------temp(x,1,z) at t=',t,'-------------'
        write(*,fm1) ((q4(i,1,k),i=1,nx),k=1,nzp)
        write(*,'(//,A,G12.4,A)') '-----------p(x,1,z) at t=',t,'-------------'
        write(*,fm1) ((q5(i,1,k),i=1,nx),k=1,nzp)
    end if
    deallocate(q1,q2,q3,q4,q5)
end subroutine printvars


subroutine write_grid()
    ! write out grid file for restarting
    use mpicom, only: myid
    use dim
    use grid
    use parameters, only: xlen, ylen, zlen, x0
    use iofiles
    ! local vars
    integer :: i, ii
    ii = iunit_grid
    if (myid==0) then
        open(ii,file=grid_file,form='unformatted',status='unknown')
        ! write out grid parameters
        write(ii) nx,ny,nzp
        write(ii) x0,xl,xlen,ylen,zlen
        write(ii) (xpts(i),i=1,nx)
        write(ii) (ypts(i),i=1,ny)
        write(ii) (zpts(i),i=1,nzp)
        write(ii) (gf(i),i=1,nzp)
        close(ii)
    end if
end subroutine write_grid

subroutine read_grid()
    ! write out grid file for restarting
    use dim
    use mpicom
    use grid
    use flags, only: debug, nstart
    use parameters, only: xlen, ylen, zlen, x0
    use iofiles
    use formatstrings
    ! local vars
    logical :: i_exist
    integer :: i,j,k, ii!, nxin, nyin, nzpin
    real :: twopi, xlenin, ylenin, zlenin, x0in
   
    inquire(file=grid_file, exist=i_exist) ! check if file already exists 
    if (i_exist .eqv. .false.) write(*,*) "** error: no grid file in folder **"

    ii = 130 + myid
    open(ii,file=grid_file,form='unformatted',status='old')
    ! write out grid parameters
    read(ii) nxin,nyin,nzpin
    read(ii) x0in,xl,xlenin,ylenin,zlenin
    ! allocate input grid
    nyplin = nyin/nproch
    if (.not. allocated(xin)) allocate( xin(nxin),yin(nyin),yin_loc(nyplin),zin(nzpin),gfin(nzpin) )
    read(ii) (xin(i),i=1,nxin)
    read(ii) (yin(j),j=1,nyin)
    read(ii) (zin(k),k=1,nzpin)
    read(ii) (gfin(k),k=1,nzpin)
    close(ii)

    ! proc specific y coord in phys space
    do j=1,nyplin
        yin_loc(j) = yin(myid*nyplin + j)
    end do

    if (debug>=1 .and. myid==0) then
        write(*,*) 'xin ='
        write(*,fm1) xin(1:nx)
        write(*,*) 'yin ='
        write(*,fm2) yin(1:ny)
        write(*,*) 'zin, gfin'
        write(*,'(2(G12.4,'',''))') (zin(k),gfin(k),k=1,nzp)
    end if

    if (nstart==1) then
    
        ! do some dimension error checking
        if ((nxin /= nx).or.(nyin /= ny).or.(nzpin /= nzp)) then
            if (myid == 0) then
                write(*,*) 'input resolution in x,y or z direction not consistent w/ grid file!'
                write(*,*) 'nx, ny, nzp in grid file = ',nxin,nyin,nzpin
                write(*,*) 'nx, ny, nzp in input = ',nx,ny,nzp
                write(*,*) 'Program is write(*,*)ped. TRY AGAIN !'
            endif
            close(ii)
            call MPI_FINALIZE(ierr)
        endif
        ! domain dimensions
!    x0in = xpts(1)
!    xlenin = xpts(nx) - xpts(1)
!    ylenin = ypts(ny) - ypts(1)
!    zlenin = zpts(1) - zpts(nzp)
        ! write(*,*) if dimensions are inconsistent
        if ((xlenin /= xlen).or.(ylenin /= ylen).or.(zlenin /= zlen)) then
            if (myid == 0) then
                write(*,*) 'Physical dimensions of computational domain are not consistent !'
                write(*,*) 'x0, xlen, ylen, zlen in grid file = ',x0in,xlenin,ylenin,zlenin
                write(*,*) 'x0, xlen, ylen, zlen in input = ',x0,xlen,ylen,zlen
            endif
!        x0 = x0in
!        xlen = xlenin
!        ylen = ylenin
!        zlen = zlenin
        endif
        ! allocate grid vars
        if (.not. allocated(xpts)) allocate(xpts(nx),ypts(ny),y(nypl),zpts(nzp),wavx(nxpl),wavy(ny),wavxx(nxpl),wavyy(ny),gf(nzp))
        ! if everything is good transfer input grid to actual grid
        xpts = xin
        ypts = yin
        y = yin_loc
        zpts = zin
        gf = gfin
    
        ! set up rest of grid parameters here    
        dx = xlen/(nx-1)
        dy = ylen/(ny-1)
        ! set up wave numbers
        twopi = 8.*atan(1.)
        alpha = twopi/xlen
        do i=1,nxhpl
            wavx(2*i-1:2*i) = (myid*nxhpl + i - 1.)*alpha
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
        
        if (debug>=1 .and. myid==0) then
            write(*,*) 'xpts ='
            write(*,fm1) xpts(1:nx)
            write(*,*) 'ypts ='
            write(*,fm2) ypts(1:ny)
            write(*,*) 'zpts, gf'
            write(*,'(2(G12.4,'',''))') (zpts(k),gf(k),k=1,nzp)
        end if

    end if

end subroutine read_grid


subroutine write_paral()
    ! write out grid file for restarting
    use dim
    use mpicom, only: myid
    use paral
    use iofiles
    ! local vars
    integer :: i, ii
    ii = iunit_paral
    if (myid==0) then ! only process 0 outputs file
        open(ii,file=paral_file,form='unformatted',status='unknown')
        ! write out parallel qties
        write(ii) nzp
        write(ii) (ubar(i),i=1,nzp)
        write(ii) (ubarp(i),i=1,nzp)
        write(ii) (ubarpp(i),i=1,nzp)
        write(ii) (wbar(i),i=1,nzp)
        write(ii) (wbarp(i),i=1,nzp)
        write(ii) (wbarpp(i),i=1,nzp)
        write(ii) (tbar(i),i=1,nzp)
        write(ii) (tbarp(i),i=1,nzp)
        write(ii) (tbarpp(i),i=1,nzp)
        close(ii)
    end if
end subroutine write_paral

subroutine read_paral()
    ! write out grid file for restarting
    use dim
    use mpicom
    use paral
    use iofiles
    ! local vars
    integer :: i, ii
    logical :: i_exist

    inquire(file=paral_file, exist=i_exist) ! check if file already exists 
    if (i_exist .eqv. .false.) write(*,*) "** error: no paral file in folder **"
    
    ii = 500 + myid ! each process opens and reads the file
    open(ii,file=paral_file,form='unformatted',status='old')
    ! write out parallel qties
    read(ii) nzpin
    ! do some dimension error checking
    if (nzpin /= nzp) then
        if (myid == 0) then
            write(*,*) 'input resolution in z direction not consistent w/ paral file!'
            write(*,*) 'nzp in grid file = ',nzpin
            write(*,*) 'nzp in input = ',nzp
            write(*,*) 'Program is write(*,*)ped. TRY AGAIN !'
        endif
        close(ii)
        call MPI_FINALIZE(ierr)
        write(*,*)
    endif
    read(ii) (ubar(i),i=1,nzp)
    read(ii) (ubarp(i),i=1,nzp)
    read(ii) (ubarpp(i),i=1,nzp)
    read(ii) (wbar(i),i=1,nzp)
    read(ii) (wbarp(i),i=1,nzp)
    read(ii) (wbarpp(i),i=1,nzp)
    read(ii) (tbar(i),i=1,nzp)
    read(ii) (tbarp(i),i=1,nzp)
    read(ii) (tbarpp(i),i=1,nzp)
    close(ii)
end subroutine read_paral



subroutine write_restart_file()
    ! Simplified version, where each processor dumps out its own file
    ! containing specific block of data.
    use mpicom
    use dim
    use core
    use grid
    use paral
    use time
    use runparam, only: iomod,nsteps
    use parameters
    use iofiles
    ! local vars
    integer :: i,j,k, prev_itn, l
    logical :: i_exist
    character (len=160) :: fout, prev_fout
    real :: timeend, timestart

    !***TAK 4-26-2012: now automatically adjust the length of the file name
    !     with extension having N_CHAR_EXTENT characters. following variables
    !     are added.
    integer, parameter :: N_CHAR_EXTENT=4
    integer :: foutlen,ip1,ip2, iunit
    !***TAK 4-26-2012

    character (len=8) :: fmt1

    ! Create output file
    ! CAREFUL: Naming scheme can handle up to 999 processors !

    if (myid == 0) write (*,*) 'DUMPING OUT RESTART FILE'
!     timestart = MPI_WTIME()
    call cpu_time(timestart)
!     fout = RESTART_FILE
    if (istep==nsteps) then
        fout = 'start.dat' !restart_file
    else
        prev_itn = int(isum/iomod)
        write(prev_fout,'(I5.5)') prev_itn
        fout = 'start'//trim(prev_fout)//'.dat'
    end if

!     write(fout,'(I5.5)') int(isum/iomod)
!     fout = 'start'//trim(fout)//'.dat'
    !    fout= trim(output_dir)//trim(fout)

    !***TAK 4-26-2012: automatically adjust the length of the file name
    !     using pointers ip1 and ip2.0
    foutlen=len_trim(fout)
    ip1=foutlen+1
    ip2=foutlen+N_CHAR_EXTENT

    if ((myid >= 0) .and. (myid < 10) )        fmt1 = '(a3,i1)'
    if ((myid >= 10) .and. (myid < 100))     fmt1 = '(a2,i2)'
    if ((myid >= 100) .and. (myid < 1000)) fmt1 = '(a1,i3)'

    if ((myid >= 0) .and. (myid < 10) ) write(unit=fout(ip1:ip2),fmt=fmt1) '_00',myid ! was fout(36:39)
    if ((myid >= 10) .and. (myid < 100)) write(unit=fout(ip1:ip2),fmt=fmt1) '_0',myid ! was fout(36:39)
    if ((myid >= 100) .and. (myid < 1000)) write(unit=fout(ip1:ip2),fmt=fmt1) '_',myid ! was fout(36:39)

    if (myid == 0) write(*,*) 'Filename: ',trim(fout)

    ! copy previous restart file to prev_start.dat_00
    inquire(file=fout, exist=i_exist)
    if (i_exist .eqv. .true.) then
        call rename(fout,'prev_'//fout)
    end if
    ! Open output file (each processor does that separately, on his own channel)
    iunit = 2000+myid
    open(iunit,file=fout,form='unformatted',status='unknown')
    !    write(*,*) 'proc & file',myid,fout

    ! Write to output file
    write(iunit) t,dt,dt2,dt3,isum
    write(iunit) nx,ny,nz,nproch,nprocv
    !
    ! velocity and temperature at (n) time level in spectral space
    !
    write(iunit) ( ((uf(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
    write(iunit) ( ((vf(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
    write(iunit) ( ((wf(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
    write(iunit) ( ((tempf(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
    !
    !
    ! non-linear velocity and temperature advection terms
    ! at (n-1) time level in spectral space
    ! (for Adams-Bashforth)
    ! unif, unmif, nlnif, nlnmif
    do l = 1, 3
        write(iunit) ( ((unif(i,j,k,l),i=1,nxpl),j=1,ny),k=1,nzpl )
        write(iunit) ( ((unmif(i,j,k,l),i=1,nxpl),j=1,ny),k=1,nzpl )
        write(iunit) ( ((nlnif(i,j,k,l),i=1,nxpl),j=1,ny),k=1,nzpl )
        write(iunit) ( ((nlnmif(i,j,k,l),i=1,nxpl),j=1,ny),k=1,nzpl )
    end do
!     write(iunit) ( ((unm(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
!     write(iunit) ( ((vnm(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
!     write(iunit) ( ((wnm(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
    write(iunit) ( ((tnf(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)

    ! Close output file
    close(iunit)

!     timeend = MPI_WTIME()
    call cpu_time(timeend)

    ! Save current time and timestep in log file
!     if (myid == 0) then
!         open(IUNIT_RESTART_LOG,file=RESTART_LOG_FILE,position='append')
!         write(IUNIT_RESTART_LOG,*)'t=',t, ', dt=',dt
!         close(IUNIT_RESTART_LOG)
!     end if
    if (myid == 0) write(*,*) 'Time spent writing restart file',timeend-timestart

    call MPI_BARRIER(comm,ierr)

end subroutine write_restart_file


subroutine read_restart_file()
    ! reads all values from a previous restart file
    use dim
    use grid
    use core
    use time
    use mpicom
    use iofiles
    use formatstrings
    use flags, only: debug
    ! local vars
    logical :: i_exist
    integer :: i,j,k,l,nx_in,ny_in,nz_in,nprochin,nprocvin,iunit
!     real :: tin,xlenin,ylenin,zlenin
!     real :: twopi
    real :: timeend, timestart
    !***TAK 4-26-2012: now automatically adjust the length of the file name
    !     with extension having N_CHAR_EXTENT characters. following variables
    !     are added.
    integer, parameter :: N_CHAR_EXTENT=4 ! TAK 4-26-2012: add
    character (len=160) :: fin ! TAK 4-26-2012: was 39
    integer :: finlen,ip1,ip2 ! TAK 4-26-2012: add
    !**TAK 4-26-2012

    character (len=8) :: fmt1

    ! read from output file
    ! CAREFUL: Naming scheme can handle up to 999 processors !

    call cpu_time(timestart)
    !***TAK 4-26-2012: automatically adjust the length of the file name
    !     using pointers ip1 and ip2.0
    fin = RESTART_FILE
    finlen=len_trim(fin)
    ip1=finlen+1
    ip2=finlen+N_CHAR_EXTENT

    if ((myid >= 0) .and. (myid < 10) )        fmt1 = '(a3,i1)'
    if ((myid >= 10) .and. (myid < 100))     fmt1 = '(a2,i2)'
    if ((myid >= 100) .and. (myid < 1000)) fmt1 = '(a1,i3)'

    if ((myid >= 0) .and. (myid < 10) ) write(unit=fin(ip1:ip2),fmt=fmt1) '_00',myid ! was fout(36:39)
    if ((myid >= 10) .and. (myid < 100)) write(unit=fin(ip1:ip2),fmt=fmt1) '_0',myid ! was fout(36:39)
    if ((myid >= 100) .and. (myid < 1000)) write(unit=fin(ip1:ip2),fmt=fmt1) '_',myid ! was fout(36:39)

    ! check to make sure the restart file exists
    inquire(file=fin, exist=i_exist) ! check if a restart file exists
    if (i_exist .neqv. .true.) then
        if (myid==0) write(*,'(A)') "** Error: no restart file "//fin//" available **"
        call mpi_barrier(comm,ierr)
        call mpi_finalize(ierr)
    end if

    ! Open output file (each processor does that separately, on his own channel)
    iunit = 2000+myid
    open(iunit,file=fin,form='unformatted',status='old')

    ! Write to output file
    read(iunit) t,dt,dt2,dt3,isum
    read(iunit) nx_in,ny_in,nz_in,nprochin,nprocvin

    ! -------------------------------------------

    if (myid == 0) write(*,*) 'PROC-0: Restart time is t=',t

    if (myid == 0) write(*,*) 'nx, ny, nz in restart file: ',nx_in, ny_in, nz_in

    if ((nx_in /= nx).or.(ny_in /= ny).or.(nz_in /= nz)) then
        if (myid == 0) then
            write(*,*) 'Total resolution in x,y or z direction not consistent with input/grid file !'
        endif
        call MPI_FINALIZE(ierr)
    endif

    if ((nproch /= nprochin).or.(nprocv /= nprocvin)) then
        if (myid == 0) then
            write(*,*) 'Number of processors in horizontal or vertical direction not consistent !'
            write(*,*) 'These restart files were generated with ',nprochin,' horizontal procs. and ', nprocvin, ' vertical procs.'
        endif
        call MPI_FINALIZE(ierr)
    endif

    ! velocity and temperature at (n) time level in spectral-physical space
    read(iunit) ( ((uf(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
    read(iunit) ( ((vf(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
    read(iunit) ( ((wf(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
    read(iunit) ( ((tempf(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)

    ! non-linear velocity and temperature advection terms
    ! at (n-1) time level in spectral space
    ! (for Adams-Bashforth)
!     read(iunit) ( ((unm(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
!     read(iunit) ( ((vnm(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
!     read(iunit) ( ((wnm(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
    do l = 1, 3
        write(iunit) ( ((unif(i,j,k,l),i=1,nxpl),j=1,ny),k=1,nzpl )
        write(iunit) ( ((unmif(i,j,k,l),i=1,nxpl),j=1,ny),k=1,nzpl )
        write(iunit) ( ((nlnif(i,j,k,l),i=1,nxpl),j=1,ny),k=1,nzpl )
        write(iunit) ( ((nlnmif(i,j,k,l),i=1,nxpl),j=1,ny),k=1,nzpl )
    end do
    read(iunit) ( ((tnf(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)


    ! Close output file
    close(iunit)

    call MPI_BARRIER(comm,ierr)

    call cpu_time(timeend)

    if (myid == 0) write(*,*) 'Time spent reading restart file',timeend-timestart

end subroutine read_restart_file



subroutine snapshot(uf,vf,wf,tempf,pf,name)
    ! gather all qties in physical space and dump into single file for later
    ! post-processing using any non-parallel post-processor
    use mpicom
    use flags, only: debug
    use dim
    use modhorfft
!     use core
    use paral
    use grid
    use runparam
    use parameters
    use formatstrings
    use time
    use iofiles
    ! input/output
    integer, intent(in) :: name
    real, dimension(nxpl,ny,nzpl), intent(inout) :: uf,vf,wf,tempf,pf
    character (len=160) :: fout, filenum
    ! other vars
    logical :: i_exist
    integer :: nq, i, j , k, l, jglob
    real, allocatable, dimension(:,:,:,:) :: q,qout
    real :: timeend, timestart
    integer :: AllocateStatus
    real, allocatable, dimension(:,:,:) :: u, v, w, temp, dynp

!     timestart = MPI_WTIME()
    call cpu_time(timestart)

    allocate( u(nxpp,nypl,nzpl), v(nxpp,nypl,nzpl), w(nxpp,nypl,nzpl), temp(nxpp,nypl,nzpl), dynp(nxpp,nypl,nzpl), stat=AllocateStatus )
    if (AllocateStatus /= 0) then
        write(*,*) "**Not Enough Memory - writephys core arrays allocation **"
    end if

    allocate(q(nxpp,ny,nzpl,5),qout(nxpp,ny,nzpl,5),stat=AllocateStatus )
    if (AllocateStatus /= 0) then
        write(*,*) "**Not Enough Memory - writephys mpi arrays allocation **"
    end if

    if (myid == 0) write (*,*) 'DUMPING OUT PHYSICAL QTY FILE'
    
    write(filenum,'(I5.5)') int(isum/iomod)
    if (name==1) then
        fout = 'qavg'//trim(filenum)//'.dat'
    else
        fout = 'q'//trim(filenum)//'.dat'
    end if

    if (myid==0) then
        inquire(file=fout, exist=i_exist) ! check if file already exists 
        if (i_exist .eqv. .true.) call rename(fout, 'prev_'//fout) ! rename instead of overwriting
    end if

    ! transfer all qties to phys space
    call horfft(u,uf,1)
    call horfft(v,vf,1)
    call horfft(w,wf,1)
    call horfft(temp,tempf,1)
    call horfft(dynp,pf,1)

!     dynp = dynp/dt - ( u**2 + v**2 + w**2 )

    ! add base flow velocity to u
    do k=1,nzp
        u(:,:,k) = u(:,:,k) + ubar(k)
    end do

    ! zero any part that has no associated value on that processor
    q=0.0
    ! form global data structure to mpi_reduce to proc0
    do j=1,nypl
        jglob = (myid*nypl)+j
        q(:,jglob,:,1) = u(:,j,:)
        q(:,jglob,:,2) = v(:,j,:)
        q(:,jglob,:,3) = w(:,j,:)
        q(:,jglob,:,4) = temp(:,j,:)
        q(:,jglob,:,5) = dynp(:,j,:)
    end do

    if (myid==1 .and. debug>1) then
    ! for debugging just print to file
        do l=1,1
            write(*,*) '****************************'
            write(*,*) 'q(x,4,z,',l,') before reduce'
            write(*,*) '****************************'
            write(*,fm1) ((q(i,4,k,l),i=1,nx),k=1,nzpl)
        end do
    end if
    ! now each proc has a version of u0 with zeros everywhere
    ! except where its global values are located

    ! size of q
    nq = nxpp*ny*nzpl*5
    ! adding together all versions of q into the proc0 version
    call mpi_reduce(q,qout,nq,mpi_double_precision,MPI_SUM,0,comm,ierr)

    if (myid==0) then
    ! now we can move on to printing all qties to file
        open(IUNIT_SOL,file=fout,form='unformatted',status='unknown')

        ! Write time and resolution related qties to file for adaptive reading
        write(IUNIT_SOL) t,dt,xlen,ylen,zlen
        write(IUNIT_SOL) nx,ny,nz,nproch,nprocv

        ! Write avg qties
        do l=1,5
            write(IUNIT_SOL) ( ((qout(i,j,k,l),i=1,nx),j=1,ny),k=1,nzpl)
        end do

        close(IUNIT_SOL)

        ! for debugging just print to file
        if (debug>1) then
            do l=1,5
                write(*,*) '***********************'
                write(*,*) 'q(x,1,z,',l,')'
                write(*,*) '***********************'
                write(*,fm1) ((qout(i,1,k,l),i=1,nx),k=1,nzpl)
            end do
        end if

    end if

    deallocate( u, v, w, temp, q, qout )

!     timeend = MPI_WTIME()
    call cpu_time(timeend)
    if (myid == 0) write(*,'(A,A,G12.4,A)') trim(fout),' written in ',timeend-timestart,'seconds'

end subroutine snapshot



subroutine savegrid()
    ! save grid points in separate .dat files
    use dim
    use grid
    use mpicom, only: myid
    use iofiles

    if ( myid==0 ) then
        open(IUNIT_XGRID,file=XGRID_FILE,status='replace',form='unformatted')
        write(IUNIT_XGRID) nx
        write(IUNIT_XGRID) xpts(1:nx)
        close(IUNIT_XGRID)
        open(IUNIT_YGRID,file=YGRID_FILE,status='replace',form='unformatted')
        write(IUNIT_YGRID) ny
        write(IUNIT_YGRID) ypts(1:ny)
        close(IUNIT_YGRID)
        open(IUNIT_ZGRID,file=ZGRID_FILE,status='replace',form='unformatted')
        write(IUNIT_ZGRID) nzp
        write(IUNIT_ZGRID) zpts(1:nzp)
        write(IUNIT_ZGRID) gf(1:nzp)
        close(IUNIT_ZGRID)
    end if

end subroutine savegrid






end module io