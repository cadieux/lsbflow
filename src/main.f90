program main
    ! 3D N-S Solver
    use mpicom
    use dim
    use flags
    use runparam
    use parameters
    use sgsinput
    use time
    use init
    use navierstokes
    use orrsomm
    use io
    use parallel_io
    use iofiles
    use core
    use tavg
    ! vars
    real :: cpu_t1, cpu_t2
    integer :: ii, nanflag

    nanflag = 0

    ! initialize mpi
    call init_mpi()

    ! Read inputs in parallel
    ii = 100+myid
    open(unit=ii,file=input_file)
    read(ii,dimen)
    read(ii,flaglist)
    read(ii,run)
    read(ii,param)
    read(ii,sgsin)
    close(unit=ii)

    n_re_d = 1; n_alphastar = 1;

    if (linstab==2) then
        n_re_d = 5; n_alphastar = 5; ! set up for linear stability test
    end if

    ! create Reynolds number and alpha* array to loop over
    allocate( re_d(n_re_d), alphastar(n_alphastar) )

    do i_re_d = 1, n_re_d ! loop over Reynolds numbers base on BL thickness
    do j_alphastar = 1, n_alphastar ! loop over different wavelength/domain size alpha*
    
    if (linstab>=1) then
        if (linstab==2) call set_linstab_param() ! set Re and alpha*
        if (.not. allocated(enr1)) allocate( enr1(nsteps), enr2(nsteps), grow1(nsteps), grow2(nsteps) )
    end if
    
    ! Call overarching init/setup subroutine
    call initialize()

    ! measure cpu time before main loop
!     cpu_t1 = mpi_wtime()
    call cpu_time(cpu_t1)

    ! main loop: timestep()
    do istep=1,nsteps

        call timestep() ! integrate Navier-Stokes in time

        ! linear stability growth rate test
        if (linstab>=1) call lintest(uf,vf,wf,enr1,enr2,grow1,grow2)
        
        ! time averaging
        if (otfavg==1 .and. its>istartavg) then
            uf_tavg = uf_tavg + uf/iomod
            vf_tavg = vf_tavg + vf/iomod
            wf_tavg = wf_tavg + wf/iomod
            tempf_tavg = tempf_tavg + tempf/iomod
            pf_tavg = pf_tavg + pf/iomod
        end if

        ! print out if iomod
        if (mod(istep,iomod)==0) then
            ! check for NaN and abort if so
            if (any(isnan(uf)) .or. any(isnan(vf)) .or. any(isnan(wf))) then
                if (myid==0) write(*,*) 'NaN detected at t=',t,'iteration #:',its,'. Abort'
!                 call cpu_time(cpu_t2)
                cpu_t2 = mpi_wtime()
                call write_endofsim_stats(cpu_t1,cpu_t2)
                call mpi_finalize(ierr)
            end if
            
            if (io_mpi==1) then
                call mpi_write_restart_file()
                if (iwrite==1) call mpi_snapshot(uf,vf,wf,tempf,pf,0)
            else
                call write_restart_file()
                if (iwrite==1) call snapshot(uf,vf,wf,tempf,pf,0)
            end if
            
            ! output avg qties
            if ( otfavg==1 .and. its>istartavg ) then
                if (io_mpi==1) then
                    call mpi_snapshot(uf_tavg,vf_tavg,wf_tavg,tempf_tavg,pf_tavg,1)
                else
                    call snapshot(uf_tavg,vf_tavg,wf_tavg,tempf_tavg,pf_tavg,1)
                end if
                ! zero avg arrays
                uf_tavg = 0.; vf_tavg = 0.; wf_tavg = 0.0
                tempf_tavg = 0.; pf_tavg = 0.;
            end if

        end if ! iomod writeout over

    end do ! main loop

    ! measure cpu time after main loop
!     cpu_t2 = mpi_wtime()
    call cpu_time(cpu_t2)
    call write_endofsim_stats(cpu_t1,cpu_t2)
    
    if (linstab>=1) call write_linstab_results(nsteps)

    end do ! j_alphastar
    end do ! i_re_d
    !     ! postproc
    !     call postproc()

    ! end gracefully
    call finalize()

end program main