module dealias

! dealiasing.f90:
! 
! this file contains subroutines needed in dealiasing nonlinear terms 
! in parallel 2d (x-y periodic). implementation of this 
! package requires the following three steps:
!
! 1. initialization
!
!         call setup_dealiasing(comm,myid,1)
!
!         input:
!             comm (integer) : mpi communicator handle
!             myid (integer) : calling proccess id
!             (the last argument must be 1)
! 
! 2. perform dealiasing of a nonlinear product
!
!         call dealiasing(a,b,ab,comm,myid,direction)
!
!         input:
!             a, b (real) : field arrays defined in physical space
!
!         output:
!             ab (real) : dealiased product of arrays a and b (a*b in 
!                                     physical space). 
!
!         note: initialization must be called only once before the first 
!                     call of dealiasing.
!
! 3. finalization
!
!         call seutp_dealiasing(myid,comm,0)
!
!         input:
!             comm (integer) : mpi communicator handle
!             myid (integer) : calling proccess id
!             (the last argument must be 0)
!
!         this deallocates all allocatable arrays that are defined in the 
!         dealiasing package.
!
! requirements:
!
!     1. 2d fft routines (parallel version):
!            horfft (in horfft.f)
!
!     2. library:
!            mpich2 (mpi parallel library)
!            fftw 3.3.2
!
!
! original: programmed by t.sakai 4-27-2012 @ cornell university



!
! this defines shared variables used in dealiasing routines
!
! revisions:
!     10-17-2012 t.sakai: add fftw plan variables.

implicit none
public
save

! array dimensions for dealiasing
integer :: nx_lg
integer :: ny_lg
integer :: nxpp_lg
integer :: nxpl_lg
integer :: nypl_lg
integer :: nxf_lg
integer :: nyf_lg

! fftw plans
integer, parameter, private :: kd=8
integer(kind=kd) :: plan_x_lg_fwd
integer(kind=kd) :: plan_y_lg_fwd
integer(kind=kd) :: plan_x_lg_backwd
integer(kind=kd) :: plan_y_lg_backwd

contains


subroutine setup_dealiasing(ioption)
    ! setup_dealiasing initializes dealiasing operations if ioption=1, 
    ! and finalize (release array memory) if ioption=0.
    !
    ! other inputs are:
    ! myid (integer) proccess id
    ! comm (integer) communicator handle
    !
    ! revisions:
    !     10-17-2012 t.sakai: add setup for fftw.
    use dim
    use mpicom
    use, intrinsic :: iso_c_binding
    include 'fftw3.f03'

    ! argument variables
    integer, intent(in) :: ioption
    ! local variables
    integer :: errorcode,istat

    !<<< if case: ioption = 1 >>>
    if ( ioption == 1 ) then

        ! horizontal mesh must be even (x-only) and divisible by number of 
        ! processors
        if ( mod(nx,2) /= 0 ) then
            if ( myid == 0 ) then
                write(*,*) 'setup_dealiasing: nx must be even.' 
            end if
            call mpi_abort(comm,errorcode,ierr)
        elseif ( mod(nx,nproch) /= 0 .or. mod(ny,nproch) /= 0 ) then
            if ( myid == 0 ) then
                write(*,*) 'setup_dealiasing: (nx,ny) must be divisible by nproch.'
            end if
            call mpi_abort(comm,errorcode,ierr)
        end if 

        ! define the size of refined fourier grid in x and y directions 
        ! in accordance with 3/2 rule (canuto 2006 i pp.134-135).
        if ( nx < 4 .or. ny < 1 ) then
            if ( myid == 0 ) then 
                write(*,*)'setup_dealiasing: must nx >= 4 and ny >= 1.'
            end if
            call mpi_abort(comm,errorcode,ierr)
        end if

        ! define expanded grid size in x-direction
        if ( mod(nx/nproch,2) == 0 .and. mod(nx,4) == 0 ) then
            ! the last condition ensures nx_lg to be even
            nx_lg = 3*nx/2
        elseif ( mod(nx/nproch,3) == 0 ) then
            nx_lg = 5*nx/3
        elseif ( mod(nx/nproch,5) == 0 ) then
            nx_lg = 8*nx/5
        else
            nx_lg = 2*nx
        end if

        ! define expanded grid size in y-direction
        if ( mod(ny/nproch,2) == 0 ) then
            ny_lg = 3*ny/2
        elseif ( mod(ny/nproch,3) == 0 ) then
            ny_lg = 5*ny/3
        elseif ( mod(ny/nproch,5) == 0 ) then
            ny_lg = 8*ny/5
        else
            ny_lg = 2*ny
        end if
             
        nxpp_lg = nx_lg + 2
        nxf_lg    = 3*nx_lg/2 + 1
        nyf_lg    = 2*ny_lg 
        nypl_lg = ny_lg/nproch
        nxpl_lg = nx_lg/nproch

        ! setup fftw plans for extended mesh
        call setup_dealiasing_fftw

        if ( myid == 0 ) then
            write(*,*) 'setup_dealiasing: fft initialization complete.'
        end if

    ! <<< if case: ioption = 0 >>>
    elseif ( ioption == 0 ) then

    ! destroy fftw plans
        call dfftw_destroy_plan(plan_x_lg_fwd)
        call dfftw_destroy_plan(plan_y_lg_fwd)
        call dfftw_destroy_plan(plan_x_lg_backwd)
        call dfftw_destroy_plan(plan_y_lg_backwd)

    ! <<< if case: error >>>
    else

        if ( myid == 0 ) then
            write(*,*)'setup_dealiasing: invalid option.'
        end if
        call mpi_abort(comm,errorcode,ierr)
    end if 
    !<<< end of if case option >>> 

end subroutine setup_dealiasing


subroutine setup_dealiasing_fftw
    !
    ! setup fftw plans for 3d dealiasing.
    !
    use dim
    use, intrinsic :: iso_c_binding
    include 'fftw3.f03'

    ! local variables
    real, dimension(nxpl_lg,ny_lg) :: arr2df
    real, dimension(nxpp_lg,nypl_lg) :: arr2d
    integer :: n_rank,n_howmany,idist,odist,istride,ostride
    integer, dimension(1) :: n_len,inembed,onembed

    ! setup plans for fftw

    ! 1d complex-to-complex backward fft in y-directions.
    n_rank = 1 
    n_len(1) = ny_lg
    n_howmany = nxpl/2
    idist = 1
    odist = 1
    istride = nxpl_lg/2
    ostride = nxpl_lg/2
    inembed(1) = ny_lg
    onembed(1) = ny_lg

    call dfftw_plan_many_dft(plan_y_lg_backwd,n_rank,n_len,n_howmany,arr2df(1,1),inembed,istride,idist,arr2df(1,1),onembed,ostride,odist,fftw_backward,FFTW_ESTIMATE) 

    ! 1d complex-to-real backward fft in x-directions.
    n_rank = 1
    n_len(1) = nx_lg
    n_howmany = nypl_lg
    idist = nx_lg/2+1
    odist = nxpp_lg
    istride = 1
    ostride = 1
    inembed(1) = nx_lg/2+1
    onembed(1) = nxpp_lg

    call dfftw_plan_many_dft_c2r(plan_x_lg_backwd,n_rank,n_len,n_howmany,arr2d(1,1),inembed,istride,idist,arr2d(1,1),onembed,ostride,odist,FFTW_ESTIMATE) 

    ! 1d real-to-complex forward fft in x-directions.
    n_rank = 1
    n_len(1) = nx_lg
    n_howmany = nypl_lg
    idist = nxpp_lg
    odist = nx_lg/2+1
    istride = 1
    ostride = 1
    inembed(1) = nxpp_lg
    onembed(1) = nx_lg/2+1

    call dfftw_plan_many_dft_r2c(plan_x_lg_fwd,n_rank,n_len,n_howmany,arr2d(1,1),inembed,istride,idist,arr2d(1,1),onembed,ostride,odist,FFTW_ESTIMATE) 

    ! 1d complex-to-complex forward fft in y-directions.
    n_rank = 1
    n_len(1) = ny_lg
    n_howmany = nxpl/2 
    idist = 1
    odist = 1
    istride = nxpl_lg/2
    ostride = nxpl_lg/2
    inembed(1) = ny_lg
    onembed(1) = ny_lg 

    call dfftw_plan_many_dft(plan_y_lg_fwd,n_rank,n_len,n_howmany,arr2df(1,1),inembed,istride,idist,arr2df(1,1),onembed,ostride,odist,fftw_forward,FFTW_ESTIMATE) 

end subroutine setup_dealiasing_fftw


subroutine dealiasing_2d_fftw(a,b,ab)

    !
    ! perform dealiasing of a product (ab) of 3d physical field varialbes 
    ! a and b (periodic in x & y and finite in z).
    !
    ! <<< dealiasing is performed in 2 (x,y) directions. >>>
    !
    ! input:
    !        a(nxpp,nypl,nzp), b(nxpp,nypl,nzp) (real)
    !        communicator handle (integer comm) and process id (integer myid)
    !
    ! output:
    !        ab(nxpp,nypl,nzp) (= dealiased a.*b) (real)
    !
    ! before the first call, the initialization routine (setup_dealiasing 
    ! with ioption=1) must be called.
    !


    use dim
    use mpicom
    use modhorfft
    use, intrinsic :: iso_c_binding
    include 'fftw3.f03'

    ! argument variables
    real, dimension(nxpp,nypl,nzp), intent(inout) :: a, b
    real, dimension(nxpp,nypl,nzp), intent(out) :: ab

    ! local parameters
    ! fwd:    forward fft (physical space to fourier coefficient space)
    ! back: backward fft (fourier coefficient space to physical space)
    ! do not change these parameter values!
    integer, parameter :: fwd = -1, back = +1

    ! local storage arrays 
    real, dimension(:,:,:), allocatable :: af, bf
    real, dimension(:,:,:), allocatable :: a_lg, b_lg
    real, dimension(:,:,:), allocatable :: af_lg, bf_lg
        
    real, dimension(:,:), allocatable :: tmp

    ! misc variables
    integer :: i,j,k,n,ks,k1,k2,i1,i2
    integer :: err,errorcode,ifftsign
    real :: scalefactor

    ! allocate local arrays
    allocate(af(nxpl,ny,nzp), bf(nxpl,ny,nzp), a_lg(nxpp_lg,nypl_lg,nzp), b_lg(nxpp_lg,nypl_lg,nzp), af_lg(nxpl_lg,ny_lg,nzp), bf_lg(nxpl_lg,ny_lg,nzp), tmp(nxpl,nypl_lg), stat=err)
    if (err /= 0) print *, "dealiasing: Allocation request denied"
    
    ! initialize array variables.
    ab         = 0.0
    a_lg     = 0.0
    b_lg     = 0.0
    af_lg    = 0.0
    bf_lg    = 0.0


    ! take forward 2d-fft on a and b to get corresponding arrays 
    ! af and bf in fourier space.
    ifftsign = fwd
    call horfft(a,af,ifftsign)
    call horfft(b,bf,ifftsign)

    ! store in extended fourier arrays (af_lg and bf_lg).
    af_lg(1:nxpl,1:ny,1:nzp) = af
    bf_lg(1:nxpl,1:ny,1:nzp) = bf
    
    ! take backward complex fft in y-directions: work on extended arrays.
    do k = 1,nzp
        call dfftw_execute_dft(plan_y_lg_backwd,af_lg(1,1,k),af_lg(1,1,k))
        call dfftw_execute_dft(plan_y_lg_backwd,bf_lg(1,1,k),bf_lg(1,1,k))
    end do

    ! data flip from y to x to get x-oriented data (a_lg and b_lg).
    call dataflip_ytox_dealias(a_lg,af_lg)
    call dataflip_ytox_dealias(b_lg,bf_lg)

    ! at this point the obtained x-oriented data contains gaps (blocks of    
    ! zero subarrays) as a result of using larger arrays.
    ! compress(shift backward) all finite subarray blocks in x-direction 
    ! and pad the rest of array space with zeros. work on a_lg and b_lg.
    if ( nproch > 1 ) then
        do k = 1,nzp
            do n = 2,nproch

                ! variable a_lg:

                ! store subarray in temporary array.
                i1 = (n-1)*nxpl_lg + 1
                i2 = (n-1)*nxpl_lg + nxpl 
                tmp(1:nxpl,1:nypl_lg) = a_lg(i1:i2,1:nypl_lg,k)
                ! paste it back to inboard location.
                i1 = (n-1)*nxpl + 1
                i2 = n*nxpl
                a_lg(i1:i2,1:nypl_lg,k) = tmp(1:nxpl,1:nypl_lg)

                ! variable b_lg:

                ! store subarray in temporary array.
                i1 = (n-1)*nxpl_lg + 1
                i2 = (n-1)*nxpl_lg + nxpl
                tmp(1:nxpl,1:nypl_lg) = b_lg(i1:i2,1:nypl_lg,k)
                ! paste it back to inboard location.
                i1 = (n-1)*nxpl + 1
                i2 = n*nxpl
                b_lg(i1:i2,1:nypl_lg,k) = tmp(1:nxpl,1:nypl_lg)

            end do
        end do
    end if

    ! zero padding the rest of array space.
    i1 = nproch*nxpl + 1
    i2 = nxpp_lg
    a_lg(i1:i2,1:nypl_lg,1:nzp) = 0.0
    b_lg(i1:i2,1:nypl_lg,1:nzp) = 0.0

    ! take backward real fft in x-directions.
    do k = 1,nzp
        call dfftw_execute_dft_c2r(plan_x_lg_backwd,a_lg(1,1,k),a_lg(1,1,k))
        call dfftw_execute_dft_c2r(plan_x_lg_backwd,b_lg(1,1,k),b_lg(1,1,k))
    end do

    ! now a_lg and b_lg are in physical space. compute their product.
    ! overwrite a_lg with the product to save memory.
    a_lg = a_lg*b_lg

    ! in the follwing we dealias the product (a_lg) and project back to 
    ! the original grid.

    ! take forward real fft in x-directions.
    do k = 1,nzp
        call dfftw_execute_dft_r2c(plan_x_lg_fwd,a_lg(1,1,k),a_lg(1,1,k))
    end do

    ! truncate upper 1/3 x-fourier modes and equally divide the remaining 
    ! modes and scatter outboard (shift forward) in x-directions.
    if ( nproch > 1 ) then
        do k = 1,nzp
            do n = nproch,2,-1
                ! store subarray in temporary array.
                i1 = (n-1)*nxpl + 1
                i2 = n*nxpl
                tmp(1:nxpl,1:nypl_lg) = a_lg(i1:i2,1:nypl_lg,k)
                ! paste the temporary array in outboard location.
                i1 = (n-1)*nxpl_lg + 1
                i2 = (n-1)*nxpl_lg + nxpl
                a_lg(i1:i2,1:nypl_lg,k) = tmp(1:nxpl,1:nypl_lg)
            end do
        end do
        ! note: zero padding the higher modes (and array gaps) is 
        ! not needed here because we set the subsequent fft to perform 
        ! only to the lowest 2/3 x-fourier modes. 
    end if

    ! data flip x to y to get y-oriented data (af_lg).
    ! here we reuse the array af_lg to save memory.
    call dataflip_xtoy_dealias(a_lg,af_lg)

    ! take forward complex fft in y-directions.
    scalefactor = 1.0/(real(nx_lg)*real(ny_lg))
    do k = 1,nzp
        call dfftw_execute_dft(plan_y_lg_fwd,af_lg(1,1,k),af_lg(1,1,k))
        ! multiply a scaling factor to obtain nonscaled fourier 
        ! coefficients because fftw does not scale.
        af_lg(:,:,k) = scalefactor*af_lg(:,:,k)
    end do

    ! truncate upper 1/3 y-fourier modes then, store in smaller 
    ! array (af). reuse the array af to save memory.
    af = af_lg(1:nxpl,1:ny,1:nzp)

    ! here we have obtained a dealiased product (af) in fourier space. 
    ! now we take fft the product back to physical space.
    ! the final output is ab.
    ifftsign = back
    call horfft(ab,af,ifftsign)
    ! not sure if it makes sense to take it back to phys space...

    if (allocated(af)) deallocate(af, bf, a_lg, b_lg, af_lg, bf_lg, tmp, stat=err)
    if (err /= 0) print *, "dealiasing: Deallocation request denied"

end subroutine dealiasing_2d_fftw

subroutine dealiasf_2d_fftw(af,bf,abf)

    !
    ! perform dealiasing of a product (ab) of 3d physical field varialbes 
    ! a and b (periodic in x & y and finite in z).
    !
    ! <<< dealiasing is performed in 2 (x,y) directions. >>>
    !
    ! input:
    !        a(nxpp,nypl,nzp), b(nxpp,nypl,nzp) (real)
    !        communicator handle (integer comm) and process id (integer myid)
    !
    ! output:
    !        ab(nxpp,nypl,nzp) (= dealiased a.*b) (real)
    !
    ! before the first call, the initialization routine (setup_dealiasing 
    ! with ioption=1) must be called.
    !
    use dim
    use mpicom
    use modhorfft
    use, intrinsic :: iso_c_binding
    include 'fftw3.f03'

    ! argument variables
    real, dimension(nxpl,ny,nzp), intent(inout) :: af, bf
    real, dimension(nxpl,ny,nzp), intent(out) :: abf

    ! local parameters
    ! fwd:    forward fft (physical space to fourier coefficient space)
    ! back: backward fft (fourier coefficient space to physical space)
    ! do not change these parameter values!
    integer, parameter :: fwd = -1, back = +1

    ! local storage arrays 
!     real, dimension(:,:,:), allocatable :: af, bf
    real, dimension(:,:,:), allocatable :: a_lg, b_lg
    real, dimension(:,:,:), allocatable :: af_lg, bf_lg
        
    real, dimension(:,:), allocatable :: tmp

    ! misc variables
    integer :: i,j,k,n,ks,k1,k2,i1,i2
    integer :: err,errorcode,ifftsign
    real :: scalefactor

    ! allocate local arrays
    allocate( a_lg(nxpp_lg,nypl_lg,nzp), b_lg(nxpp_lg,nypl_lg,nzp), af_lg(nxpl_lg,ny_lg,nzp), bf_lg(nxpl_lg,ny_lg,nzp), tmp(nxpl,nypl_lg), stat=err )
    if (err /= 0) print *, "dealiasing: Allocation request denied"
    
    ! initialize array variables.
    abf      = 0.0
    a_lg     = 0.0
    b_lg     = 0.0
    af_lg    = 0.0
    bf_lg    = 0.0

    ! store in extended fourier arrays (af_lg and bf_lg).
    af_lg(1:nxpl,1:ny,1:nzp) = af
    bf_lg(1:nxpl,1:ny,1:nzp) = bf
    
    ! take backward complex fft in y-directions: work on extended arrays.
    do k = 1,nzp
        call dfftw_execute_dft(plan_y_lg_backwd,af_lg(1,1,k),af_lg(1,1,k))
        call dfftw_execute_dft(plan_y_lg_backwd,bf_lg(1,1,k),bf_lg(1,1,k))
    end do

    ! data flip from y to x to get x-oriented data (a_lg and b_lg).
    call dataflip_ytox_dealias(a_lg,af_lg)
    call dataflip_ytox_dealias(b_lg,bf_lg)

    ! at this point the obtained x-oriented data contains gaps (blocks of    
    ! zero subarrays) as a result of using larger arrays.
    ! compress(shift backward) all finite subarray blocks in x-direction 
    ! and pad the rest of array space with zeros. work on a_lg and b_lg.
    if ( nproch > 1 ) then
        do k = 1,nzp
            do n = 2,nproch

                ! variable a_lg:

                ! store subarray in temporary array.
                i1 = (n-1)*nxpl_lg + 1
                i2 = (n-1)*nxpl_lg + nxpl 
                tmp(1:nxpl,1:nypl_lg) = a_lg(i1:i2,1:nypl_lg,k)
                ! paste it back to inboard location.
                i1 = (n-1)*nxpl + 1
                i2 = n*nxpl
                a_lg(i1:i2,1:nypl_lg,k) = tmp(1:nxpl,1:nypl_lg)

                ! variable b_lg:

                ! store subarray in temporary array.
                i1 = (n-1)*nxpl_lg + 1
                i2 = (n-1)*nxpl_lg + nxpl
                tmp(1:nxpl,1:nypl_lg) = b_lg(i1:i2,1:nypl_lg,k)
                ! paste it back to inboard location.
                i1 = (n-1)*nxpl + 1
                i2 = n*nxpl
                b_lg(i1:i2,1:nypl_lg,k) = tmp(1:nxpl,1:nypl_lg)

            end do
        end do
    end if

    ! zero padding the rest of array space.
    i1 = nproch*nxpl + 1
    i2 = nxpp_lg
    a_lg(i1:i2,1:nypl_lg,1:nzp) = 0.0
    b_lg(i1:i2,1:nypl_lg,1:nzp) = 0.0

    ! take backward real fft in x-directions.
    do k = 1,nzp
        call dfftw_execute_dft_c2r(plan_x_lg_backwd,a_lg(1,1,k),a_lg(1,1,k))
        call dfftw_execute_dft_c2r(plan_x_lg_backwd,b_lg(1,1,k),b_lg(1,1,k))
    end do

    ! now a_lg and b_lg are in physical space. compute their product.
    ! overwrite a_lg with the product to save memory.
    a_lg = a_lg*b_lg

    ! in the follwing we dealias the product (a_lg) and project back to 
    ! the original grid.

    ! take forward real fft in x-directions.
    do k = 1,nzp
        call dfftw_execute_dft_r2c(plan_x_lg_fwd,a_lg(1,1,k),a_lg(1,1,k))
    end do

    ! truncate upper 1/3 x-fourier modes and equally divide the remaining 
    ! modes and scatter outboard (shift forward) in x-directions.
    if ( nproch > 1 ) then
        do k = 1,nzp
            do n = nproch,2,-1
                ! store subarray in temporary array.
                i1 = (n-1)*nxpl + 1
                i2 = n*nxpl
                tmp(1:nxpl,1:nypl_lg) = a_lg(i1:i2,1:nypl_lg,k)
                ! paste the temporary array in outboard location.
                i1 = (n-1)*nxpl_lg + 1
                i2 = (n-1)*nxpl_lg + nxpl
                a_lg(i1:i2,1:nypl_lg,k) = tmp(1:nxpl,1:nypl_lg)
            end do
        end do
        ! note: zero padding the higher modes (and array gaps) is 
        ! not needed here because we set the subsequent fft to perform 
        ! only to the lowest 2/3 x-fourier modes. 
    end if

    ! data flip x to y to get y-oriented data (af_lg).
    ! here we reuse the array af_lg to save memory.
    call dataflip_xtoy_dealias(a_lg,af_lg)

    ! take forward complex fft in y-directions.
    scalefactor = 1.0/(real(nx_lg)*real(ny_lg))
    do k = 1,nzp
        call dfftw_execute_dft(plan_y_lg_fwd,af_lg(1,1,k),af_lg(1,1,k))
        ! multiply a scaling factor to obtain nonscaled fourier 
        ! coefficients because fftw does not scale.
        af_lg(:,:,k) = scalefactor*af_lg(:,:,k)
    end do

    ! truncate upper 1/3 y-fourier modes then, store in smaller 
    ! array (af). reuse the array af to save memory.
    abf = af_lg(1:nxpl,1:ny,1:nzp)

    if (allocated(a_lg)) deallocate( a_lg, b_lg, af_lg, bf_lg, tmp, stat=err )
    if (err /= 0) print *, "dealiasing: Deallocation request denied"

end subroutine dealiasf_2d_fftw

    
subroutine dataflip_xtoy_dealias(x,xf)

    !
    ! subroutine that reorders data from domain decomposition in
    ! x-direction to d.d. in y-direction
    !
    ! this is a modified version of dataflip_xtoy - increased array size 
    ! to use for dealiasing in which larger arrays are needed.
    !
    ! input: 
    ! x(nxpp_lg,nypl_lg,nzp): nzp 2-d slices of data partitioned normal to y-direction
    ! output:
    ! xf(nxpl_lg,ny_lg,nzp) : nzp 2-d slices of data partitioned normal to x-direction
    !
                 
    use dim
    use mpicom

    ! argument variables
    real, dimension(nxpp_lg,nypl_lg,nzp), intent(in) :: x
    real, dimension(nxpl_lg,ny_lg,nzp), intent(out)    :: xf

    ! local variables
    real, dimension(:,:,:,:), allocatable :: xtempin,xtempout
    integer :: errorcode, err
    integer :: i,j,k, size_of_block, iglob, jglob

    allocate(xtempin(nxpl_lg,nzp,nypl_lg,nproch), xtempout(nxpl_lg,nzp,nypl_lg,nproch), stat=err)
    if (err /= 0) print *, "dataflip_xtoy_dealias: Allocation request denied"
    
    ! pack input array into send buffer. 
    do hproc=1,nproch
        do k=1,nzp
            do j=1,nypl_lg
                do i=1,nxpl_lg
                    iglob = (hproc-1)*nxpl_lg + i
                    xtempin(i,k,j,hproc) = x(iglob,j,k)                
                enddo
            enddo
        enddo
    enddo

    ! size of blocks to be transmitted to other processors
    size_of_block = nxpl_lg*nypl_lg*nzp

    ! syncronize processors right before the global collective 
    ! communication. this is to ensure the portability of the code.
!     call mpi_barrier(comm,ierr)

    ! perform data transposition
    call mpi_alltoall(xtempin,size_of_block,mpi_double_precision,xtempout,size_of_block,mpi_double_precision,comm,ierr)

    if (ierr /= 0) then 
        if (myid == 0) write(*,*) 'mpi_alltoall error in dataflip_xtoy_dealias'
        call mpi_abort(comm,errorcode,ierr)
    endif

    ! now copy the receive buffer from xtempout to the output array
    do hproc=1,nproch
        do k=1,nzp
            do j=1,nypl_lg
                do i=1,nxpl_lg
                    jglob = (hproc-1)*nypl_lg + j
                    xf(i,jglob,k) = xtempout(i,k,j,hproc)
                enddo
            enddo
        enddo
    enddo
     
    ! add a barrier to be on the safe side
!     call mpi_barrier(comm,ierr)

    if (allocated(xtempin)) deallocate(xtempin, xtempout, stat=err)
    if (err /= 0) print *, "dataflip_xtoy_dealias: Deallocation request denied"
    
end subroutine dataflip_xtoy_dealias

    
subroutine dataflip_ytox_dealias(x,xf)

    !
    ! subroutine that reorders data from domain decomposition in
    ! y-direction to d.d. in x-direction
    !
    ! this is a modified version of dataflip_ytox - increase array size 
    ! to use for dealiasing in which larger arrays are needed.
    !
    ! input:
    ! xf(nxpl_lg,ny_lg,nzp) : nzp 2-d slices of data partitioned
    !                                                         normal to x-direction
    !
    ! output: 
    ! x(nxpp_lg,nypl_lg,nzp): nzp 2-d slices of data partitioned 
    !                                                         normal to y-direction
    !

    use dim
    use mpicom

    ! argument variables
    real, dimension(nxpp_lg,nypl_lg,nzp), intent(out) :: x
    real, dimension(nxpl_lg,ny_lg,nzp), intent(in)    :: xf

    ! local variables
    real, dimension(:,:,:,:), allocatable :: xtempin,xtempout
    integer :: errorcode, err
    integer :: i,j,k, size_of_block, iglob, jglob

    allocate(xtempin(nypl_lg,nzp,nxpl_lg,nproch), xtempout(nypl_lg,nzp,nxpl_lg,nproch), stat=err)
    if (err /= 0) print *, "dataflip_ytox_dealias: Allocation request denied"

    ! initialize output array x because it contains extra padding in x
    x=0.0

    ! pack input array into send buffer.
    do hproc=1,nproch         
        do k=1,nzp
            do j=1,nypl_lg
                do i=1,nxpl_lg
                     jglob = (hproc-1)*nypl_lg + j
                     xtempin(j,k,i,hproc) = xf(i,jglob,k)
                 enddo
             enddo
         enddo
    enddo

    ! size of blocks to be transmitted to other processors
    size_of_block = nxpl_lg*nypl_lg*nzp

    ! syncronize processors right before the global collective 
    ! communication. this is to ensure the portability of the code.
!     call mpi_barrier(comm,ierr)
        
    ! perform data transposition
    call mpi_alltoall(xtempin,size_of_block,mpi_double_precision,xtempout,size_of_block,mpi_double_precision,comm,ierr)

    if (ierr /= 0) then
        if (myid == 0) write(*,*) 'mpi_alltoall error in dataflip_ytox_dealias'
        call mpi_abort(comm,errorcode,ierr)
    endif

    ! now copy the receive buffer from xtempout to the output array
    do hproc=1,nproch
        do k=1,nzp
            do j=1,nypl_lg
                do i=1,nxpl_lg
                    iglob = (hproc-1)*nxpl_lg + i
                    x(iglob,j,k) = xtempout(j,k,i,hproc)
                enddo
            enddo
        enddo
    enddo

    ! add a barrier to be on the safe side
!     call mpi_barrier(comm,ierr)

    if (allocated(xtempin)) deallocate(xtempin, xtempout, stat=err)
    if (err /= 0) print *, "dataflip_ytox_dealias: Deallocation request denied"

end subroutine dataflip_ytox_dealias



end module dealias