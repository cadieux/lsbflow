!**********************************************************************
! 
! SETS UP FFT PLANS TO DRIVE FFTW ROUTINES. PLANS DEFINED IN SETUP_FFTW
! ARE USED IN HORIZONTAL 2D FFTS (horfft).
!
! USAGE:
!
!     1. TO INITIALIZE PLANS USE OPTION 1:
!
!    CALL SETUP_FFTW(1)
!
!    THIS INITIALIZATION IS DONE ONLY ONCE AND THE CREATED PLANS
!    ARE USED THROUGHOUT THE MAIN PROGRAM.
!
!     2.0 TO DESTROY PLANS USE OPTION 0:
!
!    CALL SETUP_FFTW(0)
!
!    THIS IS USUALLY DONE AT THE END OF MAIN PROGRAM.
!
! REQUIREMENTS:
!
!     1. FFTW 3.3.2
!
!     2.0 ARRAY DIMENSIONS
!
! REVISION HISTORY:
!
!     ORIGINAL: T.SAKAI 10-16-2012 @ CORNELL UNIVERSITY
!**********************************************************************

module modfftw


! THIS MODULE DEFINES FLAGS AND FFTW PLAN VARIABLES USED IN HORIZONTAL
! FFTS.
!
! DESCRIPTION OF VARIABLES:
!
! PLAN_X_FWD     : PLAN FOR 1D real-TO-COMPLEX FORWARD FFT IN X-DIR.
! PLAN_Y_FWD     : PLAN FOR 1D COMPLEX-TO-COMPLEX FORWARD FFT IN Y-DIR.
! PLAN_X_BACKWD: PLAN FOR 1D COMPLEX-TO-real FORWARD FFT IN X-DIR.
! PLAN_Y_BACKWD: PLAN FOR 1D COMPLEX-TO-COMPLEX BACKWARD FFT IN Y-DIR.
use, intrinsic :: iso_c_binding
implicit none
public
save
include 'fftw3.f03'
! include 'fftw3.f'

! integer, parameter, private :: KD=8
integer(kind=8) :: PLAN_X_FWD
integer(kind=8) :: PLAN_Y_FWD
integer(kind=8) :: PLAN_X_BACKWD
integer(kind=8) :: PLAN_Y_BACKWD

contains 

subroutine setup_fftw(ioption)
    ! THIS SUBROUTINE SETS UP FFTW PLANS FOR 2D (HORIZONTAL) FFTS. SET
    ! ioption = 1 (integer) FOR INITIALIZATION OR ioption = 0 FOR
    ! FINALIZATION(FREE ALL PLANS).
    use dim, only: nx, nxpp, ny, nypl, nxpl, nz, nzpl

    ! Argument variables
    integer, intent(IN) :: ioption

    ! Local variables
    real, dimension(nxpp,nypl)    :: arr2d
    real, dimension(nxpl,ny)    :: arr2df
    integer :: n_rank,n_howmany,idist,odist,istride,ostride
    integer, dimension(1) :: n_len,inembed,onembed

    !     integer :: nthreads, tid, iret

    if ( ioption == 1 ) then
        ! Add threading support
    !         call dfftw_init_threads(iret)
    !         print *, iret
    !         call dfftw_plan_with_nthreads(4)

        ! Plan 1D real-to-Complex Forward FFT in x-direction
        n_rank = 1
        n_len(1) = nx
        n_howmany = nypl
        idist = nxpp
        odist = nx/2+1
        istride = 1
        ostride = 1
        inembed(1) = nxpp
        onembed(1) = nx/2+1

        call dfftw_plan_many_dft_r2c(PLAN_X_FWD,n_rank,n_len,n_howmany,arr2d(1,1),inembed,istride,idist,arr2d(1,1),onembed,ostride,odist,FFTW_MEASURE)

        ! Plan 1D Complex-to-Complex Forward FFT in y-direction
        n_rank = 1
        n_len(1) = ny
        n_howmany = nxpl/2
        idist = 1
        odist = 1
        istride = nxpl/2
        ostride = nxpl/2
        inembed(1) = ny
        onembed(1) = ny

        call dfftw_plan_many_dft(PLAN_Y_FWD,n_rank,n_len,n_howmany,arr2df(1,1),inembed,istride,idist,arr2df(1,1),onembed,ostride,odist,FFTW_FORWARD,FFTW_MEASURE)

        ! Plan 1D Complex-to-Complex Backward FFT in y-direction
        n_rank = 1
        n_len(1) = ny
        n_howmany = nxpl/2
        idist = 1
        odist = 1
        istride = nxpl/2
        ostride = nxpl/2
        inembed(1) = ny
        onembed(1) = ny

        call dfftw_plan_many_dft(PLAN_Y_BACKWD,n_rank,n_len,n_howmany,arr2df(1,1),inembed,istride,idist,arr2df(1,1),onembed,ostride,odist,FFTW_BACKWARD,FFTW_MEASURE)

        ! Plan 1D Complex-to-real Backward FFT in x-direction
        n_rank = 1
        n_len(1) = nx
        n_howmany = nypl
        idist = nx/2+1
        odist = nxpp
        istride = 1
        ostride = 1
        inembed(1) = nx/2+1
        onembed(1) = nxpp

        call dfftw_plan_many_dft_c2r(PLAN_X_BACKWD,n_rank,n_len,n_howmany,arr2d(1,1),inembed,istride,idist,arr2d(1,1),onembed,ostride,odist,FFTW_MEASURE)


    elseif ( ioption == 0 ) then

        ! Destroy plans
        call dfftw_destroy_plan(PLAN_X_FWD)
        call dfftw_destroy_plan(PLAN_Y_FWD)
        call dfftw_destroy_plan(PLAN_X_BACKWD)
        call dfftw_destroy_plan(PLAN_Y_BACKWD)

    else

        write(*,*)'setup_fftw: invalid option.'

    end if

end subroutine setup_fftw

end module modfftw