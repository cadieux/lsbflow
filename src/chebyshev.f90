module chebyshev

implicit none
public
save

contains

subroutine build_cheb_zder_matrix(dor,d,d2)
    use dim, only: nzp
    use grid, only: gf, gf2
    use flags, only: debug
    use mpicom, only: myid
    use formatstrings
    ! i/o
    real, dimension(nzp,nzp), intent(out) :: dor, d, d2
    ! vars
    integer :: i,j
    ! allocate derivative matrix
!     allocate(dor(nzp,nzp),d(nzp,nzp),d2(nzp,nzp))
    ! compute first derivative matrix
    do j=1,nzp
        do i=1,nzp
            dor(i,j) = del(i,j,nzp)
            d(i,j) = del(i,j,nzp)*gf(i)
        end do
    end do
    ! compute second derivative matrix
    d2 = matmul(d,d)
!     d2 = matmul(dor,dor)
!     do i=1,nzp
!         d2(:,i) = d2(:,i)*gf2(:)
!     end do
!     do i=1,nzp
!         do j=1,nzp
!             d2(i,j) = 0.0
!             do k=1,nzp
!                 d2(i,j) = d2(i,j) + d(i,k)*d(k,j)
!             end do
!         end do
!     end do
    ! debug: print out derivative matrix for visual inspection
    if (debug>=2 .and. myid==0) then
        write(*,*) 'dor(1:nzp,1:nzp)'
        write(*,fm3) dor(1:nzp,1:nzp)
    end if

end subroutine build_cheb_zder_matrix


real function del(k,j,n)
    ! returns Chebyshev a single derivative coefficient 
    ! of d/d(xi) matrix given the desired index (k,j)
    ! where xi = cos(pi*(k-1)/nz) = [-1,1]
    use grid, only: xl
    integer, intent(in) :: k,j,n
    integer :: kr,jr,nr
    real :: fac
    kr = k - 1
    jr = j - 1
    nr = n - 1
    fac = 1.
    if( kr .eq. nr ) then
        ! take care of boundary condition k=nz (bottom row), 
        ! make it the same as first row by setting kr = 0
        kr = 0
        ! but reverse order of bottom row, 
        ! and change its sign (compared to top row)
        jr = nr - jr
        fac = -1.
    end if
    del = fac*xnum(kr,jr,nr)/den(kr,jr,nr)
    ! scale the d/d(xi) entry to reflect 
    ! xi = xl/2*cos(pi*(k-1)/nz) = [-xl/2,xl/2]
    del = del*(2.0/xl)

end function del


real function den(kr,jr,nr)
    ! returns denominator of single Chebyshev derivative coefficient
    ! of d/d(xi) matrix given the desired index (k-1,j-1)
    integer, intent(in) :: kr,jr,nr
    real :: pi
    pi = 4.*atan(1.)
    ! den = nz*sin(pi(k-1)/nz) except first and last row
    den = nr*sin(pi*kr/nr)
    ! and at first or last column, den is double
    if(jr .eq. 0 .or. jr .eq. nr) den = den*2.0

    ! this section targets first and last row of matrix
    if (kr==0) then
        ! for first and last row, den = nz/2 for even (j-1)
        den = .5*nr
        ! for odd (j-1), den = -nz/2
        if(mod(jr,2) .eq. 1) den = - den
        ! except at four corners of matrix, where den = 1
        if(jr .eq. 0 .or. jr .eq. nr) den = 1.
    end if

end function den


real function xnum(kr,jr,nr)
    ! returns numerator of single Chebyshev derivative coefficient
    ! of d/d(xi) matrix given the desired index (k-1,j-1)
    integer, intent(in) :: kr,jr,nr
    real :: pi
    pi = 4.*atan(1.)
    
    ! nominally, except at first & last row, 
    ! xnum = +/- (nz/2)/tan(pi*(k+j-2)/(2*nz)) +/- (nz/2)/tan(pi*(k-j)/(2*nz))
    xnum = ff(kr+jr,nr) + ff(kr-jr,nr)
    
    if (kr==0) then
        ! for interior points of first and last row
        xnum = .5*(nr+.5)*((nr+.5) &
            & + (1./tan(pi*jr/(2.0*nr)))**2) &
            & + 1./8. - .25/(sin(pi*jr/(2.0*nr))**2.0) &
            & - .5*nr*nr
        
        ! for the 2 right corners (jr=nz,kr=0) and (jr=nz,kr=nz)
        ! xnum = 0.5
        if(jr .eq. nr) then
            xnum = .5
            ! if nz is not divisible by 2, xnum = -0.5 for same 2 corners
            if( mod(nr,2) .eq. 1 ) xnum = - xnum
        end if
        ! for 2 left corners (kr=0,jr=0) and (kr=nz,jr=0),
        ! xnum = (nz^2)/3 + 1/6
        if(jr.eq.0) xnum = (1./3.)*nr*nr + 1./6.
            
    end if

end function xnum


real function ff(i,nr)
    ! computes half of the numerator for 
    ! a single Chebyshev derivative coefficient
    ! of d/d(xi) matrix given the desired index i = (k+j-2) or i = (k-j)
    integer, intent(in) :: i, nr
    real :: pi
    pi = 4.*atan(1.)
    ! ff = nz/2 / tan(pi*(k-j)/(2*nz)) OR ff = nz/2 / tan(pi*(k+j-2)/(2*nz))
    ff = nr*.5/tan(pi*i/(2.0*nr))
    
    ! for even (k+j-2) and (k-j), ff is negative
    if (mod(i,2) .eq. 0) ff = -ff

    ! zero diagonal component (k=j) from i=(k-j),
    ! and from (k+j-2) when k=1,j=1
    ! because tan(0) = 0 and 1/tan(0) would go to infinity otherwise
    ! hence why boundaries are handled in xnum
    if(i .eq. 0) then
        ff = 0.0
    end if 

end function ff


! subroutine dct_diff(dydx, y, half_interval, n)
!     include 'fftw3.f'
!     ! use the discrete cosine transform from fftw to calculate the derivative
!     integer, intent(in) :: n ! number of data points 
!     real, intent(in), dimension(n) :: y ! function values 
!     real, intent(out), dimension(n) :: dydx ! derivative values 
!     real, intent(in) :: half_interval ! half the interval length 
!     real, dimension(n) :: beta ! derivative coefficients
!     real, dimension(n) :: alpha ! function coefficients 
!     integer(kind=8) :: plan1, plan2
!     integer :: i, n_logical ! the logical size of the transform, size of
!     ! the corresponding real symmetric DFT

!     ! the logical size depends on the type of transform, check the docs:
!     ! http://www.fftw.org/doc/1d-real_002deven-DFTs-_0028DCTs_0029.html 
!     n_logical = 2*(n-1) 

!     ! forward DCT: 
!     call dfftw_plan_r2r_1d(plan1, n, y, alpha, FFTW_REDFT00, FFTW_ESTIMATE)
!     call dfftw_execute_r2r(plan1, y, alpha)
!     call dfftw_destroy_plan(plan1)
!     alpha = alpha / half_interval 

!     ! recurrence for the derivative coefficients: 
!     beta(n) = 0
!     beta(n-1) = 2 * real(n-1) * alpha(n)
!     do i = n-1, 2, -1
!     beta(i-1) = beta(i+1) + 2 * real(i-1) * alpha(i)
!     end do
! !     beta = - beta ! this makes it work, but not sure why!

!     ! inverse DCT: 
!     call dfftw_plan_r2r_1d(plan2, n, beta, dydx, FFTW_REDFT00, FFTW_ESTIMATE)
!     call dfftw_execute_r2r(plan2, beta, dydx)
!     call dfftw_destroy_plan(plan2)

!     ! FFTW computes the un-normalized transforms, normalize by logical size
!     dydx = dydx / real(n_logical)

! end subroutine dct_diff


subroutine dct_diff(n1, n2, n, y, dydx, p)
    ! use the discrete cosine transform from fftw to calculate the derivative
    ! taken from: variousconsequences.com/2009/05/fftw-discrete-cosine-transform.html
!     use dim
    use grid, only: xl, gf
!     include 'fftw3.f'
    use, intrinsic :: iso_c_binding
    include 'fftw3.f03'
    integer, intent(in) :: n1, n2, n, p ! p = 1 : first derivative; p = 2 : second derivative
    real, intent(in), dimension(n1,n2,n) :: y ! function values 
    real, intent(out), dimension(n1,n2,n) :: dydx ! derivative values 
    real :: half_interval ! half the interval length 
    real, allocatable, dimension(:) :: alpha, beta ! work arrays
!     real, dimension(n) :: beta ! derivative coefficients
!     real, dimension(n) :: alpha ! function coefficients 
    integer(kind=8) :: plan1, plan2
    integer :: i,j,k, l, n_logical, err ! the logical size of the transform, size of
    ! the corresponding real symmetric DFT

    half_interval = 0.5*xl

    ! the logical size depends on the type of transform, check the docs:
    ! http://www.fftw.org/doc/1d-real_002deven-DFTs-_0028DCTs_0029.html 
    n_logical = 2*(n-1)

    if (.not. allocated(alpha)) allocate(alpha(n), beta(n), stat=err)
    if (err /= 0) print *, "array: Allocation request denied"

    call dfftw_plan_r2r_1d(plan1, n, beta, alpha, FFTW_REDFT00, FFTW_ESTIMATE)
    call dfftw_plan_r2r_1d(plan2, n, beta, alpha, FFTW_REDFT00, FFTW_ESTIMATE)

    do j=1,n2
        do i=1,n1
            beta = y(i,j,:) ! copy vertical data into beta
            do l=1,p ! loop to calculate higher derivatives
                ! forward DCT: 
                call dfftw_execute_r2r(plan1, beta, alpha)
                ! divide by half interval (because domain is not from -1 to 1)
                alpha(2:n-1) = alpha(2:n-1) / half_interval

                beta = 0.0! re-zero beta
                ! recurrence for the derivative coefficients: 
                beta(n) = 0.0
                beta(n-1) = 2.0 * (n - 1.) * alpha(n) ! 2.0*n*alpha(n)
                do k = n-1, 2, -1
                    beta(k-1) = beta(k+1) + 2.0 * (k - 1.) * alpha(k) !2.0*k*alpha(k) 
                end do
                
                alpha = 0.0! re-zero alpha
                ! inverse DCT: 
                call dfftw_execute_r2r(plan2, beta, alpha)
                ! FFTW computes the un-normalized transforms
                ! normalize by logical size
                ! and apply mapping (chain rule)
                beta = alpha*gf / n_logical
            end do
            dydx(i,j,:) = beta
        end do
    end do

    call dfftw_destroy_plan(plan1)
    call dfftw_destroy_plan(plan2)

    if (allocated(alpha)) deallocate(alpha, beta, stat=err)
    if (err /= 0) print *, "array: Deallocation request denied"

end subroutine dct_diff


! subroutine cheb_fft_diff(uf, dudzf)
!     ! use the discrete cosine transform from fftw to calculate the derivative
!     ! taken from: variousconsequences.com/2009/05/fftw-discrete-cosine-transform.html
!     use dim
!     use grid, only: xl, gf
! !     include 'fftw3.f'
!     use, intrinsic :: iso_c_binding
!     include 'fftw3.f03'
! !     integer, intent(in) :: n ! number of data points 
!     real, intent(in), dimension(nxpl,ny,nl) :: uf ! function values 
!     real, intent(out), dimension(nxpl,ny,nzpl) :: dudzf ! derivative values 
!     real :: half_interval, pi ! half the interval length
! !     real, allocatable, dimension(:) :: u, uhat, what, w, t
!     complex*16, allocatable, dimension(:) :: u, uhat, what, w, t
!     integer :: plan1, plan2
!     integer :: i,j,k, n, err, wavz

!     half_interval = 0.5*xl

!     pi = 4.*atan(1.)

!     n = 2*(nzp-1)

!     allocate(u(n), uhat(n), what(n), w(n), t(nzp), stat=err)
!     if (err /= 0) print *, "cheb_fft_diff: Allocation request denied"
    
!     call dfftw_plan_dft_1d(plan1, n, u, uhat, FFTW_FORWARD, FFTW_ESTIMATE)
!     call dfftw_plan_dft_1d(plan2, n, what, w, FFTW_BACKWARD, FFTW_ESTIMATE)

!     do j=1,ny
!         do i=1,nxpl !,2
!             ! extend data to 2*(nzp-1): u_0, u_1,...,u_n, u_n-1,...,u_1 
!             u(1:nzp) = uf(i,j,:) ! + (0.0,1.)*uf(i+1,j,:)
!             do k=2,nz
!                 u(nzp+k-1) = uf(i,j,nzp-k+1) ! + (0.0,1.)*uf(i+1,j,nzp-k+1)
!             end do

!             ! forward fourier transform: 
!             call dfftw_execute(plan1, u, uhat)
!             ! normalize by number of actual grid points
! !             uhat = uhat / nzp
!             uhat = real( uhat )

!             ! calculate derivative: what(kz) = i kz uhat, where kz = -nzp+1,...,nzp
!             do k = 1,n
!                 wavz = -nzp+k
!                 what(k) = -(0.0,1.)*wavz*uhat(k)
!             end do
!             what(nzp) = 0.0! for first derivative

!             ! inverse fourier transform: 
!             call dfftw_execute(plan2, what, w)

!             w = w / real(n) !nzp

!             ! Chebyshev derivatives
!             t = 0.0
!             do k = 1,nz
! !                 t(k) = -w(k)/sqrt( 1. - cos((nz-k+1)*pi/real(nz)) )
!                 t(k) = - real( w(k) )/sqrt( 1. - cos((nz-k+1)*pi/real(nz)) )
!             end do
!             ! derivatives for end points
!             t(1) = 0.5*nzp*uhat(nzp+1)
!             t(nzp) = 0.5*nzp*uhat(nzp+1)
!             do k=1,nzp
!                 t(1) = t(1) + real(k-1)**2*uhat(k)/real(nzp)
!                 t(nzp) = t(nzp) + real((-1)**(k-1)*(k-1))**2*uhat(k)/real(nzp)
!             end do

!             ! apply mapping (chain rule)
!             dudzf(i,j,:) = t*gf
! !             dudzf(i,j,:) = real(t)*gf
! !             dudzf(i+1,j,:) = aimag(t)*gf
!         end do
!     end do

!     call dfftw_destroy_plan(plan1)
!     call dfftw_destroy_plan(plan2)

!     if (allocated(uhat)) deallocate(u, uhat, what, w, t, stat=err)
!     if (err /= 0) print *, "cheb_fft_diff: Deallocation request denied"

! end subroutine cheb_fft_diff



end module chebyshev