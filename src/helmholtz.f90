module helmholtz
! contains all variables and subroutines necessary to solve
! Helmholtz equation of the form: (D^2 - k^2 - dtnui)u = rhs
! with Dirichlet boundary conditions at zpts(1) and zpts(nzp)
! setting dtfac to zero solves the corresponding Poisson eqn

implicit none
public
save

contains

subroutine diagonalize_matrix(D2,n,lamda,e,ei)
    ! find eigenvalues and eigenvectors of d^2/dz^2 matrix D2
    ! to diagonalize it s.t. D2 = E lamda E^-1, where diag is diagonal,
    ! and E are the matrices whose columns are eigenvectors
    use flags, only: debug
    !     use dim, only: nzp,nz,nzm
    !     use mats
    use mpicom, only: myid
    use formatstrings
    integer, intent(in) :: n
    real, dimension(n,n) :: D2
    real, dimension(n), intent(out) :: lamda
    real, dimension(n,n), intent(out) :: e, ei
    ! local vars req'd for LAPACK
    character*1 :: jobvl, jobvr
    integer :: k, i, info, lwork, lda, ldb, ldvl, ldvr
    integer, allocatable, dimension(:) :: ipiv
    real, allocatable, dimension(:,:) :: A, B, vl, vr
    real, allocatable, dimension(:) :: work, alphar, alphai, beta
    !     real :: tmp

    ! setup input for LAPACK SGGEV which calculates eigenvalues, 
    ! and eigenvectors for problem A v = lamda*B v
    jobvl = 'N'
    jobvr = 'V'
    !     n = nzm
    lda = n
    ldb = n
    ldvl = n
    ldvr = n
    lwork = 10*n
    allocate( A(lda,n), B(ldb,n), vl(ldvl,n), vr(ldvr,n), work(lwork), alphar(n), alphai(n), beta(n), ipiv(n) )

    ! A = D2(2:nz,2:nz) because we always solve helmholtz and poisson eqn
    ! with homogeneous boundary conditions e.g.:
    ! [D^2 - k^2 - 2/(nu*dt)] phi = RHS + hbot(z)bc1(x,y) + h2(z)bc2(x,y)
    ! where phi(z=0,H or -1,1) = 0
    !     A = D2(ns:n+1,ns:n+1)
    A = D2

    ! B is the identity matrix in this case because we want the solution to
    ! the classical eigenvalue problem A v = lamda*v
    B = 0.0
    do k=1,n
        B(k,k) = 1.
    end do

    ! DGGEV computes for a pair of N-by-N real nonsymmetric matrices (A,B)
    !     the generalized eigenvalues, and optionally, the left and/or right
    !     generalized eigenvectors.
    !     A generalized eigenvalue for a pair of matrices (A,B) is a scalar
    !     lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
    !     singular. It is usually represented as the pair (alpha,beta), as
    !     there is a reasonable interpretation for beta=0, and even for both
    !     being zero.
    !     The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
    !     of (A,B) satisfies
    !              A * v(j) = lambda(j) * B * v(j).
    !     The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
    !     of (A,B) satisfies
    !              u(j)**H * A  = lambda(j) * u(j)**H * B .
    !     where u(j)**H is the conjugate-transpose of u(j).
    call dggev( jobvl, jobvr, n, A, lda, B, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork, info )
    ! check for errors
    if (info>=2 .and. myid==0) write(*,'((//,A),I11.4)') 'info =', info

    if (debug>=2) then
        ! scale each eigen vector by its first value for comparison with previous code results
        do k=1,n
            vr(:,k) = vr(:,k)/vr(1,k)
        end do
    end if

    ! store calculate eigenvalues in lamda, and eigenvectors in e
    lamda = 0.0
    lamda(1:n) = alphar/beta ! assuming only real eigenvalues
    e(1:n,1:n) = vr

    ! switch last two eigenvalues if they are in wrong order
    !     if (lamda(n)<lamda(n-1)) then
    !         tmp = lamda(n)
    !         lamda(n) = lamda(n-1)
    !         lamda(n-1) = tmp
    !         vl(:,n) = e(:,n)
    !         e(:,n) = e(:,n-1)
    !         e(:,n-1) = vl(:,n)
    !     end if

    ! now find E^-1 = vr^-1
    ! first get LU decomposition
    call dgetrf( n, n, vr, lda, ipiv, info )
    ! then solve for inverse
    call dgetri( n, vr, lda, ipiv, work, lwork, info )
    ! check for errors
    if (debug>=2 .and. myid==0) write(*,'((//,A),I11.4)') 'info =', info

    ! store E^-1 into ei
    ei(1:n,1:n) = vr

    if (debug>=1 .and. myid==0) then
        write(*,'((//,A),I11.4)') 'info =', info
        write(*,'(//,A)') 'eigenvalues ='
        write(*,'(2(E14.7,'',''))') (lamda(i),alphai(i),i=1,n)
    !         write(*,'(//,A)') 'right eigenvectors ='
    !         write(*,fm3) ((e(i,k),k=1,n),i=1,n)
    !         write(*,'(//,A)') 'inv(E) ='
    !         write(*,fm3) ((ei(i,k),k=1,n),i=1,n)
    end if


    ! deallocate all local vars
    deallocate( A, B, vl, vr, work, alphar, alphai, beta, ipiv )

end subroutine diagonalize_matrix



subroutine solve_helmholtz(rhs,dtnui,n,diag,smat,simat)
    ! solves the Helmholtz or Poisson eqn with homogeneous boundary conditions
    use dim
    use grid
    !     use mats
    ! i/o
    integer, intent(in) :: n
    real, dimension(nxpl,ny,nzpl), intent(inout) :: rhs
    real, intent(in) :: dtnui
    real, dimension(n), intent(in) :: diag
    real, dimension(n,n), intent(in) :: smat,simat
    ! vars
    integer :: i, j, k, ks, ke
    real, allocatable, dimension(:) :: denom
    real, allocatable, dimension(:,:) :: dfac
    real, allocatable, dimension(:,:,:) :: q, trhs
    real :: mindiag

    mindiag = 1.E-4
    ks = 2; ke = nz;
    if (n==nz) ks = 2; ke = n+1;

    ! allocate local vars
    allocate( dfac(nxpl,ny), q(nxpl,ny,n), trhs(nxpl,ny,n), denom(n) )
    ! copy rhs into q
    q = rhs(:,:,ks:ke)
    ! denom
    i = 1
    do k=ks,ke
        denom(i) = diag(k-1)
        i = i +1
    end do
    ! solve (D^2 - k^2 - dtnui)u = rhs
    ! using previously decomposed matrix D^2 = E diag E^-1 s.t.
    ! (D^2 - k^2 - dtnui)q = (E diag E^-1 - I (k^2 + dtnui))q = rhs
    ! left multiply eqn by E^-1, right multiply by E and obtain
    ! (diag - I (k^2 + dtnui))q = E^-1 rhs E
    ! q(z) = E^-1 rhs(z) E / (diag(z) - I (k^2 + dtnui)) for each i,j
    do j=1,ny
        do i=1,nxpl
            ! dfac = -k^2 -dtnui
            dfac(i,j) = -wavxx(i) - wavyy(j) - dtnui
            ! left multiply eqn by E^-1,
!              trhs(i,j,:) = matmul(simat,q(i,j,:))
          call DGEMV ( 'N', n, n, 1.0, simat, n, q(i,j,:), 1, 0.0, trhs(i,j,:), 1 )
!                 trhs(i,j,ks:ke) = matmul(rhs(i,j,ks:ke),simat)
            
            ! q(z) = E^-1 rhs(z) E / (diag(z) - I (k^2 + dtnui))
!             trhs(i,j,ks:ke) = trhs(i,j,ks:ke)/(diag(1:nzm) + dfac(i,j))
    !             trhs(i+1,j,ks:ke) = trhs(i+1,j,ks:ke)/(diag(1:nzm) + dfac(ii,j))
            trhs(i,j,:) = trhs(i,j,:)/(denom(:) + dfac(i,j))
!             do k=ks,ke
!                 if (abs(diag(k-1)+ dfac(i,j)) >= mindiag ) then
!                     trhs(i,j,k) = trhs(i,j,k)/(diag(k-1) + dfac(i,j))
!                 else
!                     trhs(i,j,k) = 0.0
!                 end if
!             end do
            ! right multiply by E
!              q(i,j,:) = matmul(smat,trhs(i,j,:))
           call DGEMV ( 'N', n, n, 1.0, smat, n, trhs(i,j,:), 1, 0.0, q(i,j,:), 1 )
    !             q(i,j,ks:ke) = matmul(trhs(i,j,ks:ke),smat)
    !             q(i,j,ks:ke) = q(i,j,ks:ke)/( diag(1:nzm) + dfac(ii,j) )
            ! copy soln into imaginary part of mode since dfac doesn't change
    !             q(i+1,j,ks:ke) = q(i,j,ks:ke)
        end do
    end do

    rhs(:,:,ks:ke) = q

    ! apply homogeneous boundary conditions
!     if (n==nzm) then
        rhs(:,:,1) = 0.0
        rhs(:,:,nzp) = 0.0
!     end if
    ! deallocate local vars
    deallocate(dfac, q, trhs, denom)

end subroutine solve_helmholtz


subroutine solve_nhom_helmholtz(rhs, bc, dtnui, diag, smat, simat)
    ! solves the Helmholtz or Poisson eqn with non-homogeneous boundary conditions
    ! bc(x,y,z=zpts(nzp)) and bc(x,y,z=zpts(1)) using Green's function method:
    ! (D^2 - k^2 - dtnui)u = rhs, let u = phi + htop(z)f1(x,y) + hbot(z)f2(x,y) s.t.
    ! (D^2 - k^2 - dtnui)phi = rhs - (D^2 - k^2 - dtnui)(htop(z)f1(x,y) + hbot(z)f2(x,y))
    ! where phi has homogeneous boundary conditions
    use dim
    use grid
    use runparam, only: mapping
    use paral, only: wbar, wbarpp
    !     use mats
    use flags, only: debug
    use formatstrings
    use mpicom, only: myid
    ! i/o
    real, dimension(nxpl,ny,nzpl), intent(inout) :: rhs
    real, dimension(nxpl,ny,nzpl), intent(in) :: bc
    real, intent(in) :: dtnui
    !     integer, intent(in) :: mapping
    real, dimension(nzpl), intent(in) :: diag
    real, dimension(nzpl,nzpl), intent(in) :: smat,simat
    ! vars
    integer :: i, j, k
    real, allocatable, dimension(:) :: hbot, htop, d2hbot, d2htop
    real :: btop, bbot, atop, abot, width, h

    ! allocate vars
    allocate( hbot(nzp), htop(nzp), d2hbot(nzp), d2htop(nzp) )
    ! based on mapping, compute what hbot and h2 should be:
    select case(mapping)
        case(0) ! channel flow
            ! linear htop, hbot
            width = zpts(1) - zpts(nzp)
            htop(1:nzp) = (zpts(1:nzp) - zpts(nzp))/width ! hbot(1) = 1, hbot(nzp) = 0
            hbot(1:nzp) = -(zpts(1:nzp) - zpts(1))/width  ! h2(1) = 0, h2(nzp) = 1
            d2htop = 0.0
            d2hbot = 0.0
        case(1) ! semi-infinite BL
            ! exponential from blayscal
            h = zpts(1) - zpts(nzp)
            btop = 0.2
            bbot = 50.0
    !         htop = wbar/wbar(1)
    !         htop(nzp) = 0.0
    !         d2htop = wbarpp/wbar(1)
            htop(1:nzp) = 0.0!exp( -btop*(h-zpts(1:nzp)) )
            d2htop = 0.0!btop**2.0*htop
!             hbot = ( exp(-bbot*zpts) - exp(bbot*(zpts-2.0*h)) )
!             hbot = hbot/( 1.-exp(-2.0*bbot*h) )
!             d2hbot = bbot**2.0*hbot
            hbot(1:nzp) = exp( -bbot*zpts(1:nzp) )
            hbot(1) = 0.0
            d2hbot(:) = bbot**2.0*hbot(:)
        case(2) ! finite BL
            ! exponential for bottom from blayscal, CTR blending-function for top
            ! to ensure dhtop/dz | z=zpts(1) = 0
            h = zpts(1) - zpts(nzp)
            btop = 8.*atan(1.)
    !         btop = .001
            bbot = 5.
            htop(1:nzp) = zpts(1:nzp)/h - sin( btop*zpts(1:nzp)/h )/btop
            d2htop = btop/h**2.0*sin( btop*zpts(1:nzp)/h )
!             htop = wbar/wbar(1)
!             htop(nzp) = 0.0
!             d2htop = wbarpp/wbar(1)
    !         htop(1:nzp) = exp( -btop*(h-zpts(1:nzp)) )
    !         d2htop = btop**2.0*htop
    !         htop = exp(btop*(zpts-h))*(1.-exp(-2.0*btop*zpts))
    !         htop = htop/( 1.-exp(-2.0*btop*h) )
    !         d2htop = btop**2.0*htop
            hbot = ( exp(-bbot*zpts) - exp(bbot*(zpts-2.0*h)) )
            hbot = hbot/( 1.-exp(-2.0*bbot*h) )
!             hbot(1:nzp) = exp( -bbot*zpts(1:nzp) )
            d2hbot = bbot**2*hbot
        case(3) ! channel flow from z = 0 to zlen
            btop = 50.0
            bbot = 50.0
            ! exponential from blayscal but equal at both ends
    !         htop(1:nzp) = exp( -btop*(zpts(1)-zpts(1:nzp)) )
    !         htop(nzp) = 0.0
    !         hbot(1:nzp) = exp( -bbot*zpts(1:nzp) )
    !         hbot(1) = 0.0
    !         d2htop = btop**2.0*htop
    !         d2hbot = bbot**2.0*hbot
            ! exponential with proper zero BC at other end
            h = zpts(1) - zpts(nzp)
            hbot = ( exp(-bbot*zpts) - exp(bbot*(zpts-2.0*h)) )
            hbot = hbot/( 1.-exp(-2.0*bbot*h) )
            d2hbot = bbot**2.0*hbot
            htop = exp(btop*(zpts-h))*(1.-exp(-2.0*btop*zpts))
            htop = htop/( 1.-exp(-2.0*btop*h) )
            d2htop = btop**2.0*htop
            ! CTR belnding function using x - sin(2pi*x/xlen)/2pi
    !         h = zpts(1) - zpts(nzp)
    !         btop = 8.*atan(1.)
    !         htop(1:nzp) = zpts(1:nzp)/h - sin( btop*zpts(1:nzp)/h )/btop
    !         d2htop = btop/h**2.0*sin( btop*zpts(1:nzp)/h )
    !         hbot = 1. - htop
    !         d2hbot = -d2htop
        case(4) ! equally space points from z = 0 to zlen
            ! linear htop, hbot
            width = zpts(1) - zpts(nzp)
            htop(1:nzp) = (zpts(1:nzp) - zpts(nzp))/width ! hbot(1) = 1, hbot(nzp) = 0
            hbot(1:nzp) = -(zpts(1:nzp) - zpts(1))/width  ! h2(1) = 0, h2(nzp) = 1
            d2htop = 0.0
            d2hbot = 0.0

        case(22) ! well-resolved Blasius BL problem
            h = zpts(1) - zpts(nzp)
            bbot = 5.
            hbot = ( exp(-bbot*zpts) - exp(bbot*(zpts-2.0*h)) )
            hbot = hbot/( 1.-exp(-2.0*bbot*h) )
            d2hbot = bbot**2*hbot
            htop = wbar/wbar(1)
            htop(nzp) = 0.
            d2htop = wbarpp/wbar(1)
            
    end select

    if (debug>=2 .and. myid==0) then
        write(*,'(//,A)') '-----------------helmholtz----------------'
        write(*,*) 'hbot(z)='
        write(*,fm3) hbot(1:nzp)
        write(*,*) 'd^2 hbot/dz^2 ='
        write(*,fm3) d2hbot(1:nzp)
        write(*,*) 'htop(z)='
        write(*,fm3) htop(1:nzp)
        write(*,*) 'd^2 htop/dz^2 ='
        write(*,fm3) d2htop(1:nzp)
        write(*,*) 'h(z)*bc = '
        write(*,fm3) hbot(1:nzp)*bc(1,1,nzp) + htop(1:nzp)*bc(1,1,1)
    end if


    ! rhs = rhs - (D^2 - k^2 - dtnui)(htop(z)f1(x,y) + hbot(z)f2(x,y))
    do k=1,nzp
        do j=1,ny
            do i=1,nxpl
                abot = ( d2hbot(k) + (-wavxx(i) - wavyy(j) - dtnui)*hbot(k) )
                rhs(i,j,k) = rhs(i,j,k) - abot*bc(i,j,nzp)

                atop = ( d2htop(k) + (-wavxx(i) - wavyy(j) - dtnui)*htop(k) )
                rhs(i,j,k) = rhs(i,j,k) - atop*bc(i,j,1)
            end do
        end do
    end do

    if (debug>=2 .and. myid==0) then
        write(*,*) 'rhs(x,z) before helmholtz='
        write(*,fm1f) (rhs(1:nxpl,1,k),k=1,nzp)
    end if

    ! now solve for phi
    call solve_helmholtz(rhs, dtnui, nzm, diag, smat, simat)

    if (debug>=2 .and. myid==0) then
        write(*,*) 'rhs(x,z) after helmholtz='
        write(*,fm1f) (rhs(1:nxpl,1,k),k=1,nzp)
    end if

    ! form complete field u = phi + htop(z)f1(x,y) + hbot(z)f2(x,y)
    do k=1,nzpl
        rhs(:,:,k) = rhs(:,:,k) + hbot(k)*bc(:,:,nzp) + htop(k)*bc(:,:,1)
    end do


    deallocate(htop,hbot,d2htop,d2hbot)

end subroutine solve_nhom_helmholtz



end module helmholtz
