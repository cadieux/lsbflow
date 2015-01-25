module derivatives

implicit none
public
save

contains

subroutine zinteg(tmp,tint,n1,n2,mapping,per,rev)
    ! integrates a given field tmp in the vertical direction
    ! using trapezoidal rule (per=0) or Gauss quadrature (per=1)
    ! an adjustement is made for mappings with singularities
    use dim
    use grid
    use mpicom, only: myid
    ! input/output
    integer, intent(in) :: n1, n2, per, mapping, rev
    real, dimension(n1,n2,nzp), intent(inout) :: tmp
    real, dimension(n1,n2,nzp), intent(out) ::tint
    ! other vars
    integer :: i, m, j, k
    real :: wts(nzm),tj
    real :: integ, pi

    if (per==1) then
        ! compute integration weights
        pi=4.*atan(1.)
        do j=1,nzm
            tj=pi*float(j)/float(nz)
            wts(j)=0.0
            do m=1,nzm
                wts(j)=wts(j)+sin(float(m)*tj)*(1.0-cos(pi*float(m)))/float(m)
            end do
            wts(j)=2.00/float(nz)*sin(tj)*wts(j)
        end do
        ! integration with change of variable
        tint = 0.0
        do j=1,n2
            do i=1,n1
                integ=0.0
                do k=nz,2,-1
                    integ=integ+wts(k-1)*0.5*xl/gf(k)*tmp(i,j,k)
                    tint(i,j,k) = integ
                end do
                if (mapping==1) tint(i,j,1) = tint(i,j,2)
            end do
        end do        
    else
        ! use good old trapeze rule
        tint = 0.0   
        do j=1,n2
            do i=1,n1
                integ=0.0
                if (rev==0) then
                    do k=nzp,2,-1
                        integ = integ + 0.5*(tmp(i,j,k-1)+tmp(i,j,k))*(zpts(k-1)-zpts(k))
                        tint(i,j,k-1) = integ
                    end do
                    if (mapping==1) tint(i,j,1) = tint(i,j,2)
                else
                    do k=2,nzp
                        integ = integ + 0.5*(tmp(i,j,k-1)+tmp(i,j,k))*(zpts(k-1)-zpts(k))
                        tint(i,j,k) = integ
                    end do
                end if
            end do
        end do
    end if

end subroutine zinteg


subroutine test_zinteg(avgerr)
    ! test the accuracy of vertical integration
    use dim
    use runparam, only: mapping
    use paral
    use mpicom, only: myid
    ! i/o
    real, intent(out) :: avgerr
    ! vars
    integer :: k, AllocateStatus
    real, allocatable, dimension(:,:,:) :: dfdz, f_exact, f

    allocate(dfdz(nxpp,nypl,nzpl), f_exact(nxpp,nypl,nzpl), f(nxpp,nypl,nzpl), stat=AllocateStatus)
    if (AllocateStatus /= 0) then
        write(*,*) '**Error: not enough memory for allocation in test_zinteg**'
    end if
    ! let 
    do k=1,nzp
        dfdz(:,:,k) = ubarp(k)
        f_exact(:,:,k) = (ubar(k)-ubar(nzp))
    end do
    call zinteg(dfdz,f,nx,nypl,mapping,0,0)
    avgerr = 0.0
    do k=1,nzp
        avgerr = avgerr + abs( f(1,1,k) - f_exact(1,1,k) )/nzp
    end do
    if (myid==0) then
        write(*,'(//,A)') '--------------test_zinteg results------------------'
        write(*,*) 'dfdz(z), f_exact(z), f_comp(z),'
        do k=1,nzp
            write(*,'(3(G12.4,'',''))') dfdz(1,1,k), f_exact(1,1,k), f(1,1,k)
        end do
        write(*,'(//,A,G12.4)') 'zinteg error =',avgerr
    end if
    deallocate(dfdz, f_exact, f)

end subroutine test_zinteg


subroutine dfdz(uf,dudzf)
    use dim
    use chebyshev
    use mats
    ! i/o
    real, dimension(nxpl,ny,nzpl), intent(in) :: uf
    real, dimension(nxpl,ny,nzpl), intent(out) :: dudzf
    ! vars
    integer :: i,j

    if (nz>32) then
        call dct_diff(nxpl,ny,nzp,uf,dudzf,1)
    else
        do j=1,ny
            do i=1,nxpl
                dudzf(i,j,:) = matmul(d,uf(i,j,:))
            end do
        end do
    end if

end subroutine dfdz


subroutine dpdz(n1,n2,n,u,dudz,p)
!     use dim
    use chebyshev
    use mats
    ! i/o
    integer, intent(in) :: n1, n2, n, p
    real, dimension(n1,n2,n), intent(in) :: u
    real, dimension(n1,n2,n), intent(out) :: dudz
    ! vars
    integer :: i,j,l

    if (n>32) then
        call dct_diff(n1,n2,n,u,dudz,p)
    else
        do j=1,n2
            do i=1,n1
                do l=1,p
                    dudz(i,j,:) = matmul(d,u(i,j,:))
                end do
            end do
        end do
    end if

end subroutine dpdz


subroutine d2fdz2(uf,dudz2f)
    use dim
    use chebyshev
    use mats
    ! i/o
    real, dimension(nxpl,ny,nzpl), intent(in) :: uf
    real, dimension(nxpl,ny,nzpl), intent(out) :: dudz2f
    ! vars
    integer :: i,j, err
    real, allocatable, dimension(:,:,:) :: dudzf

    allocate(dudzf(nxpl,ny,nzpl), stat=err)
    if (err /= 0) print *, "dudzf: Allocation request denied"
    
    if (nz>32) then
        call dct_diff(nxpl,ny,nzp,uf,dudz2f,2)
    else
        do j=1,ny
            do i=1,nxpl
                dudz2f(i,j,:) = matmul(d2,uf(i,j,:))
            end do
        end do
    end if

    if (allocated(dudzf)) deallocate(dudzf, stat=err)
    if (err /= 0) print *, "dudzf: Deallocation request denied"

end subroutine d2fdz2


subroutine test_dfdz()
    ! test accuracy of derivative
    use dim
    use grid
    use parameters, only: xlen
    use mpicom, only: myid
    use modhorfft
    use core
    use mats
    use chebyshev
    real, allocatable, dimension(:,:,:) ::  qf, dqf, q, dq, dq_exact
    real :: pi, avgerr
    integer :: i,j,k, AllocateStatus
    pi = 4.*atan(1.)
    ! allocate vars
    allocate( q(nxpp,nypl,nzpl), qf(nxpl,ny,nzpl), dq(nxpp,nypl,nzpl), dqf(nxpl,ny,nzpl), dq_exact(nxpp,nypl,nzpl), stat=AllocateStatus)
    if (AllocateStatus /= 0) then
        write(*,*) '**Error: not enough memory for allocation in test_dfdz**'
    end if
    ! test function
    qf = uf
    call horfft(q,qf,1)
    do j = 1, ny    
        do i = 1, nxpl
            dqf(i,j,:) = matmul(d2,qf(i,j,:))
        end do
    end do
    call horfft(dq_exact,dqf,1)
    dqf = 0.0
!     do k = 1, nzp
!         do i = 1, nx
!             q(i,:,k) = exp(-((xpts(i)-xpts(nx/2))/(4.*dx))**2.0)*ubarp(k)
!             dq_exact(i,:,k) = exp(-((xpts(i)-xpts(nx/2))/(4.*dx))**2.0)*ubarpp(k)
!         end do
!     end do
    ! transfer to spectral space
!     call horfft(q,qf,-1)
    ! take derivative
!     call dfdz(qf,dqf)
    call d2fdz2(qf,dqf)
!     call dfdz(dqf,qf)
    ! transfer to phys space
    call horfft(dq,dqf,1)

    ! compute error
    avgerr = 0.0
    do k=1,nzp
        do j=1,nypl
            do i=1,nx
                avgerr = avgerr + abs( dq_exact(i,j,k) - dq(i,j,k) )/nzp/nx/nypl
            end do
        end do
    end do

    if (myid==2) then
        write(*,'(//,A)') '--------------test_dfdz results------------------'
        write(*,*) 'f(z), dfdz_exact, dfdz '
        do k=1,nzpl
            write(*,'(3(G15.7,'',''))') q(3,2,k), dq_exact(3,2,k), dq(3,2,k)
        end do
        write(*,'(//,A,G12.4)') 'dfdz avg error =',avgerr
!         do k=1,nzpl
!             write(*,'(32(G15.7,'',''))') dq(1:nx,2,k)
!         end do
    end if

    deallocate(q,qf,dq,dqf,dq_exact)
end subroutine test_dfdz


subroutine dfdx(uf,dudxf)
    ! uf (any primary var) in spectral space -> dudxf is its x-derivative
    use dim
    use grid, only: wavx
    real, dimension(nxpl,ny,nzpl), intent(in) :: uf
    real, dimension(nxpl,ny,nzpl), intent(out) :: dudxf
    integer :: i, k
    ! compute x derivative in Fourier space
    ! du/dx = ikx ( uf(2*i-1,:,:) + i*uf(2*i,:,:) )
    !        = ikx uf(2*i-1,:,:) - kx uf(2*i,:,:)
    ! dudxf(2*i,:,:) = kx uf(2*i-1,:,:)
    ! dudxf(2*i-1,:,:) = - kx uf(2*i,:,:)
    do k=1,nzpl
        do i=1,nxpl,2
            dudxf(i,:,k) = - wavx(i)*uf(i+1,:,k)
            dudxf(i+1,:,k) = wavx(i)*uf(i,:,k)
        end do
    end do

end subroutine dfdx


subroutine dfdy(uf,dudyf)
    ! uf (any primary var) in spectral space -> dudyf is its y-derivative
    use dim
    use grid, only: wavy
    real, dimension(nxpl,ny,nzpl), intent(in) :: uf
    real, dimension(nxpl,ny,nzpl), intent(out) :: dudyf
    integer :: i,k
    ! compute y derivative in Fourier space
    ! du/dy = iky ( uf(2*i-1,:,:) + i*uf(2*i,:,:) )
    !        = iky uf(2*i-1,:,:) - ky uf(2*i,:,:)
    ! dudyf(2*i,:,:) = ky uf(2*i-1,:,:)
    ! dudyf(2*i-1,:,:) = - ky uf(2*i,:,:)

    do k=1,nzpl
        do i=1,nxpl,2
            dudyf(i,:,k) = - wavy(:)*uf(i+1,:,k)
            dudyf(i+1,:,k) = wavy(:)*uf(i,:,k)
        end do
    end do

    end subroutine dfdy


subroutine lapl(uf,lapluf)
    ! takes laplacian (D^2-k^2) of uf
    use dim
    use grid, only: wavxx,wavyy
    ! input/output
    real, dimension(nxpl,ny,nzpl), intent(in) :: uf
    real, dimension(nxpl,ny,nzpl), intent(out) :: lapluf
    ! local vars
    integer :: j,k
    !     real :: wav2
    real, allocatable, dimension(:,:,:) :: dudz2f !, dudxf, dudyf
    allocate( dudz2f(nxpl,ny,nzpl) ) !, dudxf(nxpl,ny,nzpl), dudyf(nxpl,ny,nzpl) )
    call d2fdz2(uf,dudz2f) ! lapluf = d^2 u /dz^2
    do k=1,nzp
        do j=1,ny
            lapluf(:,j,k) = dudz2f(:,j,k) - (wavxx(:) + wavyy(j))*uf(:,j,k)
        enddo
    enddo
    deallocate( dudz2f )!, dudxf, dudyf )

end subroutine lapl


subroutine div(uf,vf,wf,divf)
    ! divf = ikx u + iky v + D w
    use dim
    use grid, only: wavx, wavy
    real, dimension(nxpl,ny,nzpl), intent(in) :: uf,vf,wf
    real, dimension(nxpl,ny,nzpl), intent(out) :: divf
    ! local vars
    !     integer :: i,j,k
    !     complex :: uderv, vderv
    real, allocatable, dimension(:,:,:) :: dqf
    allocate( dqf(nxpl,ny,nzpl) )
    ! compute dwdz
    call dfdz(wf,divf)
    ! compute x & y derivative in Fourier space
    call divuv(uf,vf,dqf)
    divf = dqf + divf
    deallocate( dqf )
end subroutine div

subroutine divcheck(uf,vf,wf,diverr)
    ! checks how far from zero is divergence of given field
    use dim
    use time, only: t,dt
    use modhorfft
    use formatstrings
    use mpicom, only: myid
    use flags, only: debug
    ! i/o
    real, dimension(nxpl,ny,nzpl), intent(in) :: uf,vf,wf
    real, intent(out) :: diverr
    ! local vars
    integer :: k!, i,j
    real, allocatable, dimension(:,:,:) :: divf, divp, dudxf,dvdyf,dwdzf
    allocate( divf(nxpl,ny,nzpl), divp(nxpp,nypl,nzpl), dudxf(nxpl,ny,nzpl), dvdyf(nxpl,ny,nzpl), dwdzf(nxpl,ny,nzpl) )
    ! calculate divergence
    call dfdx(uf,dudxf)
    call dfdy(vf,dvdyf)
    call dfdz(wf,dwdzf)
    divf = dudxf + dvdyf + dwdzf
    !     call div(uf,vf,wf,divf)
    call horfft(divp,divf,1)
    ! divp(nxp:nxpp,:,:) = 0.0
    ! calculate avg error
    diverr = sum( abs(divp) )/(nx*nypl*nzpl)
!     do k=1,nzp
!         do j=1,nypl
!             do i =1,nx
!                 diverr = diverr + abs(divp(i,j,k))/(nx*nypl*nzpl)
!             end do            
!         end do
!     end do
    !     diverr = 0.0
    !     do k=1,nzp
    !         do j=1,nypl
    !             diverr = diverr + sum(abs(divp(1:nx,j,k)))/(nx*nypl*nzpl)
    !         end do
    !     end do

    if (myid==0) write(*,'(//,A,G14.6,A,G14.6,A,G14.6)') 'divergence avg error =', diverr, ', at t=',t,', and dt=',dt
    if (myid==0 .and. debug>=2) write(*,fm1) (divp(1:nx,1,k),k=1,nzp)

    deallocate(divf,divp,dudxf,dvdyf,dwdzf)

end subroutine divcheck

subroutine divuv(uf,vf,divf)
    ! divf = ikx u + iky v + D w
    use dim
    use grid, only: wavx, wavy
    real, dimension(nxpl,ny,nzpl), intent(in) :: uf,vf
    real, dimension(nxpl,ny,nzpl), intent(out) :: divf
    integer :: i, k
    ! compute x & y derivative in Fourier space
    do k=1,nzpl
        do i=1,nxpl,2
            divf(i,:,k) = - wavx(i)*uf(i+1,:,k) - wavy(:)*vf(i+1,:,k)
            divf(i+1,:,k) = wavx(i)*uf(i,:,k)   + wavy(:)*vf(i,:,k)
        end do
    end do

end subroutine divuv


subroutine test_dfdx()
    ! test accuracy of derivative
    use dim
    use grid
    use parameters, only: xlen
    use mpicom, only: myid
    use modhorfft
    real, allocatable, dimension(:,:,:) ::  qf, dqf, q, dq, dq_exact
    real :: pi, avgerr
    integer :: i
    pi = 4.*atan(1.)
    ! allocate vars
    allocate( q(nxpp,nypl,nzpl), qf(nxpl,ny,nzpl), dq(nxpp,nypl,nzpl), dqf(nxpl,ny,nzpl), dq_exact(nxpp,nypl,nzpl))
    ! test function
    do i=1,nx
        q(i,:,:) = exp(-((xpts(i)-xpts(nx/2))/(4.*dx))**2.0) 
        !cos(6.*pi*(xpts(i)-xpts(1))/xlen)
        dq_exact(i,:,:) = -2.0*(xpts(i)-xpts(nx/2))*exp(-((xpts(i)-xpts(nx/2))/(4.*dx))**2.0)/(4.*dx)**2.0 
        !-6.*pi/xlen*sin(6.*pi*(xpts(i)-xpts(1))/xlen)
    end do

    ! transfer to spectral space
    call horfft(q,qf,-1)
    ! take derivative
    call dfdx(qf,dqf)
    ! transfer to phys space
    call horfft(dq,dqf,1)

    ! compute error
    avgerr = 0.0
    do i=1,nx
        avgerr = avgerr + abs( dq_exact(i,1,1) - dq(i,1,1) )/nx
    end do

    if (myid==0) then
        write(*,'(//,A)') '--------------test_dfdx results------------------'
        write(*,*) 'f(x), dfdx_exact, dfdx '
        do i=1,nx
            write(*,'(3(G15.7,'',''))') q(i,2,3), dq_exact(i,2,3), dq(i,2,3)
        end do
        write(*,'(//,A,G12.4)') 'dfdx avg error =',avgerr
    end if

    deallocate(q,qf,dq,dqf,dq_exact)

end subroutine test_dfdx


subroutine test_dfdy()
    ! test accuracy of derivative
    use dim
    use grid
    use parameters, only: ylen
    use mpicom
    use modhorfft
    real, allocatable, dimension(:,:,:) ::  qf, dqf, q, dq, dq_exact
    real :: pi, avgerr, totalavgerr
    integer :: j
    pi = 4.*atan(1.)
    ! allocate vars
    allocate( q(nxpp,nypl,nzpl), qf(nxpl,ny,nzpl), dq(nxpp,nypl,nzpl), dqf(nxpl,ny,nzpl), dq_exact(nxpp,nypl,nzpl))
    ! test function
    do j=1,nypl
        q(:,j,:) = exp(-((y(j)-ypts(ny/2))/(4.*dy))**2.0) 
        !cos(6.*pi*y(j)/ylen)
        dq_exact(:,j,:) = -2.0*(y(j)-ypts(ny/2))*exp(-((y(j)-ypts(ny/2))/(4.*dy))**2.0)/(4.*dy)**2.0
        !-6.*pi/ylen*sin(6.*pi*y(j)/ylen)
    end do

    ! transfer to spectral space
    call horfft(q,qf,-1)
    ! take derivative
    call dfdy(qf,dqf)
    ! transfer to phys space
    call horfft(dq,dqf,1)

    ! compute error
    avgerr = 0.0
    do j=1,nypl
        avgerr = avgerr + abs( dq_exact(1,j,1) - dq(1,j,1) )/ny
    end do

    call mpi_reduce(avgerr,totalavgerr,1,mpi_double_precision,MPI_SUM,0,comm,ierr)

    if (myid==0) then
        write(*,'(//,A)') '--------------test_dfdy results------------------'
        write(*,*) 'f(y), dfdy_exact, dfdy '
        do j=1,nypl
            write(*,'(3(G15.7,'',''))') q(1,j,1), dq_exact(1,j,1), dq(1,j,1)
        end do
        write(*,'(//,A,G12.4)') 'dfdy avg error =',avgerr
    end if

    deallocate(q,qf,dq,dqf,dq_exact)

end subroutine test_dfdy


subroutine test_lapl()
    ! computes laplacian of known function to check accuracy
    use dim
    use mpicom
    use modhorfft
    use paral
    use grid
    ! use parameters, only: xlen, ylen
    use formatstrings
    ! local vars
    integer :: i,j,k
    real :: avgerr, totalavgerr
    real, allocatable, dimension(:,:,:) :: q, laplq, qf, laplqf, laplq_exact
    ! real :: lx, ly
    real :: a, b, x, f

    allocate( q(nxpp,nypl,nzpl), laplq(nxpp,nypl,nzpl), laplq_exact(nxpp,nypl,nzpl), qf(nxpl,ny,nzpl), laplqf(nxpl,ny,nzpl) )

    b = xpts(nx/2)
    a = 4.*dx
    ! lx = xpts(nx) - xpts(1)
    ! ly = ypts(ny) - ypts(1)
    ! pi = 4.*atan(1.)
    do k=1,nzpl
        do j=1,nypl
            do i=1,nx
                x = xpts(i)
                f = exp(-((x-b)/a)**2.0)
                q(i,j,k) = ubar(k) + 0.1*ubar(k)*f
                laplq_exact(i,j,k) = ubarpp(k) + 0.1*ubarpp(k)*f + 0.1*ubar(k)*f*(4.*(b-x)**2.0-2.0*a**2.0)/a**4.
    !             q(i,j,k) = wbar(k) + 0.1*wbar(k)*sin(4.*pi*(xpts(i)-xpts(1))/lx) !+ 0.1*wbar(k)*sin(4.*pi*y(j)/ly)
    !             laplq_exact(i,j,k) = wbarpp(k) + 0.1*wbarpp(k)*sin(4.*pi*(xpts(i)-xpts(1))/lx)  - (4.*pi/lx)*0.1*wbar(k)*sin(4.*pi*(xpts(i)-xpts(1))/lx) ! + 0.1*wbarpp(k)*sin(4.*pi*y(j)/ly) - (4.*pi/ly)*0.1*wbar(k)*sin(4.*pi*y(j)/ly)
            end do
        end do
    end do

    call horfft(q,qf,-1)
    call lapl(qf,laplqf)
    call horfft(laplq,laplqf,1)

    ! compute error
    avgerr = 0.0
    do k=1,nzpl
        do j=1,nypl
            do i=1,nx
                avgerr = avgerr + abs( laplq_exact(i,j,k) - laplq(i,j,k) ) /(nx*ny*nzp)
            end do
        end do
    end do

    call mpi_reduce(avgerr,totalavgerr,1,mpi_double_precision,MPI_SUM,0,comm,ierr)

    if (myid==0) then
        write(*,'(//,A)') '--------------test_lapl results------------------'
        write(*,'(A,G12.4)') 'lapl avg error =',totalavgerr
        write(*,*) 'q(x,1,nz) ='
        write(*,fm1) q(1:nx,2,nz/2)
        write(*,*) 'q(1,y,nz) ='
        write(*,fm2) q(4,1:nypl,nz/2)
        write(*,*) 'q(1,1,z) ='
        write(*,fm3) q(4,2,1:nzp)
        write(*,*) 'laplq_exact(x,1,nz) ='
        write(*,fm1) laplq_exact(1:nx,2,nz/2)
        write(*,*) 'laplq_exact(1,y,nz) ='
        write(*,fm2) laplq_exact(4,1:nypl,nz/2)
        write(*,*) 'laplq_exact(1,1,z) ='
        write(*,fm3) laplq_exact(4,2,1:nzp)
        write(*,*) 'laplq(x,1,nz) ='
        write(*,fm1) laplq(1:nx,2,nz/2)
        write(*,*) 'laplq(1,y,nz) ='
        write(*,fm2) laplq(4,1:nypl,nz/2)
        write(*,*) 'laplq(1,1,z) ='
        write(*,fm3) laplq(4,2,1:nzp)
    end if

    deallocate( q, laplq, qf, laplqf, laplq_exact )


end subroutine test_lapl


end module derivatives