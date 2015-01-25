module orrsomm

implicit none
public
save

integer :: n_re_d, n_alphastar, i_re_d, j_alphastar
complex*16, allocatable, dimension(:) :: uvec,vvec,wvec
real, allocatable, dimension(:) :: alphastar,re_d,enr1,enr2,grow1,grow2
real :: avgdedt,dedterr,avgerr
complex*16 :: omega

contains

subroutine set_linstab_param()
    ! sets up the parameter space for linear stability test
    use mpicom, only: myid
    use runparam, only: astar
    use parameters, only: x0, cguess, xnu, u0
    ! explore Re #s just below and above Re_c
    re_d(i_re_d) = 1805. + (i_re_d-1)*3000./(n_re_d-1)
    ! explore a* just below and above critical region
    alphastar(j_alphastar) = 0.5 + (j_alphastar-1)*.7/(n_alphastar-1)
    ! adjust relevant qties
    x0 = re_d(i_re_d)**2*xnu/(6.02**2*u0)
    astar = alphastar(j_alphastar)
    cguess = cmplx( 0.4 - (i_re_d-1)*0.1/(n_re_d-1)  , 0.0)

    if (myid==0) write(*,*) 'Re_d = ',re_d(i_re_d), ', x0 =', x0, ', astar =', astar, ', cguess =', cguess

end subroutine set_linstab_param


subroutine add_linstab_perturb(uf,vf,wf)
    ! adds more unstable mode from linear stability theory to velocities
    use dim
    use flags, only: debug
    use io
    use mpicom, only: myid
    use grid, only: wavx
    ! i/o
    real, dimension(nxpl,ny,nzpl), intent(inout) :: uf,vf,wf

    ! add most unstable mode to initial velocities
    if (wavx(1)==0.) then
        uf(3,1,:) = uf(1,1,:) + real(uvec(:))
        uf(4,1,:) = uf(2,1,:) + aimag(uvec(:))
        vf(3,1,:) = vf(1,1,:) + real(vvec(:))
        vf(4,1,:) = vf(2,1,:) + aimag(vvec(:))
        wf(3,1,:) = wf(1,1,:) + real(wvec(:))
        wf(4,1,:) = wf(2,1,:) + aimag(wvec(:))
    endif
    ! print perturbation for debugging
    if (debug>=1) then
        if (myid==0) write(*,'(//,A)') '-------- orr-sommerfeld perturbation --------'
        call printuvw(uf,vf,wf)
    end if
    deallocate( uvec, vvec, wvec )

end subroutine add_linstab_perturb


subroutine write_linstab_results(nsteps)
    ! write energy growth rate and corresponding theoretical prediction to file
    use flags, only: debug
    use runparam, only: linstab
    use parameters, only: x0, xnu, u0
    use grid, only: alpha
    use iofiles
    use time, only: tseries
    use mpicom, only: myid
    ! i/o
    integer, intent(in) :: nsteps
    ! local vars
    integer :: i
    
    avgdedt = sum(grow2)/(nsteps-1) ! avg energy growth/decay
    avgerr = abs( (avgdedt - aimag(omega))/aimag(omega) )*100.0
    
    if (myid==0) then
        if (debug>=1) then
            write(*,'(//,A)') '---------- Linear Stability Results ------------'
            write(*,'(2(A,G12.4,'',''))') 're_d = ',6.02*sqrt(x0*xnu/u0)*u0/xnu,'alpha* = ',alpha
            write(*,*) 'i, energy, dEtotal/dt, dEnomean/dt'
            do i=1,nsteps
                write(*,'(I5,'','',3(G14.6,'',''))') i, enr2(i), grow1(i), grow2(i)
            end do
            write(*,*) 'omega_i, < dE/dt >, < % error >'
            write(*,'(3(G14.6,'',''))') aimag(omega), avgdedt, avgerr
        end if

        if (linstab==2) then
            if (i_re_d==1 .and. j_alphastar==1) then
                open(unit=iunit_energy,file=energy_file,status='unknown')
                write(iunit_energy,*) 'R, alpha*, omega_i, < dE/dt >, < % error >,'
            end if
            write(iunit_energy,'(5(G14.6,'',''))') re_d(i_re_d),alphastar(j_alphastar),aimag(omega), avgdedt, avgerr

            if (i_re_d==n_re_d .and. j_alphastar==n_alphastar) then
                close(iunit_energy)
                deallocate( alphastar,re_d,enr1,enr2,grow1,grow2 )
            end if
        end if

    end if

end subroutine write_linstab_results


subroutine solve_orrsomm(uvec,vvec,wvec,eps,alphax,betax,omega,iglb)
    ! solves the Orr-Sommerfeld equation for the most unstable mode
    ! everything is non-dimensionalized by the boundary layer thickness
    ! so Reynolds # = 1/nu (assuming u0 = 1.). 
    ! This should also affect derivatives... multiply by BL thickness
    use dim
    use parameters, only: xnu
    use paral !, only: ubar
    use grid, only: zpts
    use mats !, only: dor
    use mpicom, only: myid
    use flags, only: debug
    use prng
    ! i/o
    complex*16, dimension(nzp), intent(out) :: uvec, vvec, wvec
    complex*16, intent(inout) :: omega
    real, intent(in) :: eps, alphax, betax
    integer, intent(in) :: iglb
    ! local vars
    real, allocatable, dimension(:,:) :: did, d3, d4 !,d2 
    complex*16, allocatable, dimension(:,:) :: a, b, aa, bb, c, ck, ckp1, btrans
    complex*16, allocatable, dimension(:) :: rhs_phi, rhs_gam, residualk, residualkp1
    real, allocatable, dimension(:) ::  reig, aeig !,ubarp, ubarpp
!     complex, allocatable, dimension(:,:) :: fic, fic1, fic2, fic3, ficb, ficb1, ficb2, ficb3
    integer :: i, ii, j, k, iter, iiter, AllocateStatus
    real :: c01, c02, c03, c04 ! only used for membrane bottom wall
    real :: rey, wav2, dmag, etot, xnorm, residualk_norm, residualkp1_norm, inner_rk_norm, inner_rkp1_norm
    complex*16 :: omk, omkp1
    ! local vars req'd for LAPACK
    character*1 :: balanc, jobvl, jobvr, sense, trans
    integer :: n, info, lwork, lrwork, lda, ldb, ldvl, ldvr, ilo, ihi
    integer, allocatable, dimension(:) :: iwork
    logical, allocatable, dimension(:) :: bwork
    double precision :: abnrm, bbnrm
    complex*16, allocatable, dimension(:,:) ::   vl, vr
    complex*16, allocatable, dimension(:) :: phi, gam, phikp1, gamkp1, work
    double precision, allocatable, dimension(:) :: lscale, rscale, rconde, rcondv
    real, allocatable, dimension(:) :: rwork
    integer, allocatable, dimension(:) :: ipiv
    ! rigid wall parameters 
    ! (not needed if you change dimensions of the main arrays 
    !    (A and B) to A(NZP,NZP) and B(NZP,NZP))
    ! [for special version of linear with kaplan membrane wall boundary cond.
    !    C01 = D/M       C02 = 1/M       C03 = K/M        C04 = T/M ]
    c01 = 1.e-40
    c02 = 1.e-40
    c03 = 1.e-40
    c04 = 1.e-40
    ! set size of problem for the remainder of orr-sommerfeld solver
    n = nzpp2 ! for membrane bottom wall
!     n = nzp
    ! Reynolds number and square of the wavenumber.
    ! [for 2-D problem betax is zero (wavenumber in y-direction)]
    rey = 1./xnu
    wav2 = alphax**2.0 + betax**2.0
    ! allocate and zero a and b
    allocate( did(nzp,nzp), d3(nzp,nzp), d4(nzp,nzp), a(n,n), b(n,n), c(n,n), stat=AllocateStatus )
!     allocate( ubarp(nzp), ubarpp(nzp), d2(nzp,nzp),  stat=AllocateStatus )
    if (AllocateStatus /= 0) then
        write(*,*) "**Not Enough Memory - solve_orrsomm array allocation **"
    end if
    a = 0.; b = 0;
    ! identity matrix
    did = 0.0
    do i=1,nzp
        did(i,i) = 1.
    end do
    ! calculate higher derivatives in vertical
!     d2 = matmul(dor,dor)
!     d3 = matmul(dor,d2)
    d4 = matmul(d2,d2)
!     ubarp = matmul(dor,ubar)
!     ubarpp = matmul(d2,ubar)
    ! We set the Orr-Sommerfeld equation here 
    ! in a form A*psi = omega B*psi where
    ! psi is a streamfunction and omega is a complex
    ! eigenvalue omega=c*alphax. See any book on
    ! Orr-Sommerfeld eq.
    if (iglb<=0) then
        do j=1,nzp
            do i=1,nzp
                a(i,j) = d4(i,j) - 2.0*wav2*d2(i,j) + wav2*wav2*did(i,j) - (0.0,1.)*alphax*rey*(ubar(i)*( d2(i,j) - wav2*did(i,j) ) - ubarpp(i)*did(i,j))
                b(i,j) = -(0.0,1.)*rey*( d2(i,j) - wav2*did(i,j) )
            end do
        end do
        ! boundary conditions will need to be changed for channel
        ! flow. I think they should be A(1,1)=1, A(1,J)=0, B(1,J)=0
        ! (this gives psi=0 at the upper wall); A(2,J)=D(1,J), B(2,J)=0
        ! (this gives d psi/dz =0 at the upper wall); and corresponding
        ! conditions at the lower wall: A(NZP,NZP)=1, A(NZP,J)=0, B(NZP,J)=0;
        ! and A(NZ,J)=D(NZP,J),B(NZ,J)=0.0
        ! boundary conditions for Blasius/boundary layer flow (same as channel flow)
        ! psi = phi(z)*exp(i*(alpha*x - beta*t))
        ! u' = d psi /dz, w' = d psi /dx
        ! u' = w' = du'/dz = dw'/dz = 0 as z -> infty ==> psi = 0, d psi/dz = 0
        ! u' = w' = 0 at z = 0 ==> psi = 0, d psi/dz = 0

        ! boundary conditions for rigid wall boundary layer flow
        ! psi(z=0) = 0 -> A(nzp,1) = 1, A(nzp,2:) = 0, B(nzp,:) = 0.0
        ! psi(z=infty) = 0 -> A(1,1) = 1, A(1,2:) = 0, B(1,:) = 0.0
        ! D psi(z=0) = 0 -> A(nzpp,1:nzp) = D(nzp,1:nzp), B(nzpp,:) = 0.0
        ! D psi(z=infty) -> A(nzpp2,1:nzp) = D(1,1:nzp), B(nzpp2,:) = 0.0
        ! matrix A and B are now size (nzpp2,nzp) ? overdetermined!
        ! so try channel flow BCs instead
!         a(1,:) = 0.0
!         a(1,1) = 1.
!         b(1,:) = 0.0
!         a(2,:) = d(1,:)
!         b(2,:) = 0.0
!         a(nzpp,:) = 0.0
!         a(nzpp,nzpp) = 1.
!         b(nzpp,:) = 0.0
!         a(nzpp2,1:nzp) = d(nzp,1:nzp)
!         b(nzpp2,:) = 0.0

        ! boundary conditions appropriate for kaplan membrane wall below
        do j=1,nzp
            a(n,j) = 0.0
            do k=1,nzp
                a(n,j) = a(n,j) - (0.0,1.)*c01/(alphax*rey)*d2(nzp,k)*d(k,j)
            end do
        end do
        do j=1,n
            b(1,j) = 0.0
            b(2,j) = 0.0
            b(3,j) = 0.0
            b(n-1,j) = 0.0
            b(n,j) = 0.0
            a(1,j) = 0.0
            a(3,j) = 0.0
            a(n-1,j) = 0.0
        end do
        do j=1,nzp
            a(2,j) = d(nzp,j)
        end do
        a(1,1) = 1.0
        a(3,nzp) = 1.0
        b(3,n-1) = 1.0/alphax
        a(2,n-1) = ubarp(nzp)
        a(2,n) = 0.0
        a(n-1,n) = (1., 0.)
        b(n-1,n-1) = (0.0,-1.)
        a(n,n-1) = c03 + c04*wav2 - (0.0,1.)*wav2*(c02/rey)*(ubarp(nzp)/alphax)
        a(n,n) = c01
        b(n,n) = (0.0,1.)


        if (iglb<0) then 
            ! solve for ALL eigenvalues, but only need one with the highest growth rate
            ! or lowest damping c_i
            ! setup input for LAPACK ZGGEV which calculates eigenvalues, 
            ! and eigenvectors for problem A v = lamda*B v
            balanc = 'N'
            jobvl = 'N'
            jobvr = 'V'
            sense = 'N'
            lda = n
            ldb = n
            ldvl = n
            ldvr = n
            lwork = 4*n
            lrwork = 8*n
            allocate( aa(n,n), bb(n,n), reig(n), aeig(n), phi(n), gam(n), vl(ldvl,n), vr(ldvr,n), work(lwork), rwork(lrwork), stat=AllocateStatus )
            allocate( iwork(n+2), bwork(n), rconde(n), rcondv(n), lscale(n), rscale(n),  stat=AllocateStatus )
            if (AllocateStatus /= 0) then
                write(*,*) "**Not Enough Memory - solve_orrsomm_zggevx array allocation **"
            end if
            
            aa = a
            bb = b
!             call zggev( jobvl, jobvr, n, c, lda, b, ldb, phi, gam, vl, ldvl, vr, ldvr, work, lwork, rwork, info )
            call zggevx( balanc, jobvl, jobvr, sense, n, aa, lda, bb, ldb, phi, gam, vl, ldvl, vr, ldvr,  ilo,  ihi,  lscale,  rscale, abnrm, bbnrm, rconde, rcondv, work, lwork, rwork, iwork, bwork, info )
            if (info/=0 .and. myid==0) write(*,'(//,A,I5)') 'zggev error: info=',info
            
            reig = real(phi)
            aeig = aimag(phi)
            do k=1,n
                if (abs(gam(k))>0.) then
                    reig(k) = reig(k) / abs(gam(k))
                    aeig(k) = aeig(k) / abs(gam(k))
                end if
!                 if (real(gam(k))>0.) reig(k) = reig(k)/real( gam(k) )
!                 if (aimag(gam(k))>0.)    aeig(k) = aeig(k)/aimag( gam(k) )
            end do 

            ! print eigenvalues to screen
            if (myid==0) then
                write(*,'(//,A)') '---- eigenvalues from Orr-Sommerfeld equation ----'
                write(*,*) 'number, omega real, omega imag, c real, c imag'
                do k=1,n
!                     write(*,'(I5,'','',4(G14.6,'',''))') k, reig(k), aeig(k), reig(k)/alphax, aeig(k)/alphax
                    write(*,'(I5,'','',6(G14.6,'',''))') k, phi(k), gam(k), reig(k)/alphax, aeig(k)/alphax
                end do
            end if

            ! find min value of omega and its location
            j = maxloc(aeig/alphax, 1, reig/alphax>0 .and. reig/alphax<0.42)
!             omega = phi(j) / gam(j)
            ! corresponding eigenvector
!             wvec = 0.0
!             wvec(1:nzp) = vr(:,j)

            ! print normalized eigen vector
            if (myid==0 .and. debug>=2) then
                write(*,'(//,A,2G13.6)') 'most unstable c', reig(j)/alphax,aeig(j)/alphax
                write(*,'(A)') 'eigenvector from zggevx (streamfunction) real, imag'
                write(*,'(2(G14.6,'',''))') vr(1:n,j)
            end if

            deallocate( aa, bb, reig, aeig, phi, gam, vl, vr, work, rwork )
            deallocate( iwork, bwork, rconde, rcondv, lscale, rscale )

        end if ! (iglb < 0)

        
        ! calculate most unstable mode of Orr-Sommerfeld eqn iteratively
        if (myid==0) then
            write(*,'(//,A)') '---- solving for most unstable eigenvalue & mode iteratively ----'
            write(*,*) 'iteration, omega real, omega imag, cr, ci'
        end if

        iiter = 10    ! max outer loop # of iterations
        iter = 10    ! max inner loop # of iterations
!         n = n    ! dimension for LU decomp zgetrf
        lda = n     ! leading dimension of first input matrix for zgetrf
        ldb = n     ! leading dimension of input vector for zgetrs
        lwork = n     ! size of work array for forward/backward substitution
        ! allocate local vars
        allocate( residualk(n), residualkp1(n), reig(n), aeig(n), ck(n,n), ckp1(n,n), phi(n), gam(n), phikp1(n), gamkp1(n), btrans(n,n), rhs_phi(n), rhs_gam(n), ipiv(n), work(lwork) )
        ! set initial omega guess to input omega
        omk = omega
        ! set initial guess for eigenvectors to be [1,1,1,...,1]
        phi = 1.0
        gam = 1.0
        ! set initial guess to a random vector
!         call init_random_seed(n)
!         call random_number(reig)
!         call init_random_seed(n)
!         call random_number(aeig)
!         phi = reig + (0.0,1.)*aeig
!         phi = phi / sum(abs(phi))
!         gam = conjg(phi)
!         btrans = transpose(b)
        btrans = conjg(transpose(b))
        
        do i=1,iiter ! outer loop

            ck = a - omk*b ! form LHS matrix - unchanged during inner loop
!             ck = a ! form LHS matrix - unchanged during inner loop
            residualk = matmul(ck,phi) ! residual
            residualk_norm = sum(abs(residualk)) ! norm of residualk
            c = ck ! copy of ck which will be factorized
            call zgetrf( n, n, c, lda, ipiv, info ) ! LU decomp of matrix c stored into c
            if (info/=0 .and. myid==0) write(*,'(//,A,I5)') 'zgetrf error: info=',info

            do ii=1,iter ! inner loop

                inner_rk_norm = sum( abs( matmul(ck,phi)) )

                ! solve (A - omk*B)*phi^k+1 = B*phi^k for phi^k+1
                rhs_phi = matmul(b,phi) ! RHS = B*phi
!                 rhs_phi = matmul(omk*b,phi) ! RHS = B*phi
                trans = 'N'
                call zgetrs( trans, n, 1, c, lda, ipiv, rhs_phi, ldb, info )
                if (info/=0 .and. myid==0) write(*,'(//,A,I5)') 'zgetrs error: info=',info
                
                ! solve (A- omk*B)^H*gam^k+1 = B^H*gam^k for gam^k+1
                rhs_gam = matmul(btrans,gam) ! RHS = B^H*gam
!                 rhs_gam = matmul(omk*btrans,gam) ! RHS = B^H*gam
!                 trans = 'T'
                trans = 'C'
                call zgetrs( trans, n, 1, c, lda, ipiv, rhs_gam, ldb, info )
                if (info/=0 .and. myid==0) write(*,'(//,A,I5)') 'zgetrs error: info=',info
                
                ! normalize and store solns rhs_phi & rhs_gam into phikp1 and gamkp1
!                 rhs_gam(1) = 0.0
!                 rhs_phi(1) = 0.0
                phikp1 = rhs_phi/sum(abs(rhs_phi))
                gamkp1 = rhs_gam/sum(abs(rhs_gam))
!                 phi(1) = 0.0
!                 gam(1) = 0.0
                
                inner_rkp1_norm = sum( abs( matmul(ck,phikp1) ) )

                ! convergence check
                if ( sum(abs(phikp1-phi))<5.E-5 .or. inner_rkp1_norm < 5.E-6 ) then
                    ! check if previous answer was better and if so use it
                    if (inner_rkp1_norm > inner_rk_norm) then
                        phikp1 = phi
                        gamkp1 = gam
                    end if
                    phi = phikp1
                    gam = gamkp1
                    ! print that convergence was reached
                    if (myid==0 .and. debug>=1) then
                        write(*,*) 'inner loop convergence reached at ii =',ii
                        write(*,*) '|R_inner| = ', inner_rkp1_norm
                    end if
                    exit
                end if

                ! check for bad starting guess
                if (inner_rkp1_norm >= inner_rk_norm .and. ii<iiter-1) then
                    ! start over with new starting guess
                    if (myid==0 .and. debug>=1) then
                        write(*,*) 'bad eigenvector guess detected at ii =',ii
                        write(*,'(2(A,G17.10,'',''))') '|R_inner^k| =', inner_rk_norm, '|R_inner^k+1| =', inner_rkp1_norm
                    end if
                    k = n
                    reig = 0.0; aeig = 0.0
                    call init_random_seed(k)
                    call random_number(reig)
                    call init_random_seed(k)
                    call random_number(aeig)
                    phikp1 = phi + 0.001*(reig + (0.0,1.)*aeig)
                    phikp1 = phikp1 / sum(abs(phikp1))
                    gamkp1 = gam + 0.001*(reig + (0.0,1.)*aeig)
                    gamkp1 = gamkp1 / sum(abs(gamkp1))
                end if

                ! update phi and gam
                phi = phikp1
                gam = gamkp1

             end do ! inner loop

             ! update eigenvalue omega
            rhs_phi = matmul(b,phikp1) ! B*phi^k+1
            rhs_gam = matmul(a,phikp1) ! A*phi^k+1
            ! omega^k+1 = gam^k+1 * (A * phi^k+1) / ( gam^k+1 (B * phi^k+1) )
            omkp1 = dot_product(gamkp1,rhs_gam)/dot_product(gamkp1,rhs_phi)
            
            ! write results
            if (myid==0) write(*,'(I5,'','',4(G17.10,'',''))') i*iter, omkp1, omkp1/alphax
            ! normalize phikp1, again?!
!             dmag = sum(abs(phikp1))
!             gamkp1 = gamkp1 / dmag
!             phikp1 = phi / dmag
            ! normalize and sign eigenvector based on largest abs value for printing
            dmag = maxval(abs(phikp1))
            k = maxloc(abs(phikp1),1)
            dmag = sign( dmag, real( phikp1(k) ) )
            ! print normalized eigen vector
            if (myid==0 .and. debug>=2) then
                write(*,'(//,A)') 'normalized eigenvector (streamfunction) real, imag'
                write(*,'(65(2(G14.6,'',''),/))') phikp1(1:nzp)!/dmag !conjg(phi(1:n))/dmag
            end if

            ! calculate residual^k+1
            ckp1 = a - omkp1*b
            residualkp1 = matmul(ckp1,phikp1)
            residualkp1_norm = sum(abs(residualkp1))

            ! write new residual at each outer iteration
            if (myid==0 .and. debug>=1) write(*,'(//,2(A,G17.10,'',''),//)') '|R^k| =', residualk_norm, '|R^k+1| =', residualkp1_norm

            ! convergence check
            if (abs(omkp1 - omk) < 1.E-14 .and. residualkp1_norm < 1.5E-6) then
                if (residualkp1_norm >= residualk_norm) then
                    omkp1 = omk
                    phikp1 = phi
                    gamkp1 = gam
                end if
                if (myid==0) write(*,*) 'outer convergence reached at i =',i*iter
                exit
            end if

            ! check for bad guess or bad behavior
            if (abs(omkp1 - omk) < 1.E-14 .and. residualkp1_norm > 1.5E-6) then
                if (myid==0 .and. debug>=1) write(*,*) 'bad initial omega guess detected at i =',i*iter
                k = 1 
                call init_random_seed(k)
                call random_number(reig(k))
                omkp1 = (0.3 + 0.1*reig(k))*alphax
                if (myid==0) write(*,*) 'try with new omega guess =',omkp1
            end if

            ! write(*,*) oscillations
            if (residualkp1_norm >= residualk_norm) then
                if (myid==0 .and. debug>=1) write(*,*) 'regression detected & avoided at i =',i*iter
                omkp1 = omk
                phikp1 = phi
                gamkp1 = gam
            end if

            ! update omega
            omk = omkp1
            
        end do ! outer loop

        omega = omkp1

        ! normalize eigenvector (mode)
!         dmag = 0.0
!         do k=1,nzp
!             if (abs(phi(k)) > abs(dmag)) then
!                 dmag = abs(phi(k))
!                 dmag = sign( dmag, real(phi(k)) )
!             end if
!         end do
!         ! print normalized eigen vector
!         if (myid==0 .and. debug>=1) then
!             write(*,'(//,A)') 'normalized eigenvector (streamfunction) real, imag'
!             write(*,'(2(G14.6,'',''))') phi(1:n)/dmag !conjg(phi(1:n))/dmag
!         end if

        ! wvec = streamfunction psi
        wvec = 0.0
        wvec(1:nzp) = phikp1(1:nzp)
!         wvec = phi(2:nz)
!         wvec = aimag(phi(1:nzp)) + (0.0,1.)*real(phi(1:nzp))
!         wvec(1:nzp) = phi(2:n-1)
!         wvec(2:nzp) = phi
!         wvec = phi

        deallocate( residualk, residualkp1, reig, aeig, ck, ckp1, phi, gam, phikp1, gamkp1, btrans, rhs_phi, rhs_gam, ipiv, work )
    
    end if ! (iglb <= 0)

    ! add in checks to see if eqns are really satisfied

    ! find velocities
    ! uvec = d(psi)/dz
    uvec = 0.0
    uvec = matmul(d,wvec)
!     uvec(nzp) = 0.0
!     uvec = matmul(dor,wvec)
    ! wvec = -d(psi)/dx
    wvec = -(0.0,1.)*alphax*wvec
!     wvec(1) = 0.0
!     wvec(nzp) = 0.0
    ! for 2D, vvec = 0.0
    vvec = 0.0

    ! normalize uvec, vvec, wvec so that they form a T-S wave
    ! w a small fraction of total energy of the flow
    etot = energy(uvec,zpts,nzp) + energy(vvec,zpts,nzp) + energy(wvec,zpts,nzp)
    xnorm = sqrt(eps/etot)
    uvec = uvec*xnorm
    vvec = vvec*xnorm
    wvec = wvec*xnorm

    if (myid==0 .and. debug>=1) then
        write(*,'(//,A)') 'normalized orr-somm perturbations'
        write(*,*) 'uvec, wvec,'
        do i=1,nzp
            write(*,'(4(G14.6,'',''))') uvec(i), wvec(i)
        end do
    end if

    deallocate( did, d3, d4, a, b, c )

end subroutine solve_orrsomm



subroutine lintest(uf,vf,wf,enr1,enr2,grow1,grow2)
    ! calculate energy growth to compare to linear stability theory
    use dim
    use time
    use mpicom
    ! inputs/outputs
    real, dimension(nxpl,ny,nzpl), intent(inout) :: uf, vf, wf
    real, dimension(:), intent(inout) :: enr1, enr2, grow1, grow2
    ! other req'd variables
    real, allocatable, dimension(:,:,:) :: up, vp, wp
    real :: en1,en2
    integer :: itot

    allocate( up(nxpl,ny,nzpl), vp(nxpl,ny,nzpl), wp(nxpl,ny,nzpl) )
    !     enr1= mean + perturbation
    !     en2= only perturbation
    itot=1
    ! get total energy integrated through boundary layer
    call ener(uf,vf,wf,itot,en1)
    !     print *,'myid=',myid,'en1=',en1
    ! get mean-removed energy integrated through boundary layer
    up = uf
    vp = vf
    wp = wf
    if (myid==0) then
        up(1:2,1,:)=0.0
        vp(1:2,1,:)=0.0
        wp(1:2,1,:)=0.0
    endif
    itot=0
    call ener(up,vp,wp,itot,en2)
    !     print *,'myid=',myid,'en2=',en2

    call mpi_allreduce(en1,enr1(istep),1,mpi_double_precision,MPI_SUM,comm,ierr)
    call mpi_allreduce(en2,enr2(istep),1,mpi_double_precision,MPI_SUM,comm,ierr)

    !     print *,'myid=',myid,'enr1=',enr1(istep)
    !     print *,'myid=',myid,'enr2=',enr2(istep)
    ! calculate energy growth    
    !    grow1(istep)=log(en1/enold1)/(2.0*dt)
    !    grow2(istep)=log(en2/enold2)/(2.0*dt)
    if (istep>1) then
        grow1(istep) = 1./(2.0*enr1(istep))*(enr1(istep)-enr1(istep-1))/dt
!         grow1(istep)=log(enr1(istep)/enr1(istep-1))/(2.0*dt)
        grow2(istep) = 1./(2.0*enr2(istep))*(enr2(istep)-enr2(istep-1))/dt
!         grow2(istep) = log(enr2(istep)/enr2(istep-1))/(2.0*dt)
    else
        grow1(istep)=0.0
        grow2(istep)=0.0
    endif
    !    if (abs(en1).ge.1.0e14) iabort=1
    !    if (abs(en2).ge.1.0e14) iabort=1
    !     if (myid==0) print *, 'grow1=',grow1(its),'grow2=',grow2(its)

    deallocate( up, vp, wp )

end subroutine lintest


subroutine ener(uf,vf,wf,itot,en)
    use dim
    use grid, only: zpts
    use mpicom
    use paral, only: ubar
    use modhorfft
    ! inputs/outputs
    real, dimension(nxpl,ny,nzpl), intent(inout) :: uf,vf,wf
    integer,intent(in) :: itot
    real, intent(out) :: en
    ! other reqd vars
    real, dimension(nxpp,nypl,nzpl) :: u,v,w
    real, dimension(nzpl) :: ee
    integer :: k0, i, j, k
    ! where integration write(*,*)s in z-direction
    k0 = 2 !17 !nz/4 + 1
    ! set total energy per x-y plane to zero for each z
    ee = 0.0
    ! transfer velocity to physical space
    call horfft(u,uf,1)
    call horfft(v,vf,1)
    call horfft(w,wf,1)
    if (itot==1) then
        do k=1,nzp
            u(:,:,k)=u(:,:,k)+ubar(k)
        end do
    endif
    ! cal!total energy 1/2 v^2 per x-y plane
    do j=1,nypl
        do i=1,nx
            ee = ee+0.5*(u(i,j,:)**2.0+v(i,j,:)**2.0+w(i,j,:)**2.0)
        end do
    end do
    ! integrate total energy over z-direction
    en=0.0
    do k=nz,k0,-1
!     do k=k0,nzp
!         en = en + 0.5*(ee(k-1)+ee(k))*(zpts(k-1) - zpts(k))
        en = en + 0.5*(ee(k+1)+ee(k))*(zpts(k) - zpts(k+1))
    end do

end subroutine ener


real function energy(q,zpts,nzp)
    ! integrate energy of complex q (single mode) over z
    complex*16, dimension(nzp), intent(in) :: q
    real, dimension(nzp), intent(in) :: zpts
    integer, intent(in) :: nzp
    ! local vars
    integer :: k
    ! integrate over z using trapeze rule
    energy = 0.0
    do k=2,nzp-1
        energy = energy - 0.5*( abs(q(k+1))**2 + abs(q(k))**2 )*( zpts(k+1)-zpts(k) )
    end do
    energy = energy*2.0/7.
end function energy

end module orrsomm