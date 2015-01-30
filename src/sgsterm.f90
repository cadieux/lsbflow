module sgs_models

use sgsinput
implicit none
public
save

contains

subroutine tau_sgs(ui,uiuj,duidxj,tauij)
    ! calls subgrid-scale modeling subroutines
    use dim
    use mpicom
    use flags, only: debug
    ! i/o
    real, dimension(nxpp,nypl,nzp,3), intent(in) :: ui
    real, dimension(nxpp,nypl,nzp,3, 3), intent(in) :: uiuj, duidxj
    real, dimension(nxpp,nypl,nzp,3, 3), intent(out) :: tauij
    
    ! call model based on input
    select case(sgsmodel)
        case(0) ! classic Smagorinsky model
            call static_smagorinsky(CS,duidxj,tauij)

        case(1) ! dynamic Smagorinsky model
            call dynamic_smagorinsky(ui, uiuj, duidxj, tauij)
            if(myid==0 .and. debug>=2) write(*,*) "dyn_smag complete"  

        case(2) ! sigma model
            call sigma_model(duidxj,tauij)

        case(3) ! interscale energy transfer model
            call interscale_model(ui,tauij)

    end select

end subroutine tau_sgs

subroutine static_smagorinsky(CS, duidxj, tauij)
    ! classic static Smagorinsky model
    use dim
    use mpicom
    use grid
    use parameters, only: xnu
    use wallunit
    ! i/o
    real, intent(in) :: CS
    real, dimension(nxpp,nypl,nzp,3, 3), intent(in) :: duidxj
    real, dimension(nxpp,nypl,nzp,3, 3), intent(out) :: tauij
    ! local vars
    real, dimension(:,:), allocatable :: utau
    real, dimension(:,:,:), allocatable :: snorm, nu_sgs
    real :: delta(nzp), utau_loc_avg, utau_avg!, utau_x_avg(nypl), utau_avg
    integer :: i, j, k, err, wallcount

    if (.not. allocated(zplus)) allocate( snorm(nxpp,nypl,nzp), nu_sgs(nxpp,nypl,nzp), zplus(nzp), utau(nxpp,nypl), stat=err )
    if ( err/=0 ) write(*,*) "** static_smagorinsky: allocation error **"

    ! compute friction velocity: u_tau = nu*du/dz|wall
    utau = 0.0
    wallcount = 0
    do k = 1, nzp
        if ( zpts(k)==0.0 .or. abs(zpts(k))==0.5 ) then
            utau = utau + xnu*duidxj(:,:,k,1,3)
            wallcount = wallcount + 1
        end if
    end do
    utau = utau/wallcount
    
    ! average utau in x direction
!     utau_x_avg = sum(utau(1:nx,:),1)/nx
    ! sum over all values and divide by total
    utau_loc_avg = sum(utau(1:nx,:))/(nx*ny)

    ! average utau in y direction
    call mpi_allreduce(utau_loc_avg,utau_avg,1,mpi_double_precision,mpi_sum,comm,ierr)
!     call mpi_allreduce(utau_x_avg,utau_avg,nypl,mpi_double_precision,mpi_sum,comm,ierr)
!     utau_avg = utau_avg/ny

    ! final avg friction velocity
    utau_avg = sqrt(abs(utau_avg))

    ! compute wall units for Van Driest Damping
    zplus = zpts*utau_avg/xnu

    ! compute |S| = sqrt(2SijSij)) where Sij = 0.5*(du_i/dx_j + du_j/dx_i)
    snorm = 0.0
    do j=1,3
        do i=1,3
            snorm = snorm + ( 0.5*( duidxj(:,:,:,i,j) + duidxj(:,:,:,j,i) ) )**2
        end do
    end do
    snorm = sqrt( 2.*snorm )

    ! compute delta = (dx*dy*dz)^(1/3)
    dz = abs( zpts(1) - zpts(2) )
    delta(1) = (dx*dy*dz)**(1./3.)
    do k = 2, nz
        dz = (zpts(k-1) - zpts(k+1))/2.
        delta(k) = (dx*dy*dz)**(1./3.)
    end do
    dz = abs( zpts(nzp) - zpts(nzp-1) )
    delta(1) = (dx*dy*dz)**(1./3.)

    ! compute nu_sgs = -(CS*delta*f_van_driest)^2*|S|
    do k = 1, nzp
        nu_sgs(:,:,k) = -((CS*delta(k)*(1.-exp(-zplus(k)/25.)))**2)*snorm(:,:,k)
    end do

    ! compute tauij = nu_sgs*S_ij
    do j = 1, 3
        do i = 1, 3
            tauij(:,:,:,i,j) = nu_sgs*0.5*( duidxj(:,:,:,i,j) + duidxj(:,:,:,j,i) )
        end do
    end do

    if ( allocated(zplus) ) deallocate( zplus, snorm, nu_sgs, utau, stat=err )
    if ( err/=0 ) write(*,*) "** static_smagorinsky: deallocation error **"

end subroutine static_smagorinsky


subroutine dynamic_smagorinsky(ui, uiuj, duidxj, tauij)
    ! computes   tau_ij = -2*nu_sgs*S_ij
    ! where        nu_sgs = C*delta^2*|S|   
    !              C = < M_ij L_ij > / < M_ij M_ij >
    !              L_ij = [u_i u_j] - [u_i][u_j]
    !              M_ij = 2 [delta^2 |S|S_ij] - 2 [delta]^2 [|S|][S_ij]
    use mpicom
    use dim
    use flags, only: debug
    use runparam, only: iomod
    use time, only: istep
    use parameters, only: xnu
    use grid
    use filters
    use formatstrings
    ! i/o
    real, dimension(nxpp,nypl,nzp,3), intent(in) :: ui
    real, dimension(nxpp,nypl,nzp,3,3), intent(in) :: uiuj, duidxj
    real, dimension(nxpp,nypl,nzp,3,3), intent(out) :: tauij
    ! local vars
    integer :: i, j, k, l, m, nq, AllocateStatus, interpflag
    real, dimension(:,:,:), allocatable :: snorm, snormhat, ctot, numc, denomc, cd, delta, nu_sgs
    real, dimension(:,:,:,:), allocatable :: uihat
    real, dimension(:,:,:,:,:), allocatable :: sij, sijhat, mij, lij
    ! allocate variables
    nq = 3
    if (.not. allocated(snorm)) allocate( snorm(nxpp,nypl,nzp), snormhat(nxpp,nypl,nzp), ctot(nxpp,nypl,nzp), cd(nxpp,nypl,nzp), delta(nxpp,nypl,nzp), numc(nxpp,nypl,nzp), denomc(nxpp,nypl,nzp), nu_sgs(nxpp,nypl,nzp), uihat(nxpp,nypl,nzp,nq), sij(nxpp,nypl,nzp,nq,nq), sijhat(nxpp,nypl,nzp,nq,nq), mij(nxpp,nypl,nzp,nq,nq), lij(nxpp,nypl,nzp,nq,nq), stat=AllocateStatus )
    if ( AllocateStatus/=0 ) write(*,*) "** dynamic_smagorinsky: allocation error **"

    interpflag = 1 ! needed if grid is non-uniform

    if(myid==0 .and. debug>=2) write(*,*) "array allocation successful"

    ! compute M_ij = 2 [delta^2 |S|S_ij] - 2 [delta]^2 [|S|][S_ij]

    ! compute delta = (dx*dy*dz(k))^(1/3)
    delta(:,:,1) = (dx*dy*(zpts(1)-zpts(2)))**(1./3.)
    do k = 2, nzp-1
        dz = 0.5*( zpts(k-1) - zpts(k+1) )
        delta(:,:,k) = (dx*dy*dz)**(1./3.)
    end do
    delta(:,:,nzp) = (dx*dy*(zpts(nzp-1)-zpts(nzp)))**(1./3.)
    ! obtain [delta] by filtering - in practice it is often assumed to be twice delta
!     call filter_q_3pt(delta,cd,1,1,a,b,interpflag)
    cd = 2.0*delta

    if(myid==0 .and. debug>=2) write(*,*) "delta filtering successful"   

    ! compute S_ij, |S| = sqrt(2*SijSij)
    snorm = 0.0
    do m = 1, nq
        do l = 1, nq
            sij(:,:,:,l,m) = 0.5*( duidxj(:,:,:,l,m) + duidxj(:,:,:,m,l) )
            snorm = snorm + 2.*sij(:,:,:,l,m)*sij(:,:,:,l,m)
        end do
    end do
    snorm = sqrt(snorm)

    if(myid==0 .and. debug>=2) write(*,*) "snorm and sij computed"

    ! compute [S_ij]
    call filter_q_3pt( sij, sijhat, nq, nq, a, b, interpflag )
!     call filter_q_3pt( snorm, snormhat, nq, 1, a, b, interpflag )
    
    ! compute |[S]|
    snormhat = 0.0
    do m = 1, nq
        do l = 1, nq
            snormhat = snormhat + 2.*sijhat(:,:,:,l,m)*sijhat(:,:,:,l,m)
        end do
    end do
    snormhat = sqrt(snormhat)

    if(myid==0 .and. debug>=2) write(*,*) "sij, snorm filtering successful" 
    
    do m = 1, nq
        do l = 1, nq
            ! let lij = delta^2 |S|S_ij
            lij(:,:,:,l,m) = delta*delta*snorm*sij(:,:,:,l,m)
            ! let sijhat = [delta]^2 |[S]|[S_ij]
            sijhat(:,:,:,l,m) = cd*cd*snormhat*sijhat(:,:,:,l,m)
        end do
    end do
    ! compute [ delta^2 |S|S_ij ] and store in M_ij
    call filter_q_3pt(lij, mij, nq, nq, a, b, interpflag)

    if(myid==0 .and. debug>=2) write(*,*) "deta^2 |S| S_ij filtering successful" 

    ! M_ij = 2 [delta^2 |S|S_ij] - 2 [delta]^2 [|S|][S_ij]
    mij = 2.0*( mij - sijhat )


    ! Compute L_ij = [u_i u_j] - [u_i][u_j]
    lij = 0.0 ! to be safe
    ! compute uihat, uiujhat
    call filter_q_3pt( ui, uihat, nq, 1, a, b, interpflag )
    call filter_q_3pt( uiuj, lij, nq, nq, a, b, interpflag )

    if(myid==0 .and. debug>=2) write(*,*) "uiuj filtering successful" 

    do m = 1, nq
        do l = 1, nq
            ! L_ij = [u_i u_j] - [u_i][u_j]
            lij(:,:,:,l,m) = lij(:,:,:,l,m) - uihat(:,:,:,l)*uihat(:,:,:,m)
        end do
    end do

    ! compute C = < M_ij L_ij > / < M_ij M_ij >
    numc = 0.0
    denomc = 0.0
    do m = 1, nq
        do l = 1, nq
            numc = numc + mij(:,:,:,l,m)*lij(:,:,:,l,m) ! M_ij L_ij
            denomc = denomc + mij(:,:,:,l,m)*mij(:,:,:,l,m) ! M_ij M_ij
        end do
    end do

    ! perform local average using filtering with equal weights
!     call filterbwa3pt_interp_ub(numc, ctot, a, b) ! < M_ij L_ij >
!     call filterbwa3pt_interp_ub(denomc, cd, a, b) ! < M_ij M_ij >
    call filter_q_3pt(numc, ctot, 1, 1, 1.0, 1.0, interpflag) ! < M_ij L_ij >
    call filter_q_3pt(denomc, cd, 1, 1, 1.0, 1.0, interpflag) ! < M_ij M_ij >

    if(myid==0 .and. debug>=2) write(*,*) "numc, denomc filtering successful" 

    do k = 1, nzp
        do j = 1, nypl
            do i = 1, nx
                ctot(i,j,k) = ctot(i,j,k) / cd(i,j,k)
                ! clip negative values
                if ( ctot(i,j,k)<0.0 .or. isnan( ctot(i,j,k) ) ) ctot(i,j,k) = 0.0 
                if ( ctot(i,j,k)>0.2**2 ) ctot(i,j,k) = 0.2**2 ! clip values above 0.2
            end do
        end do
    end do

    ! form interim quantity nu_sgs
    nu_sgs = ctot*delta*delta*snorm

    ! tau_ij = -2*C*delta^2*|S|*S_ij, but we let sij = delta^2*|S|*S_ij earlier
    do m = 1, nq
        do l = 1, nq
            tauij(:,:,:,l,m) = -2.0*nu_sgs*sij(:,:,:,l,m)
        end do
    end do

    if(myid==0 .and. debug>=1 .and. mod(istep,iomod)==0) then
        write(*,'(//,A,G14.6)') 'DynSmag: glob avg nu_sgs=',sum(nu_sgs)/(nx*nypl*nzp)
        write(*,'(/,A,/)') '------- DynSmag: nu_sgs(x,1,z) -------'
        write(*,fm1) (nu_sgs(1:nx,1,k),k=1,nzp)
    end if

    ! deallocate arrays
!     if (allocated(snorm)) 
    deallocate( snorm, snormhat, ctot, numc, denomc, cd, delta, nu_sgs, uihat, sij, sijhat, mij, lij, stat=AllocateStatus )
    if ( AllocateStatus/=0 ) write(*,*) "** dynamic_smagorinsky: deallocation error **"

end subroutine dynamic_smagorinsky


subroutine sigma_model(duidxj, tauij)
    ! computes the local eddy viscosity as defined by the sigma model by Nicoud et al
    use dim
    use grid, only: zpts, dx, dy
    ! to remove
    use mpicom, only: myid
    use flags, only: debug
    use runparam, only: iomod
    use time, only: istep, t
    use formatstrings
    ! i/o
    real, dimension(nxpp,nypl,nzp,3,3), intent(in) :: duidxj
    real, dimension(nxpp,nypl,nzp,3,3), intent(out) :: tauij
    ! local vars
    integer :: i,j,k,l,m, nq, AllocateStatus
    real :: pi, dz, delta, eps
    real, dimension(:,:,:), allocatable :: nu_sgs
    real :: invariant_1, invariant_2, invariant_3
    real :: alpha_1, alpha_2, alpha_3, pre_alpha_3
    real :: sigma_1, sigma_2, sigma_3
    real :: C_sigma, D_sigma
    real, dimension(3,3) :: gij, gji, G, G2 !,G_tran
    
    if (.not. allocated(nu_sgs)) allocate(nu_sgs(nxpp,nypl,nzp), stat=AllocateStatus)
    if (AllocateStatus/=0) write(*,*) "** sigma_model: allocation error **"

    nq = 3

    pi = 4.*atan(1.)

    eps = 1.E-10

    C_sigma = 1.35 ! parameter from Nicoud et al.

    tauij = 0.0
    nu_sgs = 0.0

    do k = 1, nzp
        do j = 1, nypl
            do i = 1, nx
                ! build gij and G = gki*gkj
                gij = duidxj(i,j,k,:,:)

                gji = transpose(gij)

                G = matmul(gji,gij)
                G2 = matmul(G,G)
!                 G_tran = transpose(G)
!                 G2 = matmul(G_tran,G)

                ! compute invariants
                invariant_1 = G(1,1) + G(2,2) + G(3,3)
                invariant_2 = 0.5*( invariant_1**2 - ( G2(1,1) + G2(2,2) + G2(3,3) ) )
                invariant_3 = G(1,1)*(G(2,2)*G(3,3) - G(2,3)*G(3,2)) - G(1,2)*(G(2,1)*G(3,3) - G(2,3)*G(3,1)) + G(1,3)*(G(2,1)*G(3,2) - G(2,2)*G(3,1))

                ! compute angles alpha_i
                alpha_1 = invariant_1**2/9. - invariant_2/3.
                alpha_2 = invariant_1**3/27. - invariant_1*invariant_2/6. + 0.5*invariant_3
                ! transfer angles to radians
!                 alpha_1 = alpha_1/180.*pi
!                 alpha_2 = alpha_2/180.*pi
                ! calc alpha_3
                if ( alpha_1<0.0 ) alpha_1 = abs(alpha_1) ! avoid sqrt of negative value
                pre_alpha_3 = alpha_2 / sqrt(alpha_1)**3.
                if ( pre_alpha_3>=1. .or. isnan(pre_alpha_3) ) then
                    alpha_3 = 0.0
                else if ( pre_alpha_3<=-1. ) then
                    alpha_3 = pi/3.
                else
                    alpha_3 = acos( pre_alpha_3 )/3.
                end if
!                 alpha_3 = alpha_3/180.*pi


                ! compute singular values sigma_i
                sigma_1 = sqrt( abs( invariant_1/3. + 2.*sqrt(alpha_1)*cos(alpha_3) ) )
                sigma_2 = sqrt( abs( invariant_1/3. - 2.*sqrt(alpha_1)*cos(pi/3.+alpha_3) ) )
                sigma_3 = sqrt( abs( invariant_1/3. - 2.*sqrt(alpha_1)*cos(pi/3.-alpha_3) ) )

                ! compute D_sigma
                D_sigma = sigma_3*( sigma_1 - sigma_2 )*( sigma_2 - sigma_3 )
                if (sigma_1/=0.) D_sigma = D_sigma/sigma_1**2
                

                ! compute delta
                if (k==1) then
                    dz = zpts(1) - zpts(2)
                else if (k==nzp) then
                    dz = zpts(nz) - zpts(nzp)
                else
                    dz = 0.5*(zpts(k-1)-zpts(k+1))
                end if

                delta = (dx*dy*dz)**(1./3.)      

!                 if(myid==0 .and. debug>=2 .and. mod(istep,iomod)==0) then
                if (D_sigma<-1.E-10 .or. isnan(D_sigma)) then
                    write(*,'(//,A,3(G12.4,","))') 'gij(1,:) = ', gij(1,:)
                    write(*,'(A,3(G12.4,","))') 'gij(2,:) = ', gij(2,:)
                    write(*,'(A,3(G12.4,","))') 'gij(3,:) = ', gij(3,:)
                    write(*,'(A,3(G12.4,","))') 'G(1,:) = ', G(1,:)
                    write(*,'(A,3(G12.4,","))') 'G(2,:) = ', G(2,:)
                    write(*,'(A,3(G12.4,","))') 'G(3,:) = ', G(3,:)
                    write(*,'(A)') 'invariant_1, invariant_2, invariant_3,'
                    write(*,'(3(G12.4,","))') invariant_1, invariant_2, invariant_3
                    write(*,'(A)') 'alpha_1, alpha_2, pre_alpha_3, alpha_3,'
                    write(*,'(4(G12.4,","))') alpha_1, alpha_2, pre_alpha_3, alpha_3
                    write(*,'(A)') 'sigma_1, sigma_2, sigma_3,'
                    write(*,'(3(G12.4,","))') sigma_1, sigma_2, sigma_3
                    write(*,'(A,G12.4)') 'D_sigma = ', D_sigma
                    write(*,'(A,G12.4)') 'dz = ', dz
                    write(*,'(A,G12.4)') 'delta = ', delta
                end if

                ! clip values below machine precision
                if (D_sigma<0.0 .or. isnan(D_sigma)) D_sigma = 0.0

                ! compute nu_sgs = (C*delta)^2 D_sigma
                nu_sgs(i,j,k) = -(C_sigma*delta)**2*D_sigma 

            end do
        end do
    end do

    if(myid==0 .and. debug>=1 .and. mod(istep,iomod)==0) then
        write(*,'(A,G14.6)') 'sigma-model: avg nu_sgs =',sum(nu_sgs)/(nx*nypl*nzp)
        write(*,'(A)') 'sigma-model: nu_sgs ='
        write(*,fm1) (nu_sgs(1:nx,1,k),k=1,nzp)
    end if

    ! tau_ij = -2*nu*S_ij
    do m = 1, nq
        do l = 1, nq
            tauij(:,:,:,l,m) = 2.0*nu_sgs*0.5*( duidxj(:,:,:,l,m) + duidxj(:,:,:,m,l) )
        end do
    end do

    if (allocated(nu_sgs)) deallocate(nu_sgs, stat=AllocateStatus)
    if (AllocateStatus/=0) write(*,*) "** sigma_model: deallocation error **"

end subroutine sigma_model


subroutine interscale_model(ui, tauij)
    !60-Galilean Invariant (Interscale Transfer Model): (Cancel: -T112-T122) uib_local*ujb_local-[uib*ujb]b+ub*up-[ub*up]b
    use dim
    use filters
    ! i/o
    real, dimension(nxpp,nypl,nzp,3), intent(in) :: ui
    real, dimension(nxpp,nypl,nzp,3, 3), intent(out) :: tauij
    ! local vars
    integer :: l, m, nq, AllocateStatus, interpflag
    real, dimension(:,:,:,:), allocatable :: uibar, uibar2
    real, dimension(:,:,:,:,:), allocatable ::  uibuj, uibujbar
    ! allocate variables
    nq = 3
    if (.not. allocated(uibar)) allocate( uibar(nxpp,nypl,nzp,nq), uibar2(nxpp,nypl,nzp,nq), uibuj(nxpp,nypl,nzp,nq,nq), uibujbar(nxpp,nypl,nzp,nq,nq), stat=AllocateStatus )
    if ( AllocateStatus/=0 ) write(*,*) "** interscale_model: allocation error **"

    interpflag = 1

    ! obtain [ui] and [[ui]]
    call filter_q_3pt( ui, uibar, nq, 1, a, b, interpflag )
    call filter_q_3pt( uibar, uibar2, nq, 1, a, b, interpflag )

    do m = 1, nq
        do l = 1, nq
            uibuj(:,:,:,l,m) = uibar(:,:,:,l)*ui(:,:,:,m)
        end do
    end do

    ! obtain [[ui]uj]
    call filter_q_3pt( uibuj, uibujbar, nq, nq, a, b, interpflag )

    ! tau_ij = [[u_i]u_j] - [[u_i]][[u_j]] - [u_i]( u_j - [u_j] )
    do m = 1, nq
        do l = 1, nq
            tauij(:,:,:,l,m) = uibujbar(:,:,:,l,m) - uibar2(:,:,:,l)*uibar2(:,:,:,m) - uibar(:,:,:,l)*(ui(:,:,:,m) - uibar(:,:,:,m))
        end do
    end do

    if (allocated(uibar)) deallocate( uibar, uibar2, uibuj, uibujbar, stat=AllocateStatus )
    if ( AllocateStatus/=0 ) write(*,*) "** interscale_model: deallocation error **"

end subroutine interscale_model

end module sgs_models