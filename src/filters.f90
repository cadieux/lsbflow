module filters

implicit none
public
save

contains

subroutine filterall3pt(uf,vf,wf,tempf,pf)
    ! filter all primary qties
    use dim
    use runparam, only: iomod
    use time, only: isum
    use sgsinput
    use modhorfft
    ! to remove
    use flags, only: debug
    use mpicom, only: myid
    use formatstrings
    ! i/o
    real, dimension(nxpl,ny,nzpl), intent(inout) :: uf, vf, wf, tempf, pf
    ! local vars
    integer :: err, nvel!, l, k
    real, allocatable, dimension(:,:,:,:) :: q, qfilt
    nvel = 5
    if (.not. allocated(q)) allocate(q(nxpp,nypl,nzp,nvel), qfilt(nxpp,nypl,nzp,nvel), stat=err)
    if (err /= 0) write(*,*) "filterall3pt: Allocation request denied" 
    ! transfer all qties to physical space
    q = 0.0
    call horfft(q(:,:,:,1),uf,1)
    call horfft(q(:,:,:,2),vf,1)
    call horfft(q(:,:,:,3),wf,1)
    call horfft(q(:,:,:,4),pf,1)
    if(any(abs(tempf)>1.0E-10)) call horfft(q(:,:,:,5),tempf,1)

    call modecutFC(q,qfilt,nvel,a,b)

    call horfft(qfilt(:,:,:,1),uf,-1)
    call horfft(qfilt(:,:,:,2),vf,-1)
    call horfft(qfilt(:,:,:,3),wf,-1)
    call horfft(qfilt(:,:,:,4),pf,-1)
    if(any(abs(tempf)>1.0E-10)) call horfft(qfilt(:,:,:,5),tempf,-1)
!     call horfft(dynp_filt,pf,-1)
!     if (debug>=1 .and. myid==0) write(*,*) 'filter: transfer back to fourier space completed successfully'
    ! compute energy removed
!     if (mod(isum,iomod)==0) then
!         ekin_rem = ( u**2 + v**2 + w**2 ) - ( u_filt**2 + v_filt**2 + w_filt**2 )
!         ekin_rem(:,1,:) = sum(ekin_rem,2)/real(nypl) ! avg over y
    
!         if (debug>=1 .and. myid==0) then 
!             write(*,'(//,A)') '<E_kin>_y removed by filtering operation'
!             write(*,fm1) (ekin_rem(1:nx,1,k),k=1,nzp)
!         end if
!      end if

    if (allocated(q)) deallocate(q, qfilt, stat=err)
    if (err /= 0) write(*,*) "filterall3pt: Deallocation request denied"

end subroutine filterall3pt


subroutine tns_automatic_filtering(uf,vf,wf,tempf,pf)
    ! compute energy difference between two filtering levels
    ! <I(h)/I(Delta)>(z) = < [(u - ub1)**2 + (v - vb1)**2 + (w - wb1)**2]/[(u - ub2)**2 + (v - vb2)**2 + (w - wb2)**2] >
    ! <I(h)/I(Delta)> = int_z <I(h)/I(Delta)>(z) from z=0 to z=H
    use dim
    use runparam, only: iomod
    use time, only: isum, t
    use sgsinput
    use modhorfft
    use grid, only: zpts, dz
    ! to remove
    use flags, only: debug
    use mpicom
    use formatstrings
    ! i/o
    real, dimension(nxpl,ny,nzpl), intent(inout) :: uf, vf, wf, tempf, pf
    ! local vars
    integer :: err, k, nq, nvel, l
    real, allocatable, dimension(:,:,:,:) :: q, qfilt, qfilt_2kc
    real, allocatable, dimension(:,:) :: energy
    real :: Ih, Idelta, ekin_ratio_limit
    nq = 4
    nvel = 3
    if (.not. allocated(q)) allocate(q(nxpp,nypl,nzp,nq), qfilt(nxpp,nypl,nzp,nq), qfilt_2kc(nxpp,nypl,nzp,nq), energy(nzp,4), stat=err)
    if (err /= 0) write(*,*) "tns_automatic_filtering: Allocation request denied"

    ! set appropriate energy ratio maximum
    ekin_ratio_limit = 0.007 ! as per Tantikul, Domaradzki. J turb, 2011

    ! transfer all qties to physical space
    q = 0.0
    call horfft(q(:,:,:,1),uf,1)
    call horfft(q(:,:,:,2),vf,1)
    call horfft(q(:,:,:,3),wf,1)
    call horfft(q(:,:,:,4),pf,1)
!     if(any(abs(tempf)>1.0E-10)) call horfft(q(:,:,:,5),tempf,1)

    call modecutFC(q,qfilt_2kc,nq,1.0,6.0) ! 1/8 3/4 1/8 sharper filter
!     call modecutFC(q,qfilt_2kc,nq,1.0,9.0) ! 1/8 3/4 1/8 sharper filter

    call modecutFC(q,qfilt,nq,1.0,2.0) ! 1/4 1/2 1/4 trapezoidal filter
!     call modecutFC(q,qfilt,nq,1.0,4.0) ! 1/6 2/3 1/6 simpson's filter
!     call modecutFC(q,qfilt,nq,1.0,3.0) ! 1/5 3/5 1/5 in between filter

    ! compute energy difference averaged over planes
    energy = 0.0
    do k = 1, nzp
!         energy(k,:) = 0.0
        do l = 1, nvel
            energy(k,1) = energy(k,1) + sum( (q(:,:,k,l) - qfilt(:,:,k,l))**2 )/(ny*nx)
            energy(k,2) = energy(k,2) + sum( (q(:,:,k,l) - qfilt_2kc(:,:,k,l))**2 )/(ny*nx)         
        end do
    end do
    ! add up each section in y-direction
    call mpi_allreduce(energy(:,1:2),energy(:,3:4),2*nzp,mpi_double_precision,MPI_SUM,comm,ierr)
    ! compute integral in vertical
    Ih = 0.0
    Idelta = 0.0
    do k = nz, 1, -1
        dz = (zpts(k) - zpts(k+1))
        Ih = Ih + 0.5*(energy(k,4) + energy(k+1,4))*dz
        Idelta = Idelta + 0.5*(energy(k,3) + energy(k+1,3))*dz
    end do

    if (debug>=2 .and. myid==0) then
        write(*,'(A,G12.4)') "TNS: <I(h)>          =",Ih
        write(*,'(A,G12.4)') "TNS: <I(delta)>      =",Idelta
        write(*,'(A,F12.4)') "TNS: <I(h)/I(delta)> =",Ih/Idelta
    end if


    if (Ih/Idelta > ekin_ratio_limit) then
        call horfft(qfilt_2kc(:,:,:,1),uf,-1)
        call horfft(qfilt_2kc(:,:,:,2),vf,-1)
        call horfft(qfilt_2kc(:,:,:,3),wf,-1)
        call horfft(qfilt_2kc(:,:,:,4),pf,-1)
!         call modecutFC(q,qfilt,nq,1.0,4.0) ! use sharper filter in between
!         call horfft(qfilt(:,:,:,1),uf,-1)
!         call horfft(qfilt(:,:,:,2),vf,-1)
!         call horfft(qfilt(:,:,:,3),wf,-1)
!         call horfft(qfilt(:,:,:,4),pf,-1)
!         if(any(abs(tempf)>1.0E-10)) call horfft(qfilt(:,:,:,5),tempf,-1)
        if (myid == 0 .and. debug>=1) then
            write(*,'(A,G12.4)') "TNS: <I(h)>          =",Ih
            write(*,'(A,G12.4)') "TNS: <I(delta)>      =",Idelta
            write(*,'(A,F12.4)') "TNS: <I(h)/I(delta)> =",Ih/Idelta
            write(*,'(A,F14.6,A,I14.6,/)') "TNS: filtering applied at time=",t,", it # =", isum
        end if
!     else
!         call horfft(qfilt_2kc(:,:,:,1),uf,-1)
!         call horfft(qfilt_2kc(:,:,:,2),vf,-1)
!         call horfft(qfilt_2kc(:,:,:,3),wf,-1)
!         call horfft(qfilt_2kc(:,:,:,4),pf,-1)
    end if


    if (allocated(q)) deallocate(q, qfilt, qfilt_2kc, energy, stat=err)
    if (err /= 0) write(*,*) "tns_automatic_filtering: Deallocation request denied"

end subroutine tns_automatic_filtering


subroutine filterqspec(qf,nq)
    ! filter all primary qties
    use dim
    use sgsinput, only: kc, filt_order, filt_amp
    ! i/o
    integer, intent(in) :: nq
    real, dimension(nxpl,ny,nzpl,nq), intent(inout) :: qf
    ! local vars
    integer :: err, filtery, iq
    real, allocatable, dimension(:,:,:,:) :: qfilt
!     real :: kc, p, amp

    if (.not. allocated(qfilt)) allocate(qfilt(nxpl,ny,nzpl,nq),stat=err)
    if (err /= 0) write(*,*) "filterallspec: Allocation request denied"

    filtery = 1
    ! kc: fraction of max wavenumber to set as cutoff (usually set to half = twice grid size)
    ! filt_order: order of exponential filter
    ! filt_amp:  amplitude of filter
    do iq = 1, nq
        call specfilter_exp_xy(qf(:,:,:,iq),qfilt(:,:,:,iq),kc,filt_order,filt_amp,filtery)
    end do

    qf = qfilt

    if (allocated(qfilt)) deallocate(qfilt, stat=err)
    if (err /= 0) write(*,*) "filterallspec: Deallocation request denied"

end subroutine filterqspec



subroutine filterallspec(uf,vf,wf,cutoff)
    ! filter all primary qties
    use dim
    use sgsinput, only: kc, filt_order, filt_amp
    ! i/o
    integer, intent(in) :: cutoff
    real, dimension(nxpl,ny,nzpl), intent(inout) :: uf,vf,wf
    ! local vars
    integer :: err, filtery
    real, allocatable, dimension(:,:,:) :: qfilt
!     real :: kc, p, amp

    if (.not. allocated(qfilt)) allocate(qfilt(nxpl,ny,nzpl),stat=err)
    if (err /= 0) write(*,*) "filterallspec: Allocation request denied"

    filtery = 1
    ! kc: fraction of max wavenumber to set as cutoff (usually set to half = twice grid size)
    ! filt_order: order of exponential filter
    ! filt_amp:  amplitude of filter
    
    if (cutoff==1) then
        call specfilter_cutoff_xy(uf,qfilt,kc,filtery)
        uf = qfilt
        qfilt = 0.0
        call specfilter_cutoff_xy(vf,qfilt,kc,filtery)
        wf = qfilt
        qfilt=0.0
        call specfilter_cutoff_xy(wf,qfilt,kc,filtery)
        wf = qfilt
    else
        call specfilter_exp_xy(uf,qfilt,kc,filt_order,filt_amp,filtery)
        uf = qfilt
        qfilt = 0.0
        call specfilter_exp_xy(vf,qfilt,kc,filt_order,filt_amp,filtery)
        vf = qfilt
        qfilt = 0.0
        call specfilter_exp_xy(wf,qfilt,kc,filt_order,filt_amp,filtery)
        wf = qfilt
    end if

    if (allocated(qfilt)) deallocate(qfilt, stat=err)
    if (err /= 0) write(*,*) "filterallspec: Deallocation request denied"

end subroutine filterallspec


subroutine specfilter_cutoff_xy(uf,uff,kc,filtery)
    ! Fourier cut-off filter
    use dim
    use grid, only: wavx, wavy, alpha, beta
    ! i/o
    real, dimension(nxpl,ny,nzpl), intent(in) :: uf
    real, intent(in) :: kc ! fraction of max wav number above which to filter
    real, dimension(nxpl,ny,nzpl), intent(out) :: uff
    integer, intent(in) :: filtery
    ! local vars
    integer :: i,j!,k
    real :: wavxmax, wavymax, wavxcut, wavycut
    wavxmax = nxh*alpha
    wavymax = nyh*beta
    wavxcut = wavxmax*kc
    wavycut = wavymax*kc
    
    do i=1,nxpl
        ! remove contribution above cutoff
        if ( abs(wavx(i)) > wavxcut ) then
            uff(i,:,:) = 0.0
        else
            uff(i,:,:) = uf(i,:,:)
        end if
    end do
    if (filtery==1) then
        do j=1,ny
            ! remove contribution above cutoff
            if ( abs(wavy(j)) > wavycut ) then
                uff(:,j,:) = 0.0
            else
                uff(:,j,:) = uf(:,j,:)
            end if
        end do
    end if

end subroutine specfilter_cutoff_xy    


subroutine specfilter_exp_xy(uf,uff,kc,p,amp,filtery)
    ! Fourier exponential filter
    use dim
    use grid, only: wavx, wavy, alpha, beta, wavxx, wavyy
    ! i/o
    real, dimension(nxpl,ny,nzpl), intent(in) :: uf
    real, intent(in) :: kc, p, amp
    real, dimension(nxpl,ny,nzpl), intent(out) :: uff
    integer, intent(in) :: filtery
    ! local vars
    integer :: i,j,k
    real :: fr, fi, fmag, fmagf, cosc, sinc, phase, twopi, wavxmax, wavymax,ak,akmax
    ! parameters
    twopi = 8.*atan(1.)
    wavxmax = nxh*alpha
    wavymax = nyh*beta
    ! use akmax to combine filtering
    akmax = wavxmax
    if (filtery==1) akmax = sqrt( wavxmax**2 + wavymax**2 )
    ! express complex value of fourier coeffs in terms of fmag*(cosc + i*sinc)
    do k=1,nzp
        do j=1,ny
            do i=1,nxpl,2
                fr = uf(i,j,k)
                fi = uf(i+1,j,k)
                fmag = sqrt( fr*fr + fi*fi )
                if (fmag /= 0.) then ! only filter non-zero amp modes
                    cosc = fr/fmag
                    sinc = fi/fmag
                    ! get phase
                    if (sinc >= 0.) phase = acos(cosc)
                    if (sinc <  0.) phase = -acos(cosc)
                    ! filter according to wavenumber magnitude in x
                    ak = wavx(i)
                    if (filtery==1) ak = sqrt( wavxx(i) + wavyy(j) )
                    fmagf = fmag * specfilter_fourier2d(ak,akmax,kc,p,amp)
!                     fmagf = fmag * specfilter_fourier2d(wavx(i),wavxmax,kc,p,amp)
                    ! filter according to wavenumber magnitude in y
!                     if (filtery==1) then
!                         fmagf = fmagf * specfilter_fourier2d(wavy(i),wavymax,kc,p,amp)
!                     end if
                    ! store back into a + bi form
                    fr = fmagf*cos(phase)
                    fi = fmagf*sin(phase)
                else
                    fr = 0.0
                    fi = 0.0
                end if
                uff(i,j,k) = fr
                uff(i+1,j,k) = fi
            end do
        end do
    end do

end subroutine specfilter_exp_xy


real function specfilter_fourier2d(ak,akmax,kc,p,amp)
    !-(PD:4/23/02). Exponential filter in Fourier wave# space.
    !-Arguments:
    !-a) ak = Fourier wavenumber under consideration
    !-b) akmax = Maximum Fourier wavenumber in 2-D Fourier expansion
    !-c) kc = fraction of akmax above which filter is applied
    !-d) p = filter order
    !-e) amp = filter amplitude
    ! i/o
    real, intent(in) :: ak, akmax, kc, p, amp
    ! local vars
    real :: akc, theta
    akc = akmax*kc
     
    if (abs(ak).ge.akc) then
        theta = (abs(ak)-akc)/(akmax-akc)
!    theta = 0.0
        specfilter_fourier2d = exp(-amp*theta**(2*p))
    else
        specfilter_fourier2d = 1.0
    endif

end function specfilter_fourier2d


subroutine filter_q_3pt(q,q_filt,n4,n5,a,b,interpflag)
    ! filter a number of quantities held in tensor q in x, y, z
    use dim
    use runparam, only: mapping
    use grid, only: zpts
    ! i/o
    integer, intent(in) :: n4,n5,interpflag
    real, intent(in) :: a,b
    real, dimension(nxpp,nypl,nzp,n4,n5), intent(in) :: q
    real, dimension(nxpp,nypl,nzp,n4,n5), intent(out) :: q_filt
    ! local vars
    integer :: i, k, err
    real :: dltm, dltp
    real, dimension(:,:,:,:,:), allocatable :: q_hat
    real, dimension(:), allocatable :: az, bz, cz
    ! first allocate intermediate array qhat
    if(.not. allocated(q_hat)) allocate( q_hat(nxpp,nypl,nzp,n4,n5), stat=err )
    if ( err/=0 ) write(*,*) "** filter_q_3pt allocation error **"

    ! filter all qties in x-dir
    q_hat(1,:,:,:,:) = a*q(nx,:,:,:,:) + b*q(1,:,:,:,:) + a*q(2,:,:,:,:) ! left boundary
    do i=2,nx-1
        q_hat(i,:,:,:,:) = a*q(i-1,:,:,:,:) + b*q(i,:,:,:,:) + a*q(i+1,:,:,:,:) ! interior
    end do
    q_hat(nx,:,:,:,:) = a*q(nx-1,:,:,:,:) + b*q(nx,:,:,:,:) + a*q(1,:,:,:,:) !rightboundary
    q_hat = q_hat/( 2.*a + b )

    ! filtering all qties in y-dir implies communication
    call filtery_q_3pt(q_hat,q_filt,n4,n5,a,b)

    q_hat = q_filt ! copy xy-filtered fields back into intermediate array

    ! filter all qties in z-dir
    if (.not. allocated(az)) allocate(az(nzpl), bz(nzpl), cz(nzpl), stat=err)
    if (err /= 0) write(*,*) "filter_q_3pt: Allocation request denied"
    if ( interpflag==1 ) then
        if ( mapping==1 .or. mapping==2 ) then
            ! use only one sided interpolating filter
            do k=2,nzp-1 ! do not include last point z(1) = +infty
                dltm=zpts(k)-zpts(k-1)
                dltp=zpts(k+1)-zpts(k)
                az(k)=a*2.0*dltm/(dltp+dltm)
                bz(k)=a*2.0*(1.-dltm/dltp)+b
                cz(k)=a*2.0*dltm*dltm/(dltp*(dltp+dltm))
            end do
            if ( mapping==1 ) then
                az(2)=0.0
                bz(2)=(2.0*a+b)
                cz(2)=0.0
            end if
        else
            !Interpolating filter factors from bottom wall to mid point
            do k=2,(nzp+1)/2-1
                dltm=zpts(k)-zpts(k-1)
                dltp=zpts(k+1)-zpts(k)
                az(k)=a*2.0*dltm/(dltp+dltm)
                bz(k)=a*2.0*(1.-dltm/dltp)+b
                cz(k)=a*2.0*dltm*dltm/(dltp*(dltp+dltm))
            end do

            !Interpolating filter factors from mid-point to top wall
            do k=(nzp+1)/2+1,nzp-1
                dltm=zpts(k)-zpts(k-1)
                dltp=zpts(k+1)-zpts(k)
                az(k)=a*2.0*dltp*dltp/(dltm*(dltm+dltp))
                bz(k)=a*2.0*(1.-dltp/dltm)+b
                cz(k)=a*2.0*dltp/(dltm+dltp)
            end do

            az((nzpl+1)/2)=a
            bz((nzpl+1)/2)=b
            cz((nzpl+1)/2)=a
        end if

        az(1)=0.0
        bz(1)=(2.0*a+b)
        cz(1)=0.0

        az(nzpl)=0.0
        bz(nzpl)=(2.0*a+b)
        cz(nzpl)=0.0
    else
        az = a
        bz = b
        cz = a
    end if

    do k=2,nzpl-1
        q_filt(:,:,k,:,:)=az(k)*q_hat(:,:,k-1,:,:)+bz(k)*q_hat(:,:,k,:,:)+cz(k)*q_hat(:,:,k+1,:,:)
        q_filt(:,:,k,:,:)=q_filt(:,:,k,:,:)/(2.0*a+b)
    end do

    if (allocated(q_hat)) deallocate( q_hat, az, bz, cz, stat=err )
    if ( err/=0 ) write(*,*) "** filter_q_3pt deallocation error **"

end subroutine filter_q_3pt


subroutine filtery_q_3pt(q,q_filt,n4,n5,a,b)
    ! Filter quantities in the y direction using a 3pt rule
    ! //////////// needs work and testing //////////////
    use dim
    use mpicom
    use flags, only: debug
    ! i/o
    integer, intent(in) :: n4, n5
    real, intent(in) :: a, b
    real, dimension(nxpp,nypl,nzpl,n4,n5), intent(in) :: q
    real, dimension(nxpp,nypl,nzpl,n4,n5), intent(out) :: q_filt
    ! local vars
    real, dimension(:,:,:,:),allocatable :: buf1, buf2
    real, dimension(:,:,:,:,:),allocatable :: q_hat
!     real, dimension(:),allocatable :: ubot1, utop1, utop2, ubot2
    integer :: j, bufsize, err
    integer :: tagin, tagout, source1, source2, dest1, dest2, request1, request2, status(mpi_status_size)
    ! mpi-io error checking
!     character*(MPI_MAX_ERROR_STRING) :: error_string
!     integer ::  ierror, resultlen

    if(.not. allocated(q_hat)) allocate(q_hat(nxpp,0:nyplp,nzpl,n4,n5), buf1(nx,nzp,n4,n5), buf2(nx,nzp,n4,n5), stat=err)
    if (err /= 0) print *, "filtery_q_3pt: Allocation request denied"

    bufsize = nx*nzp*n4*n5

    q_filt = 0.0 ! ensure nxp, nxpp values are zero
    q_hat(:,1:nypl,:,:,:)=q

    ! filter interior values
    if (nypl>2) then
        do j=2,nypl-1
            q_filt(:,j,:,:,:)=a*q_hat(:,j-1,:,:,:)+b*q_hat(:,j,:,:,:)+a*q_hat(:,j+1,:,:,:)
            q_filt(:,j,:,:,:)=q_filt(:,j,:,:,:)/(2.00*a+b)
        end do
    end if

    if ( debug>=2 .and. myid==0 ) write(*,*) "filtery_q_3pt y interior filter completed"

    ! communicate boundary values only if there are more than 1 proc
    if ( nproch > 1 ) then   
        if ( mod(myid,2)==0 .or. myid==0 ) then ! even processors send data
            dest1 = myid + 1 ! destination
            if ( dest1 > nproch - 1 ) dest1 = 0 ! procmax sends to proc0
            tagout = 2000 + myid
            
            buf1 = q_hat(1:nx,nypl,1:nzp,1:n4,1:n5)
            if ( debug>=2 ) write(*,*) "filtery_q_3pt just before SEND 1 from",myid,"to",dest1
            call mpi_send(buf1,bufsize,mpi_double_precision,dest1,tagout,comm,ierr)
            if ( debug>=2 .and. ierr==0 ) write(*,*) "filtery_q_3pt EVEN proc SEND 1 completed for proc = ",myid

            ! send beginning boundary data to processor just before itself
            dest2 = myid - 1
            if ( dest2 < 0 ) dest2 = nproch - 1 ! proc0 sends to procmax
            ! set some tags to keep track of things
            tagout = 1000 + myid
            
            buf2 = q_hat(1:nx,1,1:nzp,1:n4,1:n5)
!             buf2 = q_hat(1:nx,nypl,1:nzp,1:n4,1:n5)
            if ( debug>=2  ) write(*,*) "filtery_q_3pt just before EVEN SEND 2 from",myid,"to",dest2
            call mpi_send(buf2,bufsize,mpi_double_precision,dest2,tagout,comm,ierr)
            if ( debug>=2 .and. ierr==0 ) write(*,*) "filtery_q_3pt EVEN proc SEND 2 completed for proc = ",myid
            
            ! receive some data later
            ! set source of data current proc needs
            source1 = myid - 1
            if ( source1 < 0 ) source1 = nproch - 1 ! proc0 sends to procmax
            ! tags are different from those we send, because other proc's have different myid
            tagin = MPI_ANY_TAG !2000 + source1
            if ( debug>=2  ) write(*,*) "filtery_q_3pt just before RECV 1 from",source1,"to",myid
            call mpi_irecv(buf1,bufsize,mpi_double_precision,source1,tagin,comm,request1,ierr)
            call mpi_wait(request1,status,ierr)
           
            q_hat(1:nx,0,1:nzp,1:n4,1:n5) = buf1
            if ( debug>=2 .and. ierr==0 ) write(*,*) "filtery_q_3pt EVEN proc RECV 1 completed for proc = ",myid
            ! set source of data current proc needs
            source2 = myid + 1
            if ( source2 > nproch - 1 ) source2 = 0 ! procmax sends to proc0
            tagin = MPI_ANY_TAG !source2 
            ! prepare to receive one piece of data
            if ( debug>=2  ) write(*,*) "filtery_q_3pt just before RECV 2 from",source2,"to",myid
            call mpi_irecv(buf2,bufsize,mpi_double_precision,source2,tagin,comm,request2,ierr)
            call mpi_wait(request2,status,ierr)
            
            q_hat(1:nx,nyplp,1:nzp,1:n4,1:n5) = buf2
            if ( debug>=2 .and. ierr==0 ) write(*,*) "filtery_q_3pt EVEN proc RECV 2 completed for proc = ",myid


        else ! for odd processors
            ! set source of data current proc needs
            source1 = myid - 1
            if ( source1 < 0 ) source1 = nproch - 1 ! proc0 sends to procmax
            ! tags are different from those we send, because other proc's have different myid
            tagin = MPI_ANY_TAG !2000 + source1
            if ( debug>=2  ) write(*,*) "filtery_q_3pt just before RECV 1 from",source1,"to",myid
            call mpi_irecv(buf1,bufsize,mpi_double_precision,source1,tagin,comm,request1,ierr)
            call mpi_wait(request1,status,ierr)
            q_hat(1:nx,0,1:nzp,1:n4,1:n5) = buf1
            if ( debug>=2 .and. ierr==0 ) write(*,*) "filtery_q_3pt ODD proc RECV 1 completed for proc = ",myid
            ! set source of data current proc needs
            source2 = myid + 1
            if ( source2 > nproch - 1 ) source2 = 0 ! procmax sends to proc0
            tagin = MPI_ANY_TAG !source2 
            ! prepare to receive one piece of data
            if ( debug>=2  ) write(*,*) "filtery_q_3pt just before RECV 2 from",source2,"to",myid
            call mpi_irecv(buf2,bufsize,mpi_double_precision,source2,tagin,comm,request2,ierr)
            call mpi_wait(request2,status,ierr)
            q_hat(1:nx,nyplp,1:nzp,1:n4,1:n5) = buf2
            if ( debug>=2 .and. ierr==0 ) write(*,*) "filtery_q_3pt ODD proc RECV 2 completed for proc = ",myid

            ! now send data
            dest1 = myid + 1 ! destination
            if ( dest1 > nproch - 1 ) dest1 = 0 ! procmax sends to proc0
            tagout = 2000 + myid
            buf1 = q_hat(1:nx,nypl,1:nzp,1:n4,1:n5)
            if ( debug>=2  ) write(*,*) "filtery_q_3pt just before SEND 1 from",myid,"to",dest1
            call mpi_send(buf1,bufsize,mpi_double_precision,dest1,tagout,comm,ierr)
            if ( debug>=2 .and. ierr==0 ) write(*,*) "filtery_q_3pt ODD proc SEND 1 completed for proc = ",myid
            ! send beginning boundary data to processor just before itself
            dest2 = myid - 1
            if ( dest2 < 0 ) dest1 = nproch - 1 ! proc0 sends to procmax
            ! set some tags to keep track of things
            tagout = myid
            buf2 = q_hat(1:nx,1,1:nzp,1:n4,1:n5)
            if ( debug>=2  ) write(*,*) "filtery_q_3pt just before ODD SEND 2 from",myid,"to",dest2
            call mpi_send(buf2,bufsize,mpi_double_precision,dest2,tagout,comm,ierr)
            if ( debug>=2 .and. ierr==0 ) write(*,*) "filtery_q_3pt ODD proc SEND 2 completed for proc = ",myid

        end if
    
    else ! for one proc
        q_hat(:,0,:,:,:) = q(:,ny,:,:,:)
        q_hat(:,nyplp,:,:,:) = q(:,1,:,:,:)
    end if
        
    ! compute filtered boundary values (j=1, and j=nypl)
    if (nypl>2) then
        do j=1,nypl,nypl-1
            q_filt(:,j,:,:,:)=a*q_hat(:,j-1,:,:,:)+b*q_hat(:,j,:,:,:)+a*q_hat(:,j+1,:,:,:)
            q_filt(:,j,:,:,:)=q_filt(:,j,:,:,:)/(2.00*a+b)
        enddo
    else if (nypl<=2) then
        do j=1,nypl
            q_filt(:,j,:,:,:)=a*q_hat(:,j-1,:,:,:)+b*q_hat(:,j,:,:,:)+a*q_hat(:,j+1,:,:,:)
            q_filt(:,j,:,:,:)=q_filt(:,j,:,:,:)/(2.00*a+b)
        enddo
    end if


    if (allocated(q_hat)) deallocate(q_hat, buf1, buf2, stat=err)
    if (err /= 0) print *, "filtery_q_3pt: Deallocation request denied"

end subroutine filtery_q_3pt



subroutine modecutFC(q,qb,nq,a,b)
    ! This routine obtain large component of the full velocity field by special
    ! filtering.
    use dim
    ! i/o
    integer, intent(in) :: nq
    real, intent(in) :: a, b
    real, dimension(nxpp,nypl,nzpl,nq), intent(in) :: q
    real, dimension(nxpp,nypl,nzpl,nq), intent(out) :: qb
    ! local vars
    integer :: AllocateStatus, l, nfilt, nvel
    real, dimension(:,:,:,:,:), allocatable :: q_filt
    real, dimension(:), allocatable :: coef

    nvel = nq
    nfilt = 6
    allocate(q_filt(nxpp,nypl,nzpl,nvel,nfilt), coef(nfilt), stat=AllocateStatus)
    if(AllocateStatus/=0) then
        write(*,*)" ** write(*,*) - Not Enough Memory - modecutFC - "
    endif

    ! simple procedure - filter 6 times, keep each quantity in q_filt(:,:,:,l)
    qb = 0.0
    coef = (/ 6., -15., 20., -15., 6., -1.  /)
    l = 1
    call filter_q_3pt(q,q_filt(:,:,:,:,l),nvel,1,a,b,1)
    qb = qb + coef(l)*q_filt(:,:,:,:,l)
    do l=2,nfilt
        call filter_q_3pt(q_filt(:,:,:,:,l-1),q_filt(:,:,:,:,l),nvel,1,a,b,1)
        qb = qb + coef(l)*q_filt(:,:,:,:,l)
    end do

    deallocate(q_filt, coef, stat=AllocateStatus)
    if(AllocateStatus/=0) then
        write(*,*)" ** write(*,*) - Error Deallocating- modecutFC - "
    endif  

end subroutine modecutFC


end module filters
