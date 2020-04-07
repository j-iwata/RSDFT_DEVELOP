module rotorb_mbp_module

  use parallel_module
  use array_bound_module, only: MBZ_0,MBZ_1,MSP_0,MSP_1
  use wf_module, only: unk
  use cpmd_variables, only: gam,scr,wrk,tau,sig,gamn,psi_v,psi_n &
                           ,mstocck,dt,disp_switch,MBC,MBT,MB_0_CPMD,MB_1_CPMD
  use overlap_cpmd_mbp_module, only: overlap2_mbp, overlap4_mbp, overlap5_mbp
  use rsdft_mpi_module, only: rsdft_bcast, rsdft_allreduce_sum
  use scalapack_module, only: MBSIZE, NBSIZE, NPROW, NPCOL, usermap, DESCA, DESCB, DESCZ
  use rsdft_mpi_module, only: rsdft_bcast
! use rotorb_sub_module, only: rotorb_sub
  use watch_module

  implicit none
! include 'mpif.h'

  PRIVATE
  PUBLIC :: rotorb_mbp
  PUBLIC :: rotorb2_mbp

  integer,allocatable :: ir_i(:),id_i(:)
  integer :: ls,le,li,n1,n2,ML0,MB0
  real(8),allocatable :: psi_tmp(:,:)
! logical,parameter :: idebug=.true.
  logical,parameter :: idebug=.false.

CONTAINS

  subroutine rotorb_mbp
    implicit none
    integer,parameter :: maxit=100
    integer :: i,j,k,s,n,it
    real(8),parameter :: eps=1.0d-8
    real(8),parameter :: hf = 0.5d0
    real(8),parameter :: hm =-0.5d0
    real(8),parameter :: on = 1.0d0
    real(8),parameter :: zr = 0.0d0
    real(8) :: error,tmp1,tttt(2,16),ttmp(2)
    integer :: ii0, jj0

    !call write_border( 1, "rotorb(start)" )
    !call watchb( ttmp, barrier="on" ); tttt=0.0d0

    n1  = idisp(myrank)+1
    n2  = idisp(myrank)+ircnt(myrank)
    ML0 = ircnt(myrank)
    MB0 = MB_1_CPMD - MB_0_CPMD + 1
    allocate( psi_tmp(n1:n2,MB_0_CPMD:MB_1_CPMD) ) ; psi_tmp=0.0d0

    !call watchb( ttmp, tttt(:,1), barrier="on" )

    do s=MSP_0,MSP_1
    do k=MBZ_0,MBZ_1

       !call watchb( ttmp, barrier="on" )

       MBT=mstocck(k,s)

       call overlap2_mbp(s,k,ii0,jj0) ! sig

       !call watchb( ttmp, tttt(:,2), barrier="on" )

       call overlap4_mbp(s,k,ii0,jj0) ! tau

       !call watchb( ttmp, tttt(:,3), barrier="on" )

       gam(1:ii0,1:jj0) = sig(1:ii0,1:jj0)*hf      ! (I-A)/2

       !call watchb( ttmp, tttt(:,4), barrier="on" )

!      if ( ls <= le ) then
!         call DGEMM &
!              ('N','N',MBT,li,MBT,hf,tau,MBC,tau(1,ls),MBC,hf,sig(1,ls),MBC)
!      end if
       call pdgemm('n','n',MBT,MBT,MBT,hf,tau,1,1,DESCA,tau,1,1,DESCA,hf,sig,1,1,DESCA)

       !call watchb( ttmp, tttt(:,5), barrier="on" )

       do it=1,maxit

          !call watchb( ttmp, barrier="on" )

          gamn(1:ii0,1:jj0) = sig(1:ii0,1:jj0)
          wrk(1:ii0,1:jj0)  = tau(1:ii0,1:jj0) - gam(1:ii0,1:jj0)

          !call watchb( ttmp, tttt(:,6), barrier="on" )
          !call watchb( ttmp, tttt(:,7), barrier="on" )

          call pdgemm('n','n',MBT,MBT,MBT,hm,wrk,1,1,DESCA,wrk,1,1,DESCA,on,gamn,1,1,DESCA)

          !call watchb( ttmp, tttt(:,8), barrier="on" )

          error=0.0d0
          do j=1,jj0
          do i=1,ii0
             error = error + ( gamn(i,j)-gam(i,j) )**2
          end do
          end do

          call rsdft_allreduce_sum( error, MPI_COMM_WORLD )
          error=sqrt(error)/MBT
          if ( myrank == 0 .and. error > on ) write(*,*) k,s,error

          !call watchb( ttmp, tttt(:,9), barrier="on" )

          gam(1:ii0,1:jj0) = gamn(1:ii0,1:jj0)

          !call watchb( ttmp, tttt(:,10), barrier="on" )

          if ( error <= eps ) exit
          if ( it == maxit .and. myrank == 0 ) then
             write(*,*) "WARNING: iteration was not converged in rotorb"
          end if

       end do ! it

       !call watchb( ttmp, barrier="on" )
       !call watchb( ttmp, tttt(:,11), barrier="on" )

       wrk=gam
       call rotv_rotorb(ML0,MBC,MB0,unk(n1,MB_0_CPMD,k,s),psi_tmp(n1,MB_0_CPMD), 1.0d0, 0.0d0 )

       !call watchb( ttmp, tttt(:,12), barrier="on" )

!$OMP parallel workshare
       psi_n(n1:n2,MB_0_CPMD:MB_1_CPMD,k,s) = &
       psi_n(n1:n2,MB_0_CPMD:MB_1_CPMD,k,s) + psi_tmp(n1:n2,MB_0_CPMD:MB_1_CPMD)
!$OMP end parallel workshare

       !call watchb( ttmp, tttt(:,13), barrier="on" )

       tmp1=1.0d0/dt
!$OMP parallel workshare
       psi_v(n1:n2,MB_0_CPMD:MB_1_CPMD,k,s) = &
       psi_v(n1:n2,MB_0_CPMD:MB_1_CPMD,k,s) + tmp1*psi_tmp(n1:n2,MB_0_CPMD:MB_1_CPMD)
!$OMP end parallel workshare

       !call watchb( ttmp, tttt(:,14), barrier="on" )

    end do ! k
    end do ! s

    !call watchb( ttmp, barrier="on" )

!$OMP parallel workshare
    unk(:,MB_0_CPMD:MB_1_CPMD,:,:) = psi_n(:,MB_0_CPMD:MB_1_CPMD,:,:)
!$OMP end parallel workshare


    !call watchb( ttmp, tttt(:,15), barrier="on" )

!    if ( myrank == 0 ) then
!       do i=1,15
!          write(*,'("time_rotorb(",i2.2,")",2f10.5)') i,tttt(1:2,i)
!       end do
!    end if

    deallocate (psi_tmp)
    !call write_border( 1, "rotorb(end)" )

  end subroutine rotorb_mbp

  SUBROUTINE rotv_rotorb(ii,MB0,MB1,u,u2,a1,a2)
  implicit none
! integer,intent(IN) :: ii, MB0, MB1, MB_0, MB_1
  integer,intent(IN) :: ii, MB0, MB1
! real(8),intent(IN)    :: u(ii,MB0)
! real(8),intent(INOUT) ::u2(ii,MB0)
  real(8),intent(IN)    :: u(ii,MB1)
  real(8),intent(INOUT) ::u2(ii,MB1)
  real(8),intent(IN) :: a1, a2
  real(8),allocatable :: utmp(:,:),utmp2(:,:)
  real(8),parameter :: on = 1.0d0
  real(8),parameter :: zr = 0.0d0
  integer :: ns,ne,nn,ms,me,mm,i0,j0
  integer :: IPROW, IPCOL, iroot1, iroot2
  integer :: irank_g,ierr
  integer :: MBLK

!
! block-cycle, band-parallel
!________________________________________________________________________________
    integer :: ncycle, ncycle2, icycle, icycle2, irankb, nbss, irankc
    integer :: msV, meV, nsV, neV

    real(8), allocatable :: aa(:,:)
    real(8),parameter :: zero=0.0d0

    MBLK=max(MBSIZE,NBSIZE)
    ncycle  = (MB0-1)/MBLK+1
    ncycle2 = (ncycle-1)/np_band+1
    allocate ( aa(ii,MBLK) ) ; aa=zero


!   allocate( utmp(ii,MB_0:MB_1) )
    allocate( utmp(ii,MB1) )
    utmp=0.d0

       i0=0
       do icycle=1,ncycle
          ms = MBLK*(icycle-1)+1
          me = min(ms+MBLK-1,MB0)
          mm = me-ms+1

          msV = MBLK*int((icycle-1)/np_band) + 1
          meV = min(msV+MBLK-1, MB1)

          irankb = mod( icycle-1, np_band )

          IPROW=mod( (ms-1)/MBLK,NPROW )
          if (irankb == myrank_b ) aa(:,1:meV-msV+1) = u(:,msV:meV)

          call rsdft_bcast ( aa, irankb, comm_band )

          j0=0
          do icycle2=1,ncycle2
             nbss = (icycle2-1)*np_band+myrank_b + 1
             ns = MBLK*(nbss-1)+1
             ne = min(ns+MBLK-1,MB0)
             nn = ne-ns+1
             nsV = MBLK*(icycle2-1)+1
             neV = min( nsV+MBLK-1, MB1 )

             IPCOL=mod( (ns-1)/MBLK,NPCOL )

!            write(*,*) "# myrank=",myrank," ms,me=",ms,me,"  iprow=",iprow

             iroot1 = usermap(IPROW,IPCOL,1)
             iroot2 = usermap(IPROW,IPCOL,2)

             if ( mm < 1 .or. nn < 1 ) cycle

             irankc = mod( (ns-1)/MBLK, np_band )
             if ( irankc == myrank_b ) then     !(#1)

                allocate( utmp2(ms:me,ns:ne) )

                if ( iroot1 == myrank ) then
!$OMP parallel workshare
!                  utmp2(ms:me,ns:ne)=gam(i0+1:i0+mm,j0+1:j0+nn)
!                  utmp2(ms:me,ns:ne)=gam(i0+1:i0+mm,j0+1:j0+nn)
                   utmp2(ms:me,ns:ne)=wrk(i0+1:i0+mm,j0+1:j0+nn)
!$OMP end parallel workshare
                   j0=j0+nn
                end if
                call mpi_bcast(utmp2,mm*nn,MPI_REAL8,iroot2,comm_grid,ierr)
!               call dgemm('N','N',ii,nn,mm,on,u(1,ms),ii,utmp2(ms,ns),mm,on,utmp(1,ns),ii)
!               call dgemm('N','N',ii,nn,mm,on,aa,ii,utmp2(ms,ns),mm,on,utmp(1,nsV),ii)
                call dgemm('N','N',ii,nn,mm,a1,aa,ii,utmp2(ms,ns),mm,on,utmp(1,nsV),ii)
                deallocate( utmp2 )

             endif ! (#1)


          end do ! ns

          if ( any(usermap(IPROW,0:NPCOL-1,1)==myrank) ) i0=i0+mm

       end do ! ms

!      u2(:,:)=utmp(:,:)
       u2(:,:)=a2*u2(:,:)+utmp(:,:)

    deallocate( utmp )

    return
  END SUBROUTINE rotv_rotorb

!-----------------------------------------------------------------------
!     orthogonalization of wave function velocity
!----------------------------------------------------------------------- 
  subroutine rotorb2_mbp

    implicit none
    integer :: k,s,i
    real(8),parameter :: on= 1.0d0
    real(8),parameter :: om=-1.0d0
    real(8) :: ttmp(2),tttt(2,2)
    integer :: ii0, jj0

    !call write_border( 1, "rotorb2(start)" )
    !call watchb( ttmp, barrier="on" ); tttt=0.0d0
!debug
    n1  = idisp(myrank)+1
    n2  = idisp(myrank)+ircnt(myrank)
    ML0 = ircnt(myrank)
    MB0 = MB_1_CPMD - MB_0_CPMD + 1
!debug

    do s=MSP_0,MSP_1
    do k=MBZ_0,MBZ_1

       MBT = mstocck(k,s)

       !call watchb( ttmp, barrier="on" )
!
       wrk=0.d0
!
       call overlap5_mbp(s,k,ii0,jj0)

       !call watchb( ttmp, tttt(:,1), barrier="on" )

! ---(1)
!      call DGEMM('N','N',ML0,MB0,MBT,om,unk(n1,1,k,s),ML0 &
!           ,wrk(1,MB_0_CPMD),MBC,on,psi_v(n1,MB_0_CPMD,k,s),ML0)
       call rotv_rotorb(ML0,MBC,MB0,unk(n1,MB_0_CPMD,k,s),psi_v(n1,MB_0_CPMD,k,s), -1.0d0, 1.0d0)

       !call watchb( ttmp, tttt(:,2), barrier="on" )

    end do ! k
    end do ! s

    !call write_border( 1, "rotorb2(end)" )

  end subroutine rotorb2_mbp

end module rotorb_mbp_module
