MODULE overlap_cpmd_module

  use wf_module, only: unk
  use rgrid_module, only: dV
  use parallel_module
  use cpmd_variables, only: wrk,tau,sig,scr,MBC,MBT,psi_n,psi_v &
       ,ir_band_cpmd,id_band_cpmd,MB_0_CPMD,MB_1_CPMD
  use watch_module
  use calc_overlap_module
  use calc_overlap_bp_module
  use calc_overlap_bp2_module
  use rsdft_mpi_module

  implicit none

  PRIVATE
  PUBLIC :: overlap2, overlap4, overlap5

  real(8),allocatable :: w1(:),w2(:)

CONTAINS


  SUBROUTINE overlap2(s,k)
    implicit none
    integer,intent(IN) :: s,k
    integer :: i,j,l,m,n,mbn,mblocaldim
    integer :: n1,n2,ML0,nn1,nn2,irank_b,ns,ne,ms,me
    integer :: mm1,mm2,ierr,m1,m2
    real(8) :: memax,mem,nop_tot,nop_max,nop_min,nop_0
    real(8) :: nop(13),ct(13),et(13)
    real(8) :: ctime0,ctime1,etime0,etime1,ctt,ett,ttmp(2),tttt(2,9)
    real(8) :: tmp1
    real(8),parameter :: zero=0.d0
!    integer,allocatable :: ir(:),id(:)
    integer :: nbss,k1,ib,NBAND_BLK,ncycle,mrnk
!    integer,allocatable :: ir_loc(:),id_loc(:)
    integer :: ls_loc,le_loc,li_loc
    complex(8),allocatable :: ztmp(:,:)

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp ); tttt=0.0d0

    n1    = idisp(myrank)+1
    n2    = idisp(myrank)+ircnt(myrank)
    ML0   = n2-n1+1
    mrnk  = id_class(myrank,4)
    m1    = id_band_cpmd(myrank_b)+1
    m2    = id_band_cpmd(myrank_b)+ir_band_cpmd(myrank_b)

!    if ( np_band > 1 ) then
!       allocate( ir(0:np_band-1) ) ; ir=0
!       allocate( id(0:np_band-1) ) ; id=0
!       ir(0:np_band-1)=ir_band_cpmd(0:np_band-1)*ML0
!       id(0:np_band-1)=id_band_cpmd(0:np_band-1)*ML0
!       call rsdft_allgatherv(psi_n(:,m1:m2,k,s),psi_n(:,:,k,s),ir,id,comm_band)
!!       call mpi_allgatherv(psi_n(n1,MB_0_CPMD,k,s),ir(mrnk),MPI_REAL8 &
!!            ,psi_n(n1,1,k,s),ir,id,MPI_REAL8,comm_band,ierr)
!    end if

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp, tttt(:,1) )

! --- (1)
!
!    call dsyrk('L','t',MBT,ML0,-dV,psi_n(n1,1,k,s),ML0,zero,sig,MBC)
!
! --- (2)
!
    call calc_overlap_bp2( MBT, psi_n(:,m1:m2,k,s), psi_n(:,m1:m2,k,s), -dV, sig )

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp, tttt(:,2) )

    call mochikae_matrix2( sig, 2 )
!
! -------

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp, tttt(:,3) )

!    call allreduce_sub( sig )
#ifdef test
    n = (MBC*(MBC+1))/2
    if ( .not.allocated(w1) ) then
       allocate( w1(n) ) ; w1=0.0d0
    end if
    do j=1,MBC
    do i=j,MBC
       m=(j-1)*MBC-(j*(j-1))/2+i
       w1(m)=sig(i,j)
    end do
    end do
    call watchb( ttmp, tttt(:,3) )
    !call mpi_allreduce(MPI_IN_PLACE,w1,n,mpi_real8,mpi_sum,comm_grid,ierr)
    call rsdft_allreduce_sum( w1(1:n), comm_grid )
    call watchb( ttmp, tttt(:,4) )
    do j=1,MBC
    do i=j,MBC
       m=(j-1)*MBC-(j*(j-1))/2+i
       sig(i,j)=w1(m)
    end do
    end do
#endif

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp, tttt(:,5) )

!$OMP parallel do
    do j=1  ,MBC
       !do i=j+1,MBC
       !   sig(j,i) = sig(i,j)
!       do i=1,j-1
!          sig(i,j) = sig(j,i)
!       end do
       sig(j,j)=sig(j,j)+1.0d0
    end do
!$OMP end parallel do

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp, tttt(:,6) )

!    if ( np_band > 1 ) deallocate( id,ir )

!    if ( myrank == 0 ) then
!       do i=1,6
!          write(*,'(2x,"time_overlap2(",i1,")",2f10.5)') i,tttt(:,i)
!       end do
!    end if

    return

  END SUBROUTINE overlap2


  SUBROUTINE overlap4(s,k)
    implicit none
    integer,intent(IN) :: s,k
    integer :: i,j,ij,l,m,n,mbn,mblocaldim
    integer :: n1,n2,ML0,nn1,nn2,irank_b,ns,ne,ms,me
    integer :: mm1,mm2,ierr,m1,m2
    real(8) :: memax,mem,nop_tot,nop_max,nop_min,nop_0
    real(8) :: nop(13),ct(13),et(13)
    real(8) :: ctime0,ctime1,etime0,etime1,ctt,ett,ttmp(2),tttt(2,9)
    real(8) :: tmp1
    real(8),parameter :: zero=0.d0
!    integer,allocatable :: ir(:),id(:)
    integer :: nbss,k1,ib,NBAND_BLK,ncycle,mrnk
    integer,allocatable :: ir_loc(:),id_loc(:)
    integer ls_loc,le_loc,li_loc
    complex(8),allocatable :: ztmp(:,:)

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp ); tttt=0.0d0

    n1    = idisp(myrank)+1
    n2    = idisp(myrank)+ircnt(myrank)
    ML0   = n2-n1+1
    mrnk  = id_class(myrank,4)
    m1    = id_band_cpmd(myrank_b)+1
    m2    = id_band_cpmd(myrank_b)+ir_band_cpmd(myrank_b)

!    if ( np_band > 1 ) then
!       allocate( ir(0:np_band-1) ) ; ir=0
!       allocate( id(0:np_band-1) ) ; id=0
!       ir(0:np_band-1)=ir_band_cpmd(0:np_band-1)*ML0
!       id(0:np_band-1)=id_band_cpmd(0:np_band-1)*ML0
!       call rsdft_allgatherv(unk(:,m1:m2,k,s),unk(:,:,k,s),ir,id,comm_band)
!!       call mpi_allgatherv(  unk(n1,MB_0_CPMD,k,s),ir(mrnk),MPI_REAL8 &
!!            ,  unk(n1,1,k,s),ir,id,MPI_REAL8,comm_band,ierr)
!       call rsdft_allgatherv(psi_n(:,m1:m2,k,s),psi_n(:,:,k,s),ir,id,comm_band)
!!       call mpi_allgatherv(psi_n(n1,MB_0_CPMD,k,s),ir(mrnk),MPI_REAL8 &
!!            ,psi_n(n1,1,k,s),ir,id,MPI_REAL8,comm_band,ierr)
!    end if

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp, tttt(:,1) )

#ifdef _DRSDFT_
!
! --- (1)
!
!    call calc_overlap(ML0,MBT,unk(n1,1,k,s),psi_n(n1,1,k,s),-dV,tau)
!
! --- (2)
!
    call calc_overlap_bp2( MBT, unk(:,m1:m2,k,s), psi_n(:,m1:m2,k,s), -dV, tau )

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp, tttt(:,2) )

    call mochikae_matrix2( tau, 2 )
    !call rsdft_allreduce_sum( tau, comm_grid )
!
! -------
!
#endif

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp, tttt(:,3) )

!    call allreduce_sub( tau )

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp, tttt(:,4) )

!$OMP parallel do
    do j=1,MBC
       !do i=j+1,MBC
       !   tau(j,i) = tau(i,j)
!       do i=1,j-1
!          tau(i,j) = tau(j,i)
!       end do
       tau(j,j)=tau(j,j)+1.0d0
    end do
!$OMP end parallel do

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp, tttt(:,5) )

!    if ( np_band > 1 ) deallocate( id,ir )

!    if ( myrank == 0 ) then
!       do i=1,5
!          write(*,'(2x,"time_overlap4(",i1,")",2f10.5)') i,tttt(:,i)
!       end do
!    end if

    return

  END SUBROUTINE overlap4


  SUBROUTINE overlap5(s,k)
    implicit none
    integer,intent(IN) :: s,k
    integer :: i,j,ij,l,m,n,mbn,mblocaldim
    integer :: n1,n2,ML0,nn1,nn2,irank_b,ns,ne,ms,me
    integer :: mm1,mm2,ierr,m1,m2
    real(8) :: memax,mem,nop_tot,nop_max,nop_min,nop_0
    real(8) :: nop(13),ct(13),et(13)
    real(8) :: ctime0,ctime1,etime0,etime1,ctt,ett,ttmp(2),tttt(2,9)
    real(8) :: tmp1
    real(8),parameter :: zero=0.d0
!    integer,allocatable :: ir(:),id(:)
    integer :: nbss,k1,ib,NBAND_BLK,ncycle,mrnk
    integer,allocatable :: ir_loc(:),id_loc(:)
    integer ls_loc,le_loc,li_loc
    complex(8),allocatable :: ztmp(:,:)

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp ); tttt=0.0d0

    n1    = idisp(myrank)+1
    n2    = idisp(myrank)+ircnt(myrank)
    ML0   = n2-n1+1
    mrnk  = id_class(myrank,4)
    m1    = id_band_cpmd(myrank_b)+1
    m2    = id_band_cpmd(myrank_b)+ir_band_cpmd(myrank_b)

!    if ( np_band > 1 ) then
!       allocate( ir(0:np_band-1) ) ; ir=0
!       allocate( id(0:np_band-1) ) ; id=0
!       ir(0:np_band-1)=ir_band_cpmd(0:np_band-1)*ML0
!       id(0:np_band-1)=id_band_cpmd(0:np_band-1)*ML0
!       call rsdft_allgatherv(unk(:,m1:m2,k,s),unk(:,:,k,s),ir,id,comm_band)
!!       call mpi_allgatherv(  unk(n1,MB_0_CPMD,k,s),ir(mrnk),MPI_REAL8 &
!!            ,unk(n1,1,k,s),ir,id,MPI_REAL8,comm_band,ierr)
!       call rsdft_allgatherv(psi_v(:,m1:m2,k,s),psi_v(:,:,k,s),ir,id,comm_band)
!!       call mpi_allgatherv(psi_v(n1,MB_0_CPMD,k,s),ir(mrnk),MPI_REAL8 &
!!            ,psi_v(n1,1,k,s),ir,id,MPI_REAL8,comm_band,ierr)
!    end if

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp, tttt(:,1) )

#ifdef _DRSDFT_
!
! --- (1)
!
!    call calc_overlap(ML0,MBT,unk(n1,1,k,s),psi_v(n1,1,k,s),dV,wrk)
!
! --- (2)
!
    call calc_overlap_bp2( MBT, unk(:,m1:m2,k,s), psi_v(:,m1:m2,k,s), dV, wrk )

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp, tttt(:,2) )

    call mochikae_matrix2( wrk, 2 )
    !call rsdft_allreduce_sum( wrk, comm_grid )
!
! -------
#endif

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp, tttt(:,3) )

!    call allreduce_sub( wrk )

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp, tttt(:,4) )

#ifdef test
!$OMP parallel do
    do j=1,MBC
    !do i=j+1,MBC
    !   wrk(j,i) = wrk(i,j)
    do i=1,j-1
       wrk(i,j) = wrk(j,i)
    end do
    end do
!$OMP end parallel do 
#endif

    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call watchb( ttmp, tttt(:,5) )

!    if ( np_band > 1 ) deallocate( id,ir )

!    if ( myrank == 0 ) then
!       do i=1,5
!          write(*,'(2x,"time_overlap5(",i1,")",2f10.5)') i, tttt(:,i)
!       end do
!    end if

    return

  END SUBROUTINE overlap5


  subroutine allreduce_sub( sig )
    implicit none
    real(8),intent(inout) :: sig(:,:)
    integer :: n,i,j,m
    n = (MBC*(MBC+1))/2
    if ( .not.allocated(w1) ) then
       allocate( w1(n) ) ; w1=0.0d0
    end if
    do j=1,MBC
    do i=j,MBC
       m=(j-1)*MBC-(j*(j-1))/2+i
       w1(m)=sig(i,j)
    end do
    end do
    !call mpi_allreduce(MPI_IN_PLACE,w1,n,mpi_real8,mpi_sum,comm_grid,ierr)
    call rsdft_allreduce_sum( w1(1:n), comm_grid )
    do j=1,MBC
    do i=j,MBC
       m=(j-1)*MBC-(j*(j-1))/2+i
       sig(i,j)=w1(m)
    end do
    end do
  end subroutine allreduce_sub


END MODULE overlap_cpmd_module
