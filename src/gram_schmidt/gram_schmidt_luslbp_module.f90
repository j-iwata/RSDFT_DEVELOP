module gram_schmidt_luslbp_module

  use sl_tools_module, only: slinfo

  implicit none

  private
  public :: gram_schmidt_luslbp

  integer :: icontxt_sys
  integer :: LDR,LDC
  character(1) :: UPLO='U'

  logical :: has_init_done = .false.

  interface gram_schmidt_luslbp
    module procedure d_gram_schmidt_luslbp, z_gram_schmidt_luslbp
  end interface

  type(slinfo) :: sl

contains

  subroutine init_luslbp( n )
    use sl_tools_module, only: restore_param_sl
    use pinfo_module, only: get_world_rank_pinfo
    use parallel_module, only: myrank_k, myrank_s, myrank
    use io_tools_module, only: IOTools_readIntegerKeyword
    implicit none
    integer,intent(in) :: n
    integer :: i,j,k,l,ierr,itmp(4)
    integer,allocatable :: usermap(:,:)
    integer,external :: NUMROC

    if ( has_init_done ) return

    call write_border( 0, " init_luslbp(start)" )

    sl%nprow =0
    sl%npcol =0
    sl%mbsize=0
    sl%nbsize=0

    call restore_param_sl( sl )

    itmp = (/ sl%nprow, sl%npcol, sl%mbsize, sl%nbsize /)
    call IOTools_readIntegerKeyword( "SCLGS", itmp )

    if( myrank == 0 )then
      write(*,*) "sl%nprow =",sl%nprow
      write(*,*) "sl%npcol =",sl%npcol
      write(*,*) "sl%mbsize=",sl%mbsize
      write(*,*) "sl%nbsize=",sl%nbsize
    end if

! ---

    allocate( usermap(0:sl%nprow-1,0:sl%npcol-1) ); usermap=0

    k = get_world_rank_pinfo( 0,0,myrank_k,myrank_s ) - 1
    do j = 0, sl%npcol-1
      do i = 0, sl%nprow-1
        k = k + 1
        usermap(i,j) = k
      end do
    end do

    call blacs_get( 0, 0, icontxt_sys )
    sl%icontxt_a = icontxt_sys
    call blacs_gridmap( sl%icontxt_a, usermap, sl%nprow, sl%nprow, sl%npcol )

    if( any(usermap==myrank) )then
      call blacs_gridinfo( sl%icontxt_a, sl%nprow, sl%npcol, sl%myrow, sl%mycol )
      LDR = NUMROC( n, sl%mbsize, sl%myrow, 0, sl%nprow )
      LDC = NUMROC( n, sl%nbsize, sl%mycol, 0, sl%npcol )
      call descinit( sl%desca,n,n,sl%mbsize,sl%nbsize,0,0,sl%icontxt_a,LDR,ierr )
      call descinit( sl%descz,n,n,sl%mbsize,sl%nbsize,0,0,sl%icontxt_a,LDR,ierr )
    end if

    deallocate( usermap )

    has_init_done = .true.

    call write_border( 0, " init_luslbp(end)" )

  end subroutine init_luslbp


  subroutine d_gram_schmidt_luslbp( u, dV )
    use rsdft_allreduce_module, only: rsdft_allreduce
    use sl_tools_module, only: distribute_matrix, gather_matrix
    use parallel_module, only: get_range_parallel
    use watch_module, only: watchb
    use calc_overlap_gs_module, only: calc_overlap_gs
    implicit none
    real(8),intent(inout) :: u(:,:)
    real(8),intent(in) :: dV
    real(8),allocatable :: Ssub(:,:)
    real(8),allocatable :: S(:,:),utmp(:,:),S1(:,:)
    real(8),parameter :: z0=0.0d0, z1=1.0d0
    integer :: ng,nb,i,n0,n1
    character(3),parameter :: ON_OFF='off'
    real(8) :: ttmp(2),tt(2,14),tini(2)
    logical :: disp_on
    integer :: j,i0,j0

    call write_border( 1, " d_gram_schmidt_luslbp(start)" )
    call check_disp_switch( disp_on, 0 )

    tt=0.0d0; call watchb( ttmp, barrier=ON_OFF ); tini=ttmp

    ng = size( u, 1 )
    nb = size( u, 2 ) ! assuming nb is the whole band size

    call get_range_parallel( n0, n1, 'b' )

    if ( .not.has_init_done ) call init_luslbp( nb )

    allocate( Ssub(LDR,LDC) ); Ssub=z0
    allocate( S(nb,nb) ); S=z0

    call watchb( ttmp, tt(:,1), barrier=ON_OFF )

    ! --- (1)
    ! call DSYRK( UPLO, 'C', nb, ng, dV, u, ng, z0, S, nb )
    ! call watchb( ttmp, tt(:,2), barrier=ON_OFF )
    ! call rsdft_allreduce( S, 'g' )
    ! call watchb( ttmp, tt(:,3), barrier=ON_OFF )
    ! call distribute_matrix( sl, S, Ssub )
    ! call watchb( ttmp, tt(:,4), barrier=ON_OFF )
    !
    ! --- (2)
    call calc_overlap_gs( nb, u(:,n0:n1), dV, S )
    call watchb( ttmp, tt(:,2), barrier=ON_OFF )
    call distribute_matrix( sl, S, Ssub )
    call watchb( ttmp, tt(:,4), barrier=ON_OFF )
    ! ---

    deallocate( S )

    call watchb( ttmp, tt(:,5), barrier=ON_OFF )

    call PDPOTRF( UPLO, nb, Ssub, 1, 1, sl%desca, i ) ! Cholesky factorization

    call watchb( ttmp, tt(:,6), barrier=ON_OFF )

    call PDTRTRI( UPLO, 'N', nb, Ssub, 1, 1, sl%desca, i )

    call watchb( ttmp, tt(:,7), barrier=ON_OFF )

    ! ----- (1)
    !
    allocate( S(nb,nb) ); S=z0
    call watchb( ttmp, tt(:,8), barrier=ON_OFF )
    call gather_matrix( sl, Ssub, S, sl%icontxt_a )
    call watchb( ttmp, tt(:,9), barrier=ON_OFF )
    ! --- (1-1)
    !
    !(1)
    !call DTRMM( 'R', 'U', 'N', 'N', ng, nb, z1, S, nb, u, ng )
    !call watchb( ttmp, tt(:,11), barrier=ON_OFF )
    !
    !(2)
    !allocate( utmp(ng,nb) ); utmp=u
    !call watchb( ttmp, tt(:,10), barrier=ON_OFF )
    !do i = n0, n1
    !  u(:,i) = matmul( utmp(:,1:i), S(1:i,i) )
    !end do
    !call watchb( ttmp, tt(:,11), barrier=ON_OFF )
    !deallocate( utmp )
    !call watchb( ttmp, tt(:,12), barrier=ON_OFF )
    !
    !(3)
    allocate( utmp(ng,n0:n1) ); utmp=z0
    call watchb( ttmp, tt(:,10), barrier=ON_OFF )
    call DGEMM('N','N',ng,n1-n0+1,nb,z1,u,ng,S(1,n0),nb,z0,utmp,ng)
    call watchb( ttmp, tt(:,11), barrier=ON_OFF )
    u(:,n0:n1)=utmp
    deallocate( utmp )
    call watchb( ttmp, tt(:,12), barrier=ON_OFF )
    !
    ! --- (1-2)
    ! call dtrmm_bp( UPLO, sl, u, v, Ssub, S )
    ! ---
    deallocate( S )
    !
    ! ----- (2)
    ! call dtrmm_bp( UPLO, sl, u, Ssub )
    !
    ! -----

    deallocate( Ssub )

    call watchb( ttmp, tt(:,13), barrier=ON_OFF )
    tt(:,14) = ttmp - tini

    ! if ( disp_on ) then
    !   do i = 1, 13
    !     write(*,'(1x,"time_d_gs_luslbp(",i3.3,")=",2f12.5)') i,tt(:,i)
    !   end do
    !   write(*,'(1x,"time_d_gs_luslbp(tot)=",2f12.5)') tt(:,14)
    ! end if

    call write_border( 1, " d_gram_schmidt_luslbp(end)" )

  end subroutine d_gram_schmidt_luslbp


  subroutine z_gram_schmidt_luslbp( u, dV )
    use rsdft_allreduce_module, only: rsdft_allreduce
    use sl_tools_module, only: distribute_matrix, gather_matrix
    use parallel_module, only: get_range_parallel
    use watch_module, only: watchb
    use calc_overlap_gs_module, only: calc_overlap_gs
    implicit none
    complex(8),intent(inout) :: u(:,:)
    real(8),intent(in) :: dV
    complex(8),allocatable :: Ssub(:,:)
    complex(8),allocatable :: S(:,:),utmp(:,:),S1(:,:)
    complex(8),parameter :: z0=(0.0d0,0.0d0), z1=(1.0d0,0.0d0)
    complex(8) :: zdV
    integer :: ng,nb,i,n0,n1
    character(3),parameter :: ON_OFF='off'
    real(8) :: ttmp(2),tt(2,14),tini(2)
    logical :: disp_on
    integer :: j,i0,j0

    call write_border( 1, " z_gram_schmidt_luslbp(start)" )
    call check_disp_switch( disp_on, 0 )

    tt=0.0d0; call watchb( ttmp, barrier=ON_OFF ); tini=ttmp

    ng = size( u, 1 )
    nb = size( u, 2 ) ! assuming nb is the whole band size

    zdV = dV

    call get_range_parallel( n0, n1, 'b' )

    if ( .not.has_init_done ) call init_luslbp( nb )

    allocate( Ssub(LDR,LDC) ); Ssub=z0
    allocate( S(nb,nb) ); S=z0

    call watchb( ttmp, tt(:,1), barrier=ON_OFF )

    ! --- (1)
    call ZHERK( UPLO, 'C', nb, ng, zdV, u, ng, z0, S, nb )
    call watchb( ttmp, tt(:,2), barrier=ON_OFF )
    call rsdft_allreduce( S, 'g' )
    call watchb( ttmp, tt(:,3), barrier=ON_OFF )
    call distribute_matrix( sl, S, Ssub )
    call watchb( ttmp, tt(:,4), barrier=ON_OFF )
    !
    ! --- (2)
    ! call calc_overlap_gs( nb, u(:,n0:n1), dV, S )
    ! call watchb( ttmp, tt(:,2), barrier=ON_OFF )
    ! call distribute_matrix( sl, S, Ssub )
    ! call watchb( ttmp, tt(:,4), barrier=ON_OFF )
    ! ---

    deallocate( S )

    call watchb( ttmp, tt(:,5), barrier=ON_OFF )

    call PZPOTRF( UPLO, nb, Ssub, 1, 1, sl%desca, i ) ! Cholesky factorization

    call watchb( ttmp, tt(:,6), barrier=ON_OFF )

    call PZTRTRI( UPLO, 'N', nb, Ssub, 1, 1, sl%desca, i )

    call watchb( ttmp, tt(:,7), barrier=ON_OFF )

    ! ----- (1)
    !
    allocate( S(nb,nb) ); S=z0
    call watchb( ttmp, tt(:,8), barrier=ON_OFF )
    call gather_matrix( sl, Ssub, S, sl%icontxt_a )
    call watchb( ttmp, tt(:,9), barrier=ON_OFF )
    ! --- (1-1)
    !
    !(1)
    !call ZTRMM( 'R', 'U', 'N', 'N', ng, nb, z1, S, nb, u, ng )
    !call watchb( ttmp, tt(:,11), barrier=ON_OFF )
    !
    !(2)
    !allocate( utmp(ng,nb) ); utmp=u
    !call watchb( ttmp, tt(:,10), barrier=ON_OFF )
    !do i = n0, n1
    !  u(:,i) = matmul( utmp(:,1:i), S(1:i,i) )
    !end do
    !call watchb( ttmp, tt(:,11), barrier=ON_OFF )
    !deallocate( utmp )
    !call watchb( ttmp, tt(:,12), barrier=ON_OFF )
    !
    !(3)
    allocate( utmp(ng,n0:n1) ); utmp=z0
    call watchb( ttmp, tt(:,10), barrier=ON_OFF )
    call ZGEMM('N','N',ng,n1-n0+1,nb,z1,u,ng,S(1,n0),nb,z0,utmp,ng)
    call watchb( ttmp, tt(:,11), barrier=ON_OFF )
    u(:,n0:n1)=utmp
    deallocate( utmp )
    call watchb( ttmp, tt(:,12), barrier=ON_OFF )
    !
    ! --- (1-2)
    ! call ztrmm_bp( UPLO, sl, u, v, Ssub, S )
    ! ---
    deallocate( S )
    !
    ! ----- (2)
    ! call ztrmm_bp( UPLO, sl, u, Ssub )
    !
    ! -----

    deallocate( Ssub )

    call watchb( ttmp, tt(:,13), barrier=ON_OFF )
    tt(:,14) = ttmp - tini

    ! if ( disp_on ) then
    !   do i = 1, 13
    !     write(*,'(1x,"time_z_gs_luslbp(",i3.3,")=",2f12.5)') i,tt(:,i)
    !   end do
    !   write(*,'(1x,"time_z_gs_luslbp(tot)=",2f12.5)') tt(:,14)
    ! end if

    call write_border( 1, " z_gram_schmidt_luslbp(end)" )

  end subroutine z_gram_schmidt_luslbp

end module gram_schmidt_luslbp_module
