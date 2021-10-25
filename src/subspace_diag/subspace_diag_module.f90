module subspace_diag_module

  use subspace_diag_variables, only: MB_diag, mat_block, algo_sd
  use subspace_diag_la_module, only: subspace_diag_la
  use subspace_diag_sl_module, only: subspace_diag_sl
  use subspace_sdsl_module, only: subspace_sdsl
  use subspace_mate_sl2_module, only: subspace_mate_sl2
  use subspace_solv_sl2_module, only: subspace_solv_sl2
  use subspace_rotv_sl2_module, only: d_subspace_rotv_sl2, z_subspace_rotv_sl2

  implicit none

  private
  public :: subspace_diag
  public :: init_subspace_diag

  interface subspace_diag
    module procedure d_subspace_diag, z_subspace_diag
  end interface

  integer :: ialgo_sd

contains

  subroutine d_subspace_diag( k, s, unk, esp )
    use sl_variables, only: sl1
    use watch_module, only: watchb
    use wf_module, only: gather_b_wf
    implicit none
    integer,intent(in) :: k, s
    real(8),intent(inout) :: unk(:,:)
    real(8),intent(inout) :: esp(:)
    real(8) :: ttmp(2),tt(2,4)
    logical :: disp_on
    integer :: i
    character(3),parameter :: ON_OFF='off'
#ifdef _LAPACK_
    call subspace_diag_la( k, s )
#else
    select case( ialgo_sd )
    case( 0 )
      call subspace_diag_la( k, s )
    case( 1 )
      call subspace_diag_sl( k, s )
    case( 2 )
      !tt=0.0d0; call watchb( ttmp, barrier=ON_OFF )
      call gather_b_wf( k, s )
      !call watchb( ttmp, tt(:,1), barrier=ON_OFF )
      call subspace_mate_sl2( k, s, unk )
      !call watchb( ttmp, tt(:,2), barrier=ON_OFF )
      call subspace_solv_sl2( sl1, esp )
      !call watchb( ttmp, tt(:,3), barrier=ON_OFF )
      call d_subspace_rotv_sl2( unk )
      !call watchb( ttmp, tt(:,4), barrier=ON_OFF )
      !call check_disp_switch( disp_on, 0 )
      ! if ( disp_on ) then
      !   do i = 1, 4
      !     write(*,'(1x,"time_sd(",i1,")=",2f12.5)') i, tt(:,i)
      !   end do
      ! end if
    case default
      ! call subspace_sdsl( k, s, unk, esp )
    end select
#endif
  end subroutine d_subspace_diag

  subroutine z_subspace_diag( k, s, unk, esp )
    use sl_variables, only: sl1
    use watch_module, only: watchb
    use wf_module, only: gather_b_wf
    implicit none
    integer,intent(in) :: k, s
    complex(8),intent(inout) :: unk(:,:)
    real(8),intent(inout) :: esp(:)
    real(8) :: ttmp(2),tt(2,4)
    logical :: disp_on
    integer :: i
    character(3),parameter :: ON_OFF='off'
#ifdef _LAPACK_
    call subspace_diag_la( k, s )
#else
    select case( ialgo_sd )
    case( 0 )
      call subspace_diag_la( k, s )
    case( 1 )
      call subspace_diag_sl( k, s )
    case( 2 )
      !tt=0.0d0; call watchb( ttmp, barrier=ON_OFF )
      call gather_b_wf( k, s )
      !call watchb( ttmp, tt(:,1), barrier=ON_OFF )
      call subspace_mate_sl2( k, s, unk )
      !call watchb( ttmp, tt(:,2), barrier=ON_OFF )
      call subspace_solv_sl2( sl1, esp )
      !call watchb( ttmp, tt(:,3), barrier=ON_OFF )
      call z_subspace_rotv_sl2( unk )
      !call watchb( ttmp, tt(:,4), barrier=ON_OFF )
      !call check_disp_switch( disp_on, 0 )
      ! if ( disp_on ) then
      !   do i = 1, 4
      !     write(*,'(1x,"time_sd(",i1,")=",2f12.5)') i, tt(:,i)
      !   end do
      ! end if
    case default
      ! call subspace_sdsl( k, s, unk, esp )
    end select
#endif
  end subroutine z_subspace_diag


  subroutine init_subspace_diag( MB, nprocs_MB )
    use parallel_module, only: load_div_parallel
    implicit none
    integer,intent(in) :: MB, nprocs_MB
    integer :: i,j,mm,me,tot_num_of_mate
    integer,allocatable :: ircnt_MB(:), idisp_MB(:)
    logical :: disp_on

    call write_border( 0, " init_subspace_diag(start)" )
    call check_disp_switch( disp_on, 0 )

! ---

    MB_diag = MB

    ialgo_sd = algo_sd()

    if ( ialgo_sd /= 1 ) then
      call write_border( 0, " init_subspace_diag(return)" )
      return
    end if

! ---

    allocate( ircnt_MB(0:nprocs_MB-1) ); ircnt_MB=0
    allocate( idisp_MB(0:nprocs_MB-1) ); idisp_MB=0
    call load_div_parallel( ircnt_MB, idisp_MB, MB )
    i = minval( ircnt_MB )
    j = maxval( ircnt_MB )
    if ( i /= j ) then
      write(*,*) "All ircnt_MB must be equal: ircnt_MB=",ircnt_MB
      stop 'stop@init_subspace_diag(0)'
    end if

    tot_num_of_mate = ( MB*MB + MB )/2

    call parameter_check( tot_num_of_mate, MB )

    if ( .not.allocated(mat_block) ) then
      allocate( mat_block(0:nprocs_MB-1,0:4) )
    end if
    mat_block(:,:) = 0

    do i=0,nprocs_MB-1
      me=idisp_MB(i)+ircnt_MB(i)
      mm=ircnt_MB(i)
      mat_block(i,0)=(mm*(mm+1))/2
      mat_block(i,1)=MB-me
      mat_block(i,2)=mm
      mat_block(i,3)=mat_block(i,0)+mat_block(i,1)*mat_block(i,2)
    end do

    if ( sum(mat_block(:,3)) /= tot_num_of_mate ) then
      write(*,*) sum(mat_block(:,3)),tot_num_of_mate
      stop "stop@init_subspace_diag(1)"
    end if

    if ( nprocs_MB > 1 ) then

      do j = 0, (nprocs_MB-1)/2-1
        do i = nprocs_MB-1, nprocs_MB-1-((nprocs_MB-1)/2-1)+j, -1
          mm=ircnt_MB(i)
          mat_block(j,1)=mat_block(j,1)-mm
          mat_block(i,1)=mat_block(i,1)+mm
        end do
      end do
      mat_block(:,3)=mat_block(:,0)+mat_block(:,1)*mat_block(:,2)

      if ( sum(mat_block(:,3)) /= tot_num_of_mate ) then
        write(*,*) sum(mat_block(:,3)),tot_num_of_mate
        stop "stop@init_subspace_diag(2)"
      end if

    end if

    do i=0,nprocs_MB-1
      mat_block(i,4)=sum( mat_block(0:i,3) )-mat_block(i,3)
    end do

    if ( disp_on ) then
      write(*,'(1x,6a10)') "rank_b","tri","m","n","nme","idis"
      do i=0,nprocs_MB-1
        write(*,'(1x,6i10)') i,(mat_block(i,j),j=0,4)
      end do
    end if

    deallocate( idisp_MB )
    deallocate( ircnt_MB )

    call write_border( 0, " init_subspace_diag(end)" )

  end subroutine init_subspace_diag

  subroutine parameter_check(nme,MB)
    implicit none
    integer,intent(in) :: nme,MB
    real(8) :: d_nme
    d_nme = ( dble(MB)*dble(MB)+dble(MB) )/2.0d0
    if ( abs(d_nme-nme) > 1.d-10 ) then
      write(*,*) "MB,nme,d_nme=",MB,nme,d_nme
      write(*,*) "MB may be too large"
      stop "stop@init_subspace_diag"
    end if
  end subroutine parameter_check


end module subspace_diag_module
