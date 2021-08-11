module subspace_diag_module

  use subspace_diag_variables
  use subspace_diag_la_module
  use subspace_diag_sl_module
  use subspace_sdsl_module
  use subspace_diag_ncol_module
  use io_tools_module, only: IOTools_readIntegerKeyword
  use subspace_mate_sl2_module, only: subspace_mate_sl2
  use subspace_solv_sl2_module, only: subspace_solv_sl2
  use subspace_rotv_sl2_module, only: subspace_rotv_sl2

  implicit none

  private
  public :: subspace_diag
  public :: init_subspace_diag

  integer :: ialgo_sd=1

contains


  subroutine subspace_diag( k,s,ML_0,ML_1,MK_0,MS_0,unk,esp )
    use sl_variables, only: sl1
    implicit none
    integer,intent(in) :: k,s,ML_0,ML_1,MK_0,MS_0
#ifdef _DRSDFT_
    real(8),intent(inout) :: unk(:,:,:,:)
#else
    complex(8),intent(inout) :: unk(:,:,:,:)
#endif
    real(8),intent(inout) :: esp(:,:,:)
    integer :: k0,s0
    if ( flag_noncollinear ) then
      call subspace_diag_ncol( k, ML_0,ML_1, unk, esp )
    else
#ifdef _LAPACK_
      call subspace_diag_la(k,s)
#else
      k0 = k - MK_0 + 1
      s0 = s - MS_0 + 1
      select case( ialgo_sd )
      case( 0 )
        call subspace_diag_la(k,s)
      case( 1 )
        call subspace_diag_sl(k,s)
      case( 2 )
        call subspace_mate_sl2( k, s, unk(:,:,k0,s0) )
        call subspace_solv_sl2( sl1, esp(:,k0,s0) )
        call subspace_rotv_sl2( unk(:,:,k0,s0) )
      case default
        call subspace_sdsl( k, s, unk(:,:,k0,s0), esp(:,k0,s0) )
      end select
#endif
    end if
  end subroutine subspace_diag


  subroutine init_subspace_diag( MB, nprocs_MB )
    use parallel_module, only: load_div_parallel, disp_switch_parallel
    implicit none
    integer,intent(in) :: MB, nprocs_MB
    integer :: i,j,mm,me,tot_num_of_mate
    integer,allocatable :: ircnt_MB(:), idisp_MB(:)

    call write_border( 0, " init_subspace_diag(start)" )

! ---

    MB_diag = MB

    call IOTools_readIntegerKeyword( 'IALGO_SD', ialgo_sd )
    i=0; call IOTools_readIntegerKeyword( 'SCL2', i )
    if ( i /= 0 ) ialgo_sd = 2

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

    if ( disp_switch_parallel ) then
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
