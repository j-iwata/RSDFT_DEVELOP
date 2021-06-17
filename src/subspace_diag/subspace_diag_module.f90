module subspace_diag_module

  use parallel_module, only: np_band, ir_band, id_band, myrank &
                            ,disp_switch_parallel
  use subspace_diag_variables
  use subspace_diag_la_module
  use subspace_diag_sl_module
  use subspace_sdsl_module
  use subspace_diag_ncol_module
  use io_tools_module

  implicit none

  private
  public :: subspace_diag, init_subspace_diag

  integer :: ialgo_sd=1

contains


  subroutine subspace_diag( k,s,ML_0,ML_1,MK_0,MS_0,unk,esp )
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
       select case( ialgo_sd )
       case( 0 )
         call subspace_diag_la(k,s)
       case( 1 )
         call subspace_diag_sl(k,s)
       case default
         k0 = k - MK_0 + 1
         s0 = s - MS_0 + 1
         call subspace_sdsl( k, s, unk(:,:,k0,s0), esp(:,k0,s0) )
       end select
#endif
    end if
  end subroutine subspace_diag


  subroutine init_subspace_diag( MB_in )
    implicit none
    integer,intent(in) :: MB_in
    integer :: i,j,mm,ms,me,nme,ne,nn,MB,num_sq_blocks_per_bp

    call write_border( 0, " init_subspace_diag(start)" )

    MB_diag = MB_in

    call IOTools_readIntegerKeyword( 'IALGO_SD', ialgo_sd )

    if( ialgo_sd /= 1 )then
      call write_border( 0, " init_subspace_diag(return)" )
      return
    end if

    MB  = MB_diag
    nme = (MB*MB+MB)/2

    call parameter_check(nme,MB)

    if ( .not.allocated(mat_block) ) then
       allocate( mat_block(0:np_band-1,0:4) )
    end if
    mat_block(:,:) = 0

    do i=0,np_band-1
       me=id_band(i)+ir_band(i)
       mm=ir_band(i)
       mat_block(i,0)=(mm*(mm+1))/2
       mat_block(i,1)=MB-me
       mat_block(i,2)=mm
       mat_block(i,3)=mat_block(i,0)+mat_block(i,1)*mat_block(i,2)
    end do

    if ( sum(mat_block(:,3)) /= nme ) then
       write(*,*) sum(mat_block(:,3)),myrank,nme
       stop "stop@init_subspace_diag(1)"
    end if

    if ( np_band > 1 ) then

       num_sq_blocks_per_bp = ((np_band-1)*np_band)/2 / np_band

       do j=0,np_band/2-1
          do i=np_band-1,j+1,-1
             mm=ir_band(i)
             if( num_sq_blocks_per_bp*mm < mat_block(j,1) )then
                mat_block(j,1)=mat_block(j,1)-mm
                mat_block(i,1)=mat_block(i,1)+mm
             end if
          end do
       end do
       mat_block(:,3)=mat_block(:,0)+mat_block(:,1)*mat_block(:,2)

       if ( sum(mat_block(:,3))/=nme ) then
          write(*,*) sum(mat_block(:,3)),myrank,nme
          stop "stop@init_subspace_diag(2)"
       end if

    end if

    do i=0,np_band-1
       mat_block(i,4)=sum( mat_block(0:i,3) )-mat_block(i,3)
    end do

    if ( disp_switch_parallel ) then
       write(*,'(1x,6a10)') "rank_b","tri","m","n","nme","idis"
       do i=0,np_band-1
          write(*,'(1x,6i10)') i,(mat_block(i,j),j=0,4)
       end do
    end if

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
