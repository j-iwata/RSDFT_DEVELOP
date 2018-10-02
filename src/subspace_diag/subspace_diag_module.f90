MODULE subspace_diag_module

  use parallel_module, only: np_band, ir_band, id_band, myrank &
                            ,disp_switch_parallel
  use subspace_diag_variables
  use subspace_diag_la_module
  use subspace_diag_la_bp_module
  use subspace_diag_sl_module
  use subspace_diag_ncol_module
  use vector_tools_module, only: vinfo
  use scalapack_module, only: LLD_R,LLD_C,d_scatter_matrix_scalapack,d_gatherA_matrix_scalapack &
                            , MBSIZE,NBSIZE,LLD_R_e,LLD_C_e,lda_r_e,lda_c_e,idiag

  implicit none

  PRIVATE
  PUBLIC :: subspace_diag
  PUBLIC :: init_subspace_diag
  PUBLIC :: alloc_dealloc_Hsub
  PUBLIC :: alloc_dealloc_Vsub
  PUBLIC :: prep_Vsub_and_check_data
  PUBLIC :: prep_Vsub_and_check_data0
  PUBLIC :: write_Vsub
  PUBLIC :: alloc_dealloc_eigen

CONTAINS


  SUBROUTINE subspace_diag(k,s,v)
    implicit none
    integer,intent(IN) :: k,s
    type(vinfo),intent(IN) :: v(2)
    select case( idiag )
    case default
       select case( np_band )
       case(1 ); call subspace_diag_la(k,s)
       case(2:); call subspace_diag_la_bp(k,s,v)
       end select
    case( "pdsyevd", "eigen_s" )
       call subspace_diag_sl(k,s,v)
    end select
  END SUBROUTINE subspace_diag


  SUBROUTINE init_subspace_diag( MB_in, v2 )
    implicit none
    integer,intent(IN) :: MB_in
    type(vinfo),intent(IN) :: v2
    integer :: i,j,mm,ms,me,nme,ne,nn,je,MB
    integer :: np_band
    logical :: disp_sw

    call write_border( 0, " init_subspace_diag(start)" )

    np_band = v2%pinfo%np

    MB_diag = MB_in

    MB  = MB_diag
    nme = (MB*MB+MB)/2

    call parameter_check(nme,MB)

    if ( .not.allocated(mat_block) ) then
       allocate( mat_block(0:np_band-1,0:4) )
    end if
    mat_block(:,:) = 0

    do i=0,np_band-1
       me=v2%pinfo%id(i)+v2%pinfo%ir(i)
       mm=v2%pinfo%ir(i)
       mat_block(i,0)=(mm*(mm+1))/2
       mat_block(i,1)=MB-me
       mat_block(i,2)=mm
       mat_block(i,3)=mat_block(i,0)+mat_block(i,1)*mat_block(i,2)
    end do

    if ( sum(mat_block(:,3)) /= nme ) then
       write(*,*) sum(mat_block(:,3)),nme
       stop "stop@init_subspace_diag(1)"
    end if

    if ( np_band>1 ) then

       je = int( (np_band+1)*0.5 )-1
       do j=0,je
          do i=np_band-1,j+1,-1
             mm=v2%pinfo%ir(i)
             if( ((np_band-1)*np_band*0.5/np_band+1)*mm < mat_block(j,1) )then
                mat_block(j,1)=mat_block(j,1)-mm
                mat_block(i,1)=mat_block(i,1)+mm
             end if
          end do
       end do
       mat_block(:,3)=mat_block(:,0)+mat_block(:,1)*mat_block(:,2)

       if ( sum(mat_block(:,3))/=nme ) then
          write(*,*) sum(mat_block(:,3)),nme
          stop "stop@init_subspace_diag(2)"
       end if

    end if

    do i=0,np_band-1
       mat_block(i,4)=sum( mat_block(0:i,3) )-mat_block(i,3)
    end do

    call check_disp_switch( disp_sw, 0 )
    if ( disp_sw ) then
       write(*,'(1x,6a10)') "rank_b","tri","m","n","nme","idis"
       do i=0,np_band-1
          write(*,'(1x,6i10)') i,mat_block(i,0:4)
       end do
    end if

!----
    call init_nLB( MB_in, v2 )
!----

    call write_border( 0, " init_subspace_diag(end)" )

  END SUBROUTINE init_subspace_diag


  SUBROUTINE parameter_check(nme,MB)
    implicit none
    integer,intent(IN) :: nme,MB
    real(8) :: d_nme
    d_nme = ( dble(MB)*dble(MB)+dble(MB) )/2.0d0
    if ( abs(d_nme-nme) > 1.d-10 ) then
       write(*,*) "MB,nme,d_nme=",MB,nme,d_nme
       write(*,*) "MB may be too large"
       stop "stop@init_subspace_diag"
    end if
  END SUBROUTINE parameter_check


  SUBROUTINE alloc_dealloc_Hsub( ctrl )
    implicit none
    character(1),intent(IN) :: ctrl
    if ( ctrl == "a" .or. ctrl == "A" ) then
       if ( .not.allocated(Hsub) ) then
          allocate( Hsub(LLD_R,LLD_C) )
          Hsub=zero
       end if
    else if ( ctrl == "d" .or. ctrl == "D" ) then
       if ( allocated(Hsub) ) deallocate(Hsub)
    end if
  END SUBROUTINE alloc_dealloc_Hsub


  SUBROUTINE alloc_dealloc_Vsub( ctrl )
    implicit none
    character(1),intent(IN) :: ctrl
    if ( ctrl == "a" .or. ctrl == "A" ) then
       if ( .not.allocated(Vsub) ) then
          allocate( Vsub(LLD_R,LLD_C) )
          Vsub=zero
       end if
    else if ( ctrl == "d" .or. ctrl == "D" ) then
       if ( allocated(Vsub) ) deallocate(Vsub)
    end if
  END SUBROUTINE alloc_dealloc_Vsub


  SUBROUTINE alloc_dealloc_eigen( ctrl )
    implicit none
    character(1),intent(IN) :: ctrl
    if ( ctrl == "a" .or. ctrl == "A" ) then
       if ( .not.allocated(Hsub_e) ) then
          allocate( Hsub_e(LLD_R_e,LLD_C_e) )
          Hsub_e=zero
       end if
       if ( .not.allocated(Vsub_e) ) then
          allocate( Vsub_e(LLD_R_e,LLD_C_e) )
          Vsub_e=zero
       end if
    else if ( ctrl == "d" .or. ctrl == "D" ) then
       if ( allocated(Hsub_e) ) deallocate(Hsub_e)
       if ( allocated(Vsub_e) ) deallocate(Vsub_e)
    end if
  END SUBROUTINE alloc_dealloc_eigen

  SUBROUTINE prep_Vsub_and_check_data( MB, u )
    implicit none
    integer,intent(IN) :: MB
    real(8),intent(IN) :: u(:,:)
    real(8),allocatable :: Vtot(:,:)
    allocate( Vtot(MB,MB) ) ; Vtot=zero
    call random_number( Vtot )
!    call check_rotv_0( Vtot, u )
    call d_scatter_matrix_scalapack( Vtot, Vsub )
    deallocate( Vtot )
  END SUBROUTINE prep_Vsub_and_check_data


  SUBROUTINE prep_Vsub_and_check_data0( MB, u )
    implicit none
    integer,intent(IN) :: MB
    real(8),intent(IN) :: u(:,:)
    real(8),allocatable :: Vtot(:,:)
    integer :: i,j

    allocate( Vtot(MB,MB) ) ; Vtot=zero
!   call random_number( Vtot )
    call d_gatherA_matrix_scalapack( Vsub, Vtot )

!    call check_rotv_0( Vtot, u )
!   call d_scatter_matrix_scalapack( Vtot, Vsub )
    deallocate( Vtot )
  END SUBROUTINE prep_Vsub_and_check_data0


  SUBROUTINE init_nLB( MB_in, v2 )
    implicit none
    integer,intent(IN) :: MB_in
    type(vinfo),intent(IN) :: v2
    integer :: i,j,mm,ms,me,nme,ne,nn,je,MB
    integer :: np_band
    logical :: disp_sw

    integer :: ncycle, MBLK ,av_ncycle ,ncycle2
    integer :: ncount
    integer, allocatable :: nLB_data(:)

    np_band = v2%pinfo%np
    MB_diag = MB_in
    MB  = MB_diag
    nme = (MB*MB+MB)/2

    MBLK=min(NBSIZE,MBSIZE)
    ncycle    = (MB-1)/MBLK+1
!   ncycle2   = int(ncycle/2)+1
    ncycle2   = int((ncycle-1)/2)+1

    if ( .not.allocated(nLB) ) then
       allocate( nLB(ncycle) )
    end if
    nLB(:) = 0
!   return

!   av_ncycle=int(ncycle*(ncycle+1)/2/ncycle)+1
    av_ncycle=int((ncycle*(ncycle+1)/2-ncycle)/ncycle)+1

    allocate( nLB_data(ncycle) )
    do i=1,ncycle
       nLB_data(i)=ncycle-i+1
    end do

    do i=1,ncycle2
    ncount=ncycle+1
    if ( nLB_data(i) > av_ncycle ) then
      
      nn = nLB_data(i) - av_ncycle 

      if ( nn>0 ) then
         do j=ncycle,ncycle-nn+1,-1
            nLB(j)=nLB(j) -1
            ncount=ncount-1
         enddo
         nLB(i)=ncount
      endif 
    endif
    end do

    call check_disp_switch( disp_sw, 0 )
    if ( disp_sw ) then
       write(*,*) " init_nLB >>> MBLK=",MBLK
       write(*,*) " init_nLB >>> nLB=",nLB
       write(*,*) " init_nLB >>> av_ncycle=",av_ncycle
    endif

  END SUBROUTINE init_nLB

  SUBROUTINE write_Vsub
    implicit none
    logical :: disp_sw
    call check_disp_switch( disp_sw, 0 )
    if ( disp_sw ) then
     write(*,*) Vsub(1:10,1)
    endif
  END SUBROUTINE write_Vsub


END MODULE subspace_diag_module
