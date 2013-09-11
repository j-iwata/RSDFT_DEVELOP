MODULE subspace_diag_module

  use parallel_module

  implicit none

  PRIVATE
  PUBLIC :: prep_subspace_diag,MB_diag,Hsub,Vsub,mat_block &
           ,NBLK1,NBLK2,zero,one,TYPE_MAIN

  integer,allocatable :: mat_block(:,:)

#ifdef _DRSDFT_
  integer,parameter :: TYPE_MAIN = MPI_REAL8
  real(8),allocatable :: Hsub(:,:), Vsub(:,:)
  real(8),parameter :: zero=0.d0,one=1.d0
#else
  integer,parameter :: TYPE_MAIN = MPI_COMPLEX16
  complex(8),allocatable :: Hsub(:,:), Vsub(:,:)
  complex(8),parameter :: zero=(0.d0,0.d0),one=(1.d0,0.d0)
#endif
  integer :: MB_diag,NBLK1,NBLK2

!  logical :: flag_return = .false.

CONTAINS


  SUBROUTINE prep_subspace_diag(MB_in,disp_switch)
    integer,intent(IN) :: MB_in
    logical,intent(IN) :: disp_switch
    integer :: i,j,mm,ms,me,nme,ne,nn,je,MB

!    if ( flag_return ) return
!    flag_return = .true.

    MB_diag = MB_in

    MB  = MB_diag
    nme = (MB*MB+MB)/2
    
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

    if ( sum(mat_block(:,3))/=nme ) then
       write(*,*) sum(mat_block(:,3)),myrank,nme
       stop "stop@prep_subspace_diag"
    end if

    if ( np_band>1 ) then

       je = int( (np_band+1)*0.5 )-1
       do j=0,je
          do i=np_band-1,j+1,-1
             mm=ir_band(i)
             if( ((np_band-1)*np_band*0.5/np_band+1)*mm < mat_block(j,1) )then
                mat_block(j,1)=mat_block(j,1)-mm
                mat_block(i,1)=mat_block(i,1)+mm
             end if
          end do
       end do
       mat_block(:,3)=mat_block(:,0)+mat_block(:,1)*mat_block(:,2)

       if ( sum(mat_block(:,3))/=nme ) then
          write(*,*) sum(mat_block(:,3)),myrank,nme
          stop
       end if

    end if

    do i=0,np_band-1
       mat_block(i,4)=sum( mat_block(0:i,3) )-mat_block(i,3)
    end do

    if ( DISP_SWITCH ) then
       write(*,'(1x,6a10)') "rank_b","tri","m","n","nme","idis"
       do i=0,np_band-1
          write(*,'(1x,6i10)') i,mat_block(i,0:4)
       end do
    end if

  END SUBROUTINE prep_subspace_diag

END MODULE subspace_diag_module
