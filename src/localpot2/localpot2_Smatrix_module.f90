MODULE localpot2_Smatrix_module

  use localpot2_variables, only: Igrid_dense
  use localpot2_module, only: MLpot, Lpot, vloc_nl, test2_localpot2 &
                            , flag_localpot2
  use parallel_module

  implicit none

  PRIVATE
  PUBLIC :: init_localpot2_Smatrix, op_localpot2_Smatrix

  real(8),allocatable :: Smat_lpot2(:,:)
  integer :: ML, ML_0, ML_1

CONTAINS


  SUBROUTINE init_localpot2_Smatrix( ML_in, ML_0_in, ML_1_in )
    implicit none
    integer,intent(IN) :: ML_in, ML_0_in, ML_1_in
    integer :: m1,m2,m3
    real(8),allocatable :: vloc_nl_bak(:,:), vpot(:,:,:)

    if ( .not.flag_localpot2 ) return

    ML   = ML_in
    ML_0 = ML_0_in
    ML_1 = ML_1_in

    allocate( Smat_lpot2(MLpot,ML_0:ML_1) ) ; Smat_lpot2=0.0d0

    allocate( vloc_nl_bak(MLpot,ML_0:ML_1) ) ; vloc_nl_bak=0.0d0
    vloc_nl_bak(:,:) = vloc_nl(:,:)

    m1 = Igrid_dense(2,1) - Igrid_dense(1,1) + 1
    m2 = Igrid_dense(2,2) - Igrid_dense(1,2) + 1
    m3 = Igrid_dense(2,3) - Igrid_dense(1,3) + 1
    allocate( vpot(m1,m2,m3) ) ; vpot=1.0d0

    call test2_localpot2( vpot )

    deallocate( vpot )

    Smat_lpot2(:,:) = vloc_nl(:,:)

    vloc_nl(:,:) = vloc_nl_bak(:,:)

    deallocate( vloc_nl_bak )

  END SUBROUTINE init_localpot2_Smatrix


  SUBROUTINE op_localpot2_Smatrix( f, Sf )
    implicit none
#ifdef _DRSDFT_
    real(8),intent(IN)  :: f(:)
    real(8),intent(OUT) :: Sf(:)
    real(8),allocatable :: ft(:)
    real(8),parameter :: zero=0.0d0
    real(8) :: c
    integer,parameter :: TYP=MPI_REAL8
#else
    complex(8),intent(IN)  :: f(:)
    complex(8),intent(OUT) :: Sf(:)
    complex(8),allocatable :: ft(:)
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    complex(8) :: c
    integer,parameter :: TYP=MPI_COMPLEX16
#endif
    integer :: mm,i,j,ierr

    mm = ML_1 - ML_0 + 1

    allocate( ft(ML) ) ; ft=zero
    call mpi_allgatherv(f,mm,TYP,ft,ir_grid,id_grid,TYP,comm_grid,ierr)

    do i=1,mm
       c=zero
       do j=1,MLpot
          c = c + Smat_lpot2(j,i+ML_0-1)*ft(Lpot(j,i+ML_0-1))
       end do
       Sf(i) = c
    end do

    deallocate( ft )

  END SUBROUTINE op_localpot2_Smatrix


END MODULE localpot2_Smatrix_module
