MODULE kinetic_variables

  use io_tools_module

  implicit none

  PRIVATE
  PUBLIC :: coef_lap0, coef_lap, coef_nab &
           ,coef_nabk, const_k2, zcoef_kin, ggg &
           ,flag_nab, flag_n12, flag_n23, flag_n31
  PUBLIC :: wk
  PUBLIC :: read_kinetic
  PUBLIC :: coef_kin

  integer,PUBLIC :: kin_select=0
  integer,PUBLIC :: SYStype=0
  integer,PUBLIC :: Md=6

  logical :: flag_nab, flag_n12, flag_n23, flag_n31
  real(8) :: coef_lap0, ggg(6)
  real(8),allocatable :: coef_lap(:,:), coef_nab(:,:)
  real(8),allocatable :: coef_nabk(:,:,:), const_k2(:)
  real(8),allocatable :: coef_kin(:)
  complex(8),allocatable :: zcoef_kin(:,:,:)
#ifdef _DRSDFT_
  real(8),allocatable :: wk(:,:,:,:)
#else
  complex(8),allocatable :: wk(:,:,:,:)
#endif

CONTAINS

  SUBROUTINE read_kinetic
    implicit none
    call IOTools_readIntegerKeyword( "MD", Md )
!    call IOTools_readIntegerKeyword( "SYSTYPE", SYStype )
    call IOTools_readIntegerKeyword( "KINSELECT", kin_select )
  END SUBROUTINE read_kinetic

END MODULE kinetic_variables
