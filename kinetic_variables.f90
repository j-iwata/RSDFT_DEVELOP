MODULE kinetic_variables

  implicit none

  PRIVATE
  PUBLIC :: Md, coef_lap0, coef_lap, coef_nab &
           ,coef_nabk, const_k2, zcoef_kin, ggg &
           ,flag_nab, flag_n12, flag_n23, flag_n31
  PUBLIC :: wk

  integer :: Md
  logical :: flag_nab, flag_n12, flag_n23, flag_n31
  real(8) :: coef_lap0, ggg(6)
  real(8),allocatable :: coef_lap(:,:), coef_nab(:,:)
  real(8),allocatable :: coef_nabk(:,:,:), const_k2(:)
  complex(8),allocatable :: zcoef_kin(:,:,:)
#ifdef _DRSDFT_
  real(8),allocatable :: wk(:,:,:,:)
#else
  complex(8),allocatable :: wk(:,:,:,:)
#endif

END MODULE kinetic_variables
