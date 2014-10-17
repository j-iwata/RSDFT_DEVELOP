MODULE TestModule
  use parallel_module
  implicit none
CONTAINS
  SUBROUTINE export_DensityAndWF
    implicit none
  END SUBROUTINE export_DensityAndWF
!----------------------------------------------
  SUBROUTINE import_DensityAndWF
    use wf_module
    use density_module
    use localpot_module
    implicit none
    integer :: s,k,n,i
  END SUBROUTINE import_DensityAndWF
END MODULE TestModule
