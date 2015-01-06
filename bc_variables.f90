MODULE bc_variables

  implicit none

  integer :: n_neighbor(6)
  integer,allocatable :: fdinfo_send(:,:,:),fdinfo_recv(:,:,:)

END MODULE bc_variables
