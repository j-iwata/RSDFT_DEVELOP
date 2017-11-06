MODULE current_module

  use wf_module
  use aa_module, only: aa, Va
  use lattice_module, only: get_inverse_lattice
  use momentum_module
  use parallel_module, only: comm_bzsm, comm_spin, comm_band
  use rsdft_mpi_module
  use symmetry_module

  implicit none

  PRIVATE
  PUBLIC :: calc_macro_current

CONTAINS


  SUBROUTINE calc_macro_current( jav, coordinates )

    implicit none
    real(8),intent(OUT)   :: jav(3)
    character(3),optional,intent(IN) :: coordinates
    integer :: n,k,s
    real(8) :: pxyz(3),aa_inv(3,3)

    jav(:)=0.0d0
    do s=MS_0_WF,MS_1_WF
    do k=MK_0_WF,MK_1_WF
    do n=MB_0_WF,MB_1_WF
       if ( abs(occ(n,k,s)) < 1.d-10 ) cycle
       pxyz=0.0d0
#ifndef _DRSDFT_
       call calc_expectval_momentum &
            (k,ML_0_WF,ML_1_WF,1,1,unk(ML_0_WF,n,k,s),pxyz)
#endif
       jav(:)=jav(:)+occ(n,k,s)*pxyz(:)
    end do
    end do
    end do
    jav=jav/Va
    if ( present(coordinates) ) then
       if ( coordinates == "xyz" ) then
       else
          call get_inverse_lattice( aa, aa_inv )
          pxyz=jav
          jav(:)=matmul( aa_inv, pxyz(:) )
       end if
    end if

    call rsdft_allreduce_sum( jav, comm_band )
    call rsdft_allreduce_sum( jav, comm_bzsm )
    call rsdft_allreduce_sum( jav, comm_spin )

    call sym_vector_xyz( jav )

  END SUBROUTINE calc_macro_current


END MODULE current_module
