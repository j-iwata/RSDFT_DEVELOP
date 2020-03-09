module esp_calc_module

  use hamiltonian_module
  use rgrid_module, only: dV
  use parallel_module
  use wf_module, only: hunk, iflag_hunk
  use rsdft_mpi_module
  use esp_calc_ncol_module, only: flag_noncollinear, esp_calc_ncol

  implicit none

  private
  public :: esp_calc

contains

  subroutine esp_calc(k,s,n1,n2,ns,ne,wf,e)
    implicit none
    integer,intent(in) :: k,s,n1,n2,ns,ne
#ifdef _DRSDFT_
    real(8),intent(in) :: wf(:,:,:,:)
    real(8),allocatable :: hwf(:,:)
#else
    complex(8),intent(in) :: wf(:,:,:,:)
    complex(8),allocatable :: hwf(:,:)
#endif
    real(8),intent(out) :: e(:,:,:)
    integer :: n,ierr,k0,s0

    if ( flag_noncollinear ) then
       call esp_calc_ncol( k,n1,n2,wf,e )
       return
    end if

    k0 = k - id_bzsm(myrank_k)
    s0 = s - id_spin(myrank_s)

    e(:,k,s)=0.0d0

    allocate( hwf(n1:n2,1) ) ; hwf=0.0d0

    do n=ns,ne
       if ( iflag_hunk >= 1 ) then
          hwf(:,1)=hunk(:,n,k,s)
       else
          call hamiltonian(k,s,wf(:,n:n,k0,s0),hwf,n1,n2,n,n)
       end if
#ifdef _DRSDFT_
       e(n,k,s) = sum( wf(:,n,k0,s0)*hwf(:,1) )*dV
#else
       e(n,k,s) = sum( conjg(wf(:,n,k0,s0))*hwf(:,1) )*dV
#endif
    end do

    deallocate( hwf )

    call rsdft_allreduce_sum( e(ns:ne,k,s), comm_grid )

  end subroutine esp_calc

end module esp_calc_module
