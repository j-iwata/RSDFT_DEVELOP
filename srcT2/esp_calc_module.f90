MODULE esp_calc_module

  use hamiltonian_module
  use rgrid_module, only: dV
  use parallel_module
  use wf_module, only: hunk, iflag_hunk
  use rsdft_mpi_module
  use esp_calc_ncol_module, only: flag_noncollinear, esp_calc_ncol

  implicit none

  PRIVATE
  PUBLIC :: esp_calc

CONTAINS

  SUBROUTINE esp_calc(k,s,n1,n2,ns,ne,wf,e)
    implicit none
    integer,intent(IN) :: k,s,n1,n2,ns,ne
#ifdef _DRSDFT_
    real(8),intent(IN) :: wf(:,:,:,:)
    real(8),allocatable :: hwf(:,:)
#else
    complex(8),intent(IN) :: wf(:,:,:,:)
    complex(8),allocatable :: hwf(:,:)
#endif
    real(8),intent(OUT) :: e(:,:,:)
    integer :: n,ierr,k0,m

    if ( flag_noncollinear ) then
       call esp_calc_ncol( k,n1,n2,wf,e )
       return
    end if

    k0 = k - id_bzsm(myrank_k)

    e(:,k,s)=0.0d0

    allocate( hwf(n1:n2,1) ) ; hwf=0.0d0

    do n=ns,ne
       m=n-ns+1
       if ( iflag_hunk >= 1 ) then
          hwf(:,1)=hunk(:,n,k,s)
       else
          call hamiltonian(k,s,wf(:,m:m,k0,s),hwf,n1,n2,n,n)
       end if
#ifdef _DRSDFT_
       e(n,k,s) = sum( wf(:,m,k0,s)*hwf(:,1) )*dV
#else
       e(n,k,s) = sum( conjg(wf(:,m,k0,s))*hwf(:,1) )*dV
#endif
    end do

    deallocate( hwf )

    call rsdft_allreduce_sum( e(ns:ne,k,s), comm_grid )

  END SUBROUTINE esp_calc

END MODULE esp_calc_module
