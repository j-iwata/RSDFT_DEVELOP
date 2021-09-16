MODULE esp_calc_ncol_module

  use hamiltonian_module
  use hamiltonian_ncol_module
  use rgrid_module, only: dV
  use parallel_module, only: comm_grid, id_bzsm, myrank_k
  use rsdft_mpi_module
  use noncollinear_module, only: flag_noncollinear

  implicit none

  PRIVATE
  PUBLIC :: esp_calc_ncol
  PUBLIC :: flag_noncollinear

CONTAINS


  SUBROUTINE esp_calc_ncol( k, n1,n2, wf, e )

    implicit none
    integer,intent(IN) :: k,n1,n2
#ifdef _DRSDFT_
    real(8),intent(IN) :: wf(:,:,:,:)
#else
    complex(8),intent(IN) :: wf(:,:,:,:)
#endif
    real(8),intent(INOUT) :: e(:,:,:)
    complex(8),allocatable :: hwf(:,:,:)
    complex(8),allocatable :: zw1(:,:,:),zw2(:,:,:)
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    integer :: n,s,MS,MB,ierr,k0

    call write_border( 1, "esp_calc_ncol(start)" )

    MB = size( wf, 2 )
    MS = size( wf, 4 )
    k0 = k-id_bzsm(myrank_k)

    e(:,k,:) = 0.0d0

    allocate( hwf(n1:n2,1,MS) ) ; hwf=zero
    allocate( zw1(n1:n2,1,MS) ) ; zw1=zero
    allocate( zw2(n1:n2,1,MS) ) ; zw2=zero

#ifndef _DRSDFT_
    do n=1,MB
       do s=1,MS
          call hamiltonian( wf(:,n:n,k0,s), hwf(:,:,s), n,k,s )
       end do
       zw1(:,1,:)=wf(:,n,k0,:)
       zw2(:,1,:)=hwf(:,1,:)
       call hamiltonian_ncol( k, n1,n2, zw1, zw2 )
       hwf(:,1,:)=zw2(:,1,:)
       do s=1,MS
          e(n,k,1) = e(n,k,1) + sum( conjg(wf(:,n,k0,s))*hwf(:,1,s) )*dV
       end do
    end do
#endif

    deallocate( zw2 )
    deallocate( zw1 )
    deallocate( hwf )

    call rsdft_allreduce_sum( e(:,k,1), comm_grid )

    e(:,:,2) = e(:,:,1)

    call write_border( 1, "esp_calc_ncol(end)" )

  END SUBROUTINE esp_calc_ncol


END MODULE esp_calc_ncol_module
