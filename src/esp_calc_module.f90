module esp_calc_module

  implicit none

  private
  public :: esp_calc

  interface esp_calc
    module procedure d_esp_calc, z_esp_calc
  end interface

  logical :: init_done=.false.
  real(8) :: dV=1.0d0

contains

  subroutine d_esp_calc( wf, e, n_in,k,s )
    use rsdft_allreduce_module, only: rsdft_allreduce
    ! use parallel_module, only: get_range_parallel
    use hamiltonian_module, only: hamiltonian
    use wf_module, only: hunk, iflag_hunk, USE_WORKWF_AT_ESPCAL
    implicit none
    integer,intent(in) :: n_in, k, s
    real(8),intent(in) :: wf(:,n_in:)
    real(8),intent(inout) :: e(n_in:)
    real(8),allocatable :: hwf(:,:)
    integer :: m,nn,ns,ne,n

    if ( .not.init_done ) call init( sum(abs(wf(:,n_in))**2) ) 

    ! call get_range_parallel( ns, ne, 'b' )

    nn = size(e)
    ns = n_in
    ne = ns + nn - 1

    e(:)=0.0d0

    allocate( hwf(size(wf,1),1) ); hwf=0.0d0

    if ( iflag_hunk >= 1 ) then

      if ( USE_WORKWF_AT_ESPCAL ) then
        do n = ns, ne
          e(n) = sum( wf(:,n)*hunk(:,n,k,s) )*dV
        end do
      else
        do n = ns, ne
          call hamiltonian( wf(:,n:n), hwf, n,k,s )
          hunk(:,n,k,s)=hwf(:,1)
          e(n) = sum( wf(:,n)*hwf(:,1) )*dV
        end do
      end if
    
    else

      do n = ns, ne
        call hamiltonian( wf(:,n:n), hwf, n,k,s )
        e(n) = sum( wf(:,n)*hwf(:,1) )*dV
      end do
    
    end if

    deallocate( hwf )

    call rsdft_allreduce( e, 'g' )

  end subroutine d_esp_calc
  
  
  subroutine z_esp_calc( wf, e, n_in,k,s )
    use rsdft_allreduce_module, only: rsdft_allreduce
    ! use parallel_module, only: get_range_parallel
    use hamiltonian_module, only: hamiltonian
    use wf_module, only: hunk, iflag_hunk, USE_WORKWF_AT_ESPCAL
    implicit none
    integer,intent(in) :: n_in, k, s
    complex(8),intent(in) :: wf(:,n_in:)
    real(8),intent(inout) :: e(n_in:)
    complex(8),allocatable :: hwf(:,:)
    integer :: n,m,nn,ns,ne

    if ( .not.init_done ) call init( sum(abs(wf(:,n_in))**2) ) 

    ! call get_range_parallel( ns, ne, 'b' )

    nn = size(e)
    ns = n_in
    ne = ns + nn - 1

    e(:)=0.0d0

    allocate( hwf(size(wf,1),1) ); hwf=(0.0d0,0.0d0)

    if ( iflag_hunk >= 1 ) then

      if ( USE_WORKWF_AT_ESPCAL ) then
        do n = ns, ne
          e(n) = sum( conjg(wf(:,n))*hunk(:,n,k,s) )*dV
        end do
      else
        do n = ns, ne
          call hamiltonian( wf(:,n:n), hwf(:,1:1), n,k,s )
          hunk(:,n,k,s)=hwf(:,1)
          e(n) = sum( conjg(wf(:,n))*hwf(:,1) )*dV
        end do
      end if

    else
  
      do n = ns, ne
        call hamiltonian( wf(:,n:n), hwf(:,1:1), n,k,s )
        e(n) = sum( conjg(wf(:,n))*hwf(:,1) )*dV
      end do

    end if

    deallocate( hwf )

    call rsdft_allreduce( e, 'g' )
      
  end subroutine z_esp_calc

  subroutine init( c )
    use rsdft_allreduce_module, only: rsdft_allreduce
    implicit none
    real(8),intent(in) :: c
    dV=c
    call rsdft_allreduce( dV, 'g' )
    dV=1.0d0/dV
  end subroutine init

end module esp_calc_module
