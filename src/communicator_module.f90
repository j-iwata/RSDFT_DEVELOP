module communicator_module

  use parallel_module, only: comm_grid, comm_band, comm_bzsm, comm_spin

  private
  public :: set_communicator

contains

  integer function set_communicator( comm_in )
    implicit none
    character(*),intent(in) :: comm_in
    character(9) :: comm
    comm=comm_in
    call convertToLowercase( comm )
    select case( comm )
    case( 'g', 'grid', 'comm_grid' ); set_communicator = comm_grid
    case( 'b', 'band', 'comm_band' ); set_communicator = comm_band
    case( 'k', 'bzsm', 'comm_bzsm' ); set_communicator = comm_bzsm
    case( 's', 'spin', 'comm_spin' ); set_communicator = comm_spin
    case default
      write(*,*) "comm= ",comm," is not defined."
      call stop_program('set_communicator@communicator_module')
    end select
  end function set_communicator

  subroutine convertToLowercase( cbuf, ckey )
    implicit none
    character(*),intent(inout) :: cbuf
    character(*),optional,intent(out) :: ckey
    integer :: j,k,n
    n=len_trim(cbuf)
    if ( present(ckey) ) ckey=cbuf(1:n)
    do j=1,n
      k=iachar( cbuf(j:j) )
      if ( 65 <= k .and. k <= 90 ) k=k+32
      if ( present(ckey) ) then
        ckey(j:j) = achar(k)
      else
        cbuf(j:j) = achar(k)
      end if
    end do
  end subroutine convertToLowercase

end module communicator_module
