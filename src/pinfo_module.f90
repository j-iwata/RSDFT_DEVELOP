module pinfo_module

  implicit none
  private
  public :: get_world_rank_pinfo

  integer,allocatable :: node_partition(:)
  integer :: irank,irank_world
  logical :: init_done=.false.
  integer :: icheck(4)
  
contains

  subroutine init
    use io_tools_module, only: IOTools_parseInteger
    use parallel_module, only: np=>node_partition
    implicit none
    if ( .not.allocated(node_partition) ) then
      allocate( node_partition(size(np)) ); node_partition=1
    end if
    node_partition=np
    init_done = .true.
  end subroutine init
  
  integer function get_world_rank_pinfo( g, b, k, s )
    implicit none
    integer,intent(in) :: g, b, k, s
    integer,allocatable :: itmp(:)
    integer :: n
    if ( .not.init_done ) call init
    n=size(node_partition)
    allocate( itmp(n) ); itmp=0
    irank=-1
    irank_world=-1
    icheck = (/ g, b, k, s/)
    call do_loop( n, node_partition, itmp )
    get_world_rank_pinfo = irank_world
  end function get_world_rank_pinfo

  recursive subroutine do_loop( n, nnn, itmp )
    implicit none
    integer,intent(in) :: n
    integer,intent(in) :: nnn(:)
    integer,intent(inout) :: itmp(:)
    integer :: i,irank_gbks(4)
    do i = 0, nnn(n)-1
      itmp(n)=i
      if ( n > 1 ) then
        call do_loop( n-1, nnn, itmp )
      else
        irank = irank + 1
        irank_gbks(1) = itmp(1) + itmp(2)*nnn(1) + itmp(3)*product(nnn(1:2))
        irank_gbks(2:4) = itmp(4:6)
        if ( all(irank_gbks == icheck) ) irank_world = irank
      end if
    end do
  end subroutine do_loop

end module pinfo_module