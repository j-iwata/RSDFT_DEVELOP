MODULE fock_parallel_module

  use parallel_module

  implicit none

  PRIVATE
  PUBLIC :: comm_fock, myrank_f, nprocs_f, init_fock_parallel

  integer :: comm_fock, myrank_f, nprocs_f

CONTAINS


  SUBROUTINE init_fock_parallel

    implicit none
    integer :: i,i1,i2,m,ierr

    m=0
    i=0
    do i2=0,np_spin-1
    do i1=0,np_grid-1
       i=i+1
       if ( id_class(myrank,0)==i1 .and. id_class(myrank,6)==i2 ) m=i
    end do
    end do

    call mpi_comm_split(mpi_comm_world,m,myrank,comm_fock,ierr)
    call mpi_comm_rank(comm_fock,myrank_f,ierr)
    call mpi_comm_size(comm_fock,nprocs_f,ierr)

    if ( disp_switch_parallel ) then
       write(*,'(a20," init_fock_parallel")') repeat("-",20)
       write(*,*) "nprocs_f =",nprocs_f
       write(*,*) "comm_fock=",comm_fock
    end if

  END SUBROUTINE init_fock_parallel


END MODULE fock_parallel_module
