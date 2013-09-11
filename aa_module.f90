MODULE aa_module

  implicit none

  PRIVATE
  PUBLIC :: ax,aa,Va,read_aa,send_aa,construct_aa

  real(8) :: ax,Va
  real(8) :: aa(3,3)

CONTAINS

  SUBROUTINE read_aa(unit)
    integer,intent(IN) :: unit
    read(unit,*) ax
    read(unit,*) aa(1:3,1)
    read(unit,*) aa(1:3,2)
    read(unit,*) aa(1:3,3)
    write(*,*) "ax=",ax
    write(*,'(1x,"aa(1:3,1)=",3F20.15)') aa(:,1)
    write(*,'(1x,"aa(1:3,2)=",3F20.15)') aa(:,2)
    write(*,'(1x,"aa(1:3,3)=",3F20.15)') aa(:,3)
  END SUBROUTINE read_aa


  SUBROUTINE send_aa(rank)
    integer,intent(IN) :: rank
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(ax,1,MPI_REAL8,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(aa,9,MPI_REAL8,rank,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_aa


  SUBROUTINE construct_aa
    aa(:,:)=ax*aa(:,:)
    Va = aa(1,1)*aa(2,2)*aa(3,3)+aa(1,2)*aa(2,3)*aa(3,1) &
        +aa(1,3)*aa(2,1)*aa(3,2)-aa(1,3)*aa(2,2)*aa(3,1) &
        -aa(1,2)*aa(2,1)*aa(3,3)-aa(1,1)*aa(2,3)*aa(3,2)
  END SUBROUTINE construct_aa

END MODULE aa_module
