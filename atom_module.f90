MODULE atom_module

  implicit none

  PRIVATE
  PUBLIC :: Natom,Nelement,aa_atom,ki_atom,read_atom,send_atom

  integer :: Natom,Nelement
  integer,allocatable :: ki_atom(:)
  real(8),allocatable :: aa_atom(:,:)

CONTAINS

  SUBROUTINE read_atom(unit)
    integer,intent(IN) :: unit
    integer :: i
    read(unit,*) Nelement,Natom
    allocate( aa_atom(3,Natom) )
    allocate( ki_atom(Natom) )
    do i=1,Natom
       read(unit,*) ki_atom(i),aa_atom(1:3,i)
    end do
    write(*,*) "Nelment,Natom=",Nelement,Natom
!    do i=1,Natom
!       write(*,'(1x,i5,i5,3f18.12,i5)') i,ki_atom(i),aa_atom(:,i)
!    end do
    write(*,'(8x,a7,3a18)') "ki_atom","aa_atom1","aa_atom2","aa_atom3"
    do i=1,min(5,Natom)
       write(*,'(1x,i5,2x,i7,3f18.12)') i,ki_atom(i),aa_atom(:,i)
    end do
    if ( Natom>5 ) then
       write(*,'(1x,10x,".")')
       write(*,'(1x,10x,".")')
       write(*,'(1x,10x,".")')
       do i=Natom-5,Natom
          write(*,'(1x,i5,2x,i7,3f18.12)') i,ki_atom(i),aa_atom(:,i)
       end do
    end if
  END SUBROUTINE read_atom

  SUBROUTINE send_atom(myrank)
    integer,intent(IN) :: myrank
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(Natom,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Nelement,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if  ( myrank /= 0 ) then
       allocate( aa_atom(3,Natom) )
       allocate( ki_atom(Natom) )
    end if
    call mpi_bcast(ki_atom,Natom,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(aa_atom,3*Natom,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_atom

END MODULE atom_module
