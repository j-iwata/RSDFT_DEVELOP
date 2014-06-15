MODULE atom_module

  implicit none

  PRIVATE
  PUBLIC :: Natom,Nelement,aa_atom,ki_atom,read_atom

  integer :: Natom,Nelement
  integer,allocatable :: ki_atom(:)
  real(8),allocatable :: aa_atom(:,:)

CONTAINS

  SUBROUTINE read_atom(rank,unit,ax,aa)
    implicit none
    integer,intent(IN) :: rank,unit
    real(8),intent(INOUT) :: ax,aa(3,3)
    integer :: i,iflag_format
    character(3) :: cbuf,ckey
    ax=0.0d0
    aa=0.0d0
    if ( rank == 0 ) then
       iflag_format = 0
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:2) == "AX" ) then
             backspace(unit)
             read(unit,*) cbuf,ax
             iflag_format=1
          else if ( ckey(1:2) == "A1" ) then
             backspace(unit)
             read(unit,*) cbuf,aa(1:3,1)
             iflag_format=1
          else if ( ckey(1:2) == "A2" ) then
             backspace(unit)
             read(unit,*) cbuf,aa(1:3,2)
             iflag_format=1
          else if ( ckey(1:2) == "A3" ) then
             backspace(unit)
             read(unit,*) cbuf,aa(1:3,3)
             iflag_format=1
          else if ( ckey(1:3) == "XYZ" ) then
             iflag_format=2
             exit
          else if ( ckey(1:2) == "AA" ) then
             exit
          end if
       end do
999    continue
       if ( iflag_format == 0 ) then
          rewind unit
       else
          write(*,*) "iflag_format=",iflag_format
          write(*,*) "ax=",ax
          write(*,'(1x,"a1=",3f20.15)') aa(1:3,1)
          write(*,'(1x,"a2=",3f20.15)') aa(1:3,2)
          write(*,'(1x,"a3=",3f20.15)') aa(1:3,3)
       end if
       read(unit,*) Nelement,Natom
       write(*,*) "Nelment,Natom=",Nelement,Natom
    end if
    call send_atom_1(0,ax,aa)
    allocate( aa_atom(3,Natom) ) ; aa_atom=0.d0
    allocate( ki_atom(Natom)   ) ; ki_atom=0
    if ( rank == 0 ) then
       do i=1,Natom
          read(unit,*) ki_atom(i),aa_atom(1:3,i)
       end do
       write(*,'(8x,a7,3a18)') "ki_atom","aa_atom1","aa_atom2","aa_atom3"
       if ( Natom <= 11 ) then
          do i=1,Natom
             write(*,'(1x,i5,2x,i7,3f18.12)') i,ki_atom(i),aa_atom(:,i)
          end do
       else
          do i=1,min(5,Natom)
             write(*,'(1x,i5,2x,i7,3f18.12)') i,ki_atom(i),aa_atom(:,i)
          end do
          write(*,'(1x,10x,".")')
          write(*,'(1x,10x,".")')
          write(*,'(1x,10x,".")')
          do i=Natom-5,Natom
             write(*,'(1x,i5,2x,i7,3f18.12)') i,ki_atom(i),aa_atom(:,i)
          end do
       end if
    end if
    call send_atom_2(0)
  END SUBROUTINE read_atom

  SUBROUTINE send_atom_1(myrank,ax,aa)
    implicit none
    integer,intent(IN) :: myrank
    real(8),intent(IN) :: ax,aa(3,3)
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(Natom,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Nelement,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(ax,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(aa,9,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_atom_1

  SUBROUTINE send_atom_2(myrank)
    implicit none
    integer,intent(IN) :: myrank
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(ki_atom,Natom,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(aa_atom,3*Natom,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_atom_2

END MODULE atom_module
