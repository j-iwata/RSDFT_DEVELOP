PROGRAM aa2xyz

  implicit none

  real(8) :: ax,aa(3,3),aa_tmp(3,3),aa_inv(3,3)
  character(len=8) :: cbuf
  integer,parameter :: u0=5, u1=1
  integer :: natom, nelem
  real(8),allocatable :: aa_atom(:,:), xyz_atom(:,:)
  integer,allocatable :: ki_atom(:), mi_atom(:), zi_atom(:)
  integer :: i

  rewind u0

  read(u0,*) cbuf, ax
  read(u0,*) cbuf, aa(1:3,1)
  read(u0,*) cbuf, aa(1:3,2)
  read(u0,*) cbuf, aa(1:3,3)
  read(u0,*) cbuf

  read(u0,*) nelem, natom

  allocate( xyz_atom(3,natom) ) ; xyz_atom=0.0d0
  allocate( aa_atom(3,natom)  ) ; aa_atom=0.0d0
  allocate( ki_atom(natom)    ) ; ki_atom=0
  allocate( mi_atom(natom)    ) ; mi_atom=0
  allocate( zi_atom(nelem)    ) ; zi_atom=0

  backspace(u0)

  read(u0,*) nelem, natom, zi_atom(1:nelem)
  do i=1,natom
     read(u0,*) ki_atom(i), aa_atom(:,i), mi_atom(i)
  end do

  aa_tmp = ax*aa
  call get_inverse_lattice( aa_tmp, aa_inv )

  if ( cbuf == "XYZ" ) then
     xyz_atom = aa_atom
     aa_atom  = matmul( aa_inv, xyz_atom )
  else if ( cbuf == "AA" ) then
     xyz_atom = matmul( aa_tmp, aa_atom )
  else
     stop "error!"
  end if

  open(u1,file="atomic_coordinates_aa")
  write(u1,'("AX",1f20.15)') ax
  write(u1,'("A1",3f20.15)') aa(:,1)
  write(u1,'("A2",3f20.15)') aa(:,2)
  write(u1,'("A3",3f20.15)') aa(:,3)
  write(u1,'("AA")')
  write(u1,*) nelem, natom, zi_atom(1:nelem), " /"
  do i=1,natom
     write(u1,'(1x,i5,3f20.12,i4)') ki_atom(i), aa_atom(:,i), mi_atom(i)
  end do
  close(u1)

  open(u1,file="atomic_coordinates_xyz")
  write(u1,'("AX",1f20.15)') ax
  write(u1,'("A1",3f20.15)') aa(:,1)
  write(u1,'("A2",3f20.15)') aa(:,2)
  write(u1,'("A3",3f20.15)') aa(:,3)
  write(u1,'("XYZ")')
  write(u1,*) nelem, natom, zi_atom(1:nelem), " /"
  do i=1,natom
     write(u1,'(1x,i5,3f20.12,i4)') ki_atom(i), xyz_atom(:,i), mi_atom(i)
  end do
  close(u1)

END PROGRAM aa2xyz


SUBROUTINE get_inverse_lattice( a, ainv )
  implicit none
  real(8),intent(IN)  :: a(3,3)
  real(8),intent(OUT) :: ainv(3,3)
  real(8) :: b(3,3), v
  b(1,1) = a(2,2)*a(3,3) - a(3,2)*a(2,3)
  b(2,1) = a(3,2)*a(1,3) - a(1,2)*a(3,3)
  b(3,1) = a(1,2)*a(2,3) - a(2,2)*a(1,3)
  b(1,2) = a(2,3)*a(3,1) - a(3,3)*a(2,1)
  b(2,2) = a(3,3)*a(1,1) - a(1,3)*a(3,1)
  b(3,2) = a(1,3)*a(2,1) - a(2,3)*a(1,1)
  b(1,3) = a(2,1)*a(3,2) - a(3,1)*a(2,2)
  b(2,3) = a(3,1)*a(1,2) - a(1,1)*a(3,2)
  b(3,3) = a(1,1)*a(2,2) - a(2,1)*a(1,2)
  call calc_Volume( a, v )
  ainv(:,:) = transpose( b(:,:) )/v
END SUBROUTINE get_inverse_lattice


SUBROUTINE calc_Volume( aa, Va )
  implicit none
  real(8),intent(IN)  :: aa(3,3)
  real(8),intent(OUT) :: Va
  Va = aa(1,1)*aa(2,2)*aa(3,3)+aa(1,2)*aa(2,3)*aa(3,1) &
      +aa(1,3)*aa(2,1)*aa(3,2)-aa(1,3)*aa(2,2)*aa(3,1) &
      -aa(1,2)*aa(2,1)*aa(3,3)-aa(1,1)*aa(2,3)*aa(3,2)
END SUBROUTINE calc_Volume
