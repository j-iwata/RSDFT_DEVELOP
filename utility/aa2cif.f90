program aa2cif_2

  implicit none

  real(8),parameter :: aB=0.529177d0 !aB=0.529177210903d0
  real(8) :: ax,aa(3,3),aa_tmp(3,3),aa_inv(3,3),al(3)
  real(8) :: alpha, beta, gamma, c
  character(len=8) :: cbuf
  integer,parameter :: u0=5, u1=1
  integer :: natom, nelem
  real(8),allocatable :: aa_atom(:,:), xyz_atom(:,:)
  integer,allocatable :: ki_atom(:), mi_atom(:), zi_atom(:)
  integer :: i
  character(100) :: cbuf1,cbuf2
  character(2),external :: element_name

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

! ---

  al(1)=sqrt(sum(aa_tmp(:,1)**2))
  al(2)=sqrt(sum(aa_tmp(:,2)**2))
  al(3)=sqrt(sum(aa_tmp(:,3)**2))

  open(u1,file="atomic_coordinates.cif")

  write(u1,'("data_RSDFT")')

  write(u1,'("_cell_length_a   ",1x,es15.8)') al(1)*aB
  write(u1,'("_cell_length_b   ",1x,es15.8)') al(2)*aB
  write(u1,'("_cell_length_c   ",1x,es15.8)') al(3)*aB

  c = 180.0d0/acos(-1.0d0)
  alpha = acos(sum(aa_tmp(:,2)*aa_tmp(:,3) )/( al(2)*al(3)))*c
  beta  = acos(sum(aa_tmp(:,3)*aa_tmp(:,1) )/( al(3)*al(1)))*c
  gamma = acos(sum(aa_tmp(:,1)*aa_tmp(:,2) )/( al(1)*al(2)))*c

  write(u1,'("_cell_angle_alpha",1x,es15.8)') alpha
  write(u1,'("_cell_angle_beta ",1x,es15.8)') beta
  write(u1,'("_cell_angle_gamma",1x,es15.8)') gamma

  write(u1,'("loop_")')
  write(u1,'(" _symmetry_equiv_pos_site_id")')
  write(u1,'(" _symmetry_equiv_pos_as_xyz")')
  write(u1,'(a)') " 1 'x, y, z'"
  write(u1,'("loop_")')
  write(u1,'(" _atom_site_type_symbol")')
  write(u1,'(" _atom_site_label")')
  write(u1,'(" _atom_site_fract_x")')
  write(u1,'(" _atom_site_fract_y")')
  write(u1,'(" _atom_site_fract_z")')

  do i=1,natom
     cbuf1 = element_name(zi_atom(ki_atom(i)))
     write(u1,'(a2,4x,a2,3x,3f11.6)') cbuf1,cbuf1,aa_atom(:,i)
  end do

  close(u1)

end program aa2cif_2


subroutine get_inverse_lattice( a, ainv )
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
end subroutine get_inverse_lattice


subroutine calc_Volume( aa, Va )
  implicit none
  real(8),intent(IN)  :: aa(3,3)
  real(8),intent(OUT) :: Va
  Va = aa(1,1)*aa(2,2)*aa(3,3)+aa(1,2)*aa(2,3)*aa(3,1) &
      +aa(1,3)*aa(2,1)*aa(3,2)-aa(1,3)*aa(2,2)*aa(3,1) &
      -aa(1,2)*aa(2,1)*aa(3,3)-aa(1,1)*aa(2,3)*aa(3,2)
end subroutine calc_Volume


function element_name( z )

  implicit none
  integer,intent(in) :: z
  character(2) :: element_name
  character(2),save :: name(112)
!  character(2) :: a
!  integer :: i

  data name/ "H" ,"He", &
             "Li","Be","B" ,"C" ,"N" ,"O" ,"F" ,"Ne", &
             "Na","Mg","Al","Si","P" ,"S" ,"Cl","Ar", &
             "K" ,"Ca","Sc","Ti","V" ,"Cr","Mn","Fe","Co", &
             "Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr", &
             "Rb","Sr","Y" ,"Zr","Nb","Mo","Tc","Ru","Rh", &
             "Pd","Ag","Cd","In","Sn","Sb","Te","I" ,"Xe", &
             "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu", &
             "Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu", &
             "Hf","Ta","W" ,"Re","Os","Ir", &
             "Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn", &
             "Fr","Ra","Ac","Th","Pa","U" ,"Np","Pu","Am", &
             "Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf", &
             "Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn" /

  element_name = name(z)

!  z=0
!  do i=1,size(name)
!     if ( element_name(1:2) == name(i) ) then
!        z=i
!        exit
!     else
!        a=name(i)
!        if ( element_name(1:1) == a(1:len_trim(a)) ) z=i
!     end if
!  end do

end function element_name

