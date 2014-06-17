MODULE vloc_rho_mol_module

  use bberf_module
  use rgrid_mol_module
  use atom_module
  use pseudopot_module
  use bc_module
  use array_bound_module, only: ML_0,ML_1
  use fd_module

  implicit none

  PRIVATE
  PUBLIC :: calc_vloc_rho_mol

  real(8),allocatable :: vloc_b(:)
  real(8),allocatable :: rho_vloc(:)

  logical :: flag_init = .true.
  integer :: Md
  real(8),allocatable :: lap(:)

CONTAINS


  SUBROUTINE calc_vloc_rho_mol(Vion,Md_in)
    implicit none
    real(8),intent(IN) :: Vion(:)
    integer,optional,intent(IN) :: Md_in
    integer :: a,ik,m1,m2
    integer :: i,i1,i2,i3,j
    real(8) :: p1,p2,p3,p4,c
    real(8) :: Rx,Ry,Rz,x,y,z,r

    if ( flag_init ) then
       if ( present(Md_in) ) then
          Md = Md_in
          allocate( lap(-Md:Md) ) ; lap=0.0d0
          call get_coef_lapla_fd(Md,lap)
       else
          stop "stop@calc_vloc_rho_mol(1)"
       end if
       flag_init = .false.
    end if

    m1 = 1
    m2 = size(KK,2)

    allocate( vloc_b(m1:m2) ) ; vloc_b(:) = 0.0d0

    do a=1,Natom

       ik = ki_atom(a)

       p1 = -Zps(ik)*parloc(1,ik)
       p2 = sqrt( parloc(2,ik) )
       p3 = -Zps(ik)*parloc(3,ik)
       p4 = sqrt( parloc(4,ik) )

       Rx = aa_atom(1,a)
       Ry = aa_atom(2,a)
       Rz = aa_atom(3,a)

       do i=m1,m2

          i1 = KK(1,i)
          i2 = KK(2,i)
          i2 = KK(3,i)

          x = i1*Hsize - Rx
          y = i2*Hsize - Ry
          z = i3*Hsize - Rz
          r = sqrt( x*x + y*y + z*z )

          vloc_b(i) = vloc_b(i) + ( p1*bberf(p2*r) + p3*bberf(p4*r) )/r

       end do ! i

    end do ! a

! --

    www(:,:,:,:) = 0.0d0
    do i=ML_0,ML_1
       www(LL(1,i),LL(2,i),LL(3,i),1) = Vion(i)
    end do

    call bcset(1,1,Md,0)

    do i=m1,m2
       www(KK(1,i),KK(2,i),KK(3,i),1) = vloc_b(i)
    end do

    c = 3.0d0*lap(0)/Hsize**2
    rho_vloc(:) = c*Vion(:)

    do j=1,Md
       c = lap(j)/Hsize**2
       do i=ML_0,ML_1
          i1 = LL(1,i)
          i2 = LL(2,i)
          i3 = LL(3,i)
          rho_vloc(i) = rho_vloc(i) &
               + c*(  www(i1-j,i2  ,i3  ,1) + www(i1+j,i2  ,i3  ,1) &
                    + www(i1  ,i2-j,i3  ,1) + www(i1  ,i2+j,i3  ,1) &
                    + www(i1  ,i2  ,i3-j,1) + www(i1  ,i2  ,i3+j,1) )
       end do
    end do

    c = -1.0d0/( 4.0d0*acos(-1.0d0) )
    rho_vloc(:) = c*rho_vloc(:)

    www(:,:,:,:) = 0.0d0

  END SUBROUTINE calc_vloc_rho_mol


END MODULE vloc_rho_mol_module
