MODULE kinetic_mol_module

  use rgrid_variables, only: Igrid
  use rgrid_mol_module, only: LL
  use bc_module, only: www, bcset, bcset_1, bcset_3
  use kinetic_variables, only: coef_lap0, coef_lap, coef_kin, Md

  implicit none

  PRIVATE
  PUBLIC :: op_kinetic_mol

CONTAINS

  SUBROUTINE op_kinetic_mol( tpsi, htpsi )
    implicit none
#ifdef _DRSDFT_
    real(8),intent(IN)    ::  tpsi(:,:)
    real(8),intent(INOUT) :: htpsi(:,:)
    real(8) :: d
#else
    complex(8),intent(IN)    ::  tpsi(:,:)
    complex(8),intent(INOUT) :: htpsi(:,:)
    complex(8) :: d
#endif
    integer :: i,ib,i1,i2,i3,m,n,nb,ml,n1,n2
    real(8) :: c

    ml = size( tpsi, 1 )
    nb = size( tpsi, 2 )
    n1 = Igrid(1,0)
    n2 = n1+ml-1

    do ib=1,nb
!$OMP do
       do i=n1,n2
          www( LL(1,i),LL(2,i),LL(3,i),ib ) = tpsi(i-n1+1,ib)
       end do
!$OMP end do
    end do

    call bcset_3(1,nb,Md,0)
!$OMP barrier

    do ib=1,nb
!$OMP do
       do i=1,ml
          htpsi(i,ib) = htpsi(i,ib) + coef_lap0*tpsi(i,ib)
       end do
!$OMP end do
    end do

    do ib=1,nb
!$OMP do
       do i=n1,n2
          i1 = LL(1,i)
          i2 = LL(2,i)
          i3 = LL(3,i)
          d  = htpsi(i-n1+1,ib)
          do m=1,Md
             c = coef_kin(m)
             d = d + c*( www(i1-m,i2,i3,n)+www(i1+m,i2,i3,n) &
                       + www(i1,i2-m,i3,n)+www(i1,i2+m,i3,n) &
                       + www(i1,i2,i3-m,n)+www(i1,i2,i3+m,n) )
          end do
          htpsi(i-n1+1,ib) = d
       end do
!$OMP end do
    end do ! ib

  END SUBROUTINE op_kinetic_mol

END MODULE kinetic_mol_module
