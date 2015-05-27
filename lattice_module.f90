MODULE lattice_module

  implicit none

  PRIVATE
  PUBLIC :: lattice, init_lattice, get_reciprocal_lattice &
       , construct_aa_lattice

  type lattice
     real(8) :: LatticeConstant
     real(8) :: LatticeVector(3,3)
     real(8) :: Length(3)
     real(8) :: Volume
  end type lattice

  real(8) :: ax, av(3,3)

CONTAINS


  SUBROUTINE init_lattice( ax_in, av_in )
    implicit none
    real(8),intent(IN) :: ax_in, av_in(3,3)
    ax = ax_in
    av = av_in
  END SUBROUTINE init_lattice


  SUBROUTINE construct_aa_lattice( aa )

    implicit none
    type(lattice),intent(OUT) :: aa

    aa%LatticeConstant = ax
    aa%LatticeVector = av

    aa%LatticeVector(:,:) = aa%LatticeConstant * aa%LatticeVector(:,:)

    call calc_VectorLength( aa%LatticeVector, aa%Length )

    call calc_Volume( aa%LatticeVector, aa%Volume )

  END SUBROUTINE construct_aa_lattice


  SUBROUTINE calc_VectorLength( aa, aal )
    implicit none
    real(8),intent(IN)  :: aa(3,3)
    real(8),intent(OUT) :: aal(3)
    aal(1) = sqrt( sum(aa(1:3,1)**2) )
    aal(2) = sqrt( sum(aa(1:3,2)**2) )
    aal(3) = sqrt( sum(aa(1:3,3)**2) )
  END SUBROUTINE calc_VectorLength


  SUBROUTINE calc_Volume( aa, Va )
    implicit none
    real(8),intent(IN)  :: aa(3,3)
    real(8),intent(OUT) :: Va
    Va = aa(1,1)*aa(2,2)*aa(3,3)+aa(1,2)*aa(2,3)*aa(3,1) &
        +aa(1,3)*aa(2,1)*aa(3,2)-aa(1,3)*aa(2,2)*aa(3,1) &
        -aa(1,2)*aa(2,1)*aa(3,3)-aa(1,1)*aa(2,3)*aa(3,2)
  END SUBROUTINE calc_Volume


  SUBROUTINE get_reciprocal_lattice( aa, bb )
    implicit none
    type(lattice),intent(IN)  :: aa
    type(lattice),intent(OUT) :: bb
    real(8) :: a(3,3), b(3,3), pi2
    pi2 = 2.0d0*acos(-1.0d0)
    a(:,:) = aa%LatticeVector(:,:)
    b(1,1) = a(2,2)*a(3,3) - a(3,2)*a(2,3)
    b(2,1) = a(3,2)*a(1,3) - a(1,2)*a(3,3)
    b(3,1) = a(1,2)*a(2,3) - a(2,2)*a(1,3)
    b(1,2) = a(2,3)*a(3,1) - a(3,3)*a(2,1)
    b(2,2) = a(3,3)*a(1,1) - a(1,3)*a(3,1)
    b(3,2) = a(1,3)*a(2,1) - a(2,3)*a(1,1)
    b(1,3) = a(2,1)*a(3,2) - a(3,1)*a(2,2)
    b(2,3) = a(3,1)*a(1,2) - a(1,1)*a(3,2)
    b(3,3) = a(1,1)*a(2,2) - a(2,1)*a(1,2)
    b(:,:) = b(:,:)*pi2/aa%Volume
    bb%LatticeVector(:,:) = b(:,:)
    call calc_VectorLength( bb%LatticeVector, bb%Length )
    call calc_Volume( bb%LatticeVector, bb%Volume )
    bb%LatticeConstant = pi2/aa%LatticeConstant
 END SUBROUTINE get_reciprocal_lattice


END MODULE lattice_module
