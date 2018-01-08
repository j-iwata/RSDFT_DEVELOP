MODULE lattice_module

  use io_tools_module

  implicit none

  PRIVATE
  PUBLIC :: lattice
  PUBLIC :: init_lattice
  PUBLIC :: get_reciprocal_lattice
  PUBLIC :: construct_lattice
  PUBLIC :: read_lattice
  PUBLIC :: get_inverse_lattice
  PUBLIC :: backup_aa_lattice
  PUBLIC :: get_aa_lattice
  PUBLIC :: check_lattice
  PUBLIC :: calc_volume_lattice

  type lattice
     real(8) :: LatticeConstant
     real(8) :: LatticeVector(3,3)
     real(8) :: Length(3)
     real(8) :: Volume
  end type lattice

  real(8) :: ax, av(3,3)

  type(lattice) :: aa_backup

CONTAINS


  SUBROUTINE read_lattice( aa, unit )
    implicit none
    type(lattice),intent(OUT) :: aa
    integer,optional,intent(IN) :: unit
    aa%LatticeConstant=0.0d0
    aa%LatticeVector(:,:)=0.0d0
    if ( present(unit) ) then
       call IOTools_readReal8Keyword( "AX", aa%LatticeConstant, unit )
       call IOTools_readReal8Keywords( "A1", aa%LatticeVector(:,1), unit )
       call IOTools_readReal8Keywords( "A2", aa%LatticeVector(:,2), unit )
       call IOTools_readReal8Keywords( "A3", aa%LatticeVector(:,3), unit )
    else
       call IOTools_readReal8Keyword( "AX", aa%LatticeConstant )
       call IOTools_readReal8Keywords( "A1", aa%LatticeVector(:,1) )
       call IOTools_readReal8Keywords( "A2", aa%LatticeVector(:,2) )
       call IOTools_readReal8Keywords( "A3", aa%LatticeVector(:,3) )
    end if
  END SUBROUTINE read_lattice


  SUBROUTINE init_lattice( ax_in, av_in )
    implicit none
    real(8),intent(IN) :: ax_in, av_in(3,3)
    ax = ax_in
    av = av_in
  END SUBROUTINE init_lattice


  SUBROUTINE construct_lattice( aa )
    implicit none
    type(lattice),intent(INOUT) :: aa
    aa%LatticeVector(:,:) = aa%LatticeConstant * aa%LatticeVector(:,:)
    call calc_VectorLength( aa%LatticeVector, aa%Length )
    call calc_volume_lattice( aa%LatticeVector, aa%Volume )
  END SUBROUTINE construct_lattice


  SUBROUTINE calc_VectorLength( aa, aal )
    implicit none
    real(8),intent(IN)  :: aa(3,3)
    real(8),intent(OUT) :: aal(3)
    aal(1) = sqrt( sum(aa(1:3,1)**2) )
    aal(2) = sqrt( sum(aa(1:3,2)**2) )
    aal(3) = sqrt( sum(aa(1:3,3)**2) )
  END SUBROUTINE calc_VectorLength


  SUBROUTINE calc_volume_lattice( aa, Va )
    implicit none
    real(8),intent(IN)  :: aa(3,3)
    real(8),intent(OUT) :: Va
    Va = aa(1,1)*aa(2,2)*aa(3,3)+aa(1,2)*aa(2,3)*aa(3,1) &
        +aa(1,3)*aa(2,1)*aa(3,2)-aa(1,3)*aa(2,2)*aa(3,1) &
        -aa(1,2)*aa(2,1)*aa(3,3)-aa(1,1)*aa(2,3)*aa(3,2)
  END SUBROUTINE calc_volume_lattice


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
    call calc_volume_lattice( bb%LatticeVector, bb%Volume )
    bb%LatticeConstant = pi2/aa%LatticeConstant
 END SUBROUTINE get_reciprocal_lattice


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
    call calc_volume_lattice( a, v )
    ainv(:,:) = transpose( b(:,:) )/v
  END SUBROUTINE get_inverse_lattice


  SUBROUTINE backup_aa_lattice( aa )
    implicit none
    type(lattice),intent(IN) :: aa
    aa_backup = aa
  END SUBROUTINE backup_aa_lattice

  SUBROUTINE get_aa_lattice( aa )
    implicit none
    type(lattice),intent(OUT) :: aa
    aa = aa_backup
  END SUBROUTINE get_aa_lattice


  SUBROUTINE check_lattice( a, indx )
    implicit none
    real(8),intent(IN) :: a(3,3)
    character(*),optional,intent(OUT) :: indx
    real(8) :: al(3), theta23, theta31, theta12, cos23, cos31, cos12
    logical :: disp
    al(1) = sqrt( sum(a(:,1)**2) )
    al(2) = sqrt( sum(a(:,2)**2) )
    al(3) = sqrt( sum(a(:,3)**2) )
    cos23 = sum( a(:,2)*a(:,3) )/(al(2)*al(3))
    cos31 = sum( a(:,3)*a(:,1) )/(al(3)*al(1))
    cos12 = sum( a(:,1)*a(:,2) )/(al(1)*al(2))
    theta23 = acos( cos23 )*180.0d0/acos(-1.0d0)
    theta31 = acos( cos31 )*180.0d0/acos(-1.0d0)
    theta12 = acos( cos12 )*180.0d0/acos(-1.0d0)
    call check_disp_switch( disp, 0 )
    if ( disp ) then
       write(*,'(1x,"V1",3f12.5,2x,f12.5)') a(:,1), al(1)
       write(*,'(1x,"V2",3f12.5,2x,f12.5)') a(:,2), al(2)
       write(*,'(1x,"V3",3f12.5,2x,f12.5)') a(:,3), al(3)
       write(*,'(1x,"cos23  ,cos31  ,cos12  :",3f10.5)') cos23,cos31,cos12
       write(*,'(1x,"theta23,theta31,theta12:",3f10.5)') theta23,theta31,theta12
    end if
    if ( present(indx) ) then
       indx=""
       call check_cubic( al, (/cos23,cos31,cos12/), indx )
       call check_hexagonal( al, (/cos23,cos31,cos12/), indx )
       call check_fcc( al, (/cos23,cos31,cos12/), indx )
    end if
  END SUBROUTINE check_lattice

  SUBROUTINE check_cubic( al, cosines, indx )
    implicit none
    real(8),intent(IN) :: al(3), cosines(3)
    character(*),intent(INOUT) :: indx
    real(8),parameter :: eps=1.d-3
    if ( abs(al(1)-al(2)) < eps ) then
       if ( abs(al(2)-al(3)) < eps ) then
          if ( all( abs(cosines) < eps ) ) indx = "CUBIC"
       end if
    end if
  END SUBROUTINE check_cubic

  SUBROUTINE check_hexagonal( al, cosines, indx )
    implicit none
    real(8),intent(IN) :: al(3), cosines(3)
    character(*),intent(INOUT) :: indx
    real(8),parameter :: eps=1.d-3
    if ( abs(al(1)-al(2)) < eps ) then
       if ( abs(cosines(3)-0.5d0)<eps .or. abs(cosines(3)+0.5d0)<eps ) then
          if ( abs(cosines(1))<eps .and. abs(cosines(2))<eps ) indx="HEXAGONAL"
       end if
    else if ( abs(al(2)-al(3)) < eps ) then
       if ( abs(cosines(1)-0.5d0)<eps .or. abs(cosines(1)+0.5d0)<eps ) then
          if ( abs(cosines(2))<eps .and. abs(cosines(3))<eps ) indx="HEXAGONAL"
       end if
    else if ( abs(al(3)-al(1)) < eps ) then
       if ( abs(cosines(2)-0.5d0)<eps .or. abs(cosines(2)+0.5d0)<eps ) then
          if ( abs(cosines(3))<eps .and. abs(cosines(1))<eps ) indx="HEXAGONAL"
       end if
    end if
  END SUBROUTINE check_hexagonal

  SUBROUTINE check_fcc( al, cosines, indx )
    implicit none
    real(8),intent(IN) :: al(3), cosines(3)
    character(*),intent(INOUT) :: indx
    real(8),parameter :: eps=1.d-3
    if ( abs(al(1)-al(2)) < eps ) then
       if ( abs(al(2)-al(3)) < eps ) then
          if ( all( abs(cosines-0.5d0) < eps ) ) indx = "FCC"
       end if
    end if
  END SUBROUTINE check_fcc

END MODULE lattice_module
