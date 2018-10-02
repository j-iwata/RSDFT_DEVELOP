PROGRAM disturb_atomic_coordinates
!MODULE disturb_module
!
  implicit none
!
!  PRIVATE
!  PUBLIC :: disturb

  real(8) :: ratio
  real(8) :: ax,aa(3,3)
  real(8),allocatable :: aa_atom(:,:)
  integer,allocatable :: ki_atom(:),md_atom(:)
  integer,parameter :: u=970, w=971
  integer :: nelem,natom,zatom(9),i
  character(3) :: cbuf,cformat

  write(*,*) "input the ratio of disturbance [0.0~1.0]" 
  read(*,*) ratio
  write(*,*) "atomic coordinates are read from fort.970"
  read(u,*) cbuf,ax
  read(u,*) cbuf,aa(1:3,1)
  read(u,*) cbuf,aa(1:3,2)
  read(u,*) cbuf,aa(1:3,3)
  read(u,*) cformat
  read(u,*) nelem, natom, zatom(1:nelem)
  allocate( aa_atom(3,natom) ) ; aa_atom=0.0d0
  allocate( ki_atom(natom)   ) ; ki_atom=0
  allocate( md_atom(natom)   ) ; md_atom=0
  do i=1,natom
     read(u,*) ki_atom(i),aa_atom(1:3,i),md_atom(i)
  end do

  call disturb( aa_atom, ratio )

  write(*,*) "the disturbed atomic coordinates are written in fort.971"
  
  write(w,'("AX", f20.16)') ax
  write(w,'("A1",3f20.16)') aa(1:3,1)
  write(w,'("A2",3f20.16)') aa(1:3,2)
  write(w,'("A3",3f20.16)') aa(1:3,3)
  write(w,'(a3)') cformat
  write(w,'(1x,i4,10i8)') nelem,natom,zatom(1:nelem)
  do i=1,natom
     write(w,'(1x,i5,3f20.12,i4)') ki_atom(i),aa_atom(1:3,i),md_atom(i)
  end do

CONTAINS

  SUBROUTINE disturb( a, ratio )

    implicit none
    real(8),intent(INOUT) :: a(:,:)
    real(8),intent(IN) :: ratio
    real(8),allocatable :: d(:)
    real(8) :: r
    integer :: i

    if ( ratio <= 0.0d0 ) return
    r = min( ratio, 1.0d0 )

    allocate( d(size(a,1)) ) ; d=0.0d0

    do i=1,size(a,2)
       call random_number( d )
       d(:) = ( 2.0d0*d(:) - 1.0d0 )*r
       a(:,i) = a(:,i) + d(:)
!       write(*,'(1x,i4,2(5x,3f15.10))') i,a(:,i)-d(:),a(:,i)
    end do

    deallocate( d )

  END SUBROUTINE disturb

END PROGRAM disturb_atomic_coordinates
!END MODULE disturb_module
