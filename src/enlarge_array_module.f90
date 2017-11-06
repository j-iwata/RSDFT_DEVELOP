MODULE enlarge_array_module

  implicit none

  PRIVATE
  PUBLIC :: enlarge_array

  INTERFACE enlarge_array
     MODULE PROCEDURE enlarge_int_array3, enlarge_real8_array2 &
          ,enlarge_complex16_array2
  END INTERFACE

CONTAINS


  SUBROUTINE enlarge_int_array3( a, n1,n2,n3 )
    implicit none
    integer,allocatable,intent(INOUT) :: a(:,:,:)
    integer,intent(IN) :: n1,n2,n3
    integer :: m1,m2,m3
    integer,allocatable :: b(:,:,:)
    if ( .not.allocated(a) ) then
       allocate( a(n1,n2,n3) ) ; a=0
    else
       m1=size(a,1)
       m2=size(a,2)
       m3=size(a,3)
       allocate( b(m1,m2,m3) ) ; b=0
       b=a
       deallocate( a )
       allocate( a(n1,n2,n3) ) ; a=0
       a(1:m1,1:m2,1:m3)=b
       deallocate( b )
    end if
  END SUBROUTINE enlarge_int_array3


  SUBROUTINE enlarge_real8_array2( a, n1,n2 )
    implicit none
    real(8),allocatable,intent(INOUT) :: a(:,:)
    integer,intent(IN) :: n1,n2
    integer :: m1,m2
    real(8),allocatable :: b(:,:)
    if ( .not.allocated(a) ) then
       allocate( a(n1,n2) ) ; a=0.0d0
    else
       m1=size(a,1)
       m2=size(a,2)
       allocate( b(m1,m2) ) ; b=0.0d0
       b=a
       deallocate( a )
       allocate( a(n1,n2) ) ; a=0.0d0
       a(1:m1,1:m2)=b
       deallocate( b )
    end if
  END SUBROUTINE enlarge_real8_array2


  SUBROUTINE enlarge_complex16_array2( a, n1,n2 )
    implicit none
    complex(8),allocatable,intent(INOUT) :: a(:,:)
    integer,intent(IN) :: n1,n2
    integer :: m1,m2
    complex(8),allocatable :: b(:,:)
    if ( .not.allocated(a) ) then
       allocate( a(n1,n2) ) ; a=(0.0d0,0.0d0)
    else
       m1=size(a,1)
       m2=size(a,2)
       allocate( b(m1,m2) ) ; b=(0.0d0,0.0d0)
       b=a
       deallocate( a )
       allocate( a(n1,n2) ) ; a=(0.0d0,0.0d0)
       a(1:m1,1:m2)=b
       deallocate( b )
    end if
  END SUBROUTINE enlarge_complex16_array2


END MODULE enlarge_array_module
