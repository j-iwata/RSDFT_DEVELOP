!
! Check program for the output files of ldos2gp ( fort.11 fort.12 fort.13 )
!
! Usage: ./a.out < fort.11 ( or fort.12 or fort.13 )
!
! Output: fort.20 & fort.21
!
PROGRAM a

  implicit none

  integer,parameter :: max_loop=10000000
  integer,parameter :: u0=5, u1=20
  integer :: i,ne,ng,j
  real(8) :: e0, e, dx, x0, x1, d, tmp

! ---

  read(u0,*) e, x0
  read(u0,*) e, x1
  dx=x1-x0

! ---

  rewind u0

  e0=1.0d100
  ne=0
  ng=0

  do i=1,max_loop

     read(u0,*,end=9) e
     ng=ng+1

     if ( e /= e0 ) then
        e0=e
        ne=ne+1
        ng=1
     end if

  end do

9 write(*,*) "ne,ng,ndata",ne,ng,i-1
  write(*,*) "dx=",dx

! ---

  rewind u0
  rewind u1
  do i=1,ne
     tmp=0.0d0
     do j=1,ng/2
        read(u0,*) e, x0, d
        tmp=tmp+d
     end do
     do j=ng/2+1,ng
        read(u0,*) e, x0, d
     end do
     write(u1,*) e, tmp*dx
  end do

  rewind u0
  rewind u1+1
  do i=1,ne
     tmp=0.0d0
     do j=1,ng/2
        read(u0,*) e, x0, d
     end do
     do j=ng/2+1,ng
        read(u0,*) e, x0, d
        tmp=tmp+d
     end do
     write(u1+1,*) e, tmp*dx
  end do

END PROGRAM a
