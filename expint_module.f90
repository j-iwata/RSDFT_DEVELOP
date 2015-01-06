MODULE expint_module

  implicit none

  PRIVATE
  PUBLIC :: expint

CONTAINS

  function expint(n,x)
    implicit none
    integer,intent(in) :: n
    real(8),intent(in) :: x
    real(8) :: expint
    integer,parameter :: maxit=200
    real(8),parameter :: esp=1.d-12, big=huge(x)*esp
    real(8),parameter :: euler=0.577215664901532860606512d0
    integer :: i,nm1,j
    real(8) :: a,b,c,d,del,fact,h,arsum
 
    if ( .not.(n>=0.and.x>=0.d0.and.(x>0.d0.or.n>1)) ) then
       write(*,*) 'Bad arguments in expint.f'
       stop
    end if
 
    if ( n==0 ) then
       expint=exp(-x)/x
       return
    end if
    nm1=n-1
    if ( x==0.d0 ) then
       expint=1.d0/nm1
    else if ( x>1.d0 ) then
       b=x+n
       c=big
       d=1.d0/b
       h=d
       do i=1,maxit
          a=-i*(nm1+i)
          b=b+2.d0
          d=1.d0/(a*d+b)
          c=b+a/c
          del=c*d
          h=h*del
          if ( abs(del-1.d0)<=esp ) exit
       end do
       if ( i>maxit ) then
          write(*,*) 'Continued fraction failed in expint.f'
          stop
       end if
       expint=h*exp(-x)
    else
       if ( nm1/=0 ) then
          expint=1.d0/nm1
       else
          expint=-log(x)-euler
       end if
       fact=1.d0
       do i=1,maxit
          fact=-fact*x/i
          if ( i/=nm1 ) then
             del=-fact/(i-nm1)
          else
             arsum = 0.d0
             do j=1,nm1
                arsum = arsum + 1.d0/j
             end do
             del = fact*(-LOG(x)-euler+arsum)
          end if
          expint=expint+del
          if ( abs(del)<abs(expint)*esp ) exit
       end do
       if ( i>maxit ) then
          write(*,*) 'series failed in expint.f'
          stop
       end if
    end if

  end function expint

END MODULE expint_module
