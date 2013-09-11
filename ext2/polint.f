
      SUBROUTINE polint(xa,ya,n,x,y,dy)
      implicit none
      integer,intent(IN)  :: n
      real(8),intent(IN)  :: xa(n),ya(n),x
      real(8),intent(OUT) :: y,dy
      integer :: i,m
      real(8) :: ho,hp,f0av0,f1av0,f0av,f1av
      real(8),allocatable :: f0(:),f1(:)
      f0av=ya(1)
c      f1av=( ya(2)-ya(1) )/( xa(2)-xa(1) )
      allocate( f0(n) ) ; f0=ya
c      allocate( f1(n) ) ; f1=0.d0
      do m=1,n-1
         f0av0=f0av
c         f1av0=f1av
         do i=1,n-m
            ho=x-xa(i)
            hp=x-xa(i+m)
c            f1(i)=( f0(i)-f0(i+1)+hp*f1(i)-ho*f1(i+1) )/(hp-ho)
            f0(i)=( hp*f0(i)-ho*f0(i+1) )/(hp-ho)
         end do
         f0av=sum(f0(1:n-m))/(n-m)
c         f1av=sum(f1(1:n-m))/(n-m)
      end do
      y=f0(1) ; dy=f0av-f0av0
c      y=f1(1) ; dy=f1av-f1av0
      deallocate( f0 )
c      deallocate( f1 )
      return
      END SUBROUTINE polint

      SUBROUTINE dpolint(xa,ya,n,x,y,dy)
      implicit none
      integer,intent(IN)  :: n
      real(8),intent(IN)  :: xa(n),ya(n),x
      real(8),intent(OUT) :: y,dy
      integer :: i,m
      real(8) :: ho,hp,f0av0,f1av0,f0av,f1av
      real(8),allocatable :: f0(:),f1(:)
      f0av=ya(1)
      f1av=( ya(2)-ya(1) )/( xa(2)-xa(1) )
      allocate( f0(n),f1(n) )
      f0=ya
      f1=0.d0
      do m=1,n-1
         f0av0=f0av
         f1av0=f1av
         do i=1,n-m
            ho=x-xa(i)
            hp=x-xa(i+m)
            f1(i)=( f0(i)-f0(i+1)+hp*f1(i)-ho*f1(i+1) )/(hp-ho)
            f0(i)=( hp*f0(i)-ho*f0(i+1) )/(hp-ho)
         end do
         f0av=sum(f0(1:n-m))/(n-m)
         f1av=sum(f1(1:n-m))/(n-m)
c         write(*,*) f0av,f1av,f1av-f1av0
      end do
c      y=f0(1) ; dy=f0av-f0av0
      y=f1(1) ; dy=f1av-f1av0
c      write(*,*) y,dy
      deallocate( f0,f1 )
      return
      END SUBROUTINE dpolint
