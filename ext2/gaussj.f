!--------1---------2---------3---------4---------5---------6---------7--
! Gauss-Jordan
! (from Numerical Recipes)

      SUBROUTINE gaussj(a,n,np,b,m,mp)
      implicit none
      integer,intent(IN) :: n,np,m,mp
      real(8),intent(INOUT) :: a(np,np),b(np,mp)
      integer,parameter :: NMAX=5000
      integer :: indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      integer :: i,j,k,l,irow,icol,lll
      real(8) :: dum,pivinv,big
      do j=1,n
         ipiv(j)=0
      enddo
      do i=1,n
         big=0.d0
         do j=1,n
            if(ipiv(j).ne.1) then
               do k=1,n
                  if(ipiv(k).eq.0) then
                     if(abs(a(j,k)).ge.big) then
                        big=abs(a(j,k))
                        irow=j
                        icol=k
                     endif
                  else if(ipiv(k).gt.1) then
                     write(*,*) 'singular matrix in gaussj'
                     stop
                  endif
               enddo
            endif
         enddo
         ipiv(icol)=ipiv(icol)+1
         if(irow.ne.icol) then
            do l=1,n
               dum=a(irow,l)
               a(irow,l)=a(icol,l)
               a(icol,l)=dum
            enddo
            do l=1,m
               dum=b(irow,l)
               b(irow,l)=b(icol,l)
               b(icol,l)=dum
            enddo
         endif
         indxr(i)=irow
         indxc(i)=icol
         if(abs(a(icol,icol)).eq.0.d0) then
            write(*,*) 'singular matrix in gaussj'
            stop
         endif
         pivinv=1.d0/a(icol,icol)
         a(icol,icol)=1.d0
         do l=1,n
            a(icol,l)=a(icol,l)*pivinv
         enddo
         do l=1,m
            b(icol,l)=b(icol,l)*pivinv
         enddo
         do lll=1,n
            if(lll.ne.icol) then
               dum=a(lll,icol)
               a(lll,icol)=0.d0
               do l=1,n
                  a(lll,l)=a(lll,l)-a(icol,l)*dum
               enddo
               do l=1,m
                  b(lll,l)=b(lll,l)-b(icol,l)*dum
               enddo
            endif
         enddo
      enddo
      do l=n,1,-1
         if(indxr(l).ne.indxc(l)) then
            do k=1,n
               dum=a(k,indxr(l))
               a(k,indxr(l))=a(k,indxc(l))
               a(k,indxc(l))=dum
            enddo
         endif
      enddo
      return
      END SUBROUTINE gaussj
