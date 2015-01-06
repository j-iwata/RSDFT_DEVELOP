!--------1---------2---------3---------4---------5---------6---------7--
!
! Division of Area for Multipole Expansion (MEO=2)
!
      SUBROUTINE prep_hartree2_sym
      use global_variables
      implicit none
      real(8) :: ctime0,ctime1,etime0,etime1,x,y,z,r2,r2min,H
      integer :: i,j,a,n,n1,n2,ierr

      call watch(ctime0,etime0)

      if (DISP_SWITCH) write(*,'(a60," prep_hartree2_sym")') repeat("-",60)

      n1=idisp(myrank)+1
      n2=idisp(myrank)+ircnt(myrank)

      H=H1

      n=maxval( Nstar(n1:n2) )
      allocate( Ixyz(n,n1:n2) ) ; Ixyz=0

      do i=n1,n2
         do j=1,Nstar(i)
            r2min=10.d10
            do a=1,MI
               x=LL_star(1,j,i)*H-asi(1,a)
               y=LL_star(2,j,i)*H-asi(2,a)
               z=LL_star(3,j,i)*H-asi(3,a)
               r2=x*x+y*y+z*z
               if ( r2<r2min ) then
                  r2min=r2
                  Ixyz(j,i)=a
               end if
            end do
         end do
      end do

      if ( DISP_SWITCH ) then
!         do a=1,MI
!            write(*,*) a,count(Ixyz==a)
!         end do
         write(*,*) "count(Ixyz/=0)",count(Ixyz/=0),n2-n1+1
      end if

      call watch(ctime1,etime1)
      if (DISP_SWITCH) then
         write(*,*) "TIME(PREP_HARTREE2_SYM)=",ctime1-ctime0,etime1-etime0
         write(*,*) "MEO=",MEO
      end if

      return

 900  call stop_program

      END SUBROUTINE prep_hartree2_sym
