!--------1---------2---------3---------4---------5---------6---------7--

      SUBROUTINE force_IonIon(force3)
      use global_variables
      implicit none
      real(8),intent(OUT) :: force3(3,MI)
      integer :: a,b
      real(8) :: r3,c1,ctime0,ctime1,etime0,etime1

      call watch(ctime0,etime0)

      force3(:,:)=0.d0

      do a=1,MI

         do b=1,MI

            if ( a==b ) cycle

            r3=sqrt( (asi(1,a)-asi(1,b))**2+(asi(2,a)-asi(2,b))**2+(asi(3,a)-asi(3,b))**2 )**3

            c1=Zps(Kion(a))*Zps(Kion(b))/r3

            force3(1,a)=Force3(1,a)+c1*(asi(1,a)-asi(1,b))
            force3(2,a)=Force3(2,a)+c1*(asi(2,a)-asi(2,b))
            force3(3,a)=Force3(3,a)+c1*(asi(3,a)-asi(3,b))

         end do

      end do

      call watch(ctime1,etime1)

      if (DISP_SWITCH) then
         write(*,*) "Ion-Ion",sum(force3**2)
         write(*,*) ctime1-ctime0,etime1-etime0
      end if

      return
      END SUBROUTINE force_IonIon
