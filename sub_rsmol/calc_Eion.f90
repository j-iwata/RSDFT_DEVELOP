!--------1---------2---------3---------4---------5---------6---------7--

      SUBROUTINE calc_Eion
      use global_variables
      integer :: p,ia,ib
      real(8) :: r

      Ewld=0.d0

      do ia=2,MI
         do ib=1,ia-1
            r=sqrt( sum((asi(:,ia)-asi(:,ib))**2) )
            Ewld=Ewld+Zps(Kion(ia))*Zps(Kion(ib))/r
         end do
      end do

      return
      END SUBROUTINE calc_Eion
