!--------1---------2---------3---------4---------5---------6---------7--

      SUBROUTINE calc_vloc_rho
      use global_variables
      implicit none
      integer :: a,i,i1,i2,i3,ierr,ik,j,m,m1,m2,n1,n2,Mx,My,Mz
!      integer,allocatable :: KKK1(:,:,:),KKK0(:,:,:)
      real(8) :: c1,c2,p1,p2,p3,p4,r,s(2),t(2),H
      real(8) :: ctime0,ctime1,etime0,etime1,mem,memax
      real(8),allocatable :: work(:)

      INTERFACE
         FUNCTION bberf(x)
         real(8) :: bberf
         real(8),intent(IN) :: x
         END FUNCTION bberf
      END INTERFACE

      if ( isymmetry==1 ) then
         call calc_vloc_rho_sym
         return
      end if

      call watch(ctime0,etime0)

      if( DISP_SWITCH )then
         write(*,'(a60," calc_vloc_rho")') repeat("-",60)
      end if

      n1=idisp(myrank)+1
      n2=idisp(myrank)+ircnt(myrank)

      m1=idisp2(myrank)+1
      m2=idisp2(myrank)+ircnt2(myrank)

      mem=0.d0
      memax=0.d0

      H=H1
      Mx=ML1+Md
      My=ML2+Md
      Mz=ML3+Md

!- allocate ---------------------
      call gv_alloc("vloc_b")
!--------------------------------

      Vion0(n1:n2)=Vion(n1:n2)

!
! --- Perform "Vh+Vion"-calculation or not ? ---
!

      if ( ifac_vloc==0 ) return

!- allocate ----------------------------------------------------
      allocate( LL(3,n1:n2) ) ; LL=0
      allocate( KK(3,m1:m2) ) ; KK=0
      mem=mem+bsintg*(size(LL)+size(KK)) ; memax=max(mem,memax)
!---------------------------------------------------------------

      call Make_GridMap_1(LL,n1,n2)
      call Make_GridMap_2(KK,m1,m2)

      vloc_b(:)=0.d0
      do a=1,MI
         ik=Kion(a)
         p1=-Zps(ik)*parloc(1,ik) ; p2=sqrt(parloc(2,ik))
         p3=-Zps(ik)*parloc(3,ik) ; p4=sqrt(parloc(4,ik))
         do i=m1,m2
            i1=KK(1,i)
            i2=KK(2,i)
            i3=KK(3,i)
            r=sqrt(dble(i1*i1+i2*i2+i3*i3))*H
            vloc_b(i)=vloc_b(i)+( p1*bberf(p2*r)+p3*bberf(p4*r) )/r
         end do
      end do

      rho_vloc(:)=0.d0
      www(:,:,:,:)=zero
      do i=n1,n2
         www(LL(1,i),LL(2,i),LL(3,i),1)=Vion(i)
      end do
      do i=m1,m2
         www(KK(1,i),KK(2,i),KK(3,i),1)=vloc_b(i)
      end do
      call bcset_2d(1,1,Md,0)
      c1=-1.d0/(4.d0*Pi)
      c2=-2.d0*coef_lap0*c1
      do i=n1,n2
         rho_vloc(i)=c2*Vion(i)
      end do
      do j=1,Md
         c2=c1*(-2.d0*coef_lap(1,j))
         do i=n1,n2
            i1=LL(1,i) ; i2=LL(2,i) ; i3=LL(3,i)
            rho_vloc(i)=rho_vloc(i)+c2*( www(i1-j,i2,i3,1)+www(i1+j,i2,i3,1) &
&                                       +www(i1,i2-j,i3,1)+www(i1,i2+j,i3,1) &
&                                       +www(i1,i2,i3-j,1)+www(i1,i2,i3+j,1) )
         end do
      end do

      www(:,:,:,:)=zero

!- deallocate -----------------------------------
      mem=mem-bsintg*size(KK) ; deallocate( KK )
      mem=mem-bsintg*size(LL) ; deallocate( LL )
!------------------------------------------------

!
! --- normalize ---
!
      s(1)=sum(rho_vloc(n1:n2))*dV
      call mpi_allreduce(s,t,1,mpi_real8,mpi_sum,comm_grid,ierr)
      if (DISP_SWITCH) then
         write(*,*) "sum(rho_vloc)=",t(1)
      end if

      c1=sum(Ntot(1:nspin))/abs(t(1))
!      rho_vloc(:)=c1*rho_vloc(:)

      s(1)=sum(rho_vloc(n1:n2))*dV
      call mpi_allreduce(s,t,1,mpi_real8,mpi_sum,comm_grid,ierr)
      if (DISP_SWITCH) then
         write(*,*) "sum(rho_vloc)=",t(1)
      end if

!- allocate -------------------------------------------
      allocate( work(n1:n2) ) ; work=0.d0
      mem=mem+bdreal*size(work) ; memax=max(mem,memax)
!------------------------------------------------------

      Vh(n1:n2)=Vion(n1:n2)
      call Hartree_mol(n1,n2,work,Vh,0.d0,2000,ifac_vloc)
      Vion(n1:n2)=Vion0(n1:n2)-Vh(n1:n2)
      Vh(n1:n2)=0.d0

!- deallocate ---------------------------------------
      mem=mem-bdreal*size(work) ; deallocate( work )
!----------------------------------------------------

      call watch(ctime1,etime1)
      if( DISP_SWITCH )then
         write(*,*) "TIME(CALC_VLOC_RHO)=",ctime1-ctime0,etime1-etime0
      end if

      return

      END SUBROUTINE calc_vloc_rho
