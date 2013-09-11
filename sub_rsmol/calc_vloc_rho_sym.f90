!--------1---------2---------3---------4---------5---------6---------7--

      SUBROUTINE calc_vloc_rho_sym
      use global_variables
      implicit none
      integer :: a,i,i1,i2,i3,ierr,ik,j,m,m1,m2,n1,n2
      real(8) :: c1,c2,p1,p2,p3,p4,r,s(2),t(2),H
      real(8) :: ctime0,ctime1,etime0,etime1,mem,memax
      real(8),allocatable :: work(:)

      INTERFACE
         FUNCTION bberf(x)
         real(8) :: bberf
         real(8),intent(IN) :: x
         END FUNCTION bberf
      END INTERFACE

      call bwatch(ctime0,etime0)

      if( DISP_SWITCH )then
         write(*,'(a60," calc_vloc_rho_sym")') repeat("-",60)
      end if

      n1=idisp(myrank)+1
      n2=idisp(myrank)+ircnt(myrank)

      m1=idisp2(myrank)+1+ML_irreducible
      m2=m1+ircnt2(myrank)-1

      mem=0.d0
      memax=0.d0

      H=H1

!- allocate ---------------------
      call gv_alloc("vloc_b")
!--------------------------------

      Vion0(n1:n2)=Vion(n1:n2)

!
! --- Perform "Vh+Vion"-calculation or not ? ---
!

      if ( ifac_vloc==0 ) return

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

      allocate( work(0:ML_irreducible+MK_irreducible) ) ; work=0.d0
      mem=mem+bdreal*size(work) ; memax=max(mem,memax)

      call mpi_allgatherv(Vion(n1),ir_grid(myrank_g),mpi_real8,work(1),ircnt,idisp,mpi_real8,comm_grid,ierr)
      call mpi_allgatherv(vloc_b(m1),ir_grid2(myrank_g),mpi_real8 &
&                        ,work(ML_irreducible+1),ir_grid2,id_grid2,mpi_real8,comm_grid,ierr)

      rho_vloc(:)=0.d0

      c1=-1.d0/(4.d0*Pi*H*H)
      c2=3.d0*lap(0)*c1
      do i=n1,n2
         rho_vloc(i)=c2*Vion(i)
      end do
      do j=1,Md
         c2=c1*lap(j)
         do i=n1,n2
            i1=LL(1,i) ; i2=LL(2,i) ; i3=LL(3,i)
            rho_vloc(i)=rho_vloc(i)+c2*(  work( LLL(i1+j,i2,i3) ) + work( LLL(i1-j,i2,i3) ) &
&                                        +work( LLL(i1,i2+j,i3) ) + work( LLL(i1,i2-j,i3) ) &
&                                        +work( LLL(i1,i2,i3+j) ) + work( LLL(i1,i2,i3-j) ) )
         end do
      end do

      mem=mem-bdreal*size(work) ; deallocate( work )

!
! --- normalize ---
!
      s(1)=sum(rho_vloc(n1:n2)*weight(n1:n2))*dV
      call mpi_allreduce(s,t,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
      if (DISP_SWITCH) then
         write(*,*) "sum(rho_vloc)=",t(1)
      end if

      c1=sum(Ntot(1:nspin))/abs(t(1))
!      rho_vloc(:)=c1*rho_vloc(:)

      s(1)=sum(rho_vloc(n1:n2)*weight(n1:n2))*dV
      call mpi_allreduce(s,t,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
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
         write(*,*) "TIME(CALC_VLOC_RHO_SYM)=",ctime1-ctime0,etime1-etime0
      end if

      return

      END SUBROUTINE calc_vloc_rho_sym
