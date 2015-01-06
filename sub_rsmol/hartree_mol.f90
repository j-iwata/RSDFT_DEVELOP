!--------1---------2---------3---------4---------5---------6---------7--
!
! Hartree potential ( with real(8) density )
! ( The reference of this C.-G. Poisson solver is
!                  Phys. Rev. C17, 1682 (1978). )
! MEO=1,2 ; vloc+vh
!
      SUBROUTINE Hartree_mol(n1,n2,trho,tVh,err,maxiter,ifac)
      use global_variables
      implicit none
      integer,intent(IN)    :: n1,n2,maxiter,ifac
      real(8),intent(IN)    :: trho(n1:n2),err
      real(8),intent(INOUT) :: tVh(n1:n2)
      integer :: i,j,k,k1,k2,ix,iy,iz,m,n,lm,L,a,iter,m1,m2,ierr
      integer :: ML0,MK0,M_max
      real(8) :: sum0,sum1,sum2,ak,ck,b0,const,pi4,x,y,z,r,c
      real(8) :: s(4),t(4),E_0,plm(3),H
      real(8),parameter :: eps=1.d-10,eps1=1.d-16,eps2=1.d-32
      real(8),allocatable :: tk(:),zk(:),qk(:),sk(:)
      real(8),allocatable :: rholm(:,:),rholm_0(:,:)
      real(8),allocatable :: tVh0(:),zk0(:),tn(:)
      real(8),allocatable :: rholm0(:),rholm1(:)
      real(8),allocatable :: coef(:),clm(:,:)
      real(8) :: ctime0,ctime1,etime0,etime1,ct0,ct1,et0,et1
      real(8) :: mem,memax,ct(10),et(10),flops(10)
      logical :: flag_alloc(2)

      INTERFACE
         FUNCTION Ylm(x,y,z,l,m)
         real(8) :: Ylm
         real(8),intent(IN) :: x,y,z
         integer,intent(IN) :: l,m
         END FUNCTION Ylm
         FUNCTION plgndr(l,m,x)
         real(8) :: plgndr
         real(8),intent(INOUT) :: x
         integer,intent(IN) :: l,m
         END FUNCTION plgndr
      END INTERFACE

      if ( isymmetry==1 ) then
         call Hartree_mol_sym(n1,n2,trho,tVh,err,maxiter,ifac)
         return
      end if

      call watch(ctime0,etime0)

      ct(:)=0.d0 ; et(:)=0.d0 ; nop_ht(:)=0.d0 ; flops(:)=0.d0

      mem=0.d0
      memax=0.d0

      ML0 = n2-n1+1

      m1  = idisp2(myrank)+1
      m2  = idisp2(myrank)+ircnt2(myrank)
      MK0 = ircnt2(myrank)

      pi4=16.d0*atan(1.d0)

      E_0=0.d0

      H = H1

      www(:,:,:,:)=zero

!- allocate ------------------------------------------------------------
      allocate( tn(n1:n2) ) ; mem=mem+bdreal*ML0 ; memax=max(mem,memax)
!-----------------------------------------------------------------------

      if ( ifac>0 ) then
         tn(n1:n2) = trho(n1:n2) + ifac*rho_vloc(n1:n2)
         nop_ht(1)=nop_ht(1)+2.d0*ML0
      else
         tn(n1:n2) = trho(n1:n2)
      end if

      sum0=sum( tn(n1:n2)*tn(n1:n2) )*dV
      nop_ht(1)=nop_ht(1)+2.d0*ML0+1.d0

      call mpi_allreduce(sum0,b0,1,mpi_real8,mpi_sum,comm_grid,ierr)

      if( b0==zero )then
         tVh(n1:n2)=0.d0
         return
      end if

      tn(n1:n2)=trho(n1:n2)

!- allocate -------------------------------------------
      allocate( tVh0(n1:n2) ) ; tVh0(n1:n2)=tVh(n1:n2)
      allocate( tk(n1:n2),zk(n1:n2),qk(n1:n2) )
      mem=mem+bdreal*(ML0*4.d0) ; memax=max(mem,memax)
!------------------------------------------------------

      flag_alloc(:)=.false.
      if ( .not.allocated(LL) ) then
         flag_alloc(1)=.true.
         allocate( LL(3,n1:n2) ) ; LL=0
         mem=mem+bsintg*3*ML0 ; memax=max(mem,memax)
         call Make_GridMap_1(LL,n1,n2)
      end if
      if ( .not.allocated(KK) ) then
         flag_alloc(2)=.true.
         allocate( KK(3,m1:m2) ) ; KK=0
         mem=mem+bsintg*3*MK0 ; memax=max(mem,memax)
         call Make_GridMap_2(KK,m1,m2)
      end if

!
! --- Boundary condition (multipole expansion) ---
!
      call watch(ct0,et0)

      select case( MEO )
!
! Single center expansion
!
      case(1)

!- allocate ------------------------------------------------
         allocate( rholm0(lmmax_HT),rholm1(lmmax_HT) )
         mem=mem+bdreal*(lmmax_HT*2) ; memax=max(mem,memax)
!-----------------------------------------------------------

         do lm=1,lmmax_HT
            rholm0(lm)=sum( tn(n1:n2)*shf1(n1:n2,lm) )*dV
         end do
         nop_ht(1)=nop_ht(1)+2.d0*ML0+1.d0

         call mpi_allreduce(rholm0,rholm1,lmmax_HT,mpi_real8,mpi_sum,comm_grid,ierr)

         do i=m1,m2
            ix=KK(1,i)
            iy=KK(2,i)
            iz=KK(3,i)
            www(ix,iy,iz,1)=sum(shf2(1:lmmax_HT,i)*rholm1(1:lmmax_HT))
         end do
         nop_ht(1)=nop_ht(1)+2.d0*MK0

!- deallocate -----------------------------------------------------
         mem=mem-bdreal*(lmmax_HT*2) ; deallocate( rholm1,rholm0 )
!------------------------------------------------------------------

!
! Multi center expansion
!
      case(2)

         lmmax_HT=(Lmax_HT+1)**2

!- allocate ---------------------------------------------------
         allocate( rholm_0(lmmax_HT,MI),rholm(lmmax_HT,MI) )
         mem=mem+bdreal*(lmmax_HT*MI*2) ; memax=max(mem,memax)
!--------------------------------------------------------------

         rholm_0(:,:) = 0.d0
         rholm(:,:)   = 0.d0
         www(:,:,:,:)   = 0.d0
         M_max        = Lmax_HT

!- allocate -------------------------------------------------------------------------------
         allocate( sk(0:Lmax_HT) )
         allocate( coef(0:M_max),clm(0:Lmax_HT,0:M_max) )
         mem=mem+bdreal*size(sk)+bdreal*size(coef)+bdreal*size(clm) ; memax=max(mem,memax)
!------------------------------------------------------------------------------------------

! Pmm=(-1)^m*(2m-1)!!*(1-x^2)^(m/2)

         coef(:)=1.d0
         do m=1,M_max
            do i=1,2*m-1,2
               coef(m)=coef(m)*i
            end do
            coef(m)=coef(m)*(-1)**m
         end do

         Clm(:,:)=1.d0
         do L=0,Lmax_HT
         do m=0,L
            k1=l-m ; k2=l+m
            if( k1<k2 )then
               do k=k2,k1+1,-1
                  Clm(L,m)=Clm(L,m)*k
               end do
               Clm(L,m)=1.d0/Clm(L,m)
            else if( k1>k2 )then
               do k=k1,k2+1,-1
                  Clm(L,m)=Clm(L,m)*k
               end do
            end if
            Clm(L,m)=sqrt(Clm(L,m)*(2*l+1)/pi)*0.5d0
            if ( m/=0 ) clm(l,m)=clm(l,m)*sqrt(2.d0)
         end do
         end do

         do n=1,NMadv
            a=adv(n)
            do j=1,Mdv(n)
               i=Ixyz(j,n)
               x=LL(1,i)*H-asi(1,a)
               y=LL(2,i)*H-asi(2,a)
               z=LL(3,i)*H-asi(3,a)
               r=sqrt(x*x+y*y+z*z)
               if ( r==0.d0 ) then
                  rholm_0(1,a)=rholm_0(1,a)+tn(i)*Clm(0,0)
                  cycle
               end if
               sk(0)=tn(i)
               do L=1,Lmax_HT
                  sk(L)=tn(i)*r**L
               end do
               ck=z/r ; if ( abs(ck)>1.d0 ) ck=sign(1.d0,ck)
               if( abs(x)<1.d-10 )then
                  ak=0.5d0*pi
                  if ( y<0.d0 ) ak=ak+pi
               else
                  ak=atan(y/x)
                  if ( x<0.d0 ) ak=ak+pi
               end if
               lm=0
               do m=0,0
                  plm(1)=coef(m) !*(1.d0-ck*ck)**(0.5d0*m)
                  lm=lm+1
                  rholm_0(lm,a)=rholm_0(lm,a)+plm(1)*Clm(0,0)*sk(0)
                  plm(2)=ck*(2*m+1)*plm(1)
                  lm=lm+1
                  rholm_0(lm,a)=rholm_0(lm,a)+plm(2)*Clm(1,0)*sk(1)
                  do L=m+2,Lmax_HT
                     plm(3)=( ck*(2*L-1)*plm(2)-(L+m-1)*plm(1) )/(L-m)
                     lm=lm+1
                     rholm_0(lm,a)=rholm_0(lm,a)+plm(3)*Clm(L,0)*sk(L)
                     plm(1)=plm(2)
                     plm(2)=plm(3)
                  end do
               end do
               do m=1,M_max-1
                  plm(1)=coef(m)*(1.d0-ck*ck)**(0.5d0*m)
                  lm=lm+1
                  rholm_0(lm,a)=rholm_0(lm,a)+plm(1)*Clm(m,m)*sk(m)*cos(ak*m)
                  lm=lm+1
                  rholm_0(lm,a)=rholm_0(lm,a)-plm(1)*Clm(m,m)*sk(m)*sin(ak*m)*(-1)**m
                  plm(2)=ck*(2*m+1)*plm(1)
                  lm=lm+1
                  rholm_0(lm,a)=rholm_0(lm,a)+plm(2)*Clm(m+1,m)*sk(m+1)*cos(ak*m)
                  lm=lm+1
                  rholm_0(lm,a)=rholm_0(lm,a)-plm(2)*Clm(m+1,m)*sk(m+1)*sin(ak*m)*(-1)**m
                  do L=m+2,Lmax_HT
                     plm(3)=( ck*(2*L-1)*plm(2)-(L+m-1)*plm(1) )/(L-m)
                     lm=lm+1
                     rholm_0(lm,a)=rholm_0(lm,a)+plm(3)*Clm(L,m)*sk(L)*cos(ak*m)
                     lm=lm+1
                     rholm_0(lm,a)=rholm_0(lm,a)-plm(3)*Clm(L,m)*sk(L)*sin(ak*m)*(-1)**m
                     plm(1)=plm(2)
                     plm(2)=plm(3)
                  end do
               end do
               do m=M_max,M_max
                  plm(1)=coef(m)*(1.d0-ck*ck)**(0.5d0*m)
                  lm=lm+1
                  rholm_0(lm,a)=rholm_0(lm,a)+plm(1)*Clm(m,m)*sk(m)*cos(ak*m)
                  lm=lm+1
                  rholm_0(lm,a)=rholm_0(lm,a)-plm(1)*Clm(m,m)*sk(m)*sin(ak*m)*(-1)**m
               end do
            end do ! j
         end do ! n
         rholm_0(:,:)=rholm_0(:,:)*dV
         call mpi_allreduce(rholm_0,rholm,lmmax_HT*MI,mpi_real8,mpi_sum,comm_grid,ierr)

         do a=1,MI
            do i=m1,m2
               ix=KK(1,i) ; iy=KK(2,i) ; iz=KK(3,i)
               x=ix*H-asi(1,a)
               y=iy*H-asi(2,a)
               z=iz*H-asi(3,a)
               r=sqrt(x*x+y*y+z*z)
               do L=0,Lmax_HT
                  sk(L)=pi4/((2.d0*L+1.d0)*r**(L+1))
               end do
               ck=z/r
               if( abs(x)<1.d-10 )then
                  ak=0.5d0*pi
                  if ( y<0.d0 ) ak=ak+pi
               else
                  ak=atan(y/x)
                  if ( x<0.d0 ) ak=ak+pi
               end if
               lm=0
               do m=0,0
                  plm(1)=coef(m) !*(1.d0-ck*ck)**(0.5d0*m)
                  lm=lm+1
                  www(ix,iy,iz,1)=www(ix,iy,iz,1)+plm(1)*Clm(0,0)*sk(0)*rholm(lm,a)
                  plm(2)=ck*(2*m+1)*plm(1)
                  lm=lm+1
                  www(ix,iy,iz,1)=www(ix,iy,iz,1)+plm(2)*Clm(1,0)*sk(1)*rholm(lm,a)
                  do L=m+2,Lmax_HT
                     plm(3)=( ck*(2*L-1)*plm(2)-(L+m-1)*plm(1) )/(L-m)
                     lm=lm+1
                     www(ix,iy,iz,1)=www(ix,iy,iz,1)+plm(3)*Clm(L,0)*sk(L)*rholm(lm,a)
                     plm(1)=plm(2)
                     plm(2)=plm(3)
                  end do
               end do
               do m=1,M_max-1
                  plm(1)=coef(m)*(1.d0-ck*ck)**(0.5d0*m)
                  lm=lm+1
                  www(ix,iy,iz,1)=www(ix,iy,iz,1)+plm(1)*Clm(m,m)*sk(m)*cos(ak*m)*rholm(lm,a)
                  lm=lm+1
                  www(ix,iy,iz,1)=www(ix,iy,iz,1)-plm(1)*Clm(m,m)*sk(m)*sin(ak*m)*(-1)**m*rholm(lm,a)

                  plm(2)=ck*(2*m+1)*plm(1)
                  lm=lm+1
                  www(ix,iy,iz,1)=www(ix,iy,iz,1)+plm(2)*Clm(m+1,m)*sk(m+1)*cos(ak*m)*rholm(lm,a)
                  lm=lm+1
                  www(ix,iy,iz,1)=www(ix,iy,iz,1)-plm(2)*Clm(m+1,m)*sk(m+1)*sin(ak*m)*(-1)**m*rholm(lm,a)

                  do L=m+2,Lmax_HT
                     plm(3)=( ck*(2*L-1)*plm(2)-(L+m-1)*plm(1) )/(L-m)
                     lm=lm+1
                     www(ix,iy,iz,1)=www(ix,iy,iz,1)+plm(3)*Clm(L,m)*sk(L)*cos(ak*m)*rholm(lm,a)
                     lm=lm+1
                     www(ix,iy,iz,1)=www(ix,iy,iz,1)-plm(3)*Clm(L,m)*sk(L)*sin(ak*m)*(-1)**m*rholm(lm,a)
                     plm(1)=plm(2)
                     plm(2)=plm(3)
                  end do
               end do
               do m=M_max,M_max
                  plm(1)=coef(m)*(1.d0-ck*ck)**(0.5d0*m)
                  lm=lm+1
                  www(ix,iy,iz,1)=www(ix,iy,iz,1)+plm(1)*Clm(m,m)*sk(m)*cos(ak*m)*rholm(lm,a)
                  lm=lm+1
                  www(ix,iy,iz,1)=www(ix,iy,iz,1)-plm(1)*Clm(m,m)*sk(m)*sin(ak*m)*(-1)**m*rholm(lm,a)
               end do
 10            continue
            end do ! j
         end do ! a
 4       continue

!- deallocate ----------------------------------------------------------------------
         mem=mem-bdreal*(size(clm)+size(coef)+size(sk)) ; deallocate( clm,coef,sk )
         mem=mem-bdreal*(size(rholm)+size(rholm_0)) ; deallocate( rholm,rholm_0 )
!-----------------------------------------------------------------------------------

      end select

      if ( ifac>0 ) then
         do i=m1,m2
            ix=KK(1,i)
            iy=KK(2,i)
            iz=KK(3,i)
            www(ix,iy,iz,1)=www(ix,iy,iz,1)+ifac*vloc_b(i)
         end do
         tn(n1:n2) = trho(n1:n2) + ifac*rho_vloc(n1:n2)
      end if

      call watch(ct1,et1) ; ct(1)=ct1-ct0; et(1)=et1-et0

!
! --- C.-G. minimization ---
!

      do i=n1,n2
         ix=LL(1,i) ; iy=LL(2,i) ; iz=LL(3,i)
         www(ix,iy,iz,1)=tVh(i)
      end do

      call bcset_2d(1,1,Md,0)

      const=3.d0*lap(0)/H**2
      zk(n1:n2)=-const*tVh(n1:n2)-pi4*tn(n1:n2)
      nop_ht(1)=nop_ht(1)+ML0*2.d0+ML0*2.d0

      do j=1,Md
         c=lap(j)/H**2
      do i=n1,n2
         ix=LL(1,i) ; iy=LL(2,i) ; iz=LL(3,i)
         zk(i)=zk(i)-c*( www(ix-j,iy,iz,1)+www(ix+j,iy,iz,1) &
     &                  +www(ix,iy-j,iz,1)+www(ix,iy+j,iz,1) &
     &                  +www(ix,iy,iz-j,1)+www(ix,iy,iz+j,1) )
      end do
      end do

      nop_ht(1)=nop_ht(1)+7.d0*Md*ML0

      www(:,:,:,1)=0.d0
      qk(:)=zk(:)
      sum0=sum(zk(n1:n2)*zk(n1:n2))*dV
      nop_ht(1)=nop_ht(1)+2.d0*ML0+1.d0

      call mpi_allreduce(sum0,sum1,1,mpi_real8,mpi_sum,comm_grid,ierr)

!
! --- Iteration ---
!
      Iteration : do iter=1,maxiter

         do i=n1,n2
            ix=LL(1,i) ; iy=LL(2,i) ; iz=LL(3,i)
            www(ix,iy,iz,1)=qk(i)
         end do

         call bcset_2d(1,1,Md,0)

         tk(n1:n2)=const*qk(n1:n2)
         nop_ht(1)=nop_ht(1)+ML0

         do j=1,Md
            c=lap(j)/H**2
         do i=n1,n2
            ix=LL(1,i) ; iy=LL(2,i) ; iz=LL(3,i)
            tk(i)=tk(i)+c*( www(ix-j,iy,iz,1)+www(ix+j,iy,iz,1) &
     &                     +www(ix,iy-j,iz,1)+www(ix,iy+j,iz,1) &
     &                     +www(ix,iy,iz-j,1)+www(ix,iy,iz+j,1) )
         end do
         end do
         nop_ht(1)=nop_ht(1)+7.d0*ML0*Md

         sum0=sum(zk(n1:n2)*tk(n1:n2))*dV
         nop_ht(1)=nop_ht(1)+2.d0*ML0+1.d0

         call mpi_allreduce(sum0,sum2,1,mpi_real8,mpi_sum,comm_grid,ierr)

         tVh0(n1:n2)=tVh(n1:n2)

         ak=sum1/sum2
         tVh(n1:n2) = tVh(n1:n2) + ak*qk(n1:n2)
         zk(:)=zk(:)-ak*tk(:)

         nop_ht(1)=nop_ht(1)+(ML0+ML0)*2.d0

!
! Conv. Check
!
         E_0=E_Hartree
         sum0=sum(zk(n1:n2)*zk(n1:n2))*dV
         s(1)=sum0
!         s(2)=0.5d0*sum( tVh(n1:n2)*tn(n1:n2) )*dV
         s(2)=0.5d0*sum( tVh(n1:n2)*trho(n1:n2) )*dV
         s(3)=sum((tVh(n1:n2)-tVh0(n1:n2))**2)/ML
         s(4)=maxval(abs(tVh(n1:n2)-tVh0(n1:n2)))

         nop_ht(1)=nop_ht(1)+ML0+1.d0+ML0

         call mpi_allreduce(s,t,4,mpi_real8,mpi_sum,comm_grid,ierr)

         sum2=t(1)
         E_Hartree=t(2)
         if ( t(4)<=eps1 ) exit

         ck=sum2/sum1 ; sum1=sum2

         qk(n1:n2)=ck*qk(n1:n2)+zk(n1:n2)

         nop_ht(1)=nop_ht(1)+ML0+ML0

      end do Iteration

!- deallocate ----------------------------------------------------------------
      if ( flag_alloc(1) ) then
         mem=mem-bsintg*size(LL) ; deallocate(LL)
      end if
      if ( flag_alloc(2) ) then
         mem=mem-bsintg*size(KK) ; deallocate(KK)
      end if
      mem=mem-bdreal*size(tn) ; deallocate( tn )
      mem=mem-bdreal*( size(tk)+size(zk)+size(qk) ) ; deallocate( tk,zk,qk )
      mem=mem-bdreal*size(tVh0) ; deallocate( tVh0 )
!----------------------------------------------------------------------------

      call watch(ctime1,etime1)
      nop_ht(10)=sum( nop_ht(1:9) )
      ct(10)=ctime1-ctime0
      et(10)=etime1-etime0

      where(ct>0.d0)
         flops=nop_ht/ct*1.d-6
      end where

      if (DISP_SWITCH) then
         write(*,*) "TIME(HARTREE_MOL)=",ctime1-ctime0,etime1-etime0
         write(*,*) " (BC) :",ct(1),et(1),ML
         write(*,'(1x,"iter/maxiter, conv =",i4"/",i4,3g16.7)') iter,maxiter,sum2/b0,t(3:4)
         write(*,*) "E_Hartree,dE=",E_Hartree,E_Hartree-E_0
         write(*,'(1x,"(TOT)",3f20.10)') ct(10),et(10),flops(10)
         write(*,*) "MEM(MB)",memax*B2MB,mem
      end if

      return

      END SUBROUTINE Hartree_mol
