!--------1---------2---------3---------4---------5---------6---------7--

      SUBROUTINE force_local_allel
      use global_variables
      implicit none
      integer :: a,i,ik,ierr,n1,n2,ic1,ic2,ic3,id1,id2,id3,m
      integer :: i1,i2,i3,ii1,ii2,ii3,iii1,iii2,iii3,ll_0,j1,j2,j3
      real(8) :: cnsta,cnstb,cnst1,cnst0,tmp1,tmp2,tmp3
      real(8) :: x,y,z,r,r2,r3,br,fx,fy,fz,a2,b2,rps2
      real(8) :: Rx,Ry,Rz
      real(8) :: mem,memax
      real(8),allocatable :: ftmp(:,:),rwork2(:,:)
      real(8),allocatable :: rho3(:,:,:),rwork(:,:,:)

      INTERFACE
         FUNCTION bberf(x)
         real(8) :: bberf,x
         END FUNCTION bberf
      END INTERFACE

      mem   = 0.d0
      memax = 0.d0
      n1    = idisp(myrank)+1
      n2    = idisp(myrank)+ircnt(myrank)
      cnst0 = 2.d0/sqrt(Pi)

      allocate( ftmp(3,MI) ) ; ftmp=0.d0

      do a=1,MI

         ik=Kion(a)

         cnstb = Zps(ik)*cnst0*bps(ik)
         cnst1 = cnstb*bps(ik)*bps(ik)*2.d0/3.d0

         do i=n1,n2

            x  = LL(1,i)*H - Rion(1,a)
            y  = LL(2,i)*H - Rion(2,a)
            z  = LL(3,i)*H - Rion(3,a)
            r2 = x*x+y*y+z*z
            r  = sqrt(r2)
            br = bps(ik)*r

            if ( r>1.d-7 ) then
               tmp1 = Zps(ik)*bberf(br)
               tmp2 = cnstb*exp(-br*br)
               fx = tmp1*x/(r*r2)-tmp2*(x/r2)
               fy = tmp1*y/(r*r2)-tmp2*(y/r2)
               fz = tmp1*z/(r*r2)-tmp2*(z/r2)
            else
               fx = -cnst1*x
               fy = -cnst1*y
               fz = -cnst1*z
            end if

            ftmp(1,a)=ftmp(1,a)+fx*rho(i)
            ftmp(2,a)=ftmp(2,a)+fy*rho(i)
            ftmp(3,a)=ftmp(3,a)+fz*rho(i)

         end do

      end do

      ftmp(1:3,1:MI) = ftmp(1:3,1:MI)*dV

      call mpi_allreduce(ftmp,force,3*MI,mpi_real8,mpi_sum
     &     ,mpi_comm_world,ierr)

!
! --- OH ---
!

      if ( Ndense<=0 .or. Nintp<=0 ) goto 800

      ll_0 = min(0,-Nintp+1)

!- allocate -------------------------------------------
      allocate( rho3(-Mx:Mx,-My:My,-Mz:Mz) )
      rho3=0.d0
      mem=mem+bdreal*size(rho3) ; memax=max(mem,memax)
      allocate( rwork(-Mx:Mx,-My:My,-Mz:Mz) )
      rwork=0.d0
      mem=mem+bdreal*size(rwork) ; memax=max(mem,memax)
!------------------------------------------------------

      do i=n1,n2
         i1=LL(1,i)
         i2=LL(2,i)
         i3=LL(3,i)
         rwork(i1,i2,i3)=rho(i)
      end do
      m=size(rho3)
      call mpi_allreduce(rwork,rho3,m,mpi_real8,mpi_sum
     &     ,mpi_comm_world,ierr)

!- deallocate -------------------
      mem=mem-bdreal*size(rwork)
      deallocate( rwork )
!--------------------------------

      ftmp(:,:) = 0.d0
      cnst1     = 2.d0/3.d0

      do a=1,MI

         ik    = Kion(a)
         rps2  = Rps(1,ik)*Rps(1,ik)
         a2    = aps(ik)*aps(ik)
         b2    = bps(ik)*bps(ik)
         cnsta = cnst0*aps(ik)
         cnstb = cnst0*bps(ik)

         Rx=Rion(1,a)
         Ry=Rion(2,a)
         Rz=Rion(3,a)

         ic1=nint(Rx/H)
         ic2=nint(Ry/H)
         ic3=nint(Rz/H)

         do i=1,M_grid_ion

            i1 = map_grid_ion(1,i)
            i2 = map_grid_ion(2,i)
            i3 = map_grid_ion(3,i)

            id1=i1+ic1
            id2=i2+ic2
            id3=i3+ic3

            if ( id1<Mx0 .or. Mx1<id1 .or.
     &           id2<My0 .or. My1<id2 .or.
     &           id3<Mz0 .or. Mz1<id3      ) cycle

            do j3=0,Ndense-1
               z=id3*H+j3*H3d-Rz
            do j2=0,Ndense-1
               y=id2*H+j2*H2d-Ry
            do j1=0,Ndense-1
               x=id1*H+j1*H1d-Rx

               r2=x*x+y*y+z*z

               if ( r2>rps2+1.d-10 ) cycle

               r  = sqrt(r2)
               r3 = r*r2

! ERF + ERF
               if ( r>1.d-7 ) then
                  tmp1 = bberf(aps(ik)*r) - bberf(bps(ik)*r)
                  tmp2 =-cnsta*exp(-a2*r2)+cnstb*exp(-b2*r2)
                  tmp3 = Zps(ik)*( tmp1 + r*tmp2 )
                  fx = tmp3*(x/r3)
                  fy = tmp3*(y/r3)
                  fz = tmp3*(z/r3)
               else
                  tmp3 = Zps(ik)*( cnsta*a2-cnstb*b2 )*cnst1
                  fx   = x*tmp3
                  fy   = y*tmp3
                  fz   = z*tmp3
               end if

               do ii3=ll_0,Nintp
                  iii3=id3+ii3 ; if (iii3<-Mz .or. Mz<iii3) cycle
                  tmp3=Clag3(ii3,j3)*dVd
               do ii2=ll_0,Nintp
                  iii2=id2+ii2 ; if (iii2<-My .or. My<iii2) cycle
                  tmp2=Clag2(ii2,j2)*tmp3
               do ii1=ll_0,Nintp
                  iii1=id1+ii1 ; if (iii1<-Mx .or. Mx<iii1) cycle
                  tmp1=Clag1(ii1,j1)*tmp2
                  ftmp(1,a)=ftmp(1,a)+fx*tmp1*rho3(iii1,iii2,iii3)
                  ftmp(2,a)=ftmp(2,a)+fy*tmp1*rho3(iii1,iii2,iii3)
                  ftmp(3,a)=ftmp(3,a)+fz*tmp1*rho3(iii1,iii2,iii3)
               end do
               end do
               end do

            end do
            end do
            end do

         end do ! i

      end do ! a

      mem=mem-bdreal*size(rho3)
      deallocate( rho3 )

      allocate( rwork2(3,MI) ) ; rwork2=0.d0
      call mpi_allreduce(ftmp,rwork2,3*MI,mpi_real8
     &     ,mpi_sum,mpi_comm_world,ierr)

      force(1:3,1:MI)=force(1:3,1:MI)+rwork2(1:3,1:MI)*dV

      deallocate( rwork2 )

 800  continue

      deallocate( ftmp )

      return
      END SUBROUTINE force_local_allel
