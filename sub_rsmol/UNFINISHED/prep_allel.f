!--------1---------2---------3---------4---------5---------6---------7--

      SUBROUTINE prep_allel
      use global_variables
      implicit none
      integer :: a,i,i1,i2,i3,ik,j,MG_0,MG_1,n1,n2,ierr
      integer :: mm1,mm2,mm3,mm4,ML0,ll_0,l1,l2,l3
      integer :: j1,j2,j3,ic1,ic2,ic3,io3,k1,k2,k3,id1,id2,id3
      integer :: ii1,ii2,ii3,m1,m2,m3,i1_0,i2_0,i3_0
      integer :: iii1,iii2,iii3,m
      real(8) :: a1,a2,a3,ct(3),et(3),mem,memax,Gr
      real(8) :: cnst0,cnst1,cnst2,cnsta,cnstb
      real(8) :: ctime0,ctime1,etime0,etime1,ct0,ct1,et0,et1
      real(8) :: b,c,r,g,r2,rps2,x,y,z,d1,d2,d3,b2
      real(8) :: tmp1,tmp2,tmp3,Rx,Ry,Rz,v
      real(8),allocatable :: rwork(:,:,:,:)
      logical,save :: flag_allel=.true.

      INTERFACE
         FUNCTION bberf(x)
         real(8) :: bberf,x
         END FUNCTION bberf
      END INTERFACE

      call watch(ctime0,etime0)

      if( DISP_SWITCH )then
         write(*,'(a60," prep_allel")') repeat("-",60)
      end if

      if ( pselect/=0 ) then
         write(*,*) "pselect=",pselect,myrank
         call stop_program
      end if

      ct(:) = 0.d0
      et(:) = 0.d0
      mem   = 0.d0
      memax = 0.d0
      n1    = idisp(myrank)+1
      n2    = idisp(myrank)+ircnt(myrank)
      ML0   = ircnt(myrank)
      cnst0 = 2.d0/sqrt(Pi)

      if ( flag_allel ) then

         Rnucl(:)=0.d0
         Rnucl(1:MKI)=0.5d0*H

         Rps(:,:)=0.d0

         if (DISP_SWITCH) then
            write(*,'(1x,a3,2x,4a10)')
     &           "ik","Zps(ik)","bps(ik)","aps(ik)","Rps"
         end if

         do ik=1,MKI
            do i=1,10000
               r=i*0.01d0
               Rps(1,ik)=r
               a1=Zps(ik)*(1.d0-bberf(bps(ik)*r))/r
               if ( a1<1.d-15 ) exit
            end do
         end do

         if(DISP_SWITCH)then
            do ik=1,MKI
               write(*,'(1x,i3,2x,4f10.5)')
     &              ik,Zps(ik),bps(ik),aps(ik),Rps(1,ik)
            end do
         end if

         ZNel=0.d0
         do a=1,MI
            ZNel=ZNel+Zps(Kion(a))
         end do

         flag_allel=.false.

      end if

!
! --- local potential ---
!

      call gv_alloc("prep_ps_1")

      Vion(:)=0.d0

      do a=1,MI

         ik=Kion(a)

! USC
         cnst1=Zps(ik)*1.5d0/Rnucl(ik)
         cnst2=Zps(ik)*0.5d0/Rnucl(ik)**3
! ERF
c         cnst1=Zps(ik)*cnst0*bps(ik)
c         cnst2=cnst1*bps(ik)*bps(ik)/3.d0

         do i=n1,n2

            x  = LL(1,i)*H - Rion(1,a)
            y  = LL(2,i)*H - Rion(2,a)
            z  = LL(3,i)*H - Rion(3,a)
            r2 = x*x+y*y+z*z
            r  = sqrt(r2)

! USC
            if ( r>Rnucl(ik) ) then
               Vion(i) = Vion(i) - Zps(ik)/r
            else
               Vion(i) = Vion(i) - cnst1 + cnst2*r2
            end if
! ERF
c            if ( r>1.d-7 ) then
c               Vion(i) = Vion(i) - Zps(ik)*bberf(bps(ik)*r)/r
c            else
c               Vion(i) = Vion(i) - cnst1 + cnst2*r2
c            end if

         end do

      end do

!
! --- OH Method --------------------------------------------------------
!

      if ( Ndense<=0 .or. Nintp<=0 ) then
         if ( DISP_SWITCH ) then
            write(*,*) "pselect=",pselect
            write(*,*) "Ndense =",Ndense
            write(*,*) "Nintp =",Nintp
         end if
         goto 800
      end if

      if ( Nintp>Md ) then
         if (DISP_SWITCH) then
            write(*,*) "Nintp=",Nintp
            write(*,*) "Nintp is replaced by",Md
         end if
         Nintp=Md
      end if

      H1d  = H/Ndense
      H2d  = H/Ndense
      H3d  = H/Ndense
      dVd  = (H1d*H2d*H3d)/(H*H*H)

      if(DISP_SWITCH)then
         write(*,'(1x,"H =",3f20.12)') H
         write(*,'(1x,"H1d,H2d,H3d=",3f20.12)') H1d,H2d,H3d
         write(*,*) "Ndense =",Ndense
         write(*,*) "Nintp  =",Nintp
      end if

      ll_0 = min(0,-Nintp+1)

      r = maxval( Rps(1,1:MKI) )

      call Make_MinimalBox(r,mm1,mm2,mm3,mm4)

c      mm1_ps = maxval( abs(mcube_grid_ion(:,1)) )
c      mm2_ps = maxval( abs(mcube_grid_ion(:,2)) )
c      mm3_ps = maxval( abs(mcube_grid_ion(:,3)) )

!--- Coeficients of Lagrange Interporation ---

      if ( .not.allocated(Clag1) ) then
         allocate( Clag1(ll_0:Nintp,0:Ndense-1) ) ; Clag1=0.d0
         allocate( Clag2(ll_0:Nintp,0:Ndense-1) ) ; Clag2=0.d0
         allocate( Clag3(ll_0:Nintp,0:Ndense-1) ) ; Clag3=0.d0
         Max_mem=Max_mem+bdreal*Ndense*(Nintp-ll_0+1)*3
         do j1=0,Ndense-1
         do i1=ll_0,Nintp
            Clag1(i1,j1)=1.d0
            do i2=ll_0,Nintp
               if ( i2==i1 ) cycle
               x=dble(j1-i2*Ndense)
               y=dble(i1-i2)*Ndense
               Clag1(i1,j1)=Clag1(i1,j1)*(x/y)
            end do
         end do
         end do
         Clag2(:,:)=Clag1(:,:)
         Clag3(:,:)=Clag1(:,:)
      end if

!- allocate -------------------------------------------
      allocate( rwork(-Mx:Mx,-My:My,-Mz:Mz,2) )
      rwork=0.d0
      mem=mem+bdreal*size(rwork)
      memax=max(mem,memax)
!------------------------------------------------------

      do a=1,MI

         ik    = Kion(a)
         rps2  = Rps(1,ik)*Rps(1,ik)
         a2    = aps(ik)*aps(ik)/3.d0
         b2    = bps(ik)*bps(ik)/3.d0
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

               r = sqrt(r2)
! ERF + ERF
               if ( r>1.d-7 ) then
                  v = -Zps(ik)*(bberf(aps(ik)*r)-bberf(bps(ik)*r))/r
               else
                  v = -Zps(ik)*( cnsta*(1.d0-a2*r2)
     &                          -cnstb*(1.d0-b2*r2) )
               end if

               do ii3=ll_0,Nintp
                  iii3=id3+ii3 ; if (iii3<-Mz .or. Mz<iii3) cycle
                  tmp3=Clag3(ii3,j3)*v*dVd
               do ii2=ll_0,Nintp
                  iii2=id2+ii2 ; if (iii2<-My .or. My<iii2) cycle
                  tmp2=Clag2(ii2,j2)*tmp3
               do ii1=ll_0,Nintp
                  iii1=id1+ii1 ; if (iii1<-Mx .or. Mx<iii1) cycle
                  tmp1=Clag1(ii1,j1)*tmp2
                  rwork(iii1,iii2,iii3,1)
     &                 =rwork(iii1,iii2,iii3,1)+tmp1
               end do
               end do
               end do

            end do
            end do
            end do

         end do ! i

      end do ! a

      m=size(rwork)/2
      call mpi_allreduce(rwork(-Mx,-My,-Mz,1),rwork(-Mx,-My,-Mz,2)
     &     ,m,mpi_real8,mpi_sum,mpi_comm_world,ierr)

      do i=n1,n2
         i1=LL(1,i)
         i2=LL(2,i)
         i3=LL(3,i)
         Vion(i)=Vion(i)+rwork(i1,i2,i3,2)
      end do

c      rewind 10
c      do i=-Mz,Mz
c         write(10,*) i*H,rwork(0,0,i,2)
c      end do
c      goto 900

      mem=mem-bdreal*size(rwork)
      deallocate( rwork )

 800  Max_mem_inst = max(Max_mem_inst,Max_mem+memax)

      call watch(ctime1,etime1)

      if (DISP_SWITCH) then
         write(*,*) "TIME(PREP_ALLEL)="
     &        ,ctime1-ctime0,etime1-etime0
         write(*,'(1x,"MEM(MB)",3f15.5)')
     &        mem,memax*B2MB,Max_mem_inst*B2MB
      end if

      return

 900  call stop_program

      END SUBROUTINE prep_allel
