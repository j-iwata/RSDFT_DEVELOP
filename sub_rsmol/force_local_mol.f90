!--------1---------2---------3---------4---------5---------6---------7--
! using density gradients instead of potential gradients (ref: OH book ?)

      SUBROUTINE force_local_mol(force1)
      use global_variables
      implicit none

      real(8),intent(INOUT) :: force1(3,MI)
      integer :: nn1,nn2,ML0,irank,ierr,M_irad,NRc,ir0
      integer :: s,a,ik,i,j,n,lm0,lm1,i1,i2,i3,j1,j2,j3,m,ir,MMr
      integer :: l1,l2,l3,m1,m2,m3,id1,id2,id3,mm1,mm2,mm3
      integer :: ic1,ic2,ic3,ii1,ii2,ii3,k1,k2,k3,ll_0
      integer :: i1_0,i2_0,i3_0,mm,iii1,iii2,iii3,Mx,My,Mz
      integer,allocatable :: idis(:),icnt(:),irad(:,:),LLtmp(:,:)
      real(8),parameter :: ep=1.d-8
      real(8) :: p1,p2,p3,p4,err0,dvdr0,dvdr
      real(8) :: a1,a2,a3,Gr,Gx,Gy,Gz,sum1,sum2,sum3
      real(8) :: cnst0,cnst1,err,dedr0,dedr,maxerr
      real(8) :: tmp0,tmp1,tmp2,tmp3,vx,vy,vz,x,y,z,d1,d2,d3
      real(8) :: r2,r,rps2,Rx,Ry,Rz
      real(8) :: ct0,ct1,et0,et1,ct(2),et(2)
      real(8) :: ctime0,ctime1,etime0,etime1,mem,memax
      real(8),allocatable :: rho3(:,:,:),Vxc3(:,:,:),rho1(:)
      real(8),allocatable :: rwork(:,:,:,:),rtmp(:),dtmp(:,:)
      real(8),allocatable :: ftmp(:,:)

      INTERFACE
         FUNCTION bberf(x)
         real(8) :: bberf
         real(8),intent(IN) :: x
         END FUNCTION bberf
      END INTERFACE

      call watch(ctime0,etime0)
      if (DISP_SWITCH) then
         write(*,*) "--- force_local_mol ---"
      end if

      nn1     = idisp(myrank)+1
      nn2     = idisp(myrank)+ircnt(myrank)
      ML0     = ircnt(myrank)
      ct(:)   = 0.d0
      et(:)   = 0.d0
      mem     = 0.d0
      memax   = 0.d0
      Mx      = ML1
      My      = ML2
      Mz      = ML3

      mm1 = Mz + Nintp + 1
      mm2 = My + Nintp + 1
      mm3 = Mx + Nintp + 1

      force1(1:3,1:MI) = 0.d0

      www(:,:,:,:)=zero
      select case(nspin)
      case default
         do i=nn1,nn2
            www(LL(1,i),LL(2,i),LL(3,i),1)=rho(i,1)
         end do
      case(2)   
         do i=nn1,nn2
            www(LL(1,i),LL(2,i),LL(3,i),1)=rho(i,1)+rho(i,nspin)
         end do
      end select
      call bcset_2d(1,1,Md,0)
      allocate( dtmp(nn1:nn2,3) ) ; dtmp=0.d0 ; mem=mem+bdreal*size(dtmp) ; memax=max(mem,memax)
      do m=1,Md
         d1=nab(m)/H1
         d2=nab(m)/H2
         d3=nab(m)/H3
         do i=nn1,nn2
            i1=LL(1,i)
            i2=LL(2,i)
            i3=LL(3,i)
            dtmp(i,1)=dtmp(i,1)-d1*( www(i1-m,i2,i3,1) - www(i1+m,i2,i3,1) )
            dtmp(i,2)=dtmp(i,2)-d2*( www(i1,i2-m,i3,1) - www(i1,i2+m,i3,1) )
            dtmp(i,3)=dtmp(i,3)-d3*( www(i1,i2,i3-m,1) - www(i1,i2,i3+m,1) )
         end do
      end do
      
      allocate( rtmp(ML) ) ; mem=mem+bdreal*size(rtmp) ; memax=max(mem,memax)

!      rtmp(:)=0.d0
!      do s=1,nspin
!         rtmp(nn1:nn2)=rtmp(nn1:nn2)+rho(nn1:nn2,s)
!      end do
!      call mpi_allgatherv(rtmp(nn1),ML0,mpi_real8,rtmp,ir_grid,id_grid,mpi_real8,comm_grid,ierr)

      allocate( LLtmp(3,ML) ) ; mem=mem+bsintg*size(LLtmp) ; memax=max(mem,memax)

      allocate( idis(0:np_grid-1),icnt(0:np_grid-1) ) ; mem=mem+bsintg*np_grid*2 ; memax=max(mem,memax)
      idis(0:np_grid-1)=3*id_grid(0:np_grid-1)
      icnt(0:np_grid-1)=3*ir_grid(0:np_grid-1)
      call mpi_allgatherv(LL,icnt(myrank_g),mpi_integer,LLtmp,icnt,idis,mpi_integer,comm_grid,ierr)
      mem=mem-bsintg*np_grid*2 ; deallocate(icnt,idis)

!      allocate( rho3(-mm1:mm1,-mm2:mm2,-mm3:mm3) ) ; mem=mem+bdreal*size(rho3) ; memax=max(mem,memax)
!      rho3=0.d0
!      do i=1,ML
!         rho3(LLtmp(1,i),LLtmp(2,i),LLtmp(3,i))=rtmp(i)
!      end do

      allocate( rwork(-mm1:mm1,-mm2:mm2,-mm3:mm3,3) ) ; mem=mem+bdreal*size(rwork) ; memax=max(mem,memax)
      rwork=0.d0
      do j=1,3
         rtmp(:)=0.d0
         rtmp(nn1:nn2)=dtmp(nn1:nn2,j)
         call mpi_allgatherv(rtmp(nn1),ML0,mpi_real8,rtmp,ir_grid,id_grid,mpi_real8,comm_grid,ierr)
         do i=1,ML
            rwork(LLtmp(1,i),LLtmp(2,i),LLtmp(3,i),j)=rtmp(i)
         end do
      end do

      mem=mem-bsintg*size(LLtmp) ; deallocate( LLtmp )

!      allocate( rho1(nn1:nn2) ) ; mem=mem+bdreal*size(rho1) ; memax=max(mem,memax)
!      rho1(nn1:nn2)=rtmp(nn1:nn2)

      mem=mem-bdreal*size(rtmp) ; deallocate( rtmp )

      allocate( ftmp(3,MI) ) ; mem=mem+bdreal*size(ftmp) ; memax=max(mem,memax)
      ftmp=0.d0

      select case(pselect)
      case(1,2)

         ll_0 = min(0,-Nintp+1)

!- allocate ----------------------------------------------------------
         allocate( irad(0:3000,MKI) ) ; irad=0
         mem=mem+bsintg*3001*MKI ; memax=max(mem,memax)
!---------------------------------------------------------------------

         M_irad=0
         do ik=1,MKI
            NRc=min( 3000, NRcloc(ik) )
            m=0
            irad(0,ik)=1
            do ir=1,NRc
               m=int(100.d0*rad(ir,ik))+1
               irad( m,ik )=ir
            end do
            ir=irad(0,ik)
            do i=1,m
               if ( irad(i,ik)==0 ) then
                  irad(i,ik)=ir
                  cycle
               end if
               ir=irad(i,ik)
            end do
            irad(m+1:,ik)=ir
            M_irad=max(M_irad,m)
         end do

         cnst0 = 2.d0/sqrt(Pi)
         cnst1 = 2.d0/3.d0*cnst0

         maxerr = 0.d0

         do a=1,MI

            ik    = Kion(a)
            rps2  = Rcloc(ik)*Rcloc(ik)
            NRc   = NRcloc(ik)

            p1 =-Zps(ik)*parloc(1,ik)
            p2 = sqrt(parloc(2,ik))
            p3 =-Zps(ik)*parloc(3,ik)
            p4 = sqrt(parloc(4,ik))

            Rx = asi(1,a)
            Ry = asi(2,a)
            Rz = asi(3,a)

            ic1 = nint(Rx/H1)
            ic2 = nint(Ry/H2)
            ic3 = nint(Rz/H3)

            do i=nn1,nn2
               i1=LL(1,i)
               i2=LL(2,i)
               i3=LL(3,i)
               x=i1*H1-Rx
               y=i2*H2-Ry
               z=i3*H3-Rz
               r=sqrt(x*x+y*y+z*z)
               if ( r<1.d-9 ) then
                  tmp0=-cnst0*(p1*p2+p3*p4)
                  tmp3=-cnst1*(p1*p2*p2*p2+p3*p4*p4*p4)
                  vx=x*tmp3
                  vy=y*tmp3
                  vz=z*tmp3
               else
                  tmp0=-(p1*bberf(p2*r)+p3*bberf(p4*r))/r
                  tmp3=-p1*bberf(p2*r)-p3*bberf(p4*r)+r*cnst0*( p1*p2*exp(-p2*p2*r*r)+p3*p4*exp(-p4*p4*r*r) )
                  vx=tmp3*(x/(r*r*r))
                  vy=tmp3*(y/(r*r*r))
                  vz=tmp3*(z/(r*r*r))
               end if
!               ftmp(1,a)=ftmp(1,a)+rho1(i)*vx
!               ftmp(2,a)=ftmp(2,a)+rho1(i)*vy
!               ftmp(3,a)=ftmp(3,a)+rho1(i)*vz
               ftmp(1,a)=ftmp(1,a)+dtmp(i,1)*tmp0
               ftmp(2,a)=ftmp(2,a)+dtmp(i,2)*tmp0
               ftmp(3,a)=ftmp(3,a)+dtmp(i,3)*tmp0
            end do

            do i=1,M_grid_ion_lc

               i1 = map_grid_ion_lc(1,i)
               i2 = map_grid_ion_lc(2,i)
               i3 = map_grid_ion_lc(3,i)

               id1 = i1 + ic1
               id2 = i2 + ic2
               id3 = i3 + ic3

               if ( PPP(id1,id2,id3)==myrank_g ) then

                  do j3=0,Ndense-1
                     z=id3*H3+j3*H3d-Rz
                  do j2=0,Ndense-1
                     y=id2*H2+j2*H2d-Ry
                  do j1=0,Ndense-1
                     x=id1*H1+j1*H1d-Rx

                     r2 = x*x+y*y+z*z

                     if ( r2>rps2+1.d-10 ) cycle

                     r  = sqrt(r2)
                     tmp0 = 0.d0

                     ir0=irad( int(100.d0*r),ik )
                     do ir=ir0,NRc
                        if ( r<rad(ir,ik) ) exit
                     end do

!                     if( ir<=2 )then
!                        vx=dvqls(2,ik)*(x/rad(2,ik))
!                        vy=dvqls(2,ik)*(y/rad(2,ik))
!                        vz=dvqls(2,ik)*(z/rad(2,ik))
!                        err0=0.d0
!                     else
                        err0=1.d10
                        do m=1,20
                           m1=max(1,ir-m)
                           m2=min(ir+m,NRc)
                           call polint(rad(m1,ik),vqls(m1,ik),m2-m1+1,r,dvdr,err)
!                           call polint(rad(m1,ik),dvqls(m1,ik),m2-m1+1,r,dvdr,err)
                           if ( abs(err)<err0 ) then
                              dvdr0=dvdr
                              tmp0=-dvdr
                              err0=abs(err)
                              if ( err0<ep ) exit
                           end if
                        end do
!                        vx=dvdr0*(x/r)
!                        vy=dvdr0*(y/r)
!                        vz=dvdr0*(z/r)
!                     end if
                     maxerr=max(maxerr,err0)

                     do ii3=ll_0,Nintp
                        iii3=id3+ii3
!                        tmp3=Clag3(ii3,j3)*dVd
                        tmp3=Clag3(ii3,j3)*dVd*tmp0
                     do ii2=ll_0,Nintp
                        iii2=id2+ii2
                        tmp2=Clag2(ii2,j2)*tmp3
                     do ii1=ll_0,Nintp
                        iii1=id1+ii1
!                        tmp1=Clag1(ii1,j1)*tmp2*rho3(iii1,iii2,iii3)
!                        ftmp(1,a)=ftmp(1,a)+tmp1*vx
!                        ftmp(2,a)=ftmp(2,a)+tmp1*vy
!                        ftmp(3,a)=ftmp(3,a)+tmp1*vz
                        ftmp(1,a)=ftmp(1,a)+Clag1(ii1,j1)*tmp2*rwork(iii1,iii2,iii3,1)
                        ftmp(2,a)=ftmp(2,a)+Clag1(ii1,j1)*tmp2*rwork(iii1,iii2,iii3,2)
                        ftmp(3,a)=ftmp(3,a)+Clag1(ii1,j1)*tmp2*rwork(iii1,iii2,iii3,3)
                     end do
                     end do
                     end do

                  end do
                  end do
                  end do

               end if ! PPP

            end do ! i M_grid_ion_lc

         end do ! a

         ftmp(1:3,1:MI)=ftmp(1:3,1:MI)*dV

         call mpi_allreduce(ftmp(1,1),force1(1,1),3*MI,mpi_real8,mpi_sum,comm_grid,ierr)

         mem=mem-bsintg*size(irad)  ; deallocate( irad )
         
      end select

      mem=mem-bdreal*size(rwork) ; deallocate( rwork )
      mem=mem-bdreal*size(ftmp) ; deallocate( ftmp )
!      mem=mem-bdreal*size(rho1) ; deallocate( rho1 )
!      mem=mem-bdreal*size(rho3) ; deallocate( rho3 )
      mem=mem-bdreal*size(dtmp) ; deallocate( dtmp )

 800  Max_mem_inst = max( Max_mem_inst, Max_mem+memax )

      call watch(ctime1,etime1)
      if(DISP_SWITCH)then
         write(*,*) "mem(MB)=",mem,memax*B2MB,Max_mem_inst*B2MB
         write(*,*) "time=",ctime1-ctime0,etime1-etime0
      end if

      return

 900  call stop_program

      END SUBROUTINE force_local_mol
