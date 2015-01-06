!--------1---------2---------3---------4---------5---------6---------7--

      SUBROUTINE prep_ps_mol_sym(iswitch_pp)
      use global_variables
      implicit none

      integer,intent(IN) :: iswitch_pp
      logical,save :: iswitch_prep_oh=.true.
      integer,allocatable :: LLL0(:,:,:),PPP0(:,:,:),LL0(:,:)
      integer,allocatable :: iloc(:),Jtmp(:,:,:)
      integer,allocatable :: ichk0(:,:),ichk1(:,:)
      integer,allocatable :: itmp(:,:),itmp1(:),itmp2(:,:),itmp3(:,:)
      integer,allocatable :: icheck_tmp1(:),icheck_tmp2(:)
      integer,allocatable :: icheck_tmp3(:,:,:),icheck_tmp4(:,:,:)
      integer,allocatable :: ireq(:)
      integer,allocatable :: recvmap0(:,:),sendmap0(:,:)
      integer,allocatable :: Jxyz2(:,:,:)
      integer,allocatable :: mmap_tmp(:),amap_tmp(:),lmap_tmp(:)
      integer,allocatable :: iorbmap_tmp(:)
      integer,allocatable :: iuV_tmp(:),JJ_tmp(:,:),Mps_tmp(:)
      integer,allocatable :: ircntd(:),LLd(:,:)
      integer,allocatable :: irad(:,:)
      integer,allocatable :: KKK(:,:,:),KKK0(:,:,:)
      integer :: nreqs,nreqr,istatus(MPI_STATUS_SIZE)
      integer :: n1,n2,M_irad
      integer :: a,i,k,ik,ir,m,n,m1,m2,m3,m4,dm,nx,ny,nz,iorb,morb,mm
      integer :: i1,i2,i3,j1,j2,j3,j,L,lm,NRc,lma,MMr,ierr,mtmp(3,2)
      integer :: l1,l2,l3,k1,k2,k3,ix,iy,iz,jx,jy,jz,i0,ii
      integer :: ii1,ii2,ii3,ic1,ic2,ic3,ll_0,ir0,MLK
      integer :: io1,io2,io3,jo1,jo2,jo3
      integer :: minLL1,minLL2,minLL3,maxLL1,maxLL2,maxLL3
      integer :: lma1,irank,iloc1(1)
      integer :: isum0(6),isum1(6)
      real(8),parameter :: ep=1.d-8, eps=1.d-10
      real(8),allocatable :: work(:),uVtmp(:)
      real(8),allocatable :: ww2(:,:,:,:),ww3(:,:,:,:)
      real(8),allocatable :: uVk_tmp(:)
      real(8) :: x,y,z,r,v,err,rps2,rps2_0,r2,c1,c2
      real(8) :: xk,xl,tmp1,tmp2,tmp3,xd,yd,zd,rd
      real(8) :: uVrl,maxerr,v0,err0,errmax,errav
      real(8) :: ctime0,ctime1,etime0,etime1
      real(8) :: ct0,ct1,et0,et1,ct(9),et(9)
      real(8) :: mem,memax
      logical,allocatable :: mask(:)
      logical,allocatable :: lcheck_tmp1(:),lcheck_tmp2(:)
      logical :: iflag

      integer :: id1,id2,id3,iii1,iii2,iii3,ML0,irlma,nreq
!      integer :: a1b,a2b,a3b,b1b,b2b,b3b
      integer :: mm1,mm2,mm3,mm4,m8,m9,m0,ibuf(2,3)
      integer :: MMJJ_0,lma0,nzlma_0,isym,Mx,My,Mz
      integer,allocatable :: jtmp3(:,:,:),mtmp3(:)
      integer,allocatable :: nl_rank_map_tmp(:),maps_tmp(:,:)
      integer,allocatable :: lma_nsend_tmp(:)
      real(8) :: rc2,Rx,Ry,Rz,p1,p2,p3,p4,tmp0,H
      real(8),allocatable :: rwork(:,:,:,:),rwork1(:,:,:)
      real(8),allocatable :: rtmp3(:,:),rrr(:,:,:,:),sss(:,:,:,:)

      integer :: LL_tmp(3)
!      integer,allocatable :: PPP(:,:,:)

      INTERFACE
         FUNCTION bberf(x)
         real(8) :: bberf
         real(8),intent(IN) :: x
         END FUNCTION bberf
         FUNCTION Ylm(x,y,z,l,m)
         real(8) :: Ylm
         real(8),intent(IN) :: x,y,z
         integer,intent(IN) :: l,m
         END FUNCTION Ylm
      END INTERFACE

      call bwatch(ctime0,etime0)

      if (DISP_SWITCH) then
         write(*,'(a60," prep_ps_mol_sym")') repeat("-",60)
      end if

      n1    = idisp(myrank)+1
      n2    = idisp(myrank)+ircnt(myrank)
      ML0   = ircnt(myrank)
      mem   = 0.d0
      memax = 0.d0

      Mx=ML1+Md
      My=ML2+Md
      Mz=ML3+Md

      H=H1

!
! --- Parameters for OH method ---
!

      if ( iswitch_prep_oh ) then

         iswitch_prep_oh = .false.

         if ( pselect==3 .or. pselect==2 .or. Ndense<=0  .or. Nintp<=0 ) then
            if ( DISP_SWITCH ) then
               write(*,*) "pselect=",pselect
               write(*,*) "Ndense =",Ndense
               write(*,*) "(Ndense,Nintp) is set to (1,0)."
            end if
            Ndense=1
            Nintp =0
         end if

         if ( Nintp>Md ) then
            write(*,*) "Nintp=",Nintp
            write(*,*) "Nintp is replaced by",Md
            Nintp=Md
         end if

         H1d  = H1/Ndense
         H2d  = H2/Ndense
         H3d  = H3/Ndense
         dVd  = (H1d*H2d*H3d)/(H1*H2*H3)
         ll_0 = min(0,-Nintp+1)

         if(DISP_SWITCH)then
            write(*,'(1x,"H1 ,H2 ,H3 =",3f20.12)') H1,H2,H3
            write(*,'(1x,"H1d,H2d,H3d=",3f20.12)') H1d,H2d,H3d
            write(*,*) "Ndense       =",Ndense
            write(*,*) "Nintp, ll_0  =",Nintp,ll_0
         end if

! --- Coeficients of Lagrange Interporation ---

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

      end if ! iswitch_prep_oh

! ----------------------------------------------------------------------
! --- local potential --------------------------------------------------
! ----------------------------------------------------------------------

      if (DISP_SWITCH) write(*,*) "---local potential---"

      call gv_alloc("Vion")

      Vion(:) = 0.d0

      r=maxval( Rcloc(1:MKI) )

      call Make_MinimalBox(r,mm1,mm2,mm3,mm4)

      M_grid_ion_lc = M_grid_ion

!- allocate -----------------------------
      call gv_alloc("mcube_grid_ion_lc")
!----------------------------------------

      mcube_grid_ion_lc(:,:) = mcube_grid_ion(:,:)
      map_grid_ion_lc(:,:)   = map_grid_ion(:,:)

      mm1_lc = maxval( abs(mcube_grid_ion_lc(:,1)) )
      mm2_lc = maxval( abs(mcube_grid_ion_lc(:,2)) )
      mm3_lc = maxval( abs(mcube_grid_ion_lc(:,3)) )

      if (DISP_SWITCH) then
         write(*,*) "mm1_lc=",mm1_lc
         write(*,*) "mm2_lc=",mm2_lc
         write(*,*) "mm3_lc=",mm3_lc
         write(*,*) "M_grid_ion_lc=",M_grid_ion_lc
      end if

      mm1 = mm1_lc + Nintp + 1
      mm2 = mm2_lc + Nintp + 1
      mm3 = mm3_lc + Nintp + 1

!- allocate -------------------------------------------------------
      allocate( rwork(-Mx:Mx,-My:My,-Mz:Mz,2)   ) ; rwork =0.d0
      allocate( irad(0:3000,MKI) ) ; irad=0
      mem=mem+bdreal*(2*Mx+1)*(2*My+1)*(2*Mz+1)*2
      mem=mem+bsintg*3001*MKI
      memax=max(mem,memax)
!-----------------------------------------------------------------

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

      m0     = size(rwork)/2
      c1     = 2.d0/sqrt(Pi)
      maxerr = 0.d0

      call bwatch(ct0,et0)

      do a=1,MI

         ik  = Kion(a)
         NRc = NRcloc(ik)
         Rc2 = Rcloc(ik)**2

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

         do i=n1,n2
            i1=LL(1,i)
            i2=LL(2,i)
            i3=LL(3,i)
            x=i1*H-Rx
            y=i2*H-Ry
            z=i3*H-Rz
            r=sqrt(x*x+y*y+z*z)
            if ( r<1.d-9 ) then
               Vion(i)=Vion(i)+c1*(p1*p2+p3*p4)
            else
               Vion(i)=Vion(i)+( p1*bberf(p2*r)+p3*bberf(p4*r) )/r
            end if
         end do

         do i=1,M_grid_ion_lc

            i1 = map_grid_ion_lc(1,i)
            i2 = map_grid_ion_lc(2,i)
            i3 = map_grid_ion_lc(3,i)

            id1=i1+ic1
            id2=i2+ic2
            id3=i3+ic3

            if ( abs(id1)>Mx .or. abs(id2)>My .or. abs(id3)>Mz ) cycle

            if ( n1<=LLL(id1,id2,id3) .and. LLL(id1,id2,id3)<=n2 ) then

               do j3=0,Ndense-1
                  z=id3*H3+j3*H3d-Rz
               do j2=0,Ndense-1
                  y=id2*H2+j2*H2d-Ry
               do j1=0,Ndense-1
                  x=id1*H1+j1*H1d-Rx

                  r2=x*x+y*y+z*z

                  if ( r2>Rc2+1.d-10 ) cycle

                  r  = sqrt(r2)
                  v0 = 0.d0

                  ir0=irad( int(100.d0*r),ik )
                  do ir=ir0,NRc
                     if ( r<rad(ir,ik) ) exit
                  end do

                  err0=1.d10
                  do mm=1,20
                     m8=max(1,ir-mm)
                     m9=min(ir+mm,NRc)
                     call polint(rad(m8,ik),vqls(m8,ik),m9-m8+1,r,v,err)
                     if ( abs(err)<err0 ) then
                        v0=v
                        err0=abs(err)
                        if ( err0<ep ) exit
                     end if
                  end do
                  maxerr=max(maxerr,err0)

                  do ii3=ll_0,Nintp
                     iii3=id3+ii3
                     tmp3=Clag3(ii3,j3)*v0*dVd
                  do ii2=ll_0,Nintp
                     iii2=id2+ii2
                     tmp2=Clag2(ii2,j2)*tmp3
                  do ii1=ll_0,Nintp
                     iii1=id1+ii1
                     rwork(iii1,iii2,iii3,1)=rwork(iii1,iii2,iii3,1)+Clag1(ii1,j1)*tmp2
                  end do
                  end do
                  end do

               end do
               end do
               end do

            end if ! LLL

         end do ! i

      end do ! a

      call mpi_allreduce(rwork(:,:,:,1),rwork(:,:,:,2),m0,mpi_real8,mpi_sum,comm_grid,ierr)

      call bwatch(ct1,et1)

      if (DISP_SWITCH) then
         m=count(rwork(:,:,:,2)/=0.d0)
         write(*,*) "count(rwork/=0),ML=",m,ML
         if ( m>ML ) write(*,*) "*** WARNING!!! ***"
         write(*,*) "time(local)=",ct1-ct0,et1-et0
      end if

      do i=n1,n2
         i1=LL(1,i)
         i2=LL(2,i)
         i3=LL(3,i)
         Vion(i)=Vion(i)+rwork(i1,i2,i3,2)
      end do

!- deallocate ----------------------------------
      mem=mem-bdreal*size(rwork)
      deallocate( rwork )
!-----------------------------------------------


! ----------------------------------------------------------------------
! --- Nonlocal ---------------------------------------------------------
! ----------------------------------------------------------------------

      if (DISP_SWITCH) write(*,*) "---nonlocal potential---"

      r=maxval( Rps(:,:) )

      call Make_MinimalBox(r,mm1,mm2,mm3,mm4)

      mm1_ps = maxval( abs(mcube_grid_ion(:,1)) )
      mm2_ps = maxval( abs(mcube_grid_ion(:,2)) )
      mm3_ps = maxval( abs(mcube_grid_ion(:,3)) )

      mm1 = mm1_ps + Nintp + 1
      mm2 = mm2_ps + Nintp + 1
      mm3 = mm3_ps + Nintp + 1

      if (DISP_SWITCH) then
         write(*,*) "mm1_ps=",mm1_ps
         write(*,*) "mm2_ps=",mm2_ps
         write(*,*) "mm3_ps=",mm3_ps
         write(*,*) "M_grid_ion=",M_grid_ion
         write(*,*) "mm1=",mm1
         write(*,*) "mm2=",mm2
         write(*,*) "mm3=",mm3
         write(*,*) "(2*mm1+1)*(2*mm2+1)*(2*mm3+1)=",(2*mm1+1)*(2*mm2+1)*(2*mm3+1)
      end if

      if ( Mlma<np_grid ) then
         nzlma_0 = Mlma
      else
         nzlma_0 = min(Mlma/np_grid*27,Mlma)
      end if
      if ( DISP_SWITCH ) then
         write(*,*) "nzlma(guess),Mlma =",nzlma_0,Mlma
      end if

      allocate( icheck_tmp3(-mm1:mm1,-mm2:mm2,-mm3:mm3) )
      mem=mem+bsintg*(2*mm1+1)*(2*mm2+1)*(2*mm3+1) ; memax=max(mem,memax)
      icheck_tmp3=0
      do i=1,M_grid_ion
         i1=map_grid_ion(1,i)
         i2=map_grid_ion(2,i)
         i3=map_grid_ion(3,i)
         do j3=ll_0-1,Nintp+1
         do j2=ll_0-1,Nintp+1
         do j1=ll_0-1,Nintp+1
            icheck_tmp3(i1+j1,i2+j2,i3+j3)=1
         end do
         end do
         end do
      end do
      j=sum(icheck_tmp3)
      mem=mem-bsintg*size(icheck_tmp3) ; deallocate( icheck_tmp3 )

      MMJJ_0 = j

      if ( DISP_SWITCH ) then
         write(*,*) "MMJJ_0,ML0,ML =",MMJJ_0,ML0,ML
      end if

!- allocate ------------------------------------------------------------
         L=maxval(lo)
         n=maxval(norb)
         m=(2*mm1+1)*(2*mm2+1)*(2*mm3+1)
         allocate( icheck_tmp3(MI,n,2*L+1)           ) ; icheck_tmp3=0
         allocate( icheck_tmp1(0:np_grid-1)         ) ; icheck_tmp1=0
         allocate( rtmp3(MMJJ_0,nzlma_0) ) ; rtmp3=0.d0
         allocate( rrr(-mm1:mm1,-mm2:mm2,-mm3:mm3,2*L+1) ) ; rrr=0.d0
         allocate( sss(-mm1:mm1,-mm2:mm2,-mm3:mm3,2*L+1) ) ; sss=0.d0
         allocate( uVtmp(-L:L) )               ; uVtmp=0.d0
         allocate( jtmp3(6,MMJJ_0,nzlma_0) ) ; jtmp3=-100000
         allocate( mtmp3(nzlma_0) ) ; mtmp3=0
         mem=mem+bsintg*MI*n*(2*L+1)
         mem=mem+bsintg*np_grid
         mem=mem+bdreal*MMJJ_0*nzlma_0
         mem=mem+bdreal*m*(2*L+1)*2
         mem=mem+bdreal*(2*L+1)
         mem=mem+bsintg*6*MMJJ_0*nzlma_0
         mem=mem+bsintg*nzlma_0
         allocate( PPP(-Mx:Mx,-My:My,-Mz:Mz) ) ; PPP=-1
         mem=mem+bsintg*(2*Mx+1)*(2*My+1)*(2*Mz+1)
         memax=max(mem,memax)
!-----------------------------------------------------------------------

      M_irad=0
      do ik=1,MKI
         NRc=maxval( NRps(:,ik) )
         NRc=min( 3000, NRc )
         m=0
         irad(0,ik)=1
         do ir=1,NRc
            m=int(100.d0*rad1(ir,ik))+1
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

      do irank=0,np_grid-1
         do i=id_grid(irank)+1,id_grid(irank)+ir_grid(irank)
            do isym=1,nsym
               LL_tmp(1:3)=matmul( rga(1:3,1:3,isym),LL(1:3,i) )
               PPP( LL_tmp(1),LL_tmp(2),LL_tmp(3) )=irank
            end do
         end do
      end do

      maxerr             = 0
      lma                = 0
      icheck_tmp3(:,:,:) = 0
      MMJJ               = 0
      nzlma              = 0
      lma0               = 0
      m0                 = size(sss)

      call bwatch(ct0,et0)

      do a=1,MI

         Rx = asi(1,a)
         Ry = asi(2,a)
         Rz = asi(3,a)

         ic1 = nint(Rx/H1)
         ic2 = nint(Ry/H2)
         ic3 = nint(Rz/H3)

         ik = Kion(a)

         do iorb=1,norb(ik)

            Rc2            = Rps(iorb,ik)**2
            NRc            = NRps(iorb,ik)
            L              = lo(iorb,ik)
            rrr(:,:,:,:)   = zero
            sss(:,:,:,:)   = zero
            icheck_tmp1(:) = 0

            do i=1,M_grid_ion

               i1 = map_grid_ion(1,i)
               i2 = map_grid_ion(2,i)
               i3 = map_grid_ion(3,i)

               id1 = ic1 + i1
               id2 = ic2 + i2
               id3 = ic3 + i3

               do ii3=ll_0,Nintp
                  iii3=id3+ii3
               do ii2=ll_0,Nintp
                  iii2=id2+ii2
               do ii1=ll_0,Nintp
                  iii1=id1+ii1
                  irank=PPP(iii1,iii2,iii3)
                  if ( irank<0 .or. np_grid<=irank ) cycle
                  icheck_tmp1(irank)=icheck_tmp1(irank)+1
               end do
               end do
               end do

               if ( n1<=LLL(id1,id2,id3) .and. LLL(id1,id2,id3)<=n2 ) then

                  do j3=0,Ndense-1
                     z=id3*H3+j3*H3d-Rz
                  do j2=0,Ndense-1
                     y=id2*H2+j2*H2d-Ry
                  do j1=0,Ndense-1
                     x=id1*H1+j1*H1d-Rx

                     r2=x*x+y*y+z*z

                     if ( r2>Rc2+1.d-10 ) cycle

                     r        = sqrt(r2)
                     v0       = 0.d0
                     uVtmp(:) = 0.d0
                     err0     = 0.d0

                     if ( abs(x)>1.d-14 .or. abs(y)>1.d-14 .or. abs(z)>1.d-14 .or. L==0 ) then
                        ir0=irad( int(100.d0*r),ik )
                        do ir=ir0,NRc
                           if ( r<rad1(ir,ik) ) exit
                        end do
                        if ( ir<=2 ) then
                           if ( pselect==1 ) then
                              v0=viod(2,iorb,ik)/rad1(2,ik)
                           else
                              v0=viod(ir,iorb,ik)
                           end if
                           if ( ir<1 ) then
                              write(*,*) "stop : prep_ps"
                              goto 900
                           end if
                        else if ( ir<NRc ) then
                           err0=1.d10
                           do mm=1,20
                              m1=max(1,ir-mm)
                              m2=min(ir+mm,NRc)
                              call polint(rad1(m1,ik),viod(m1,iorb,ik),m2-m1+1,r,v,err)
                              if ( abs(err)<err0 ) then
                                 v0=v
                                 err0=abs(err)
                                 if ( err0<ep ) exit
                              end if
                           end do
                           if ( pselect==1 ) v0=v0/r
                        end if
                        maxerr=max(maxerr,err0)
                        do m=-L,L
                           uVtmp(m)=v0*Ylm(x,y,z,L,m)
                        end do
                     end if

                     if ( v0==0.d0 ) cycle

                     do m=1,2*L+1
                        if (abs(uVtmp(m-1-L))<1.d-10) cycle
                        tmp0=uVtmp(m-1-L)*dVd
                        do ii3=ll_0,Nintp
                           iii3=i3+ii3
                           tmp3=Clag3(ii3,j3)*tmp0
                        do ii2=ll_0,Nintp
                           iii2=i2+ii2
                           tmp2=tmp3*Clag2(ii2,j2)
                        do ii1=ll_0,Nintp
                           iii1=i1+ii1
                           sss(iii1,iii2,iii3,m)=sss(iii1,iii2,iii3,m)+tmp2*Clag1(ii1,j1)
                        end do
                        end do
                        end do
                     end do

                     end do ! j3
                     end do ! j2
                     end do ! j1

                  end if ! LLL

               end do ! i ( 1 - M_grid_ion )


               if ( icheck_tmp1(myrank_g)>0 ) then

                  call mpi_allreduce_nlpp( sss,rrr,m0,icheck_tmp1,ierr,"R" )

                  j1=0
                  do i3=-mm3,mm3
                     id3=i3+ic3
                  do i2=-mm2,mm2
                     id2=i2+ic2
                  do i1=-mm1,mm1
                     id1=i1+ic1
                     if ( PPP(id1,id2,id3)==myrank_g ) then
                        lma = lma0
                        do m=1,2*L+1
                           if ( abs(rrr(i1,i2,i3,m))<1.d-10 ) cycle
                           j1=j1+1
                           do mm=1,2*L+1
                              lma=lma+1
                              rtmp3(j1,lma)=rrr(i1,i2,i3,mm)
                              jtmp3(1,j1,lma)=id1
                              jtmp3(2,j1,lma)=id2
                              jtmp3(3,j1,lma)=id3
                              jtmp3(4,j1,lma)=i1
                              jtmp3(5,j1,lma)=i2
                              jtmp3(6,j1,lma)=i3
                              mtmp3(lma)=j1
                           end do
                           exit
                        end do
                     end if
                  end do
                  end do
                  end do

                  if ( j1>0 ) then
                  
                  MMJJ  = max( MMJJ, j1 )
                  nzlma = lma0+2*L+1

                  if ( MMJJ>MMJJ_0 .or. nzlma>nzlma_0 ) then
                     if(DISP_SWITCH)then
                        write(*,*) "MMJJ,MMJJ_0=",MMJJ,MMJJ_0
                        write(*,*) "nzlma,nzlma_0=",nzlma,nzlma_0
                        write(*,*) "MMJJ_0 and/or nzlma_0 should be set larger"
                     end if
                     goto 900
                  end if

                  lma=lma0
                  do m=1,2*L+1
                     lma=lma+1
                     icheck_tmp3(a,iorb,m)=lma
                  end do

                  lma0=lma

                  end if

               end if ! icheck_tmp1(myrank_g)>0
                  
            end do ! iorb
         end do ! a

         call bwatch(ct1,et1)

         ibuf(1,1)=MMJJ
         ibuf(2,1)=nzlma
         call mpi_allreduce(ibuf(1,1),ibuf(1,2),2,mpi_integer,mpi_min,comm_grid,ierr)
         call mpi_allreduce(ibuf(1,1),ibuf(1,3),2,mpi_integer,mpi_max,comm_grid,ierr)

         if (DISP_SWITCH) then
            write(*,*) "lma=",lma
            write(*,*) "myrank,nzlma,MMJJ",myrank_g,nzlma,MMJJ
            write(*,*) "MMJJ(min,max)  =",ibuf(1,2),ibuf(1,3)
            write(*,*) "nzlma(min,max) =",ibuf(2,2),ibuf(2,3)
            write(*,*) "time(nloc)=",ct1-ct0,et1-et0
         end if

!ji-tmp
         goto 1000
         do a=1,MI
            ik=Kion(a)
            do iorb=1,norb(ik)
               L=lo(iorb,ik)
               do m=1,2*L+1
                  lma=icheck_tmp3(a,iorb,m)
                  if(lma/=0)then
                     write(*,'(1x,8i8)') a,L,m,lma,mtmp3(lma),myrank
                  end if
               end do
            end do
         end do
         x=0.d0
         y=1.d10
         do i=n1,n2
            r=H*H*(LL(1,i)**2+LL(2,i)**2+LL(3,i)**2)
            x=max(x,r)
            y=min(y,r)
         end do
         rewind 100+myrank
         do i=n1,n2
            write(100+myrank,*) LL(1,i),LL(2,i),LL(3,i)
            write(100+myrank,*)
         end do
         write(*,*) sqrt(x),sqrt(y),myrank
         goto 900
 1000    continue
!ji-tmp

!- deallocate --------------------------------------------------------
         mem=mem-bsintg*size(PPP)
         mem=mem-bsintg*size(irad)
         mem=mem-bdreal*size(uVtmp)
         mem=mem-bsintg*np_grid
         mem=mem-bdreal*size(rrr)
         mem=mem-bdreal*size(sss)
         deallocate( PPP )
         deallocate( irad )
         deallocate( uVtmp )
         deallocate( icheck_tmp1 )
         deallocate( sss )
         deallocate( rrr )
!---------------------------------------------------------------------

         call mpi_allreduce(nzlma,nzlma_0,1,mpi_integer,mpi_max,comm_grid,ierr)

         nzlma_0 = min(nzlma_0*2,Mlma)
         nrlma   = 0

!- allocate ----------------------------------------------------------
         allocate( lma_nsend_tmp(0:np_grid-1) )
         allocate( icheck_tmp1(0:np_grid-1)   )
         allocate( icheck_tmp2(0:np_grid-1)   )
         allocate( nl_rank_map_tmp(0:np_grid) )
         allocate( maps_tmp(nzlma_0,6) )
         allocate( sendmap0(nzlma_0,0:np_grid-1) )
         allocate( recvmap0(nzlma_0,0:np_grid-1) )
         n=max(np1,np2,np3)
         allocate( itmp(n,3) ) ; itmp=0
         mem=mem+bsintg*(np_grid*3+np_grid+1+nzlma_0*6)
         mem=mem+bsintg*(nzlma_0*np_grid*2+n*3)
         memax=max(mem,memax)
!---------------------------------------------------------------------

         maps_tmp(:,:)      = 0
         sendmap0(:,:)      = 0
         recvmap0(:,:)      = 0
         icheck_tmp1(:)     = 0
         icheck_tmp2(:)     = 0
         lma_nsend_tmp(:)   = 0
         nl_rank_map_tmp(:) =-1

         do a=1,MI
            ik=Kion(a)
         do iorb=1,norb(ik)
            L=lo(iorb,ik)
         do m=-L,L

!            if ( icheck_tmp3(a,iorb,m+L+1)/=0 ) then
!            write(*,*) icheck_tmp3(a,iorb,m+L+1),a,L,m,myrank
!            end if

            icheck_tmp1(:)=0
            call mpi_allgather(icheck_tmp3(a,iorb,m+L+1),1,mpi_integer &
&                             ,icheck_tmp1,1,mpi_integer,comm_grid,ierr)

            do n=0,np_grid-1
               if ( icheck_tmp1(n)/=0 ) then
                  icheck_tmp2(n)=icheck_tmp2(n)+1
               end if
            end do

            if ( icheck_tmp1(myrank_g)/=0 ) then
               maps_tmp(icheck_tmp2(myrank_g),1)=icheck_tmp1(myrank_g)
               maps_tmp(icheck_tmp2(myrank_g),2)=inorm(iorb,ik)
               maps_tmp(icheck_tmp2(myrank_g),3)=a
               maps_tmp(icheck_tmp2(myrank_g),4)=L
               maps_tmp(icheck_tmp2(myrank_g),5)=m
               maps_tmp(icheck_tmp2(myrank_g),6)=iorb
               do n=0,np_grid-1
                  if ( n==myrank_g .or. icheck_tmp1(n)==0 ) cycle
                  lma_nsend_tmp(n)=lma_nsend_tmp(n)+1
                  if ( lma_nsend_tmp(n)>nzlma_0 ) then
                     write(*,*) "nzlma_0 is small!"
                     write(*,*) nzlma_0,lma_nsend_tmp(n)
                     goto 900
                  end if
                  sendmap0(lma_nsend_tmp(n),n)=icheck_tmp2(myrank_g)
                  recvmap0(lma_nsend_tmp(n),n)=icheck_tmp2(n)
                  if ( any(nl_rank_map_tmp(0:nrlma)==n) ) cycle
                  nrlma=nrlma+1
                  nl_rank_map_tmp(nrlma)=n
               end do
            end if

         end do ! m
         end do ! L
         end do ! a

         nzlma = icheck_tmp2(myrank_g)

         MAXMJJ=0
         do lma=1,nzlma
            l=maps_tmp(lma,1) ; if (l==0) cycle
            MAXMJJ=max(MAXMJJ,mtmp3(l))
         end do

         MMJJ=MAXMJJ

         ibuf(1,1)=MAXMJJ
         call mpi_allreduce(ibuf(1,1),ibuf(1,2),1,mpi_integer,mpi_min,comm_grid,ierr)
         call mpi_allreduce(ibuf(1,1),ibuf(1,3),1,mpi_integer,mpi_max,comm_grid,ierr)
         if (DISP_SWITCH) then
            write(*,*) "MAXMJJ(min,max)=",ibuf(1,2:3)
         end if

!- allocate -----------------
         call gv_alloc("uV")
!----------------------------

         do lma=1,nzlma
            l=maps_tmp(lma,1) ; if (l==0) cycle
            MJJ(lma)     = mtmp3(l)
            iuV(lma)     = maps_tmp(lma,2)
            amap(lma)    = maps_tmp(lma,3)
            lmap(lma)    = maps_tmp(lma,4)
            mmap(lma)    = maps_tmp(lma,5)
            iorbmap(lma) = maps_tmp(lma,6)
         end do

         do i=1,nrlma
            nl_rank_map(i)=nl_rank_map_tmp(i)
         end do

!- allocate ------------------------------------------
         call gv_alloc("uVk_p1p2")
!-----------------------------------------------------

         if ( allocated(JJP3) ) deallocate(JJP3)
         allocate( JJP3(3,MAXMJJ,nzlma) ) ; JJP3=0

         do lma=1,nzlma
            l=maps_tmp(lma,1) ; if (l==0) cycle
            do j=1,MJJ(lma)
               uVk(j,lma,MBZ_0) = rtmp3(j,l)
               JJP3(1,j,lma)    = jtmp3(1,j,l)
               JJP3(2,j,lma)    = jtmp3(2,j,l)
               JJP3(3,j,lma)    = jtmp3(3,j,l)
            end do ! j
         end do ! lma

!- deallocate -----------------------
         mem=mem-bsintg*size(mtmp3)
         mem=mem-bsintg*size(jtmp3)
         mem=mem-bdreal*size(rtmp3)
         deallocate( mtmp3 )
         deallocate( jtmp3 )
         deallocate( rtmp3 )
!------------------------------------
!- deallocate -------------------------------
         mem=mem-bsintg*size(itmp)
         mem=mem-bsintg*size(nl_rank_map_tmp)
         mem=mem-bsintg*size(icheck_tmp2)
         mem=mem-bsintg*size(icheck_tmp1)
         mem=mem-bsintg*size(maps_tmp)
         mem=mem-bsintg*size(icheck_tmp3)
         deallocate( itmp )
         deallocate( nl_rank_map_tmp )
         deallocate( icheck_tmp2 )
         deallocate( icheck_tmp1 )
         deallocate( maps_tmp )
         deallocate( icheck_tmp3 )
!--------------------------------------------

         allocate( LLL0(-Mx:Mx,-My:My,-Mz:Mz) )
         LLL0=0 ; mem=mem+bsintg*size(LLL0) ; memax=max(mem,memax)
         do i=n1,n2
            LLL0(LL(1,i),LL(2,i),LL(3,i))=i
         end do
         allocate( uVk_tmp(MAXMJJ) )
         mem=mem+bdreal*size(uVk_tmp) ; memax=max(mem,memax)
         allocate( itmp1(n1:n2) )
         mem=mem+bsintg*size(itmp1) ; memax=max(mem,memax)
         allocate( itmp2(3,MAXMJJ) )
         mem=mem+bsintg*size(itmp2) ; memax=max(mem,memax)
         do lma=1,nzlma
            itmp1(:)=0
            itmp2(:,:)=0
            do j=1,MJJ(lma)
               i=LLL0(JJP3(1,j,lma),JJP3(2,j,lma),JJP3(3,j,lma))
               if ( i/=0 ) itmp1(i)=j
            end do ! j
            m=0
            uVk_tmp(:)=0.d0
            do i=n1,n2
               j=itmp1(i)
               if ( j/=0 ) then
                  m=m+1 ; if( m>MJJ(lma) ) goto 900
                  uVk_tmp(m)=uVk(j,lma,MBZ_0)
                  itmp2(1:3,m)=JJP3(1:3,j,lma)
               end if
            end do
            MJJ0(lma)=m
            do j=1,MJJ(lma)
               i=LLL0(JJP3(1,j,lma),JJP3(2,j,lma),JJP3(3,j,lma))
               if ( i==0 ) then
                  m=m+1
                  uVk_tmp(m)=uVk(j,lma,MBZ_0)
                  itmp2(1:3,m)=JJP3(1:3,j,lma)
               end if
            end do
            if ( m/=MJJ(lma) ) goto 900
            do j=1,MJJ(lma)
               uVk(j,lma,MBZ_0)=uVk_tmp(j)
               JJP3(1:3,j,lma)=itmp2(1:3,j)
               JJP(j,lma)=LLL(JJP3(1,j,lma),JJP3(2,j,lma),JJP3(3,j,lma))
            end do
         end do ! lma

         mem=mem-bsintg*size(itmp2) ; deallocate(itmp2)
         mem=mem-bsintg*size(itmp1) ; deallocate(itmp1)
         mem=mem-bdreal*size(uVk_tmp) ; deallocate(uVk_tmp)
         mem=mem-bsintg*size(LLL0) ; deallocate(LLL0)

         nl_max_send = maxval(lma_nsend_tmp)

         allocate( itmp1(0:np_grid-1) ) ; itmp1=0
         itmp1(myrank_g)=nl_max_send
         call mpi_allgather(nl_max_send,1,mpi_integer,itmp1,1,mpi_integer,comm_grid,ierr)
         if (DISP_SWITCH) then
            write(*,*) "nl_max_send(max)=",maxval(itmp1),maxloc(itmp1)
            write(*,*) "nl_max_send(min)=",minval(itmp1),minloc(itmp1)
         end if
         deallocate( itmp1 )

!- allocate -------------------------------
         call gv_alloc("lma_nsend")
!------------------------------------------

         do n=0,np_grid-1
            sendmap(1:nl_max_send,n) = sendmap0(1:nl_max_send,n)
            lma_nsend(n) = lma_nsend_tmp(n)
         end do

!- allocate -------------------------------------------------
         allocate( ireq(2*np_grid) )
         mem=mem+bsintg*2*np_grid ; memax=max(mem,memax)
!------------------------------------------------------------

         nreq=0
         do n=0,np_grid-1
            if ( lma_nsend(n)<=0 .or. n==myrank_g ) cycle
            nreq=nreq+1
            call mpi_isend(recvmap0(1,n),lma_nsend(n),mpi_integer,n,1,comm_grid,ireq(nreq),ierr)
            nreq=nreq+1
            call mpi_irecv(recvmap(1,n) ,lma_nsend(n),mpi_integer,n,1,comm_grid,ireq(nreq),ierr)
         end do
         call mpi_waitall(nreq,ireq,istatus,ierr)

!- deallocate -----------------------------------------------------
         deallocate( ireq ) ; mem=mem-bsintg*2*np_grid
!------------------------------------------------------------------

!- deallocate -----------------------------------------------------
         deallocate( recvmap0,sendmap0,lma_nsend_tmp )
         mem=mem-bsintg*(nzlma_0*np_grid*2 + np_grid)
!------------------------------------------------------------------

!- allocate ---------------------
      call gv_alloc("sbufnl")
!--------------------------------

      Max_mem_inst = max( Max_mem+memax, Max_mem_inst )

      call bwatch(ctime1,etime1)
      if (DISP_SWITCH) then
         write(*,*) "TIME(PREP_PS_MOL_SYM)=",ctime1-ctime0,etime1-etime0
         write(*,*) "MEM(MB)=",mem,memax*B2MB
      end if

      return

 900  call stop_program1("prep_ps_mol_sym",1)

      END SUBROUTINE prep_ps_mol_sym
