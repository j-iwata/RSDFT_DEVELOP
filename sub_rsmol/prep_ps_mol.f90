!--------1---------2---------3---------4---------5---------6---------7--

      SUBROUTINE prep_ps_mol(iswitch_pp)
      use global_variables
      implicit none

      integer,intent(IN) :: iswitch_pp
      logical,save :: iswitch_prep_oh=.true.
      integer :: Mx,My,Mz
      integer,allocatable :: LLL0(:,:,:),PPP0(:,:,:),LL0(:,:)
      integer,allocatable :: iloc(:),Jtmp(:,:,:)
      integer,allocatable :: ichk0(:,:),ichk1(:,:)
      integer,allocatable :: itmp(:,:),itmp1(:),itmp2(:),itmp3(:,:)
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
      real(8),allocatable :: uV_tmp(:,:)
      real(8) :: x,y,z,r,v,err,rps2,rps2_0,r2,c1,c2,H
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
      integer :: MMJJ_0,lma0,nzlma_0
      integer,allocatable :: jtmp3(:,:,:),mtmp3(:)
      integer,allocatable :: nl_rank_map_tmp(:),maps_tmp(:,:)
      integer,allocatable :: lma_nsend_tmp(:)
      real(8) :: rc2,Rx,Ry,Rz,p1,p2,p3,p4,tmp0
      real(8),allocatable :: rwork(:,:,:,:),rwork1(:,:,:)
      real(8),allocatable :: rtmp3(:,:),sss(:,:,:,:),rrr(:,:,:,:)

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


      if ( isymmetry==1 ) then
         call prep_ps_mol_sym(iswitch_pp)
         return
      end if

      call bwatch(ctime0,etime0)

      if (DISP_SWITCH) then
         write(*,'(a60," prep_ps")') repeat("-",60)
      end if

      n1    = idisp(myrank)+1
      n2    = idisp(myrank)+ircnt(myrank)
      ML0   = ircnt(myrank)
      mem   = 0.d0
      memax = 0.d0

      Mx = ML1
      My = ML2
      Mz = ML3

      H = H1

      if ( iswitch_prep_oh ) then

         iswitch_prep_oh=.false.
!
! --- Parameters for OH method ---
!
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

!
! --- Coeficients of Lagrange Interporation ---
!

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

!- allocate ---------------------
      call gv_alloc("Vion")
!--------------------------------

!
! --- Grid Map ---
!

!- allocate -------------------------------------------------------------
      allocate( LLp(3,0:np_grid-1) ) ; LLp=-1
      allocate( LLLp(np1,np2,np3) ) ; LLLp=-1
      mem=mem+bsintg*3*np_grid+bsintg*np_grid
      allocate( LL(3,n1:n2) ) ; LL=0
      mem=mem+bsintg*3*(n2-n1+1)
      allocate( PPP(-Mx:Mx,-My:My,-Mz:Mz)  ) ; PPP=0
      allocate( LLL(-Mx:Mx,-My:My,-Mz:Mz)  ) ; LLL=0
      allocate( PPP0(-Mx:Mx,-My:My,-Mz:Mz) ) ; PPP0=0
      mem=mem+bsintg*(2*Mx+1)*(2*My+1)*(2*Mz+1)*3 ; memax=max(mem,memax)
!------------------------------------------------------------------------

      n=-1
      do i3=1,np3
      do i2=1,np2
      do i1=1,np1
         n=n+1
         LLp(1,n)=i1
         LLp(2,n)=i2
         LLp(3,n)=i3
         LLLp(i1,i2,i3)=n
      end do
      end do
      end do

      Rc2=Rsize**2
      i=n1-1
      do i3=a3b,b3b
      do i2=a2b,b2b
      do i1=a1b,b1b
         select case(SYStype)
         case(1)
            r2=H*H*(i1*i1+i2*i2+i3*i3)
            if ( r2<Rc2+eps ) then
               i=i+1
               LL(1,i)=i1
               LL(2,i)=i2
               LL(3,i)=i3
            end if
         case(2)
            r2=H*H*(i1*i1+i2*i2)
            z=abs(i3*H)
            if ( r2<Rc2+eps .and. z<Zsize+eps ) then
               i=i+1
               LL(1,i)=i1
               LL(2,i)=i2
               LL(3,i)=i3
            end if
         end select
      end do
      end do
      end do

      m0=size(PPP0)
      PPP0(:,:,:)=0
      do i=n1,n2
         i1=LL(1,i) ; i2=LL(2,i) ; i3=LL(3,i)
         PPP0(i1,i2,i3)=myrank_g-MPI_PROC_NULL
      end do
      call mpi_allreduce(PPP0,PPP,m0,mpi_integer,mpi_sum,comm_grid,ierr)
      PPP(:,:,:)=PPP(:,:,:)+MPI_PROC_NULL
      PPP0(:,:,:)=0
      do i=n1,n2
         i1=LL(1,i) ; i2=LL(2,i) ; i3=LL(3,i)
         PPP0(i1,i2,i3)=i
      end do
      call mpi_allreduce(PPP0,LLL,m0,mpi_integer,mpi_sum,comm_grid,ierr)

      if (DISP_SWITCH) then
         write(*,*) "ML,MK,ML+MK",ML,MK,ML+MK
         write(*,*) "count(LLL)",count(LLL/=0),size(LLL)
         write(*,*) "count(PPP)",count(PPP>=0),size(PPP)
      end if

! - deallocate ---------------------------------------
      mem=mem-bsintg*size(PPP0) ; deallocate( PPP0 )
!-----------------------------------------------------


! ----------------------------------------------------------------------
! --- local potential --------------------------------------------------
! ----------------------------------------------------------------------

      if (DISP_SWITCH) write(*,*) "---local potential---"

      Vion(:) = 0.d0

      r=maxval( Rcloc(1:MKI) )

      ll_0 = min(0,-Nintp+1)

!
! --- Construct minimal cubic box including PP spheres --- 
!

      call Make_MinimalBox_Mol(r,mm1,mm2,mm3,mm4)

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

!      mm1 = mm1_lc + Nintp + 1
!      mm2 = mm2_lc + Nintp + 1
!      mm3 = mm3_lc + Nintp + 1
      mm1 = Mx + Nintp + 1
      mm2 = My + Nintp + 1
      mm3 = Mz + Nintp + 1

!- allocate -------------------------------------------------------
!      allocate( rwork(-Mx:Mx,-My:My,-Mz:Mz,2) ) ; rwork =0.d0
!      allocate( rwork(-Mx-Md:Mx+Md,-My-Md:My+Md,-Mz-Md:Mz+Md,2) ) ; rwork =0.d0
      allocate( rwork(-mm1:mm1,-mm2:mm2,-mm3:mm3,2) ) ; rwork =0.d0
      mem=mem+bdreal*size(rwork)
      allocate( irad(0:3000,MKI) ) ; irad=0
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
         rc2 = Rcloc(ik)**2

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

            if ( PPP(id1,id2,id3)==myrank_g ) then

               do j3=0,Ndense-1
                  z=id3*H3+j3*H3d-Rz
               do j2=0,Ndense-1
                  y=id2*H2+j2*H2d-Ry
               do j1=0,Ndense-1
                  x=id1*H1+j1*H1d-Rx

                  r2=x*x+y*y+z*z

                  if ( r2>rc2+1.d-10 ) cycle

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

            end if ! PPP

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

!- deallocate -----------------------------------------
      mem=mem-bdreal*size(rwork) ; deallocate( rwork )
!------------------------------------------------------

! ----------------------------------------------------------------------
! --- Nonlocal ---------------------------------------------------------
! ----------------------------------------------------------------------

      if (DISP_SWITCH) write(*,*) "---nonlocal potential---"

      r=maxval( Rps(:,:) )

      ll_0 = min(0,-Nintp+1)

!
! --- Construct minimal cubic box including PP spheres --- 
!

      call Make_MinimalBox_Mol(r,mm1,mm2,mm3,mm4)

      mm1_ps = maxval( abs(mcube_grid_ion(:,1)) )
      mm2_ps = maxval( abs(mcube_grid_ion(:,2)) )
      mm3_ps = maxval( abs(mcube_grid_ion(:,3)) )

      if (DISP_SWITCH) then
         write(*,*) "mm1_ps=",mm1_ps
         write(*,*) "mm2_ps=",mm2_ps
         write(*,*) "mm3_ps=",mm3_ps
         write(*,*) "M_grid_ion=",M_grid_ion
      end if

!
!  --- array size guess ( nzlma_0, MMJJ_0 ) ---
!

      if ( Mlma<np_grid ) then
         nzlma_0 = Mlma
      else
         nzlma_0 = min(Mlma*729/np_grid,Mlma)
         nzlma_0 = max(1,nzlma_0)
      end if
      if ( DISP_SWITCH ) then
         write(*,*) "nzlma(guess),Mlma =",nzlma_0,Mlma
      end if

      mm1 = mm1_ps + Nintp + 1
      mm2 = mm2_ps + Nintp + 1
      mm3 = mm3_ps + Nintp + 1

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
         allocate( icheck_tmp3(MI,n,2*L+1)  ) ; icheck_tmp3=0
         allocate( icheck_tmp1(0:np_grid-1) ) ; icheck_tmp1=0
!         allocate( JJ_tmp_0(0:9,MMJJ_0) ) ; JJ_tmp_0=0
         allocate( rtmp3(MMJJ_0,nzlma_0)    ) ; rtmp3=0.d0
         allocate( sss(-mm1:mm1,-mm2:mm2,-mm3:mm3,2*L+1) ) ; sss=0.d0
         allocate( rrr(-mm1:mm1,-mm2:mm2,-mm3:mm3,2*L+1) ) ; rrr=0.d0
         allocate( uVtmp(-L:L) )              ; uVtmp=0.d0
         allocate( jtmp3(6,MMJJ_0,nzlma_0)  ) ; jtmp3=-100000
         allocate( mtmp3(nzlma_0)           ) ; mtmp3=0
         mem=mem+bsintg*MI*n*(2*L+1)
         mem=mem+bsintg*np_grid
!         mem=mem+bsintg*MMJJ_0*10
         mem=mem+bdreal*MMJJ_0*nzlma_0
         mem=mem+bdreal*m*(2*L+1)*2
         mem=mem+bdreal*(2*L+1)
         mem=mem+bsintg*6*MMJJ_0*nzlma_0
         mem=mem+bsintg*nzlma_0
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

      maxerr             = 0
      lma                = 0
      icheck_tmp3(:,:,:) = 0
      MMJJ               = 0
      nzlma              = 0
      lma0               = 0
      m0                 = size(rrr)

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
            sss(:,:,:,:) = zero
            rrr(:,:,:,:) = zero
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

               if ( PPP(id1,id2,id3)==myrank_g ) then

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

                  end if ! PPP

               end do ! i ( 1 - M_grid_ion )

               if ( icheck_tmp1(myrank_g)>0 ) then

                  call mpi_allreduce_nlpp( sss,rrr,m0,icheck_tmp1,ierr,"R" )

                  j1=0
                  j2=0
                  do i3=-mm3,mm3
                     id3=i3+ic3
                  do i2=-mm2,mm2
                     id2=i2+ic2
                  do i1=-mm1,mm1
                     id1=i1+ic1
                     if ( PPP(id1,id2,id3)==myrank_g ) then
                        j1    = j1+1
                        lma   = lma0
                        do m=1,2*L+1
                           if ( abs(rrr(i1,i2,i3,m))<1.d-10 ) cycle
                           j2=j2+1
                           do mm=1,2*L+1
                              lma=lma+1
                              rtmp3(j2,lma)=rrr(i1,i2,i3,mm)
                              jtmp3(1,j2,lma)=id1
                              jtmp3(2,j2,lma)=id2
                              jtmp3(3,j2,lma)=id3
                              jtmp3(4,j2,lma)=i1
                              jtmp3(5,j2,lma)=i2
                              jtmp3(6,j2,lma)=i3
                              mtmp3(lma)=j2
                           end do
                           exit
                        end do
                     end if
                  end do
                  end do
                  end do

                  MMJJ  = max( MMJJ, j2 )
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

               end if ! icheck_tmp1(myrank_g)>0

            end do ! iorb
         end do ! a

         ibuf(1,1)=MMJJ
         ibuf(2,1)=nzlma
         call mpi_allreduce(ibuf(1,1),ibuf(1,2),2,mpi_integer,mpi_min,comm_grid,ierr)
         call mpi_allreduce(ibuf(1,1),ibuf(1,3),2,mpi_integer,mpi_max,comm_grid,ierr)

         call bwatch(ct1,et1)

         if (DISP_SWITCH) then
            write(*,*) "lma=",lma
            write(*,*) "myrank_g,nzlma,MMJJ",myrank_g,nzlma,MMJJ
            write(*,*) "MMJJ(min,max)  =",ibuf(1,2),ibuf(1,3)
            write(*,*) "nzlma(min,max) =",ibuf(2,2),ibuf(2,3)
            write(*,*) "time(nloc)=",ct1-ct0,et1-et0
         end if

!- deallocate --------------------------------------------------------
         mem=mem-bsintg*size(irad)
         mem=mem-bdreal*size(uVtmp)
         mem=mem-bsintg*np_grid
         mem=mem-bdreal*size(rrr)
         mem=mem-bdreal*size(sss)
         deallocate( irad )
         deallocate( uVtmp )
         deallocate( icheck_tmp1 )
         deallocate( rrr )
         deallocate( sss )
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

            icheck_tmp1(:)=0
            call mpi_allgather(icheck_tmp3(a,iorb,m+L+1),1,mpi_integer,icheck_tmp1,1,mpi_integer,comm_grid,ierr)

            itmp(:,:)=0
            do n=0,np_grid-1
               if ( icheck_tmp1(n)==0 ) cycle
               itmp(LLp(1,n),1)=1 !LLp(1,n)
               itmp(LLp(2,n),2)=1 !LLp(2,n)
               itmp(LLp(3,n),3)=1 !LLp(3,n)
            end do

            k1=count( itmp(:,1)>0 )
            k2=count( itmp(:,2)>0 )
            k3=count( itmp(:,3)>0 )

            ic1=0
            id1=np1
            do i=1,np1
               if ( ic1==0 .and. itmp(i,1)/=0 ) then
                  ic1=i
               else if ( ic1/=0 .and. itmp(i,1)==0 ) then
                  id1=i-1
                  exit
               end if
            end do
            if ( id1-ic1+1/=k1 ) then
               i1=0
               j1=np1
               do i=id1+1,np1
                  if ( i1==0 .and. itmp(i,1)/=0 ) then
                     i1=i
                  else if ( i1/=0 .and. itmp(i,1)==0 ) then
                     j1=i-1
                     exit
                  end if
               end do
               if ( id1-ic1+1+j1-i1+1/=k1 .or. j1/=np1 .or. ic1/=1 ) then
                  write(*,'(1x,"a",5i4)') ic1,id1,i1,j1,k1
                  goto 900
               end if
               i1=i1-np1
               j1=j1-np1
               ic1=i1
            end if
            ic2=0
            id2=np2
            do i=1,np2
               if ( ic2==0 .and. itmp(i,2)/=0 ) then
                  ic2=i
               else if ( ic2/=0 .and. itmp(i,2)==0 ) then
                  id2=i-1
                  exit
               end if
            end do
            if ( id2-ic2+1/=k2 ) then
               i2=0
               j2=np2
               do i=id2+1,np2
                  if ( i2==0 .and. itmp(i,2)/=0 ) then
                     i2=i
                  else if ( i2/=0 .and. itmp(i,2)==0 ) then
                     j2=i-1
                     exit
                  end if
               end do
               if ( id2-ic2+1+j2-i2+1/=k2 .or. j2/=np2 .or. ic2/=1 ) then
                  write(*,*) ic2,id2,i2,j2,k2
                  goto 900
               end if
               i2=i2-np2
               j2=j2-np2
               ic2=i2
            end if
            ic3=0
            id3=np3
            do i=1,np3
               if ( ic3==0 .and. itmp(i,3)/=0 ) then
                  ic3=i
               else if ( ic3/=0 .and. itmp(i,3)==0 ) then
                  id3=i-1
                  exit
               end if
            end do
            if ( id3-ic3+1/=k3 ) then
               i3=0
               j3=np3
               do i=id3+1,np3
                  if ( i3==0 .and. itmp(i,3)/=0 ) then
                     i3=i
                  else if ( i3/=0 .and. itmp(i,3)==0 ) then
                     j3=i-1
                     exit
                  end if
               end do
               if ( id3-ic3+1+j3-i3+1/=k3 .or. j3/=np3 .or. ic3/=1 ) then
                  write(*,*) ic3,id3,i3,j3,k3
                  goto 900
               end if
               i3=i3-np3
               j3=j3-np3
               ic3=i3
            end if

            if ( id1-ic1+1/=k1 .or. id2-ic2+1/=k2 .or. id3-ic3+1/=k3 ) then
               write(*,*) ic1,id1,k1,ic2,id2,k2,ic3,id3,k3
               goto 900
            end if

            do j3=ic3,id3
            do j2=ic2,id2
            do j1=ic1,id1
               k=LLLp(j1,j2,j3)
               if ( icheck_tmp1(k)==0 ) icheck_tmp1(k)=-1
            end do
            end do
            end do

            do n=0,np_grid-1
               if ( icheck_tmp1(n)/=0 ) then
                  icheck_tmp2(n)=icheck_tmp2(n)+1
               end if
            end do

            if ( icheck_tmp1(myrank_g)/=0 ) then
               if ( icheck_tmp1(myrank_g)>0 ) then
                  maps_tmp(icheck_tmp2(myrank_g),1)=icheck_tmp1(myrank_g)
               end if
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

         end do
         end do
         end do

         nzlma = icheck_tmp2(myrank_g)


!- allocate --------------------------------------------------------------
         allocate( icheck_tmp4(a1b:b1b,a2b:b2b,a3b:b3b) ) ; icheck_tmp4=0
         mem=mem+bsintg*size(icheck_tmp4) ; memax=max(mem,memax)
!-------------------------------------------------------------------------

         MAXMJJ=0
         do lma=1,nzlma
            l=maps_tmp(lma,1) ; if (l==0) cycle
            j=0
            icheck_tmp4(:,:,:)=0
            do i=1,mtmp3(l)
               i1=jtmp3(1,i,l)
               i2=jtmp3(2,i,l)
               i3=jtmp3(3,i,l)
               if ( icheck_tmp4(i1,i2,i3)==0 ) then
                  j=j+1
                  icheck_tmp4(i1,i2,i3)=j
               end if
            end do
            MAXMJJ=max(MAXMJJ,j)
         end do

         MMJJ=MAXMJJ

!- allocate ---------------------
         call gv_alloc("uV")
!--------------------------------

         do lma=1,nzlma
            l=maps_tmp(lma,1) ; if (l==0) cycle
            j=0
            icheck_tmp4(:,:,:)=0
            do i=1,mtmp3(l)
               i1=jtmp3(1,i,l)
               i2=jtmp3(2,i,l)
               i3=jtmp3(3,i,l)
               if ( icheck_tmp4(i1,i2,i3)==0 ) then
                  j=j+1
                  icheck_tmp4(i1,i2,i3)=j
               end if
            end do
            MJJ(lma)=j
         end do

         ibuf(1,1)=MAXMJJ
         call mpi_allreduce(ibuf(1,1),ibuf(1,2),1,mpi_integer,mpi_min,comm_grid,ierr)
         call mpi_allreduce(ibuf(1,1),ibuf(1,3),1,mpi_integer,mpi_max,comm_grid,ierr)
         if (DISP_SWITCH) then
            write(*,*) "MAXMJJ(min,max)=",ibuf(1,2:3)
         end if

         do lma=1,nzlma
            l=maps_tmp(lma,1) ; if ( l==0 ) cycle
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

         do lma=1,nzlma
            l=maps_tmp(lma,1) ; if (l==0) cycle
            j=0
            icheck_tmp4(:,:,:)=0
            do i=1,mtmp3(l)
               i1=jtmp3(1,i,l)
               i2=jtmp3(2,i,l)
               i3=jtmp3(3,i,l)
               j1=icheck_tmp4(i1,i2,i3)
               if ( j1==0 ) then
                  j=j+1
                  icheck_tmp4(i1,i2,i3)=j
                  uVk(j,lma,MBZ_0)  = rtmp3(i,l)
                  JJP(j,lma) = LLL(i1,i2,i3)
               else 
                  uVk(j1,lma,MBZ_0) = uVk(j1,lma,MBZ_0) + rtmp3(i,l)
               end if
            end do
         end do

! - deallocate -------------------------------------------------------
         mem=mem-bsintg*size(icheck_tmp4) ; deallocate( icheck_tmp4 )
!---------------------------------------------------------------------

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
!- allocate -------------------------------------------------------
      allocate( itmp(3,nrlma) ) ; itmp=0
      allocate( itmp1(nrlma), work(nrlma) )
      allocate( itmp2(nrlma),itmp3(3,nrlma) )
      mem=mem+bsintg*8*nrlma+bdreal*nrlma ; memax=max(mem,memax)
!------------------------------------------------------------------

      do irlma=1,nrlma
         n=nl_rank_map(irlma)
         itmp(1,irlma)=LLp(1,n)-LLp(1,myrank_g)
         itmp(2,irlma)=LLp(2,n)-LLp(2,myrank_g)
         itmp(3,irlma)=LLp(3,n)-LLp(3,myrank_g)
      end do

      nrlma_xyz(1:6)=0

      m=0
      n=0
      do i=1,nrlma
         if( itmp(2,i)==0 .and. itmp(3,i)==0 .and. itmp(1,i)>0 )then
            n=n+1
            work(n)=itmp(1,i)
            itmp2(n)=i
         end if
      end do
      if ( n>0 ) then
         call indexx(n,work,itmp1)
         do i=1,n
            j=itmp2( itmp1(i) )
            itmp3(:,m+i)=itmp(:,j)
         end do
      end if
      m=m+n
      nrlma_xyz(1)=nrlma_xyz(1)+n
      n=0
      do i=1,nrlma
         if( itmp(2,i)==0 .and. itmp(3,i)==0 .and. itmp(1,i)<0 )then
            n=n+1
            work(n)=itmp(1,i)
            itmp2(n)=i
         end if
      end do
      if ( n>0 ) then
         call indexx(n,work,itmp1)
         do i=1,n
            j=itmp2(itmp1(i))
            itmp3(:,m+n-i+1)=itmp(:,j)
         end do
      end if
      m=m+n
      nrlma_xyz(2)=nrlma_xyz(2)+n

      n=0
      do i=1,nrlma
         if( itmp(1,i)==0 .and. itmp(3,i)==0 .and. itmp(2,i)>0 )then
            n=n+1
            work(n)=itmp(2,i)
            itmp2(n)=i
         end if
      end do
      if ( n>0 ) then
         call indexx(n,work,itmp1)
         do i=1,n
            j=itmp2( itmp1(i) )
            itmp3(:,m+i)=itmp(:,j)
         end do
      end if
      m=m+n
      nrlma_xyz(3)=nrlma_xyz(3)+n
      n=0
      do i=1,nrlma
         if( itmp(1,i)==0 .and. itmp(3,i)==0 .and. itmp(2,i)<0 )then
            n=n+1
            work(n)=itmp(2,i)
            itmp2(n)=i
         end if
      end do
      if ( n>0 ) then
         call indexx(n,work,itmp1)
         do i=1,n
            j=itmp2(itmp1(i))
            itmp3(:,m+n-i+1)=itmp(:,j)
         end do
      end if
      m=m+n
      nrlma_xyz(4)=nrlma_xyz(4)+n

      n=0
      do i=1,nrlma
         if( itmp(1,i)==0 .and. itmp(2,i)==0 .and. itmp(3,i)>0 )then
            n=n+1
            work(n)=itmp(3,i)
            itmp2(n)=i
         end if
      end do
      if ( n>0 ) then
         call indexx(n,work,itmp1)
         do i=1,n
            j=itmp2( itmp1(i) )
            itmp3(:,m+i)=itmp(:,j)
         end do
      end if
      m=m+n
      nrlma_xyz(5)=nrlma_xyz(5)+n
      n=0
      do i=1,nrlma
         if( itmp(1,i)==0 .and. itmp(2,i)==0 .and. itmp(3,i)<0 )then
            n=n+1
            work(n)=itmp(3,i)
            itmp2(n)=i
         end if
      end do
      if ( n>0 ) then
         call indexx(n,work,itmp1)
         do i=1,n
            j=itmp2(itmp1(i))
            itmp3(:,m+n-i+1)=itmp(:,j)
         end do
      end if
      m=m+n
      nrlma_xyz(6)=nrlma_xyz(6)+n

!- allocate ----------------------
      call gv_alloc("num_2_rank")
!---------------------------------

      m=0
      do i=1,6
         do j=1,nrlma_xyz(i)
            m=m+1
            i1=itmp3(1,m)+LLp(1,myrank_g)
            i2=itmp3(2,m)+LLp(2,myrank_g)
            i3=itmp3(3,m)+LLp(3,myrank_g)
            k=LLLp(i1,i2,i3)
            num_2_rank(j,i)=k
         end do
      end do

!- deallocate ---------------------------------------------
      mem=mem-bsintg*(nrlma*8)-bdreal*nrlma
      deallocate( itmp,itmp1,itmp2,itmp3,work )
!----------------------------------------------------------

      do i=1,5,2
         n=max( nrlma_xyz(i),nrlma_xyz(i+1) )
         nrlma_xyz(i)=n
         nrlma_xyz(i+1)=n
      end do

!- allocate ---------------------
      call gv_alloc("sbufnl")
!--------------------------------



!- deallocate ---------------------------------------
      mem=mem-bsintg*size(PPP)  ; deallocate( PPP )
      mem=mem-bsintg*size(LLL)  ; deallocate( LLL )
      mem=mem-bsintg*size(LL)   ; deallocate( LL )
      mem=mem-bsintg*size(LLLp) ; deallocate( LLLp )
      mem=mem-bsintg*size(LLp)  ; deallocate( LLp )
!----------------------------------------------------

      Max_mem_inst = max( Max_mem+memax, Max_mem_inst )

      call bwatch(ctime1,etime1)
      if (DISP_SWITCH) then
         write(*,*) "TIME(PREP_PS_MOL)=",ctime1-ctime0,etime1-etime0
         write(*,*) "MEM(MB)=",mem,memax*B2MB,Max_mem_inst*B2MB
      end if

      return

 900  call stop_program1("prep_ps_mol",1)

      END SUBROUTINE prep_ps_mol
