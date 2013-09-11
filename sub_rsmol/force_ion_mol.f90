!--------1---------2---------3---------4---------5---------6---------7--

      SUBROUTINE force_ion_mol
      use global_variables
      implicit none
      integer :: i1_0,i2_0,i3_0,iii1,iii2,iii3,ii1,ii2,ii3,s
      integer :: i1,i2,i3,m0,ML0,i1_anear,i2_anear,i3_anear
      integer :: i,j,k,ir,iorb,L,n,L1,L1z,NRc,l2,l3,irank,jrank
      integer :: ireq(1024),istatus(MPI_STATUS_SIZE,123)
      integer :: nreq,mrnk,mx,my,mz,mtmp,Mlma_np,id1,id2,id3
      integer :: a,ik,lm,m,lm0,lm1,lma,lma0,lma1,im,m1,m2
      integer :: nn1,nn2,ierr,M_irad,ir0,j1,j2,j3,ii_0
      integer,allocatable :: irad(:,:),idis(:),icnt(:),a_rank(:)
      integer,allocatable :: icheck_tmp3(:),JJ_tmp(:,:),PPP0(:,:,:)
      real(8),parameter :: ep=1.d-8
      real(8) :: Y1(0:2,-2:2,0:3,-3:3)
      real(8) :: Y2(0:2,-2:2,0:3,-3:3)
      real(8) :: Y3(0:2,-2:2,0:3,-3:3)
      real(8) :: err,err0,maxerr,rx,ry,rz
      real(8) :: a1,a2,a3,c1,c2,c3,d1,d2,d3,const,Gx,Gy,Gz
      real(8) :: x,y,z,r,rps2,r2,Gr,kr,pi2,k1,k2,k3,kx,ky,kz
      real(8) :: tmp,tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
      real(8) :: ct(9),et(9),ct0,et0,ct1,et1
      real(8) :: ctime0,ctime1,etime0,etime1,mem,memax
      real(8),allocatable :: force1(:,:),force2(:,:),force3(:,:)
      real(8),allocatable :: yy1(:),yy2(:),yy3(:),rr(:,:),work(:),work2(:,:)
      real(8),allocatable :: rf1(:,:,:,:,:),rf2(:,:,:,:,:)
      complex(8),parameter :: z0=(0.d0,0.d0),zi=(0.d0,1.d0)
      complex(8) :: ztmp,phase
      complex(8),allocatable :: zwork(:,:,:),zwork0(:),fftwork(:)
      logical :: DISP_SWITCH_LOC,flag_alloc

      INTERFACE
         FUNCTION Ylm(x,y,z,l,m)
         real(8) :: Ylm
         real(8),intent(IN) :: x,y,z
         integer,intent(IN) :: l,m
         END FUNCTION Ylm
      END INTERFACE

      call watch(ctime0,etime0)
      if (DISP_SWITCH) write(*,*) repeat("-",60)," d_force_ion_mol"

      nn1 = idisp(myrank)+1
      nn2 = idisp(myrank)+ircnt(myrank)
      ML0 = ircnt(myrank)
      pi2 = 2.d0*Pi

      ct(:) = 0.d0
      et(:) = 0.d0
      mem   = 0.d0
      memax = 0.d0

      allocate( force1(3,MI),force2(3,MI),force3(3,MI) )
      force1=0.d0
      force2=0.d0
      force3=0.d0
      mem=mem+bdreal*3*3*MI ; memax=max(mem,memax)

!
! --- Grid Map ---
!
      flag_alloc=.false.
      if ( .not.allocated(LL) ) then
         allocate( LL(3,nn1:nn2) ) ; LL=0 ; mem=mem+bsintg*size(LL) ; memax=max(mem,memax)
         flag_alloc=.true.
      end if
      call Make_GridMap_1(LL,nn1,nn2)
      allocate( LLL(-ML1:ML1,-ML2:ML2,-ML3:ML3) ) ; LLL=0
      allocate( PPP(-ML1:ML1,-ML2:ML2,-ML3:ML3) ) ; PPP=0
      allocate( PPP0(-ML1:ML1,-ML2:ML2,-ML3:ML3) ) ; PPP0=0
      mem=mem+bsintg*(size(PPP)+size(PPP0)+size(LLL)) ; memax=max(mem,memax)
      PPP0(:,:,:)=0
      do i=nn1,nn2
         i1=LL(1,i) ; i2=LL(2,i) ; i3=LL(3,i)
         PPP0(i1,i2,i3)=myrank_g-MPI_PROC_NULL
      end do
      m0=size(PPP0)
      call mpi_allreduce(PPP0,PPP,m0,mpi_integer,mpi_sum,comm_grid,ierr)
      PPP(:,:,:)=PPP(:,:,:)+MPI_PROC_NULL
      PPP0(:,:,:)=0
      do i=nn1,nn2
         i1=LL(1,i) ; i2=LL(2,i) ; i3=LL(3,i)
         PPP0(i1,i2,i3)=i
      end do
      call mpi_allreduce(PPP0,LLL,m0,mpi_integer,mpi_sum,comm_grid,ierr)

      if (DISP_SWITCH) then
         write(*,*) "ML,MK,ML+MK,ML",MK,ML+MK
         write(*,*) "count(LLL)",count(LLL/=0),size(LLL)
         write(*,*) "count(PPP)",count(PPP>=0),size(PPP)
      end if

! - deallocate ---------------------------------------
      mem=mem-bsintg*size(PPP0) ; deallocate( PPP0 )
!-----------------------------------------------------

!
! --- Local part ---
!

      call force_local_mol(force1)

!
! --- Nonlocal ---
!

      call watch(ct0,et0)

      force2(:,:)=0.d0

      select case(pselect)

      case default

      case(1,2)

      call calcY

      mx = mm1_ps+Nintp+1
      my = mm2_ps+Nintp+1
      mz = mm3_ps+Nintp+1

! - Memory Check -
      m = (2*mx+1)*(2*my+1)*(2*mz+1)
      L = 2*maxval(lo)+1
      mem=mem+bdreal*3*m*L
      mem=mem+bdreal*3*m*L
      mem=mem+bdmain*4*nzlma*(MB_1-MB_0+1)*N_MBZ*N_MSP
      mem=mem+bdmain*4*nzlma
      mem=mem+bsintg*np_grid
      mem=mem+bdreal*L*3
      mem=mem+bsintg*3001*MKI
      mem=mem+bsintg*m*4
      mem=mem+bdreal*m*3
      mem=mem+bsintg*MI
      memax=max(mem,memax)
      if(DISP_SWITCH)then
         write(*,*) "mem(MB)=",mem*B2MB,memax*B2MB
      end if

      allocate( rf1(3,-mx:mx,-my:my,-mz:mz,L) ) ; rf1=0.d0
      allocate( rf2(3,-mx:mx,-my:my,-mz:mz,L) ) ; rf2=0.d0
      allocate( wtmp5(0:3,nzlma,MB_0:MB_1,MBZ_0:MBZ_1,MSP_0:MSP_1) ) ; wtmp5=zero
      allocate( vtmp2(0:3,nzlma) ) ; vtmp2=zero
      allocate( icheck_tmp3(0:np_grid-1) ) ; icheck_tmp3=0
      allocate( yy1(L),yy2(L),yy3(L) )
      allocate( irad(0:3000,MKI) ) ; irad=0
      allocate( JJ_tmp(0:3,m) ) ; JJ_tmp=0
      allocate( rr(3,m) ) ; rr=0.d0
      allocate( a_rank(MI) ) ; a_rank=0

      mrnk   = id_class(myrank,0)
      maxerr = 0.d0
      m0     = size(rf1)
      ii_0   = min(0,-Nintp+1)
      c1     = 1.d0/ML1
      c2     = 1.d0/ML2
      c3     = 1.d0/ML3
      M_irad = 0

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

      do a=1,MI

         ik  = Kion(a)
         lm0 = 0

         i1_anear = nint(asi(1,a)*aaL(1)/H1)
         i2_anear = nint(asi(2,a)*aaL(2)/H2)
         i3_anear = nint(asi(3,a)*aaL(3)/H3)

!         l1=i1_anear/ML1 ; if ( i1_anear<0 ) l1=(i1_anear+1)/ML1-1
!         l2=i2_anear/ML2 ; if ( i2_anear<0 ) l2=(i2_anear+1)/ML2-1
!         l3=i3_anear/ML3 ; if ( i3_anear<0 ) l3=(i3_anear+1)/ML3-1
!         i1=i1_anear-l1*ML1
!         i2=i2_anear-l2*ML2
!         i3=i3_anear-l3*ML3
         a_rank(a)=PPP(i1_anear,i2_anear,i3_anear)

         Rx=aa(1,1)*asi(1,a)+aa(1,2)*asi(2,a)+aa(1,3)*asi(3,a)
         Ry=aa(2,1)*asi(1,a)+aa(2,2)*asi(2,a)+aa(2,3)*asi(3,a)
         Rz=aa(3,1)*asi(1,a)+aa(3,2)*asi(2,a)+aa(3,3)*asi(3,a)

         do iorb=1,norb(ik)

            rps2             = Rps(iorb,ik)**2
            NRc              = NRps(iorb,ik)
            L                = lo(iorb,ik)
            rf1(:,:,:,:,:) = 0.d0
            rf2(:,:,:,:,:) = 0.d0
            icheck_tmp3(:)   = 0

            do i=1,M_grid_ion

               i1 = map_grid_ion(1,i)
               i2 = map_grid_ion(2,i)
               i3 = map_grid_ion(3,i)

               id1 = i1 + i1_anear
               id2 = i2 + i2_anear
               id3 = i3 + i3_anear

               do ii3=ii_0,Nintp
                  iii3=id3+ii3
               do ii2=ii_0,Nintp
                  iii2=id2+ii2
               do ii1=ii_0,Nintp
                  iii1=id1+ii1
                  irank=PPP(iii1,iii2,iii3)
                  if ( irank<0 .or. np_grid<irank ) cycle
                  icheck_tmp3(irank)=icheck_tmp3(irank)+1
               end do
               end do
               end do

               if ( PPP(id1,id2,id3)==myrank_g ) then

                  do j3=0,Ndense-1
                     d3=id3*H3+j3*H3d
                  do j2=0,Ndense-1
                     d2=id2*H2+j2*H2d
                  do j1=0,Ndense-1
                     d1=id1*H1+j1*H1d

                     x  = d1*ee(1,1)+d2*ee(1,2)+d3*ee(1,3)-Rx
                     y  = d1*ee(2,1)+d2*ee(2,2)+d3*ee(2,3)-Ry
                     z  = d1*ee(3,1)+d2*ee(3,2)+d3*ee(3,3)-Rz
                     r2 = x*x+y*y+z*z

                     if ( r2>rps2+1.d-10 ) cycle

                     yy1(:) = 0.d0
                     yy2(:) = 0.d0
                     yy3(:) = 0.d0
                     lm1    = lm0
                     r      = sqrt(r2)

                     ir0=irad( int(100.d0*r),ik )
                     do ir=ir0,NRc
                        if ( r<rad1(ir,ik) ) exit
                     end do

                     do L1=abs(L-1),L+1
                        lm1  = lm1 + 1
                        tmp0 = 0.d0
                        err0 = 0.d0
                        if ( abs(x)>1.d-14 .or. abs(y)>1.d-14 .or. abs(z)>1.d-14 .or. L1==0 ) then
                           if ( ir<=2 ) then
                              tmp0=dviod(2,lm1,ik)/(rad1(2,ik)**2)
                              if ( ir<1 ) then
                                 write(*,*) "stop : d_force_ion"
                                 goto 900
                              end if
                           else if ( ir<NRc ) then
                              err0=1.d10
                              do im=1,20
                                 m1=max(1,ir-im)
                                 m2=min(ir+im,NRc)
                                 call polint(rad1(m1,ik),dviod(m1,lm1,ik),m2-m1+1,r,tmp3,err)
                                 if ( abs(err)<err0 ) then
                                    tmp0=tmp3
                                    err0=abs(err)
                                    if ( err0<ep ) exit
                                 end if
                              end do
                              tmp0=tmp0/(r*r)
                           end if
                           maxerr=max(maxerr,err0)
                           do L1z=-L1,L1
                              tmp=tmp0*Ylm(x,y,z,L1,L1z)
                              do m=1,2*L+1
                                 yy1(m)=yy1(m)+tmp*Y1(L,m-L-1,L1,L1z)
                                 yy2(m)=yy2(m)+tmp*Y2(L,m-L-1,L1,L1z)
                                 yy3(m)=yy3(m)+tmp*Y3(L,m-L-1,L1,L1z)
                              end do
                           end do
                        end if
                     end do ! L1

                     if ( any(yy1(1:2*L+1)/=0.d0) .or. any(yy2(1:2*L+1)/=0.d0) .or. any(yy3(1:2*L+1)/=0.d0) ) then

                        do m=1,2*L+1
                           tmp4= yy1(m)*dVd
                           tmp5= yy2(m)*dVd
                           tmp6=-yy3(m)*dVd
                           do ii3=ii_0,Nintp
                              iii3=i3+ii3
                              tmp3=Clag3(ii3,j3)
                           do ii2=ii_0,Nintp
                              iii2=i2+ii2
                              tmp2=Clag2(ii2,j2)*tmp3
                           do ii1=ii_0,Nintp
                              iii1=i1+ii1
                              tmp1=Clag1(ii1,j1)*tmp2
                              rf1(1,iii1,iii2,iii3,m)=rf1(1,iii1,iii2,iii3,m)+tmp1*tmp4
                              rf1(2,iii1,iii2,iii3,m)=rf1(2,iii1,iii2,iii3,m)+tmp1*tmp5
                              rf1(3,iii1,iii2,iii3,m)=rf1(3,iii1,iii2,iii3,m)+tmp1*tmp6
                           end do
                           end do
                           end do
                        end do ! m
                     end if

                  end do ! j3
                  end do ! j2
                  end do ! j1

               end if ! PPP==myrank_g

            end do ! i

            if ( icheck_tmp3(mrnk)>0 ) then

               call mpi_allreduce_nlpp(rf1,rf2,m0,icheck_tmp3,ierr,"R")

               JJ_tmp(:,:) = 0
               rr(:,:)     = 0.d0
               m           = 2*L+1
               mtmp        = 0

               do i3=-mz,mz
                  j3=i3+i3_anear
               do i2=-my,my
                  j2=i2+i2_anear
               do i1=-mx,mx
                  j1=i1+i1_anear
                  if ( PPP(j1,j2,j3)/=myrank_g ) cycle
                  if( any(abs(rf2(1:3,i1,i2,i3,1:m))>1.d-10) )then
                     i=LLL(j1,j2,j3)
                     if ( i<nn1 .or. nn2<i ) cycle
                     mtmp=mtmp+1
                     JJ_tmp(0,mtmp)=i
                     JJ_tmp(1,mtmp)=i1
                     JJ_tmp(2,mtmp)=i2
                     JJ_tmp(3,mtmp)=i3
!                     rr(1,mtmp)=tmp1
!                     rr(2,mtmp)=tmp2
!                     rr(3,mtmp)=tmp3
                  end if
               end do
               end do
               end do

               do s=MSP_0,MSP_1
               do k=MBZ_0,MBZ_1
               do n=MB_0,MB_1

                  if ( occ(n,k,s)<1.d-10 ) cycle

                  const=-2.d0*occ(n,k,s)*dV*dV

                  do m=1,2*L+1

                     do lma=1,nzlma
                        if ( amap(lma)==a .and. lmap(lma)==L .and. mmap(lma)==m-L-1 ) exit
                     end do
                     if ( lma>nzlma ) then
                        write(*,*) "WARNING(d_force_ion)",myrank
                        goto 900
                     end if

                     do j=1,mtmp
                        i =JJ_tmp(0,j)
                        i1=JJ_tmp(1,j)
                        i2=JJ_tmp(2,j)
                        i3=JJ_tmp(3,j)
                        wtmp5(1,lma,n,k,s)=wtmp5(1,lma,n,k,s)+rf2(1,i1,i2,i3,m)*unk(i,n,k,s)
                        wtmp5(2,lma,n,k,s)=wtmp5(2,lma,n,k,s)+rf2(2,i1,i2,i3,m)*unk(i,n,k,s)
                        wtmp5(3,lma,n,k,s)=wtmp5(3,lma,n,k,s)+rf2(3,i1,i2,i3,m)*unk(i,n,k,s)
                     end do

                     do j=1,MJJ(lma)
                        i=JJP(j,lma)
                        ztmp=unk(i,n,k,s)
                        wtmp5(0,lma,n,k,s)=wtmp5(0,lma,n,k,s)+uVk(j,lma,k)*conjg(ztmp)
                     end do
                     wtmp5(0,lma,n,k,s)=wtmp5(0,lma,n,k,s)*const*iuV(lma)

                  end do ! m

               end do ! n
               end do ! k
               end do ! s

            end if

            lm0=lm0+2*L+1

         end do ! iorb   
      end do ! a

      mem=mem-bsintg*size(irad)
      mem=mem-bdreal*size(yy1)
      mem=mem-bsintg*size(icheck_tmp3)
      mem=mem-bsintg*size(JJ_tmp)
      mem=mem-bdreal*size(rr)
      mem=mem-bdreal*size(rf1)
      mem=mem-bdreal*size(rf2)
      deallocate( JJ_tmp )
      deallocate( rr )
      deallocate( rf1 )
      deallocate( rf2 )
      deallocate( icheck_tmp3 )
      deallocate( yy3,yy2,yy1 )
      deallocate( irad )

      do s=MSP_0,MSP_1
      do k=MBZ_0,MBZ_1
      do n=MB_0,MB_1

         if ( occ(n,k,s)<1.d-10 ) cycle

         do i=1,6
            select case(i)
            case(1,3,5)
               j=i+1
               vtmp2(:,:)=wtmp5(:,:,n,k,s)
            case(2,4,6)
               j=i-1
            end select
            do m=1,nrlma_xyz(i)
               nreq=0
               irank=num_2_rank(m,i)
               jrank=num_2_rank(m,j)
               if( irank>=0 )then
                  i1=0
                  do i2=1,lma_nsend(irank)
                  do i3=0,3
                     i1=i1+1
                     sbufnl(i1,irank)=vtmp2(i3,sendmap(i2,irank))
                  end do
                  end do
                  nreq=nreq+1
                  call mpi_isend(sbufnl(1,irank),4*lma_nsend(irank),TYPE_MAIN,irank,1,comm_grid,ireq(nreq),ierr)
               end if
               if( jrank>=0 )then
                  nreq=nreq+1
                  call mpi_irecv(rbufnl(1,jrank),4*lma_nsend(jrank),TYPE_MAIN,jrank,1,comm_grid,ireq(nreq),ierr)
               end if
               call mpi_waitall(nreq,ireq,istatus,ierr)
               if( jrank>=0 )then
                  i1=0
                  do i2=1,lma_nsend(jrank)
                  do i3=0,3
                     i1=i1+1
                     wtmp5(i3,recvmap(i2,jrank),n,k,s)=wtmp5(i3,recvmap(i2,jrank),n,k,s)+rbufnl(i1,jrank)
                  end do
                  end do
               end if
            end do ! m
         end do ! i

         do lma=1,nzlma
            a=amap(lma)
            if ( a<=0 ) cycle
            if ( a_rank(a)/=mrnk ) cycle
            force2(1,a)=force2(1,a)+dble(wtmp5(0,lma,n,k,s)*wtmp5(1,lma,n,k,s))
            force2(2,a)=force2(2,a)+dble(wtmp5(0,lma,n,k,s)*wtmp5(2,lma,n,k,s))
            force2(3,a)=force2(3,a)+dble(wtmp5(0,lma,n,k,s)*wtmp5(3,lma,n,k,s))
         end do

      end do ! n
      end do ! k
      end do ! s

      force3(:,:)=force2(:,:)
      call mpi_allreduce(force3,force2,3*MI,mpi_real8,mpi_sum,mpi_comm_world,ierr)

      mem=mem-bsintg*MI
      mem=mem-bdmain*size(wtmp5)
      mem=mem-bdmain*size(vtmp2)
      deallocate( a_rank )
      deallocate( vtmp2 )
      deallocate( wtmp5 )

      end select

!
! --- Ion-Ion ---
!

      call force_IonIon(force3)

!
! --- Total force ---
!

      force(1:3,1:MI)=force1(1:3,1:MI)+force2(1:3,1:MI)+force3(1:3,1:MI)

!
! --- symmetrization ---
!

!      call sym_force

!
! --- Constraint ---
!

      do a=1,MI
         force(1:3,a)=matmul( tim(1:3,1:3,mkd(a)), force(1:3,a) )
      end do

!
! --- Results ---
!
!      goto 2
      DISP_SWITCH_LOC=.true.
      if ( DISP_SWITCH .and. DISP_SWITCH_LOC ) then
         write(*,*) "force(local)    ",sum( force1**2 )
         do a=1,MI
            write(*,'(1x,i4,1x,4g18.10)') a,force1(:,a),sqrt(sum(force1(:,a)**2))
         end do
         write(*,*) "force(non-local)",sum( force2**2 )
         do a=1,MI
            write(*,'(1x,i4,1x,4g18.10)') a,force2(:,a),sqrt(sum(force2(:,a)**2))
         end do
         write(*,*) "force(ewald)    ",sum( force3**2 )
         do a=1,MI
            write(*,'(1x,i4,1x,4g18.10)') a,force3(:,a),sqrt(sum(force3(:,a)**2))
         end do
         write(*,*) "force(total)    ",sum( force**2 )
         do a=1,MI
            write(*,'(1x,i4,1x,4g18.10)') a,force(:,a),sqrt(sum(force(:,a)**2))
         end do
      end if
 2    continue

!- deallocate -----------------------------------------------------
      mem=mem-bsintg*(size(PPP)+size(LLL)) ; deallocate(PPP,LLL)
      if ( flag_alloc ) then
         mem=mem-bsintg*size(LL) ; deallocate(LL)
      end if
      deallocate( force1,force2,force3 ) ; mem=mem-bdreal*(3*MI)*3
!------------------------------------------------------------------

      Max_mem_inst=max(Max_mem+memax,Max_mem_inst)

      call watch(ctime1,etime1)
      if (DISP_SWITCH) then
         write(*,*) "TIME(D_FORCE_ION_MOL)=",ctime1-ctime0,etime1-etime0
         write(*,*) "mem(MB)=",mem,memax*B2MB,Max_mem_inst*B2MB
      end if

      return

 900  call stop_program

      CONTAINS

      subroutine calcY
      implicit none
      Y1=0.d0 ; Y2=0.d0 ; Y3=0.d0
      Y1( 0, 0, 1, 1) =  0.282094791773878d0
      Y1( 1,-1, 2,-2) = -0.218509686118416d0
      Y1( 1, 0, 2, 1) =  0.218509686118416d0
      Y1( 1, 1, 0, 0) =  0.282094791773878d0
      Y1( 1, 1, 2, 0) = -0.126156626101008d0
      Y1( 1, 1, 2, 2) =  0.218509686118416d0
      Y1( 2,-2, 1,-1) = -0.218509686118416d0
      Y1( 2,-2, 3,-3) = -0.226179013159540d0
      Y1( 2,-2, 3,-1) =  0.058399170081902d0
      Y1( 2,-1, 3,-2) = -0.184674390922372d0
      Y1( 2, 0, 1, 1) = -0.126156626101008d0
      Y1( 2, 0, 3, 1) =  0.202300659403421d0
      Y1( 2, 1, 1, 0) =  0.218509686118416d0
      Y1( 2, 1, 3, 0) = -0.143048168102669d0
      Y1( 2, 1, 3, 2) =  0.184674390922372d0
      Y1( 2, 2, 1, 1) =  0.218509686118416d0
      Y1( 2, 2, 3, 1) = -0.058399170081902d0
      Y1( 2, 2, 3, 3) =  0.226179013159540d0

      Y2( 0, 0, 1,-1) =  0.282094791773878d0
      Y2( 1,-1, 0, 0) =  0.282094791773878d0
      Y2( 1,-1, 2, 0) = -0.126156626101008d0
      Y2( 1,-1, 2, 2) = -0.218509686118416d0
      Y2( 1, 0, 2,-1) =  0.218509686118416d0
      Y2( 1, 1, 2,-2) = -0.218509686118416d0
      Y2( 2,-2, 1, 1) = -0.218509686118416d0
      Y2( 2,-2, 3, 1) =  0.058399170081902d0
      Y2( 2,-2, 3, 3) =  0.226179013159540d0
      Y2( 2,-1, 1, 0) =  0.218509686118416d0
      Y2( 2,-1, 3, 0) = -0.143048168102669d0
      Y2( 2,-1, 3, 2) = -0.184674390922372d0
      Y2( 2, 0, 1,-1) = -0.126156626101008d0
      Y2( 2, 0, 3,-1) =  0.202300659403421d0
      Y2( 2, 1, 3,-2) = -0.184674390922372d0
      Y2( 2, 2, 1,-1) = -0.218509686118416d0
      Y2( 2, 2, 3,-3) =  0.226179013159540d0
      Y2( 2, 2, 3,-1) =  0.058399170081902d0


      Y3( 0, 0, 1, 0) =  0.282094791773878d0
      Y3( 1,-1, 2,-1) =  0.218509686118416d0
      Y3( 1, 0, 0, 0) =  0.282094791773878d0
      Y3( 1, 0, 2, 0) =  0.252313252202016d0
      Y3( 1, 1, 2, 1) =  0.218509686118416d0
      Y3( 2,-2, 3,-2) =  0.184674390922372d0
      Y3( 2,-1, 1,-1) =  0.218509686118416d0
      Y3( 2,-1, 3,-1) =  0.233596680327607d0
      Y3( 2, 0, 1, 0) =  0.252313252202016d0
      Y3( 2, 0, 3, 0) =  0.247766695083476d0
      Y3( 2, 1, 1, 1) =  0.218509686118416d0
      Y3( 2, 1, 3, 1) =  0.233596680327607d0
      Y3( 2, 2, 3, 2) =  0.184674390922372d0
      return
      end subroutine calcY

      END SUBROUTINE force_ion_mol
