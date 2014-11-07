MODULE ForcePSnonLoc2
  use atom_module, only: ki_atom,aa_atom
  use aa_module, only: aa
  use VarPSMember
  use VarPSMemberG
  use rgrid_module
  use array_bound_module
  use parallel_module
  use localpot_module, only: Vloc
  use ParaRGridComm, only: threeWayComm
  use ForceSub
  implicit none
  PRIVATE
  PUBLIC :: calcForcePSnonLoc2
  real(8),allocatable :: y2a(:,:,:),y2b(:,:,:)

CONTAINS

  SUBROUTINE calcForcePSnonLoc2(MI,force2)
    use bz_module
    use wf_module
    use watch_module
  use ps_nloc2_variables
    implicit none
    integer,intent(IN) :: MI
    real(8),intent(OUT) :: force2(3,MI)
!    integer :: lma1,lma2
    integer :: ib1,ib2
    integer :: i,j,k,s,n,ir,L1,L1z,NRc,irank,jrank
    integer :: nreq,max_nreq
    integer :: a0,lm0,lm1,lma
    integer :: a,ik,m,L
    integer :: ierr,ir0
    integer,allocatable :: ireq(:),istatus(:,:),ilm1(:,:,:)
    real(8) :: a1,a2,a3
    real(8) :: kr,c
    real(8) :: tmp
    real(8) :: ctt(0:4),ett(0:4)
    real(8) :: yy1,yy2,yy3
    real(8),allocatable :: work2(:,:),duVdR(:,:,:)
!    real(8) :: Dij_f,const_f
    real(8) :: Dij_f(MB_0:MB_1),const_f(MB_0:MB_1)
#ifdef _DRSDFT_
    real(8) :: ztmp
    real(8),allocatable :: wtmp5(:,:,:,:,:)
    real(8),allocatable :: uVunk_tmp(:,:,:,:)
#else
    complex(8) :: ztmp
    complex(8),allocatable :: wtmp5(:,:,:,:,:)
    complex(8),allocatable :: uVunk_tmp(:,:,:,:)
#endif
    integer :: i0,iorb0
    real(8) :: forceQ(3,MI)

    real(8),parameter :: pi2=2.d0*acos(-1.d0)
! d1,d2 need to be removed from this module, define a different variable name
    real(8) :: d1,d2,d3
    integer :: i1,i2,i3

    INTERFACE
      FUNCTION Ylm(x,y,z,l,m)
        real(8) :: Ylm
        real(8),intent(IN) :: x,y,z
        integer,intent(IN) :: l,m
      END FUNCTION Ylm
    END INTERFACE
    
    call getSHY
    call getCijLM

    force2(:,:) = 0.0d0

    if ( Mlma <= 0 ) return

    maxerr=0.d0
    ctt(:)=0.d0
    ett(:)=0.d0

    call setLocalIndexForBoundary(Igrid)
    !OUT: a1b,a2b,a3b,ab1,ab2,ab3
    call setConstGridWidth(Ngrid)
    !OUT: ML1,ML2,ML3,c1,c2,c3

    if ( .not.allocated(ilm1) ) then
      L1=maxval(lo)+1
      n=maxval(norb)
      allocate( ilm1(0:L1,n,Nelement_) ) ; ilm1=0
      do ik=1,Nelement_
        lm1=0
        do iorb=1,norb(ik)
          L=lo(iorb,ik)
          do L1=abs(L-1),L+1
            lm1=lm1+1
            ilm1(L1,iorb,ik)=lm1
          end do
        end do
      end do
    end if

    if ( .not.allocated(y2b) ) then
      lm1=maxval(ilm1)
      NRc=maxval(NRps)
      allocate( y2b(NRc,lm1,Nelement_) )
      y2b=0.d0
      do ik=1,Nelement_
        do iorb=1,norb(ik)
          L=lo(iorb,ik)
          do L1=abs(L-1),L+1
            lm1=ilm1(L1,iorb,ik)
            d1=0.d0
            d2=0.d0
            call spline(rad1(1,ik),dviod(1,lm1,ik),NRps(iorb,ik),d1,d2,y2b(1,lm1,ik))
          end do
        end do
      end do
    end if

    allocate( wtmp5(0:3,nzlma,MB_0:MB_1,MBZ_0:MBZ_1,MSP_0:MSP_1) )
    allocate( duVdR(3,MMJJ,nzlma) )

    call setLocalAtoms(MI,aa_atom,Igrid)
    !OUT: isAtomInThisNode
do a=1,MI
!write(3200+myrank,*) a,isAtomInThisNode(a)
enddo

!$OMP parallel

!$OMP workshare
    wtmp5=zero
    duVdR=0.d0
!$OMP end workshare

#ifndef _SPLINE_
 call setIndexForAtomCenteredGrid(Nelement_,NRps,rad1)
 !OUT: irad(0:3000,1:Nelement_)
#endif

!$OMP master
    call watch(ctt(0),ett(0))
!$OMP end master

!!$OMP do schedule(dynamic) firstprivate( maxerr ) &
!!$OMP    private( a,L,m,iorb,ik,Rx,Ry,Rz,NRc,d1,d2,d3,x,y,z,r  &
!!$OMP            ,ir,ir0,yy1,yy2,yy3,err0,err,tmp0,tmp1,m1,m2  &
!!$OMP            ,lma,j,L1,L1z,lm1,im )
    do lma=1,nzlma
      call getAtomInfoFrom_lma(lma)
      !OUT: a,L,m,iorb,ik,notAtom
!write(3000+myrank,'(" lma=",I7," noAtomHere=",L7,I4)') lma,noAtomHere,ikind
      if (noAtomHere) cycle
      call getAtomPosition_real
      !OUT: Rx,Ry,Rz
!write(3000+myrank,'(2A7)') 'iorb','ikind'
!write(3000+myrank,'(2I7)') iorb,ikind
!write(3000+myrank,'(4A7)') 'iatom','iL','im','ikind'
!write(3000+myrank,'(4I7)') iatom,iL,im,ikind
      NRc=NRps(iorb,ikind)
!!$OMP parallel do firstprivate( maxerr ) &
!!$OMP             private( d1,d2,d3,x,y,z,r,ir0,yy1,yy2,yy3,lm1,err,err0 &
!!$OMP                     ,tmp0,tmp1,m1,m2,j,L1,im,L1z )
!write(3100+myrank,'(" lma=",I6," MJJ_MAP(lma)=",I6,A40)') lma,MJJ_MAP(lma),repeat("-",30)
!write(3300+myrank,'(2A6,4A20)') 'lma','j','x','y','z','r'
      do j=1,MJJ_MAP(lma)
        call getAtomCenteredPositionFrom_lma(lma,j)
        !OUT: x,y,z,r
!write(3300+myrank,'(2I6,4g20.7)') lma,j,x,y,z,r
#ifndef _SPLINE_
        ir0=irad( int(100.d0*r),ikind )
!write(3100+myrank,'(2A6,3A20)') 'ir0','NRc','r','rad1(ir0)','rad(NRc)'
!write(3100+myrank,'(2I6,3G20.7)') ir0,NRc,r,rad1(ir0,ikind),rad1(NRc,ikind)
        do ir=ir0,NRc
          if ( r<rad1(ir,ikind) ) exit
        end do
!write(3100+myrank,'(" j=",I6," ir=",I5)') j,ir
#endif
        yy1=0.d0
        yy2=0.d0
        yy3=0.d0
        do L1=abs(iL-1),iL+1
          lm1=ilm1(L1,iorb,ikind)
          if ( abs(x)>1.d-14 .or. abs(y)>1.d-14 .or. abs(z)>1.d-14 .or. L1==0 ) then
!write(3100+myrank,'(" j=",I6," ir=",I5)') j,ir
            call interpolate_dviod(lm1,ir,NRc)
            !OUT: tmp0
            do L1z=-L1,L1
              tmp=tmp0*Ylm(x,y,z,L1,L1z)
              yy1=yy1+tmp*SH_Y1(iL,im,L1,L1z)
              yy2=yy2+tmp*SH_Y2(iL,im,L1,L1z)
              yy3=yy3+tmp*SH_Y3(iL,im,L1,L1z)
            end do
          end if
        end do ! L1

        duVdR(1,j,lma)= yy1
        duVdR(2,j,lma)= yy2
        duVdR(3,j,lma)=-yy3

      end do ! j
!!$OMP end parallel do
    end do ! lma
!!$OMP end do

!$OMP master
    call watch(ctt(1),ett(1))
!$OMP end master

!#ifndef _SPLINE_
!!$OMP single
!    deallocate( irad )
!!$OMP end single
!#endif

    do s=MSP_0,MSP_1
      do k=MBZ_0,MBZ_1
!$OMP do schedule(dynamic) private( c,i,d1,d2,d3,kr,ztmp,i1,i2,i3 )
        do n=MB_0,MB_1
          if ( occ(n,k,s) < 1.d-10 ) cycle
!          c=-2.d0*occ(n,k,s)*dV*dV
          c=1.d0
          do lma=1,nzlma
            if ( MJJ_MAP(lma) == MJJ(lma) ) then
!write(3500+myrank,'(6A6)') 's','k','n','lma','j','i','----1'
#ifdef _DRSDFT_
              do j=1,MJJ(lma)
                i=JJP(j,lma)
                wtmp5(0,lma,n,k,s)=wtmp5(0,lma,n,k,s)+uVk(j,lma,k)*unk(i,n,k,s)
                wtmp5(1,lma,n,k,s)=wtmp5(1,lma,n,k,s)+duVdR(1,j,lma)*unk(i,n,k,s)
                wtmp5(2,lma,n,k,s)=wtmp5(2,lma,n,k,s)+duVdR(2,j,lma)*unk(i,n,k,s)
                wtmp5(3,lma,n,k,s)=wtmp5(3,lma,n,k,s)+duVdR(3,j,lma)*unk(i,n,k,s)
              end do
              wtmp5(0,lma,n,k,s)=c*wtmp5(0,lma,n,k,s)
#else
              do j=1,MJJ(lma)
                i=JJP(j,lma)
                wtmp5(0,lma,n,k,s)=wtmp5(0,lma,n,k,s)+uVk(j,lma,k)*conjg(unk(i,n,k,s))
              end do
              wtmp5(0,lma,n,k,s)=c*wtmp5(0,lma,n,k,s)
              do j=1,MJJ_MAP(lma)
                i=JJP(j,lma)
!write(3500+myrank,'(6I6)') s,k,n,lma,j,i
                d1=c1*JJ_MAP(1,j,lma)+JJ_MAP(4,j,lma)
                d2=c2*JJ_MAP(2,j,lma)+JJ_MAP(5,j,lma)
                d3=c3*JJ_MAP(3,j,lma)+JJ_MAP(6,j,lma)
                kr=pi2*(kbb(1,k)*d1+kbb(2,k)*d2+kbb(3,k)*d3)
                ztmp=dcmplx(cos(kr),sin(kr))*unk(i,n,k,s)
                wtmp5(1,lma,n,k,s)=wtmp5(1,lma,n,k,s)+duVdR(1,j,lma)*ztmp
                wtmp5(2,lma,n,k,s)=wtmp5(2,lma,n,k,s)+duVdR(2,j,lma)*ztmp
                wtmp5(3,lma,n,k,s)=wtmp5(3,lma,n,k,s)+duVdR(3,j,lma)*ztmp
              end do
#endif
            else ! --- MJJ(lma) /= MJJ_MAP(lma) ---
!write(3500+myrank,'(6A6)') 's','k','n','lma','j','i','----2'
#ifdef _DRSDFT_
              do j=1,MJJ(lma)
                i=JJP(j,lma)
                wtmp5(0,lma,n,k,s)=wtmp5(0,lma,n,k,s)+uVk(j,lma,k)*unk(i,n,k,s)
              end do
              wtmp5(0,lma,n,k,s)=c*wtmp5(0,lma,n,k,s)
              do j=1,MJJ_MAP(lma)
                i1=JJ_MAP(1,j,lma)
                i2=JJ_MAP(2,j,lma)
                i3=JJ_MAP(3,j,lma)
                i = i1-a1b + (i2-a2b)*ab1 + (i3-a3b)*ab1*ab2 + ML_0
                wtmp5(1,lma,n,k,s)=wtmp5(1,lma,n,k,s)+duVdR(1,j,lma)*unk(i,n,k,s)
                wtmp5(2,lma,n,k,s)=wtmp5(2,lma,n,k,s)+duVdR(2,j,lma)*unk(i,n,k,s)
                wtmp5(3,lma,n,k,s)=wtmp5(3,lma,n,k,s)+duVdR(3,j,lma)*unk(i,n,k,s)
              end do
#else
              do j=1,MJJ(lma)
                i=JJP(j,lma)
!write(3500+myrank,'(6I6)') s,k,n,lma,j,i
                wtmp5(0,lma,n,k,s)=wtmp5(0,lma,n,k,s)+uVk(j,lma,k)*conjg(unk(i,n,k,s))
              end do
              wtmp5(0,lma,n,k,s)=c*wtmp5(0,lma,n,k,s)
              do j=1,MJJ_MAP(lma)
                i1=JJ_MAP(1,j,lma)
                i2=JJ_MAP(2,j,lma)
                i3=JJ_MAP(3,j,lma)
                i = i1-a1b + (i2-a2b)*ab1 + (i3-a3b)*ab1*ab2 + ML_0
!write(3500+myrank,'(5I6,6I4)') s,k,n,lma,j,JJ_MAP(1:6,j,lma)
!write(3500+myrank,'(6I6)') s,k,n,lma,j,i
                d1=c1*i1+JJ_MAP(4,j,lma)
                d2=c2*i2+JJ_MAP(5,j,lma)
                d3=c3*i3+JJ_MAP(6,j,lma)
                kr=pi2*(kbb(1,k)*d1+kbb(2,k)*d2+kbb(3,k)*d3)
                ztmp=dcmplx(cos(kr),sin(kr))*unk(i,n,k,s)
                wtmp5(1,lma,n,k,s)=wtmp5(1,lma,n,k,s)+duVdR(1,j,lma)*ztmp
                wtmp5(2,lma,n,k,s)=wtmp5(2,lma,n,k,s)+duVdR(2,j,lma)*ztmp
                wtmp5(3,lma,n,k,s)=wtmp5(3,lma,n,k,s)+duVdR(3,j,lma)*ztmp
              end do
#endif
            end if
          end do ! lma
        end do ! n
!$OMP end do
      end do ! k
    end do ! s

!$OMP master
    call watch(ctt(2),ett(2))
!$OMP end master

!$OMP single
    deallocate( duVdR )
!$OMP end single

!$OMP workshare
    force2(:,:)=0.d0
!$OMP end workshare

!$OMP single
    do s=MSP_0,MSP_1
      do k=MBZ_0,MBZ_1
        do n=MB_0,MB_1,MB_d
          ib1=n
          ib2=min(ib1+MB_d-1,MB_1)
          if ( occ(n,k,s) < 1.d-10 ) cycle
          call threeWayComm(nrlma_xyz,num_2_rank,sendmap,recvmap,lma_nsend,sbufnl,rbufnl,nzlma,ib1,ib2,wtmp5(0,1,ib1,k,s),3)
        end do ! n
      end do ! k
    end do ! s
    do s=MSP_0,MSP_1
      do k=MBZ_0,MBZ_1
        do m=1,N_nzqr
          lma1=nzqr_pair(m,1)
          lma2=nzqr_pair(m,2)
          a=amap(lma1)
          if ( a <= 0 ) cycle
!if (myrank==0) write(150,*) a,isAtomInThisNode(a)
          if ( isAtomInThisNode(a) ) then
            do n=MB_0,MB_1
              Dij_f(n)=Dij(m,s)-esp(n,k,s)*qij_f(m)
              const_f(n)=Dij_f(n)*(-2.d0)*occ(n,k,s)*dV*dV
            enddo
!if (myrank==0) write(150,'(4a5,6a20)') 's','k','n','m','dV','const_f','Dij_f','Dij(m,s)','esp(n,k,s)','qij_f(m)'
do n=MB_0,MB_1
!if (myrank==0) write(150,'(4i5,6g20.7)') s,k,n,m,dV,const_f(n),Dij_f(n),Dij(m,s),esp(n,k,s),qij_f(m)
enddo
!            Dij_f(MB_0:MB_1)=Dij(m,s)-esp(MB_0:MB_1,k,s)*qij_f(m)
!            const_f(MB_0:MB_1)=Dij_f(MB_0:MB_1)*(-2.d0)*occ(MB_0:MB_1,k,s)*dV*dV
            if (lma1<lma2) stop 'Nzqr_pair is strange'
            if (lma1==lma2) then
              force2(1,a)=force2(1,a)+sum(const_f(MB_0:MB_1)*real(wtmp5(0,lma1,MB_0:MB_1,k,s)*wtmp5(1,lma2,MB_0:MB_1,k,s),8))
              force2(2,a)=force2(2,a)+sum(const_f(MB_0:MB_1)*real(wtmp5(0,lma1,MB_0:MB_1,k,s)*wtmp5(2,lma2,MB_0:MB_1,k,s),8))
              force2(3,a)=force2(3,a)+sum(const_f(MB_0:MB_1)*real(wtmp5(0,lma1,MB_0:MB_1,k,s)*wtmp5(3,lma2,MB_0:MB_1,k,s),8))
            else
              force2(1,a)=force2(1,a)+sum(const_f(MB_0:MB_1)*real(wtmp5(0,lma1,MB_0:MB_1,k,s)*wtmp5(1,lma2,MB_0:MB_1,k,s),8))
              force2(2,a)=force2(2,a)+sum(const_f(MB_0:MB_1)*real(wtmp5(0,lma1,MB_0:MB_1,k,s)*wtmp5(2,lma2,MB_0:MB_1,k,s),8))
              force2(3,a)=force2(3,a)+sum(const_f(MB_0:MB_1)*real(wtmp5(0,lma1,MB_0:MB_1,k,s)*wtmp5(3,lma2,MB_0:MB_1,k,s),8))
              force2(1,a)=force2(1,a)+sum(const_f(MB_0:MB_1)*real(wtmp5(0,lma2,MB_0:MB_1,k,s)*wtmp5(1,lma1,MB_0:MB_1,k,s),8))
              force2(2,a)=force2(2,a)+sum(const_f(MB_0:MB_1)*real(wtmp5(0,lma2,MB_0:MB_1,k,s)*wtmp5(2,lma1,MB_0:MB_1,k,s),8))
              force2(3,a)=force2(3,a)+sum(const_f(MB_0:MB_1)*real(wtmp5(0,lma2,MB_0:MB_1,k,s)*wtmp5(3,lma1,MB_0:MB_1,k,s),8))
            end if
          endif
!if (myrank==0) write(150,'(4a5,6a20)') 's','k','n','m','dV','const_f','Dij_f','Dij(m,s)','esp(n,k,s)','qij_f(m)'
!if (myrank==0) write(151,'(5a5,a20)') 's','k','n','m','lma1','wtmp5(0:3,lma1,n,k,s)'
!if (myrank==0) write(152,'(5a5,a20)') 's','k','n','m','lma2','wtmp5(0:3,lma2,n,k,s)'
do n=MB_0,MB_1
!if (myrank==0) write(150,'(4i5,6g20.7)') s,k,n,m,dV,const_f(n),Dij_f(n),Dij(m,s),esp(n,k,s),qij_f(m)
!if (myrank==0) write(151,'(5i5,8g20.7)') s,k,n,m,lma1,wtmp5(0:3,lma1,n,k,s)
!if (myrank==0) write(152,'(5i5,8g20.7)') s,k,n,m,lma2,wtmp5(0:3,lma2,n,k,s)
enddo
        end do ! m
      end do ! k
    end do ! s
! new version
! not working ---------------------------------------------------------------------------------------
!    do s=MSP_0,MSP_1
!      do k=MBZ_0,MBZ_1
!        do n=MB_0,MB_1,MB_d
!          ib1=n
!          ib2=min(ib1+MB_d-1,MB_1)
!          if ( occ(n,k,s) == 0.d0 ) cycle
!          call threeWayComm(nrlma_xyz,num_2_rank,sendmap,recvmap,lma_nsend,sbufnl,rbufnl,nzlma,ib1,ib2,wtmp5(0,1,ib1,k,s),3)
!          
!if (myrank==0) write(150,'(4a5,6a20)') 's','k','n','m','dV','const_f','Dij_f','Dij(m,s)','esp(n,k,s)','qij_f(m)'
!          do m=1,N_nzqr
!            lma1=nzqr_pair(m,1)
!            lma2=nzqr_pair(m,2)
!            a=amap(lma1)
!            if ( a <= 0 ) cycle
!            if ( a_rank(a) ) then
!              Dij_f=Dij(m,s)-esp(n,k,s)*qij_f(m)
!              const_f=Dij_f*(-2.d0)*occ(n,k,s)*dV*dV
!if (myrank==0) write(150,'(4i5,6g20.7)') s,k,n,m,dV,const_f,Dij_f,Dij(m,s),esp(n,k,s),qij_f(m)
!              if (lma1<lma2) stop 'Nzqr_pair is strange'
!              if (lma1==lma2) then
!                force2(1,a)=force2(1,a)+const_f*real(wtmp5(0,lma1,n,k,s)*wtmp5(1,lma2,n,k,s),8)
!                force2(2,a)=force2(2,a)+const_f*real(wtmp5(0,lma1,n,k,s)*wtmp5(2,lma2,n,k,s),8)
!                force2(3,a)=force2(3,a)+const_f*real(wtmp5(0,lma1,n,k,s)*wtmp5(3,lma2,n,k,s),8)
!              else
!                force2(1,a)=force2(1,a)+const_f*real(wtmp5(0,lma1,n,k,s)*wtmp5(1,lma2,n,k,s),8)
!                force2(2,a)=force2(2,a)+const_f*real(wtmp5(0,lma1,n,k,s)*wtmp5(2,lma2,n,k,s),8)
!                force2(3,a)=force2(3,a)+const_f*real(wtmp5(0,lma1,n,k,s)*wtmp5(3,lma2,n,k,s),8)
!                force2(1,a)=force2(1,a)+const_f*real(wtmp5(0,lma2,n,k,s)*wtmp5(1,lma1,n,k,s),8)
!                force2(2,a)=force2(2,a)+const_f*real(wtmp5(0,lma2,n,k,s)*wtmp5(2,lma1,n,k,s),8)
!                force2(3,a)=force2(3,a)+const_f*real(wtmp5(0,lma2,n,k,s)*wtmp5(3,lma1,n,k,s),8)
!
!              end if
!            endif
!          end do ! m
!
!        end do ! n
!      end do ! k
!    end do ! s
! not working =======================================================================================

    allocate( uVunk_tmp(nzlma,MB_0:MB_1,MBZ_0:MBZ_1,MSP_0:MSP_1) )
    uVunk_tmp=zero
#ifdef _DRSDFT_
    uVunk_tmp(:,:,:,:)=wtmp5(0,:,:,:,:)*dV
#else
    uVunk_tmp(:,:,:,:)=conjg(wtmp5(0,:,:,:,:)*dV)
#endif
    call mpi_allreduce(MPI_IN_PLACE,force2,3*MI,mpi_real8,mpi_sum,mpi_comm_world,ierr)
    deallocate( wtmp5 )
!$OMP end single

!$OMP master
    call watch(ctt(3),ett(3))
!$OMP end master
do a=1,MI
if (myrank==0) write(200,'(I5,A9,3g20.7)') a,'Dij',force2(1:3,a)
enddo

    call calcForceQ(uVunk_tmp,MI,forceQ)
do a=1,MI
if (myrank==0) write(200,'(I5,A9,3g20.7)') a,'Q',forceQ(1:3,a)
enddo

    do a=1,MI
      force2(1,a)=force2(1,a)+forceQ(1,a)
      force2(2,a)=force2(2,a)+forceQ(2,a)
      force2(3,a)=force2(3,a)+forceQ(3,a)
    enddo

!$OMP master
    call watch(ctt(3),ett(3))
!$OMP end master

!$OMP end parallel

    if ( disp_switch_parallel ) then
       write(*,*) "time(force_nloc_uspp_1)",ctt(1)-ctt(0),ett(1)-ett(0)
       write(*,*) "time(force_nloc_uspp_2)",ctt(2)-ctt(1),ett(2)-ett(1)
       write(*,*) "time(force_nloc_uspp_3)",ctt(3)-ctt(2),ett(3)-ett(2)
       write(*,*) "time(force_nloc_uspp_4)",ctt(4)-ctt(3),ett(4)-ett(3)
    end if

    return
  END SUBROUTINE calcForcePSnonLoc2

#ifdef _USPP_
  SUBROUTINE calcForceQ(uVunk_tmp,MI,forceQ)
    use bz_module
    use wf_module
    use watch_module
  use VarParaPSnonLocG
  use ps_nloc2_variables, only: lmap,amap,mmap,nzlma,iorbmap
    implicit none
    integer,intent(IN) :: MI
#ifdef _DRSDFT_
    real(8),intent(IN) :: uVunk_tmp(nzlma,MB_0:MB_1,MBZ_0:MBZ_1,MSP_0:MSP_1)
#else
    complex(8),intent(IN) :: uVunk_tmp(nzlma,MB_0:MB_1,MBZ_0:MBZ_1,MSP_0:MSP_1)
#endif
    real(8),intent(INOUT) :: forceQ(3,MI)
    integer :: i1,i2,i3
    integer :: i,j,k,s,n,ir,L,L1,L1z,NRc,irank,jrank
    integer :: iqr,ll3,cJ
    integer :: ib1,ib2
    real(8) :: yq1,yq2,yq3,tmp2,tmp3
    integer :: nreq,max_nreq
    integer :: a,a0,ik,m,lm0,lm1,lma,m1,m2
    integer :: ierr,ir0
    integer,allocatable :: ireq(:),istatus(:,:)
    real(8),parameter :: ep=1.d-8
    real(8) :: err,err0
    real(8) :: a1,a2,a3
    real(8) :: kr,c
    real(8),parameter :: pi2=2.d0*acos(-1.d0)
    real(8) :: tmp,tmp1
    real(8) :: yy1,yy2,yy3
    real(8),allocatable :: work2(:,:),dQY(:,:,:)
    real(8) :: dQY_tmp(1:3)
#ifdef _DRSDFT_
    real(8) :: ztmp
    real(8),allocatable :: rtmp5(:,:,:),rtmp2(:,:)
#else
    complex(8) :: ztmp
    complex(8),allocatable :: rtmp5(:,:,:),rtmp2(:,:)
#endif
    integer :: i0,iorb0
    integer :: k1,k2,k3
    integer :: d1,d2,d3

    INTERFACE
      FUNCTION Ylm(x,y,z,l,m)
        real(8) :: Ylm
        real(8),intent(IN) :: x,y,z
        integer,intent(IN) :: l,m
      END FUNCTION Ylm
    END INTERFACE
    

    forceQ(:,:) = 0.0d0

    maxerr=0.d0

    allocate( rtmp5(0:2,N_nzqr,MSP_0:MSP_1) )
    allocate( dQY(3,MMJJ_Q,N_nzqr) )

!$OMP parallel

!$OMP workshare
    rtmp5=zero
    dQY=0.d0
!$OMP end workshare

!if (myrank==0) write(220,'(4a5,a20)') 'iqr','j','ll3','L1','tmp0'
!!$OMP do schedule(dynamic) firstprivate( maxerr ) &
!!$OMP    private( a,L,m,ik,Rx,Ry,Rz,NRc,d1,d2,d3,x,y,z,r  &
!!$OMP            ,ir,ir0,yy1,yy2,yy3,tmp1,m1,m2  &
!!$OMP            ,lma,j,L1,L1z,lm1 )
    do iqr=1,N_nzqr
      if (iqr>N_nzqr) then
        write(200,*) 'myrank= ',myrank,' - memory leakage suspected'
        exit
      endif
      call getAtomInfoFrom_iqr(iqr)
      !OUT: iatom,ikind,ik1,ikk1,ik2
      !     lma1,lma2,l_1,l_2,m_1,m_2,iorb1,iorb2,noAtomHere
      if ( noAtomHere ) cycle
!      if ( ikk1==0    ) cycle
      call getAtomPosition_real
      !OUT: Rx,Ry,Rz
      NRc=max(NRps(iorb1,ikind),NRps(iorb2,ikind))
!write(4300+myrank,'(2A6,A20)') 'j','iqr','dQY(1:3,j,iqr)'
!write(4600+myrank,'(2A6,4A20)') 'iqr','j','x','y','z','r'
      do j=1,MJJ_MAP_Q(iqr)
        call getAtomCenteredPositionFrom_iqr(iqr,j,x,y,z,r)
        !OUT: x,y,z,r
!write(4600+myrank,'(2I6,4g20.7)') iqr,j,x,y,z,r
#ifndef _SPLINE_
        ir0=irad( int(100.d0*r),ikind )
        do ir=ir0,NRc
          if ( r<rad1(ir,ikind) ) exit
        end do
#endif
!if (myrank==0) write(220,'(6a5,a20)') 'cJ','ik','k2','ll3','m1','L1','tmp0'
        dQY_tmp=0.d0
        do ll3=1,nl3v(ik2,ikind)
        ! ll3: info of l_1,m_1,l_2,m_2
          L=l3v(ll3,ik2,ikind)-1
          ! L: {s:0,p:1,d:2,...}
          cJ=0
          ! cJ is just the index to run L1 which do not start from 1
          do L1=abs(L-1),L+1
            cJ=cJ+1
            if ( abs(x)>1.d-14 .or. abs(y)>1.d-14 .or. abs(z)>1.d-14 .or. L1==0 ) then
              call interpolate_dqrL(ll3,cJ,ir,NRc)
              !OUT: tmp0,maxerr
              yy1=0.d0
              yy2=0.d0
              yy3=0.d0
              do L1z=-L1,L1
              ! L1z: M
              ! M summation can be done, for Q depend on L, {l1,m1,l2,m2}
                yq1=0.d0
                yq2=0.d0
                yq3=0.d0
                do M=-L,L
                  yq1=yq1+C_ijLM(L,M,l_1,m_1,l_2,m_2)*SH_Y1(L,M,L1,L1z)
                  yq2=yq2+C_ijLM(L,M,l_1,m_1,l_2,m_2)*SH_Y2(L,M,L1,L1z)
                  yq3=yq3+C_ijLM(L,M,l_1,m_1,l_2,m_2)*SH_Y3(L,M,L1,L1z)
                enddo
                tmp3=Ylm(x,y,z,L1,L1z)
                yy1=yy1+yq1*tmp3
                yy2=yy2+yq2*tmp3
                yy3=yy3+yq3*tmp3
              end do ! L1z
              dQY_tmp(1)=dQY_tmp(1)+tmp0*yy1
              dQY_tmp(2)=dQY_tmp(2)+tmp0*yy2
              dQY_tmp(3)=dQY_tmp(3)+tmp0*yy3
            end if
          end do ! L1
        enddo ! ll3
        dQY(1,j,iqr)=dQY_tmp(1)
        dQY(2,j,iqr)=dQY_tmp(2)
        dQY(3,j,iqr)=dQY_tmp(3)
!write(4300+myrank,'(2I6,3g20.7)') j,iqr,dQY(1:3,j,iqr)
      end do ! j
    end do ! iqr
!!$OMP end do

!if (myrank==0) write(230,'(5a5,a20)') 's','k','n','iqr','j','ztmp'
!write(4000+myrank,'(2A8)') 'myrank','N_nzqr'
!write(4000+myrank,'(2I8)') myrank,N_nzqr
    do s=MSP_0,MSP_1
!!$OMP do schedule(dynamic) private( i,kr,ztmp )
      do iqr=1,N_nzqr
        if (iqr>N_nzqr) then
          write(200,*) 'myrank= ',myrank,' - memory leakage suspected'
          exit
        endif
!write(4000+myrank,'(4A6)') 's','iqr','j','i'
!write(4100+myrank,'(4A6,2A20)') 's','iqr','j','i','Vloc(i,s)','dQY(1:3,j,iqr)'
        do j=1,MJJ_MAP_Q(iqr)
          i1=JJ_MAP_Q(1,j,iqr)
          i2=JJ_MAP_Q(2,j,iqr)
          i3=JJ_MAP_Q(3,j,iqr)
          i=i1-a1b+(i2-a2b)*ab1+(i3-a3b)*ab1*ab2+ML_0
!write(4000+myrank,'(4I6)') s,iqr,j,i
!write(4100+myrank,'(4I6,4g20.7)') s,iqr,j,i,Vloc(i,s),dQY(1:3,j,iqr)
          ztmp=Vloc(i,s)
          rtmp5(0,iqr,s)=rtmp5(0,iqr,s)+dQY(1,j,iqr)*ztmp
          rtmp5(1,iqr,s)=rtmp5(1,iqr,s)+dQY(2,j,iqr)*ztmp
          rtmp5(2,iqr,s)=rtmp5(2,iqr,s)-dQY(3,j,iqr)*ztmp
!if (myrank==0) write(230,'(5I5,g20.7)') s,k,n,iqr,j,Vloc(i,s)
        end do ! j
      end do ! iqr
!!$OMP end do
    end do ! s
    
    do s=MSP_0,MSP_1
      do iqr=1,N_nzqr
        do i=0,2
!          write(4200+myrank+i*20,'(2g20.7)') s,iqr,rtmp5(i,iqr,s)
        enddo
      enddo
    enddo


!$OMP single
    deallocate( dQY )
!$OMP end single


!$OMP single
    do s=MSP_0,MSP_1
      call threeWayComm(nrqr_xyz,num_2_rank_Q,sendmap_Q,recvmap_Q,qr_nsend,sbufnl_Q,rbufnl_Q,N_nzqr,1,1,rtmp5(0,1,1),2)
    enddo
!if (myrank==0) write(240,'(4a5,7a20)') 's','k','n','m','tmp1','','rtmp5(0,m,n,k,s)','occ(n,k,s)','dV','','uVunk_tmp(lma1,n,k,s)'
!if (myrank==0) write(250,'(4a5,a20)') 's','k','n','m','forceQ(1,a)'
    do s=MSP_0,MSP_1
      do k=MBZ_0,MBZ_1
        do m=1,N_nzqr
          lma1=nzqr_pair(m,1)
          lma2=nzqr_pair(m,2)
          a=amap(lma1)
          if ( a <= 0 ) cycle
          if ( isAtomInThisNode(a) ) then
            if (lma1<lma2) stop 'Nzqr_pair is strange'
            tmp1=-dV
#ifdef _DRSDFT_
            if (lma1==lma2) then
              forceQ(1,a)=forceQ(1,a)+sum(occ(MB_0:MB_1,k,s)*real(uVunk_tmp(lma1,MB_0:MB_1,k,s)*uVunk_tmp(lma2,MB_0:MB_1,k,s),8)*rtmp5(0,m,s))*tmp1
              forceQ(2,a)=forceQ(2,a)+sum(occ(MB_0:MB_1,k,s)*real(uVunk_tmp(lma1,MB_0:MB_1,k,s)*uVunk_tmp(lma2,MB_0:MB_1,k,s),8)*rtmp5(1,m,s))*tmp1
              forceQ(3,a)=forceQ(3,a)+sum(occ(MB_0:MB_1,k,s)*real(uVunk_tmp(lma1,MB_0:MB_1,k,s)*uVunk_tmp(lma2,MB_0:MB_1,k,s),8)*rtmp5(2,m,s))*tmp1
            else
              forceQ(1,a)=forceQ(1,a)+sum(occ(MB_0:MB_1,k,s)*real(uVunk_tmp(lma1,MB_0:MB_1,k,s)*uVunk_tmp(lma2,MB_0:MB_1,k,s),8)*rtmp5(0,m,s))*tmp1
              forceQ(2,a)=forceQ(2,a)+sum(occ(MB_0:MB_1,k,s)*real(uVunk_tmp(lma1,MB_0:MB_1,k,s)*uVunk_tmp(lma2,MB_0:MB_1,k,s),8)*rtmp5(1,m,s))*tmp1
              forceQ(3,a)=forceQ(3,a)+sum(occ(MB_0:MB_1,k,s)*real(uVunk_tmp(lma1,MB_0:MB_1,k,s)*uVunk_tmp(lma2,MB_0:MB_1,k,s),8)*rtmp5(2,m,s))*tmp1
              forceQ(1,a)=forceQ(1,a)+sum(occ(MB_0:MB_1,k,s)*real(uVunk_tmp(lma2,MB_0:MB_1,k,s)*uVunk_tmp(lma1,MB_0:MB_1,k,s),8)*rtmp5(0,m,s))*tmp1
              forceQ(2,a)=forceQ(2,a)+sum(occ(MB_0:MB_1,k,s)*real(uVunk_tmp(lma2,MB_0:MB_1,k,s)*uVunk_tmp(lma1,MB_0:MB_1,k,s),8)*rtmp5(1,m,s))*tmp1
              forceQ(3,a)=forceQ(3,a)+sum(occ(MB_0:MB_1,k,s)*real(uVunk_tmp(lma2,MB_0:MB_1,k,s)*uVunk_tmp(lma1,MB_0:MB_1,k,s),8)*rtmp5(2,m,s))*tmp1
            endif
#else
            if (lma1==lma2) then
              forceQ(1,a)=forceQ(1,a)+sum(occ(MB_0:MB_1,k,s)*real(conjg(uVunk_tmp(lma1,MB_0:MB_1,k,s))*uVunk_tmp(lma2,MB_0:MB_1,k,s),8)*rtmp5(0,m,s))*tmp1
              forceQ(2,a)=forceQ(2,a)+sum(occ(MB_0:MB_1,k,s)*real(conjg(uVunk_tmp(lma1,MB_0:MB_1,k,s))*uVunk_tmp(lma2,MB_0:MB_1,k,s),8)*rtmp5(1,m,s))*tmp1
              forceQ(3,a)=forceQ(3,a)+sum(occ(MB_0:MB_1,k,s)*real(conjg(uVunk_tmp(lma1,MB_0:MB_1,k,s))*uVunk_tmp(lma2,MB_0:MB_1,k,s),8)*rtmp5(2,m,s))*tmp1
            else
              forceQ(1,a)=forceQ(1,a)+sum(occ(MB_0:MB_1,k,s)*real(conjg(uVunk_tmp(lma1,MB_0:MB_1,k,s))*uVunk_tmp(lma2,MB_0:MB_1,k,s),8)*rtmp5(0,m,s))*tmp1
              forceQ(2,a)=forceQ(2,a)+sum(occ(MB_0:MB_1,k,s)*real(conjg(uVunk_tmp(lma1,MB_0:MB_1,k,s))*uVunk_tmp(lma2,MB_0:MB_1,k,s),8)*rtmp5(1,m,s))*tmp1
              forceQ(3,a)=forceQ(3,a)+sum(occ(MB_0:MB_1,k,s)*real(conjg(uVunk_tmp(lma1,MB_0:MB_1,k,s))*uVunk_tmp(lma2,MB_0:MB_1,k,s),8)*rtmp5(2,m,s))*tmp1
              forceQ(1,a)=forceQ(1,a)+sum(occ(MB_0:MB_1,k,s)*real(conjg(uVunk_tmp(lma2,MB_0:MB_1,k,s))*uVunk_tmp(lma1,MB_0:MB_1,k,s),8)*rtmp5(0,m,s))*tmp1
              forceQ(2,a)=forceQ(2,a)+sum(occ(MB_0:MB_1,k,s)*real(conjg(uVunk_tmp(lma2,MB_0:MB_1,k,s))*uVunk_tmp(lma1,MB_0:MB_1,k,s),8)*rtmp5(1,m,s))*tmp1
              forceQ(3,a)=forceQ(3,a)+sum(occ(MB_0:MB_1,k,s)*real(conjg(uVunk_tmp(lma2,MB_0:MB_1,k,s))*uVunk_tmp(lma1,MB_0:MB_1,k,s),8)*rtmp5(2,m,s))*tmp1
            endif
#endif
!if (myrank==0) write(240,'(4I5,7g20.7)') s,k,n,m,tmp1,rtmp5(0,m,n,k,s),occ(n,k,s),dV,uVunk_tmp(lma1,n,k,s)
!if (myrank==0) write(250,'(4I5,g20.7)') s,k,n,m,forceQ(1,a)
          endif
        end do ! m
      end do ! k
    end do ! s

    allocate( work2(3,MI) )
!$OMP end single

!$OMP workshare
    work2(:,:)=forceQ(:,:)
!$OMP end workshare

!$OMP single
    call mpi_allreduce(work2,forceQ,3*MI,mpi_real8,mpi_sum,mpi_comm_world,ierr)
    deallocate( work2 )
    deallocate( rtmp5 )
!$OMP end single

!$OMP end parallel

    return
  END SUBROUTINE calcForceQ
#else
  SUBROUTINE calcForceQ(uVunk_tmp,MI,forceQ)
    use bz_module
    use wf_module
    use watch_module
  use VarParaPSnonLocG
  use ps_nloc2_variables, only: lmap,amap,mmap,nzlma,iorbmap
    implicit none
    integer,intent(IN) :: MI
#ifdef _DRSDFT_
    real(8),intent(IN) :: uVunk_tmp(nzlma,MB_0:MB_1,MBZ_0:MBZ_1,MSP_0:MSP_1)
#else
    complex(8),intent(IN) :: uVunk_tmp(nzlma,MB_0:MB_1,MBZ_0:MBZ_1,MSP_0:MSP_1)
#endif
    real(8),intent(INOUT) :: forceQ(3,MI)
    return
  END SUBROUTINE calcForceQ
#endif

END MODULE ForcePSnonLoc2
