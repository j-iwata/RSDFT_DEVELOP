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
  implicit none
  PRIVATE
  PUBLIC :: calcForcePSnonLoc2
  real(8),allocatable :: SH_Y1(:,:,:,:),SH_Y2(:,:,:,:),SH_Y3(:,:,:,:)
  real(8),allocatable :: C_ijLM(:,:,:,:,:,:)
  logical :: isFirstSHY=.true.,isFirstC=.true.
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
    integer :: i1,i2,i3,lma1,lma2
    integer :: ib1,ib2
    integer :: i,j,k,s,n,ir,iorb,L,L1,L1z,NRc,irank,jrank
    integer :: nreq,max_nreq
    integer :: a,a0,ik,m,lm0,lm1,lma,im,m1,m2
    integer :: ierr,M_irad,ir0
    integer,allocatable :: ireq(:),istatus(:,:),irad(:,:),ilm1(:,:,:)
    real(8),parameter :: ep=1.d-8
    real(8) :: err,err0,maxerr,Rx,Ry,Rz
    real(8) :: a1,a2,a3,c1,c2,c3,d1,d2,d3
    real(8) :: x,y,z,r,kr,pi2,c
    real(8) :: tmp,tmp0,tmp1
    real(8) :: ctt(0:4),ett(0:4)
    real(8) :: yy1,yy2,yy3
    real(8),allocatable :: work2(:,:),duVdR(:,:,:)
    real(8) :: Dij_f,const_f
#ifdef _DRSDFT_
    real(8) :: ztmp
    real(8),allocatable :: wtmp5(:,:,:,:,:)
    real(8),allocatable :: uVunk_tmp(:,:,:,:)
#else
    complex(8) :: ztmp
    complex(8),allocatable :: wtmp5(:,:,:,:,:)
    complex(8),allocatable :: uVunk_tmp(:,:,:,:)
#endif
    logical,allocatable :: a_rank(:)
    integer :: ML1,ML2,ML3,i0,iorb0
    integer :: k1,k2,k3,a1b,a2b,a3b,ab1,ab2,ab3
    real(8) :: forceQ(3,MI)

    INTERFACE
      FUNCTION Ylm(x,y,z,l,m)
        real(8) :: Ylm
        real(8),intent(IN) :: x,y,z
        integer,intent(IN) :: l,m
      END FUNCTION Ylm
    END INTERFACE
    
    if (isFirstSHY) then
      call getSHY
      isFirstSHY=.false.
    endif
    if (isFirstC) then
      call getCijLM
      isFirstC=.false.
    endif

    force2(:,:) = 0.0d0

    if ( Mlma <= 0 ) return

    pi2 = 2.d0*acos(-1.d0)

    maxerr=0.d0
    ctt(:)=0.d0
    ett(:)=0.d0

    a1b=Igrid(1,1)
    a2b=Igrid(1,2)
    a3b=Igrid(1,3)
    ab1=Igrid(2,1)-Igrid(1,1)+1
    ab2=Igrid(2,2)-Igrid(1,2)+1
    ab3=Igrid(2,3)-Igrid(1,3)+1

    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

    c1=1.d0/ML1
    c2=1.d0/ML2
    c3=1.d0/ML3

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
    allocate( a_rank(MI) )
    allocate( duVdR(3,MMJJ,nzlma) )

!$OMP parallel

!$OMP workshare
    wtmp5=zero
    a_rank(:)=.false.
    duVdR=0.d0
!$OMP end workshare

!$OMP do private( i1,i2,i3,k1,k2,k3 )
    do a=1,MI
      i1 = nint( aa_atom(1,a)*ML1 )
      i2 = nint( aa_atom(2,a)*ML2 )
      i3 = nint( aa_atom(3,a)*ML3 )
      k1 = i1/ML1 ; if ( i1<0 ) k1=(i1+1)/ML1-1
      k2 = i2/ML2 ; if ( i2<0 ) k2=(i2+1)/ML2-1
      k3 = i3/ML3 ; if ( i3<0 ) k3=(i3+1)/ML3-1
      i1 = i1 - k1*ML1
      i2 = i2 - k2*ML2
      i3 = i3 - k3*ML3
      if ( Igrid(1,1) <= i1 .and. i1 <= Igrid(2,1) .and. &
          Igrid(1,2) <= i2 .and. i2 <= Igrid(2,2) .and. &
          Igrid(1,3) <= i3 .and. i3 <= Igrid(2,3) ) then
        a_rank(a)=.true.
      end if
    end do
!$OMP end do


#ifndef _SPLINE_
!$OMP single
    allocate( irad(0:3000,Nelement_) )
    irad=0
    M_irad=0
    do ik=1,Nelement_
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
!$OMP end single
#endif

!$OMP master
    call watch(ctt(0),ett(0))
!$OMP end master

!$OMP do schedule(dynamic) firstprivate( maxerr ) &
!$OMP    private( a,L,m,iorb,ik,Rx,Ry,Rz,NRc,d1,d2,d3,x,y,z,r  &
!$OMP            ,ir,ir0,yy1,yy2,yy3,err0,err,tmp0,tmp1,m1,m2  &
!$OMP            ,lma,j,L1,L1z,lm1,im )
    do lma=1,nzlma
      a    = amap(lma)
      if ( a <= 0 ) cycle
      L    = lmap(lma)
      m    = mmap(lma)
      iorb = iorbmap(lma)
      ik   = ki_atom(a)
      Rx=aa(1,1)*aa_atom(1,a)+aa(1,2)*aa_atom(2,a)+aa(1,3)*aa_atom(3,a)
      Ry=aa(2,1)*aa_atom(1,a)+aa(2,2)*aa_atom(2,a)+aa(2,3)*aa_atom(3,a)
      Rz=aa(3,1)*aa_atom(1,a)+aa(3,2)*aa_atom(2,a)+aa(3,3)*aa_atom(3,a)
      NRc=NRps(iorb,ik)
!!$OMP parallel do firstprivate( maxerr ) &
!!$OMP             private( d1,d2,d3,x,y,z,r,ir0,yy1,yy2,yy3,lm1,err,err0 &
!!$OMP                     ,tmp0,tmp1,m1,m2,j,L1,im,L1z )
      do j=1,MJJ_MAP(lma)
        d1=c1*JJ_MAP(1,j,lma)+JJ_MAP(4,j,lma)
        d2=c2*JJ_MAP(2,j,lma)+JJ_MAP(5,j,lma)
        d3=c3*JJ_MAP(3,j,lma)+JJ_MAP(6,j,lma)
        x = aa(1,1)*d1+aa(1,2)*d2+aa(1,3)*d3-Rx
        y = aa(2,1)*d1+aa(2,2)*d2+aa(2,3)*d3-Ry
        z = aa(3,1)*d1+aa(3,2)*d2+aa(3,3)*d3-Rz
        r = sqrt(x*x+y*y+z*z)
#ifndef _SPLINE_
        ir0=irad( int(100.d0*r),ik )
        do ir=ir0,NRc-1
          if ( r<rad1(ir,ik) ) exit
        end do
#endif
        yy1=0.d0
        yy2=0.d0
        yy3=0.d0
        do L1=abs(L-1),L+1
          lm1=ilm1(L1,iorb,ik)
          if ( abs(x)>1.d-14 .or. abs(y)>1.d-14 &
              .or. abs(z)>1.d-14 .or. L1==0 ) then
#ifdef _SPLINE_
            if ( r < rad1(2,ik) ) then
              tmp0=dviod(2,lm1,ik)/(rad1(2,ik)**2)
            else
              call splint(rad1(1,ik),dviod(1,lm1,ik),y2b,NRc,r,tmp0)
              tmp0=tmp0/(r*r)
            end if
#else
            if ( ir <= 2 ) then
              err0=0.d0
              tmp0=dviod(2,lm1,ik)/(rad1(2,ik)**2)
              if ( ir < 1 ) stop "calc_force_ps_nloc_uspp"
            else if ( ir <= NRc ) then
              err0=1.d10
              do im=1,20
                m1=max(1,ir-im)
                m2=min(ir+im,NRc)
                call polint(rad1(m1,ik),dviod(m1,lm1,ik) &
                ,m2-m1+1,r,tmp1,err)
                if ( abs(err)<err0 ) then
                  tmp0=tmp1
                  err0=abs(err)
                  if ( err0<ep ) exit
                end if
              end do
              tmp0=tmp0/(r*r)
            else
              write(*,*) "force_ps_nloc_uspp",ir,NRc
              stop
            end if
            maxerr=max(maxerr,err0)
#endif
            do L1z=-L1,L1
              tmp1=tmp0*Ylm(x,y,z,L1,L1z)
              yy1=yy1+tmp1*SH_Y1(L,m,L1,L1z)
              yy2=yy2+tmp1*SH_Y2(L,m,L1,L1z)
              yy3=yy3-tmp1*SH_Y3(L,m,L1,L1z)
            end do
          end if
        end do ! L1

        duVdR(1,j,lma)=yy1
        duVdR(2,j,lma)=yy2
        duVdR(3,j,lma)=yy3

      end do ! j
!!$OMP end parallel do
    end do ! lma
!$OMP end do

!$OMP master
    call watch(ctt(1),ett(1))
!$OMP end master

#ifndef _SPLINE_
!$OMP single
    deallocate( irad )
!$OMP end single
#endif


    do s=MSP_0,MSP_1
      do k=MBZ_0,MBZ_1
!$OMP do schedule(dynamic) private( c,i,d1,d2,d3,kr,ztmp,i1,i2,i3 )
        do n=MB_0,MB_1
          if ( occ(n,k,s) == 0.d0 ) cycle
          c=-2.d0*occ(n,k,s)*dV*dV
          do lma=1,nzlma
            if ( MJJ_MAP(lma) == MJJ(lma) ) then
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
                wtmp5(0,lma,n,k,s)=wtmp5(0,lma,n,k,s) &
                     +uVk(j,lma,k)*conjg(unk(i,n,k,s))
              end do
              wtmp5(0,lma,n,k,s)=c*wtmp5(0,lma,n,k,s)
              do j=1,MJJ_MAP(lma)
                i=JJP(j,lma)
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
                wtmp5(0,lma,n,k,s)=wtmp5(0,lma,n,k,s) &
                     +uVk(j,lma,k)*conjg(unk(i,n,k,s))
              end do
              wtmp5(0,lma,n,k,s)=c*wtmp5(0,lma,n,k,s)
              do j=1,MJJ_MAP(lma)
                i1=JJ_MAP(1,j,lma)
                i2=JJ_MAP(2,j,lma)
                i3=JJ_MAP(3,j,lma)
                i = i1-a1b + (i2-a2b)*ab1 + (i3-a3b)*ab1*ab2 + ML_0
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
    max_nreq=2*maxval( nrlma_xyz )
    allocate( ireq(max_nreq) )
    allocate( istatus(MPI_STATUS_SIZE,max_nreq) )
    allocate( uVunk_tmp(nzlma,MB_0:MB_1,MBZ_0:MBZ_1,MSP_0:MSP_1) )
    uVunk_tmp=zero
    do s=MSP_0,MSP_1
      do k=MBZ_0,MBZ_1
#ifdef _DRSDFT_
        uVunk_tmp(:,:,k,s)=wtmp5(0,:,:,k,s)*dV
#else
        uVunk_tmp(:,:,k,s)=conjg(wtmp5(0,:,:,k,s)*dV)
#endif
      enddo
    enddo
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

          if ( occ(n,k,s) == 0.d0 ) cycle
          call threeWayComm(nrlma_xyz,num_2_rank,sendmap,recvmap,lma_nsend,sbufnl,rbufnl,nzlma,ib1,ib2,wtmp5(0,1,ib1,k,s),3)
          
          do m=1,N_nzqr
            lma1=nzqr_pair(m,1)
            lma2=nzqr_pair(m,2)
            a=amap(lma1)
            if ( a <= 0 ) cycle
            if ( a_rank(a) ) then
              Dij_f=Dij(m,s)-esp(n,k,s)*qij_f(m)
              const_f=Dij_f*(-2.d0)*occ(n,k,s)*dV*dV
              if (lma1<lma2) stop 'Nzqr_pair is strange'
              if (lma1==lma2) then
                force2(1,a)=force2(1,a)+const_f*real(wtmp5(0,lma1,n,k,s)*wtmp5(1,lma2,n,k,s),8)
                force2(2,a)=force2(2,a)+const_f*real(wtmp5(0,lma1,n,k,s)*wtmp5(2,lma2,n,k,s),8)
                force2(3,a)=force2(3,a)+const_f*real(wtmp5(0,lma1,n,k,s)*wtmp5(3,lma2,n,k,s),8)
              else
                force2(1,a)=force2(1,a)+const_f*real(wtmp5(0,lma1,n,k,s)*wtmp5(1,lma2,n,k,s),8)
                force2(2,a)=force2(2,a)+const_f*real(wtmp5(0,lma1,n,k,s)*wtmp5(2,lma2,n,k,s),8)
                force2(3,a)=force2(3,a)+const_f*real(wtmp5(0,lma1,n,k,s)*wtmp5(3,lma2,n,k,s),8)
                force2(1,a)=force2(1,a)+const_f*real(wtmp5(0,lma2,n,k,s)*wtmp5(1,lma1,n,k,s),8)
                force2(2,a)=force2(2,a)+const_f*real(wtmp5(0,lma2,n,k,s)*wtmp5(2,lma1,n,k,s),8)
                force2(3,a)=force2(3,a)+const_f*real(wtmp5(0,lma2,n,k,s)*wtmp5(3,lma1,n,k,s),8)

              end if
            endif
          end do ! m

        end do ! n
      end do ! k
    end do ! s

    deallocate( istatus )
    deallocate( ireq )

    allocate( work2(3,MI) )
!$OMP end single

!$OMP workshare
    work2(:,:)=force2(:,:)
!$OMP end workshare

!$OMP single
    call mpi_allreduce(work2,force2,3*MI &
         ,mpi_real8,mpi_sum,mpi_comm_world,ierr)
    deallocate( work2 )
    deallocate( a_rank )
    deallocate( wtmp5 )
!$OMP end single

!$OMP master
    call watch(ctt(3),ett(3))
!$OMP end master

    call calcForceQ(uVunk_tmp,MI,forceQ)
do a=1,MI
!if (disp_switch_parallel) write(200,'(I5,A9,3g20.7)') a,'Q',forceQ(1:3,a)
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
    integer :: i1,i2,i3,lma1,lma2
    integer :: i,j,k,s,n,ir,iorb,L,L1,L1z,NRc,irank,jrank
    integer :: iorb1,iorb2
    integer :: iqr,ll3,cJ
    integer :: l_1,l_2,m_1,m_2
    integer :: ib1,ib2
    real(8) :: yq1,yq2,yq3,tmp2,tmp3
    integer :: nreq,max_nreq
    integer :: a,a0,ik,m,lm0,lm1,lma,im,m1,m2
    integer :: ierr,M_irad,ir0
    integer,allocatable :: ireq(:),istatus(:,:),irad(:,:)
    real(8),parameter :: ep=1.d-8
    real(8) :: err,err0,maxerr,Rx,Ry,Rz
    real(8) :: a1,a2,a3,c1,c2,c3,d1,d2,d3
    real(8) :: x,y,z,r,kr,pi2,c
    real(8) :: tmp,tmp0,tmp1
    real(8) :: yy1,yy2,yy3
    real(8),allocatable :: work2(:,:),dQY(:,:,:)
#ifdef _DRSDFT_
    real(8) :: ztmp
    real(8),allocatable :: rtmp5(:,:,:,:,:),rtmp2(:,:)
#else
    complex(8) :: ztmp
    complex(8),allocatable :: rtmp5(:,:,:,:,:),rtmp2(:,:)
#endif
    logical,allocatable :: a_rank(:)
    integer :: ML1,ML2,ML3,i0,iorb0
    integer :: k1,k2,k3,a1b,a2b,a3b,ab1,ab2,ab3

    INTERFACE
      FUNCTION Ylm(x,y,z,l,m)
        real(8) :: Ylm
        real(8),intent(IN) :: x,y,z
        integer,intent(IN) :: l,m
      END FUNCTION Ylm
    END INTERFACE
    

    forceQ(:,:) = 0.0d0

    pi2 = 2.d0*acos(-1.d0)

    maxerr=0.d0

    a1b=Igrid(1,1)
    a2b=Igrid(1,2)
    a3b=Igrid(1,3)
    ab1=Igrid(2,1)-Igrid(1,1)+1
    ab2=Igrid(2,2)-Igrid(1,2)+1
    ab3=Igrid(2,3)-Igrid(1,3)+1

    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

    c1=1.d0/ML1
    c2=1.d0/ML2
    c3=1.d0/ML3

    allocate( rtmp5(0:2,c_nzqr,MB_0:MB_1,MBZ_0:MBZ_1,MSP_0:MSP_1) )
    allocate( rtmp2(1:3,c_nzqr) )
    allocate( a_rank(MI) )
    allocate( dQY(3,MMJJ_Q,c_nzqr) )

!$OMP parallel

!$OMP workshare
    rtmp5=zero
    a_rank(:)=.false.
    dQY=0.d0
!$OMP end workshare

!$OMP do private( i1,i2,i3,k1,k2,k3 )
    do a=1,MI
      i1 = nint( aa_atom(1,a)*ML1 )
      i2 = nint( aa_atom(2,a)*ML2 )
      i3 = nint( aa_atom(3,a)*ML3 )
      k1 = i1/ML1 ; if ( i1<0 ) k1=(i1+1)/ML1-1
      k2 = i2/ML2 ; if ( i2<0 ) k2=(i2+1)/ML2-1
      k3 = i3/ML3 ; if ( i3<0 ) k3=(i3+1)/ML3-1
      i1 = i1 - k1*ML1
      i2 = i2 - k2*ML2
      i3 = i3 - k3*ML3
      if ( Igrid(1,1) <= i1 .and. i1 <= Igrid(2,1) .and. &
          Igrid(1,2) <= i2 .and. i2 <= Igrid(2,2) .and. &
          Igrid(1,3) <= i3 .and. i3 <= Igrid(2,3) ) then
        a_rank(a)=.true.
      end if
    end do
!$OMP end do


#ifndef _SPLINE_
!$OMP single
    allocate( irad(0:3000,Nelement_) )
    irad=0
    M_irad=0
    do ik=1,Nelement_
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
!$OMP end single
#endif


if (myrank==0) write(220,'(4a5,a20)') 'iqr','j','ll3','L1','tmp0'
!$OMP do schedule(dynamic) firstprivate( maxerr ) &
!$OMP    private( a,L,m,iorb,ik,Rx,Ry,Rz,NRc,d1,d2,d3,x,y,z,r  &
!$OMP            ,ir,ir0,yy1,yy2,yy3,err0,err,tmp0,tmp1,m1,m2  &
!$OMP            ,lma,j,L1,L1z,lm1,im )
    do iqr=1,c_nzqr
      a=amap_Q(iqr)
      if ( a <= 0 ) cycle
      ik=ki_atom(a)
      k1=k1map_Q(iqr)
      k2=k1_to_k2(k1,ik)
      lma1=lmamap_Q(iqr,1)
      lma2=lmamap_Q(iqr,2)
      l_1=lmap(lma1)
      m_1=mmap(lma1)
      l_2=lmap(lma2)
      m_2=mmap(lma2)
      iorb1=iorbmap(lma1)
      iorb2=iorbmap(lma2)
      Rx=aa(1,1)*aa_atom(1,a)+aa(1,2)*aa_atom(2,a)+aa(1,3)*aa_atom(3,a)
      Ry=aa(2,1)*aa_atom(1,a)+aa(2,2)*aa_atom(2,a)+aa(2,3)*aa_atom(3,a)
      Rz=aa(3,1)*aa_atom(1,a)+aa(3,2)*aa_atom(2,a)+aa(3,3)*aa_atom(3,a)
      NRc=max(NRps(iorb1,ik),NRps(iorb2,ik))
!!$OMP parallel do firstprivate( maxerr ) &
!!$OMP             private( d1,d2,d3,x,y,z,r,ir0,yy1,yy2,yy3,lm1,err,err0 &
!!$OMP                     ,tmp0,tmp1,m1,m2,j,L1,im,L1z )
      do j=1,MJJ_MAP_Q(iqr)
        d1=c1*JJ_MAP_Q(1,j,iqr)+JJ_MAP_Q(4,j,iqr)
        d2=c2*JJ_MAP_Q(2,j,iqr)+JJ_MAP_Q(5,j,iqr)
        d3=c3*JJ_MAP_Q(3,j,iqr)+JJ_MAP_Q(6,j,iqr)
        x = aa(1,1)*d1+aa(1,2)*d2+aa(1,3)*d3-Rx
        y = aa(2,1)*d1+aa(2,2)*d2+aa(2,3)*d3-Ry
        z = aa(3,1)*d1+aa(3,2)*d2+aa(3,3)*d3-Rz
        r = sqrt(x*x+y*y+z*z)
#ifndef _SPLINE_
        ir0=irad( int(100.d0*r),ik )
        do ir=ir0,NRc-1
          if ( r<rad1(ir,ik) ) exit
        end do
#endif
        do ll3=1,nl3v(k2,ik)
          L=l3v(ll3,k2,ik)-1
          cJ=0
          do L1=abs(L-1),L+1
            cJ=cJ+1
            tmp0=0.d0
            err0=0.d0
            if ( abs(x)>1.d-14 .or. abs(y)>1.d-14 &
              .or. abs(z)>1.d-14 .or. L1==0 ) then
#ifdef _SPLINE_
              if ( r < rad1(2,ik) ) then
                tmp0=dqrL(2,ll3,k2,ik,cJ)/(rad1(2,ik)**3)
              else
                call splint(rad1(1,ik),dqrL(1,ll3,k2,ik,cJ),y2b,NRc,r,tmp0)
                tmp0=tmp0/(r*r)
              end if
#else
              if ( ir <= 2 ) then
                err0=0.d0
                tmp0=dqrL(2,ll3,k2,ik,cJ)/(rad1(2,ik)**3)
                if ( ir < 1 ) stop "calc_force_ps_Q"
              else if ( ir <= NRc ) then
                err0=1.d10
                do im=1,20
                  m1=max(1,ir-im)
                  m2=min(ir+im,NRc)
                  call polint(rad1(m1,ik),dqrL(m1,ll3,k2,ik,cJ) &
                  ,m2-m1+1,r,tmp1,err)
                  if ( abs(err)<err0 ) then
                    tmp0=tmp1
                    err0=abs(err)
                    if ( err0<ep ) exit
                  end if
                end do
                tmp0=tmp0/(r*r)
              else
                write(*,*) "force_ps_Q",ir,NRc
                stop
              end if
if (myrank==0) write(220,'(4I5,g20.7)') iqr,j,ll3,L1,tmp0
              maxerr=max(maxerr,err0)
#endif
              yy1=0.d0
              yy2=0.d0
              yy3=0.d0
              do L1z=-L1,L1
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
            end if
          end do ! L1
          dQY(1,j,iqr)=dQY(1,j,iqr)+tmp0*yy1
          dQY(2,j,iqr)=dQY(2,j,iqr)+tmp0*yy2
          dQY(3,j,iqr)=dQY(3,j,iqr)+tmp0*yy3
        enddo ! ll3
      end do ! j
!!$OMP end parallel do
    end do ! iqr
!$OMP end do


#ifndef _SPLINE_
!$OMP single
    deallocate( irad )
!$OMP end single
#endif


if (myrank==0) write(230,'(5a5,a20)') 's','k','n','iqr','j','ztmp'
    do s=MSP_0,MSP_1
      do k=MBZ_0,MBZ_1
!$OMP do schedule(dynamic) private( c,i,d1,d2,d3,kr,ztmp,i1,i2,i3 )
        do n=MB_0,MB_1
          if ( occ(n,k,s) == 0.d0 ) cycle
          c=-2.d0*occ(n,k,s)*dV*dV
          do iqr=1,c_nzqr
            do j=1,MJJ_MAP_Q(iqr)
              i1=JJ_MAP_Q(1,j,iqr)
              i2=JJ_MAP_Q(2,j,iqr)
              i3=JJ_MAP_Q(3,j,iqr)
              i = i1-a1b + (i2-a2b)*ab1 + (i3-a3b)*ab1*ab2 + ML_0
              ztmp=Vloc(i,s)
              rtmp5(0,iqr,n,k,s)=rtmp5(0,iqr,n,k,s)+dQY(1,j,iqr)*ztmp
              rtmp5(1,iqr,n,k,s)=rtmp5(1,iqr,n,k,s)+dQY(2,j,iqr)*ztmp
              rtmp5(2,iqr,n,k,s)=rtmp5(2,iqr,n,k,s)-dQY(3,j,iqr)*ztmp
if (myrank==0) write(230,'(5I5,g20.7)') s,k,n,iqr,j,Vloc(i,s)
            end do
          end do ! lma
        end do ! n
!$OMP end do
      end do ! k
    end do ! s


!$OMP single
    deallocate( dQY )
!$OMP end single


!$OMP single
if (myrank==0) write(240,'(4a5,7a20)') 's','k','n','m','tmp1','','rtmp5(0,m,n,k,s)','occ(n,k,s)','dV','','uVunk_tmp(lma1,n,k,s)'
if (myrank==0) write(250,'(4a5,a20)') 's','k','n','m','forceQ(1,a)'
    do s=MSP_0,MSP_1
      do k=MBZ_0,MBZ_1
        do n=MB_0,MB_1,MB_d
          ib1=n
          ib2=min(ib1+MB_d-1,MB_1)

          if ( occ(n,k,s) == 0.d0 ) cycle
          call threeWayComm(nrqr_xyz,num_2_rank_Q,sendmap_Q,recvmap_Q,qr_nsend,sbufnl_Q,rbufnl_Q,c_nzqr,ib1,ib2,rtmp5(0,1,ib1,k,s),2)
          
          do m=1,N_nzqr
            lma1=nzqr_pair(m,1)
            lma2=nzqr_pair(m,2)
            a=amap(lma1)
            if ( a <= 0 ) cycle
            if ( a_rank(a) ) then
              if (lma1<lma2) stop 'Nzqr_pair is strange'
              tmp1=-rtmp5(0,m,n,k,s)*occ(n,k,s)*dV
              tmp2=-rtmp5(1,m,n,k,s)*occ(n,k,s)*dV
              tmp3=-rtmp5(2,m,n,k,s)*occ(n,k,s)*dV
#ifdef _DRSDFT_
              if (lma1==lma2) then
                forceQ(1,a)=forceQ(1,a)+real(uVunk_tmp(lma1,n,k,s)*uVunk_tmp(lma2,n,k,s),8)*tmp1
                forceQ(2,a)=forceQ(2,a)+real(uVunk_tmp(lma1,n,k,s)*uVunk_tmp(lma2,n,k,s),8)*tmp2
                forceQ(3,a)=forceQ(3,a)+real(uVunk_tmp(lma1,n,k,s)*uVunk_tmp(lma2,n,k,s),8)*tmp3
              else
                forceQ(1,a)=forceQ(1,a)+real(uVunk_tmp(lma1,n,k,s)*uVunk_tmp(lma2,n,k,s),8)*tmp1
                forceQ(2,a)=forceQ(2,a)+real(uVunk_tmp(lma1,n,k,s)*uVunk_tmp(lma2,n,k,s),8)*tmp2
                forceQ(3,a)=forceQ(3,a)+real(uVunk_tmp(lma1,n,k,s)*uVunk_tmp(lma2,n,k,s),8)*tmp3
                forceQ(1,a)=forceQ(1,a)+real(uVunk_tmp(lma2,n,k,s)*uVunk_tmp(lma1,n,k,s),8)*tmp1
                forceQ(2,a)=forceQ(2,a)+real(uVunk_tmp(lma2,n,k,s)*uVunk_tmp(lma1,n,k,s),8)*tmp2
                forceQ(3,a)=forceQ(3,a)+real(uVunk_tmp(lma2,n,k,s)*uVunk_tmp(lma1,n,k,s),8)*tmp3
              endif
#else
              if (lma1==lma2) then
                forceQ(1,a)=forceQ(1,a)+real(conjg(uVunk_tmp(lma1,n,k,s))*uVunk_tmp(lma2,n,k,s),8)*tmp1
                forceQ(2,a)=forceQ(2,a)+real(conjg(uVunk_tmp(lma1,n,k,s))*uVunk_tmp(lma2,n,k,s),8)*tmp2
                forceQ(3,a)=forceQ(3,a)+real(conjg(uVunk_tmp(lma1,n,k,s))*uVunk_tmp(lma2,n,k,s),8)*tmp3
              else
                forceQ(1,a)=forceQ(1,a)+real(conjg(uVunk_tmp(lma1,n,k,s))*uVunk_tmp(lma2,n,k,s),8)*tmp1
                forceQ(2,a)=forceQ(2,a)+real(conjg(uVunk_tmp(lma1,n,k,s))*uVunk_tmp(lma2,n,k,s),8)*tmp2
                forceQ(3,a)=forceQ(3,a)+real(conjg(uVunk_tmp(lma1,n,k,s))*uVunk_tmp(lma2,n,k,s),8)*tmp3
                forceQ(1,a)=forceQ(1,a)+real(conjg(uVunk_tmp(lma2,n,k,s))*uVunk_tmp(lma1,n,k,s),8)*tmp1
                forceQ(2,a)=forceQ(2,a)+real(conjg(uVunk_tmp(lma2,n,k,s))*uVunk_tmp(lma1,n,k,s),8)*tmp2
                forceQ(3,a)=forceQ(3,a)+real(conjg(uVunk_tmp(lma2,n,k,s))*uVunk_tmp(lma1,n,k,s),8)*tmp3
              endif
#endif
if (myrank==0) write(240,'(4I5,7g20.7)') s,k,n,m,tmp1,rtmp5(0,m,n,k,s),occ(n,k,s),dV,uVunk_tmp(lma1,n,k,s)
if (myrank==0) write(250,'(4I5,g20.7)') s,k,n,m,forceQ(1,a)
            endif
          end do ! m
        end do ! n
      end do ! k
    end do ! s

    allocate( work2(3,MI) )
!$OMP end single

!$OMP workshare
    work2(:,:)=forceQ(:,:)
!$OMP end workshare

!$OMP single
    call mpi_allreduce(work2,forceQ,3*MI &
         ,mpi_real8,mpi_sum,mpi_comm_world,ierr)
    deallocate( work2 )
    deallocate( a_rank )
    deallocate( rtmp5 )
!$OMP end single

!$OMP end parallel
do a=1,MI
if (myrank==0) write(200,'(I5,A9,3g20.7)') a,'Q',forceQ(1:3,a)
enddo

    return
  END SUBROUTINE calcForceQ

  SUBROUTINE getSHY
    implicit none
    integer :: L,Lp1
    L=2*max_Lref
    Lp1=L+1
    allocate(SH_Y1(0:L,-L:L,-Lp1:Lp1,-Lp1:Lp1)) ; SH_Y1=0.d0
    allocate(SH_Y2(0:L,-L:L,-Lp1:Lp1,-Lp1:Lp1)) ; SH_Y2=0.d0
    allocate(SH_Y3(0:L,-L:L,-Lp1:Lp1,-Lp1:Lp1)) ; SH_Y3=0.d0

    SH_Y1(  0,  0,  1,  1)=   0.282094791773878d0

    SH_Y2(  0,  0,  1, -1)=   0.282094791773878d0

    SH_Y3(  0,  0,  1,  0)=   0.282094791773878d0

    if (L>=2) then
      SH_Y1(  1,  1,  0,  0)=   0.282094791773878d0
      SH_Y1(  1, -1,  2, -2)=  -0.218509686118416d0
      SH_Y1(  1,  1,  2,  0)=  -0.126156626101008d0
      SH_Y1(  1,  0,  2,  1)=   0.218509686118416d0
      SH_Y1(  1,  1,  2,  2)=   0.218509686118416d0
      SH_Y1(  2, -2,  1, -1)=  -0.218509686118416d0
      SH_Y1(  2,  1,  1,  0)=   0.218509686118416d0
      SH_Y1(  2,  0,  1,  1)=  -0.126156626101008d0
      SH_Y1(  2,  2,  1,  1)=   0.218509686118416d0
      SH_Y1(  2, -2,  3, -3)=  -0.226179013159540d0
      SH_Y1(  2, -1,  3, -2)=  -0.184674390922372d0
      SH_Y1(  2, -2,  3, -1)=   0.058399170081902d0
      SH_Y1(  2,  1,  3,  0)=  -0.143048168102669d0
      SH_Y1(  2,  0,  3,  1)=   0.202300659403421d0
      SH_Y1(  2,  2,  3,  1)=  -0.058399170081902d0
      SH_Y1(  2,  1,  3,  2)=   0.184674390922372d0
      SH_Y1(  2,  2,  3,  3)=   0.226179013159540d0

      SH_Y2(  1, -1,  0,  0)=   0.282094791773878d0
      SH_Y2(  1,  1,  2, -2)=  -0.218509686118416d0
      SH_Y2(  1,  0,  2, -1)=   0.218509686118416d0
      SH_Y2(  1, -1,  2,  0)=  -0.126156626101008d0
      SH_Y2(  1, -1,  2,  2)=  -0.218509686118416d0
      SH_Y2(  2,  0,  1, -1)=  -0.126156626101008d0
      SH_Y2(  2,  2,  1, -1)=  -0.218509686118416d0
      SH_Y2(  2, -1,  1,  0)=   0.218509686118416d0
      SH_Y2(  2, -2,  1,  1)=  -0.218509686118416d0
      SH_Y2(  2,  2,  3, -3)=   0.226179013159540d0
      SH_Y2(  2,  1,  3, -2)=  -0.184674390922372d0
      SH_Y2(  2,  0,  3, -1)=   0.202300659403421d0
      SH_Y2(  2,  2,  3, -1)=   0.058399170081902d0
      SH_Y2(  2, -1,  3,  0)=  -0.143048168102669d0
      SH_Y2(  2, -2,  3,  1)=   0.058399170081902d0
      SH_Y2(  2, -1,  3,  2)=  -0.184674390922372d0
      SH_Y2(  2, -2,  3,  3)=   0.226179013159540d0

      SH_Y3(  1,  0,  0,  0)=   0.282094791773878d0
      SH_Y3(  1, -1,  2, -1)=   0.218509686118416d0
      SH_Y3(  1,  0,  2,  0)=   0.252313252202016d0
      SH_Y3(  1,  1,  2,  1)=   0.218509686118416d0
      SH_Y3(  2, -1,  1, -1)=   0.218509686118416d0
      SH_Y3(  2,  0,  1,  0)=   0.252313252202016d0
      SH_Y3(  2,  1,  1,  1)=   0.218509686118416d0
      SH_Y3(  2, -2,  3, -2)=   0.184674390922372d0
      SH_Y3(  2, -1,  3, -1)=   0.233596680327607d0
      SH_Y3(  2,  0,  3,  0)=   0.247766695083476d0
      SH_Y3(  2,  1,  3,  1)=   0.233596680327607d0
      SH_Y3(  2,  2,  3,  2)=   0.184674390922372d0
    endif

    if (L>=4) then
      SH_Y1(  3, -3,  2, -2)=  -0.226179013159540d0
      SH_Y1(  3, -1,  2, -2)=   0.058399170081902d0
      SH_Y1(  3, -2,  2, -1)=  -0.184674390922372d0
      SH_Y1(  3,  1,  2,  0)=   0.202300659403421d0
      SH_Y1(  3,  0,  2,  1)=  -0.143048168102669d0
      SH_Y1(  3,  2,  2,  1)=   0.184674390922372d0
      SH_Y1(  3,  1,  2,  2)=  -0.058399170081902d0
      SH_Y1(  3,  3,  2,  2)=   0.226179013159540d0
      SH_Y1(  3, -3,  4, -4)=  -0.230329432980890d0
      SH_Y1(  3, -2,  4, -3)=  -0.199471140200716d0
      SH_Y1(  3, -3,  4, -2)=   0.043528171377568d0
      SH_Y1(  3, -1,  4, -2)=  -0.168583882836184d0
      SH_Y1(  3, -2,  4, -1)=   0.075393004386513d0
      SH_Y1(  3,  1,  4,  0)=  -0.150786008773027d0
      SH_Y1(  3,  0,  4,  1)=   0.194663900273006d0
      SH_Y1(  3,  2,  4,  1)=  -0.075393004386513d0
      SH_Y1(  3,  1,  4,  2)=   0.168583882836184d0
      SH_Y1(  3,  3,  4,  2)=  -0.043528171377568d0
      SH_Y1(  3,  2,  4,  3)=   0.199471140200716d0
      SH_Y1(  3,  3,  4,  4)=   0.230329432980890d0
      SH_Y1(  4, -4,  3, -3)=  -0.230329432980890d0
      SH_Y1(  4, -2,  3, -3)=   0.043528171377568d0
      SH_Y1(  4, -3,  3, -2)=  -0.199471140200716d0
      SH_Y1(  4, -1,  3, -2)=   0.075393004386513d0
      SH_Y1(  4, -2,  3, -1)=  -0.168583882836184d0
      SH_Y1(  4,  1,  3,  0)=   0.194663900273006d0
      SH_Y1(  4,  0,  3,  1)=  -0.150786008773027d0
      SH_Y1(  4,  2,  3,  1)=   0.168583882836184d0
      SH_Y1(  4,  1,  3,  2)=  -0.075393004386513d0
      SH_Y1(  4,  3,  3,  2)=   0.199471140200716d0
      SH_Y1(  4,  2,  3,  3)=  -0.043528171377568d0
      SH_Y1(  4,  4,  3,  3)=   0.230329432980890d0
      SH_Y1(  4, -4,  5, -5)=  -0.232932108055429d0
      SH_Y1(  4, -3,  5, -4)=  -0.208340811101706d0
      SH_Y1(  4, -4,  5, -3)=   0.034723468516951d0
      SH_Y1(  4, -2,  5, -3)=  -0.183739324706867d0
      SH_Y1(  4, -3,  5, -2)=   0.060142811686378d0
      SH_Y1(  4, -1,  5, -2)=  -0.159122922870344d0
      SH_Y1(  4, -2,  5, -1)=   0.085054779966126d0
      SH_Y1(  4,  1,  5,  0)=  -0.155288072036953d0
      SH_Y1(  4,  0,  5,  1)=   0.190188269815546d0
      SH_Y1(  4,  2,  5,  1)=  -0.085054779966126d0
      SH_Y1(  4,  1,  5,  2)=   0.159122922870344d0
      SH_Y1(  4,  3,  5,  2)=  -0.060142811686378d0
      SH_Y1(  4,  2,  5,  3)=   0.183739324706867d0
      SH_Y1(  4,  4,  5,  3)=  -0.034723468516951d0
      SH_Y1(  4,  3,  5,  4)=   0.208340811101706d0
      SH_Y1(  4,  4,  5,  5)=   0.232932108055429d0

      SH_Y2(  3,  1,  2, -2)=   0.058399170081902d0
      SH_Y2(  3,  3,  2, -2)=   0.226179013159540d0
      SH_Y2(  3,  0,  2, -1)=  -0.143048168102669d0
      SH_Y2(  3,  2,  2, -1)=  -0.184674390922372d0
      SH_Y2(  3, -1,  2,  0)=   0.202300659403421d0
      SH_Y2(  3, -2,  2,  1)=  -0.184674390922372d0
      SH_Y2(  3, -3,  2,  2)=   0.226179013159540d0
      SH_Y2(  3, -1,  2,  2)=   0.058399170081902d0
      SH_Y2(  3,  3,  4, -4)=  -0.230329432980890d0
      SH_Y2(  3,  2,  4, -3)=   0.199471140200716d0
      SH_Y2(  3,  1,  4, -2)=  -0.168583882836184d0
      SH_Y2(  3,  3,  4, -2)=  -0.043528171377568d0
      SH_Y2(  3,  0,  4, -1)=   0.194663900273006d0
      SH_Y2(  3,  2,  4, -1)=   0.075393004386513d0
      SH_Y2(  3, -1,  4,  0)=  -0.150786008773027d0
      SH_Y2(  3, -2,  4,  1)=   0.075393004386513d0
      SH_Y2(  3, -3,  4,  2)=  -0.043528171377568d0
      SH_Y2(  3, -1,  4,  2)=  -0.168583882836184d0
      SH_Y2(  3, -2,  4,  3)=   0.199471140200716d0
      SH_Y2(  3, -3,  4,  4)=  -0.230329432980890d0
      SH_Y2(  4,  2,  3, -3)=  -0.043528171377568d0
      SH_Y2(  4,  4,  3, -3)=  -0.230329432980890d0
      SH_Y2(  4,  1,  3, -2)=   0.075393004386513d0
      SH_Y2(  4,  3,  3, -2)=   0.199471140200716d0
      SH_Y2(  4,  0,  3, -1)=  -0.150786008773027d0
      SH_Y2(  4,  2,  3, -1)=  -0.168583882836184d0
      SH_Y2(  4, -1,  3,  0)=   0.194663900273006d0
      SH_Y2(  4, -2,  3,  1)=  -0.168583882836184d0
      SH_Y2(  4, -3,  3,  2)=   0.199471140200716d0
      SH_Y2(  4, -1,  3,  2)=   0.075393004386513d0
      SH_Y2(  4, -4,  3,  3)=  -0.230329432980890d0
      SH_Y2(  4, -2,  3,  3)=  -0.043528171377568d0
      SH_Y2(  4,  4,  5, -5)=   0.232932108055429d0
      SH_Y2(  4,  3,  5, -4)=  -0.208340811101706d0
      SH_Y2(  4,  2,  5, -3)=   0.183739324706867d0
      SH_Y2(  4,  4,  5, -3)=   0.034723468516951d0
      SH_Y2(  4,  1,  5, -2)=  -0.159122922870344d0
      SH_Y2(  4,  3,  5, -2)=  -0.060142811686378d0
      SH_Y2(  4,  0,  5, -1)=   0.190188269815546d0
      SH_Y2(  4,  2,  5, -1)=   0.085054779966126d0
      SH_Y2(  4, -1,  5,  0)=  -0.155288072036953d0
      SH_Y2(  4, -2,  5,  1)=   0.085054779966126d0
      SH_Y2(  4, -3,  5,  2)=  -0.060142811686378d0
      SH_Y2(  4, -1,  5,  2)=  -0.159122922870344d0
      SH_Y2(  4, -4,  5,  3)=   0.034723468516951d0
      SH_Y2(  4, -2,  5,  3)=   0.183739324706867d0
      SH_Y2(  4, -3,  5,  4)=  -0.208340811101706d0
      SH_Y2(  4, -4,  5,  5)=   0.232932108055429d0

      SH_Y3(  3, -2,  2, -2)=   0.184674390922372d0
      SH_Y3(  3, -1,  2, -1)=   0.233596680327607d0
      SH_Y3(  3,  0,  2,  0)=   0.247766695083476d0
      SH_Y3(  3,  1,  2,  1)=   0.233596680327607d0
      SH_Y3(  3,  2,  2,  2)=   0.184674390922372d0
      SH_Y3(  3, -3,  4, -3)=   0.162867503967640d0
      SH_Y3(  3, -2,  4, -2)=   0.213243618622923d0
      SH_Y3(  3, -1,  4, -1)=   0.238413613504448d0
      SH_Y3(  3,  0,  4,  0)=   0.246232521229829d0
      SH_Y3(  3,  1,  4,  1)=   0.238413613504448d0
      SH_Y3(  3,  2,  4,  2)=   0.213243618622923d0
      SH_Y3(  3,  3,  4,  3)=   0.162867503967640d0
      SH_Y3(  4, -3,  3, -3)=   0.162867503967640d0
      SH_Y3(  4, -2,  3, -2)=   0.213243618622923d0
      SH_Y3(  4, -1,  3, -1)=   0.238413613504448d0
      SH_Y3(  4,  0,  3,  0)=   0.246232521229829d0
      SH_Y3(  4,  1,  3,  1)=   0.238413613504448d0
      SH_Y3(  4,  2,  3,  2)=   0.213243618622923d0
      SH_Y3(  4,  3,  3,  3)=   0.162867503967640d0
      SH_Y3(  4, -4,  5, -4)=   0.147319200327922d0
      SH_Y3(  4, -3,  5, -3)=   0.196425600437230d0
      SH_Y3(  4, -2,  5, -2)=   0.225033795607689d0
      SH_Y3(  4, -1,  5, -1)=   0.240571246745510d0
      SH_Y3(  4,  0,  5,  0)=   0.245532000546537d0
      SH_Y3(  4,  1,  5,  1)=   0.240571246745510d0
      SH_Y3(  4,  2,  5,  2)=   0.225033795607689d0
      SH_Y3(  4,  3,  5,  3)=   0.196425600437230d0
      SH_Y3(  4,  4,  5,  4)=   0.147319200327922d0
    endif

    return
  END SUBROUTINE getSHY


  SUBROUTINE getCijLM
    implicit none
    allocate(C_ijLM(0:4,-4:4,0:2,-2:2,0:2,-2:2))
    C_ijLM=0.d0

    C_ijLM(  0,  0,  0,  0,  0,  0  )=    0.282094791773878d0
    C_ijLM(  1, -1,  1, -1,  0,  0  )=    0.282094791773878d0
    C_ijLM(  1,  0,  1,  0,  0,  0  )=    0.282094791773878d0
    C_ijLM(  1,  1,  1,  1,  0,  0  )=    0.282094791773878d0
    C_ijLM(  2, -2,  2, -2,  0,  0  )=    0.282094791773878d0
    C_ijLM(  2, -1,  2, -1,  0,  0  )=    0.282094791773878d0
    C_ijLM(  2,  0,  2,  0,  0,  0  )=    0.282094791773878d0
    C_ijLM(  2,  1,  2,  1,  0,  0  )=    0.282094791773878d0
    C_ijLM(  2,  2,  2,  2,  0,  0  )=    0.282094791773878d0
    C_ijLM(  1, -1,  0,  0,  1, -1  )=    0.282094791773878d0
    C_ijLM(  0,  0,  1, -1,  1, -1  )=    0.282094791773878d0
    C_ijLM(  2, -2,  1,  1,  1, -1  )=   -0.218509686118416d0
    C_ijLM(  2, -1,  1,  0,  1, -1  )=    0.218509686118416d0
    C_ijLM(  2,  0,  1, -1,  1, -1  )=   -0.126156626101008d0
    C_ijLM(  2,  2,  1, -1,  1, -1  )=   -0.218509686118416d0
    C_ijLM(  1, -1,  2,  0,  1, -1  )=   -0.126156626101008d0
    C_ijLM(  1, -1,  2,  2,  1, -1  )=   -0.218509686118416d0
    C_ijLM(  1,  0,  2, -1,  1, -1  )=    0.218509686118416d0
    C_ijLM(  1,  1,  2, -2,  1, -1  )=   -0.218509686118416d0
    C_ijLM(  3, -3,  2,  2,  1, -1  )=    0.226179013159540d0
    C_ijLM(  3, -2,  2,  1,  1, -1  )=   -0.184674390922372d0
    C_ijLM(  3, -1,  2,  0,  1, -1  )=    0.202300659403421d0
    C_ijLM(  3, -1,  2,  2,  1, -1  )=    0.058399170081902d0
    C_ijLM(  3,  0,  2, -1,  1, -1  )=   -0.143048168102669d0
    C_ijLM(  3,  1,  2, -2,  1, -1  )=    0.058399170081902d0
    C_ijLM(  3,  2,  2, -1,  1, -1  )=   -0.184674390922372d0
    C_ijLM(  3,  3,  2, -2,  1, -1  )=    0.226179013159540d0
    C_ijLM(  1,  0,  0,  0,  1,  0  )=    0.282094791773878d0
    C_ijLM(  0,  0,  1,  0,  1,  0  )=    0.282094791773878d0
    C_ijLM(  2, -1,  1, -1,  1,  0  )=    0.218509686118416d0
    C_ijLM(  2,  0,  1,  0,  1,  0  )=    0.252313252202016d0
    C_ijLM(  2,  1,  1,  1,  1,  0  )=    0.218509686118416d0
    C_ijLM(  1, -1,  2, -1,  1,  0  )=    0.218509686118416d0
    C_ijLM(  1,  0,  2,  0,  1,  0  )=    0.252313252202016d0
    C_ijLM(  1,  1,  2,  1,  1,  0  )=    0.218509686118416d0
    C_ijLM(  3, -2,  2, -2,  1,  0  )=    0.184674390922372d0
    C_ijLM(  3, -1,  2, -1,  1,  0  )=    0.233596680327607d0
    C_ijLM(  3,  0,  2,  0,  1,  0  )=    0.247766695083476d0
    C_ijLM(  3,  1,  2,  1,  1,  0  )=    0.233596680327607d0
    C_ijLM(  3,  2,  2,  2,  1,  0  )=    0.184674390922372d0
    C_ijLM(  1,  1,  0,  0,  1,  1  )=    0.282094791773878d0
    C_ijLM(  0,  0,  1,  1,  1,  1  )=    0.282094791773878d0
    C_ijLM(  2, -2,  1, -1,  1,  1  )=   -0.218509686118416d0
    C_ijLM(  2,  0,  1,  1,  1,  1  )=   -0.126156626101008d0
    C_ijLM(  2,  1,  1,  0,  1,  1  )=    0.218509686118416d0
    C_ijLM(  2,  2,  1,  1,  1,  1  )=    0.218509686118416d0
    C_ijLM(  1, -1,  2, -2,  1,  1  )=   -0.218509686118416d0
    C_ijLM(  1,  0,  2,  1,  1,  1  )=    0.218509686118416d0
    C_ijLM(  1,  1,  2,  0,  1,  1  )=   -0.126156626101008d0
    C_ijLM(  1,  1,  2,  2,  1,  1  )=    0.218509686118416d0
    C_ijLM(  3, -3,  2, -2,  1,  1  )=   -0.226179013159540d0
    C_ijLM(  3, -2,  2, -1,  1,  1  )=   -0.184674390922372d0
    C_ijLM(  3, -1,  2, -2,  1,  1  )=    0.058399170081902d0
    C_ijLM(  3,  0,  2,  1,  1,  1  )=   -0.143048168102669d0
    C_ijLM(  3,  1,  2,  0,  1,  1  )=    0.202300659403421d0
    C_ijLM(  3,  1,  2,  2,  1,  1  )=   -0.058399170081902d0
    C_ijLM(  3,  2,  2,  1,  1,  1  )=    0.184674390922372d0
    C_ijLM(  3,  3,  2,  2,  1,  1  )=    0.226179013159540d0
    C_ijLM(  2, -2,  0,  0,  2, -2  )=    0.282094791773878d0
    C_ijLM(  1, -1,  1,  1,  2, -2  )=   -0.218509686118416d0
    C_ijLM(  1,  1,  1, -1,  2, -2  )=   -0.218509686118416d0
    C_ijLM(  3, -3,  1,  1,  2, -2  )=   -0.226179013159540d0
    C_ijLM(  3, -2,  1,  0,  2, -2  )=    0.184674390922372d0
    C_ijLM(  3, -1,  1,  1,  2, -2  )=    0.058399170081902d0
    C_ijLM(  3,  1,  1, -1,  2, -2  )=    0.058399170081902d0
    C_ijLM(  3,  3,  1, -1,  2, -2  )=    0.226179013159540d0
    C_ijLM(  0,  0,  2, -2,  2, -2  )=    0.282094791773878d0
    C_ijLM(  2, -2,  2,  0,  2, -2  )=   -0.180223751572869d0
    C_ijLM(  2, -1,  2,  1,  2, -2  )=   -0.156078347227440d0
    C_ijLM(  2,  0,  2, -2,  2, -2  )=   -0.180223751572869d0
    C_ijLM(  2,  1,  2, -1,  2, -2  )=   -0.156078347227440d0
    C_ijLM(  4, -4,  2,  2,  2, -2  )=    0.238413613504448d0
    C_ijLM(  4, -3,  2,  1,  2, -2  )=   -0.168583882836184d0
    C_ijLM(  4, -2,  2,  0,  2, -2  )=    0.156078347227440d0
    C_ijLM(  4, -1,  2,  1,  2, -2  )=    0.063718718434028d0
    C_ijLM(  4,  0,  2, -2,  2, -2  )=    0.040299255967697d0
    C_ijLM(  4,  1,  2, -1,  2, -2  )=    0.063718718434028d0
    C_ijLM(  4,  3,  2, -1,  2, -2  )=    0.168583882836184d0
    C_ijLM(  4,  4,  2, -2,  2, -2  )=   -0.238413613504448d0
    C_ijLM(  2, -1,  0,  0,  2, -1  )=    0.282094791773878d0
    C_ijLM(  1, -1,  1,  0,  2, -1  )=    0.218509686118416d0
    C_ijLM(  1,  0,  1, -1,  2, -1  )=    0.218509686118416d0
    C_ijLM(  3, -2,  1,  1,  2, -1  )=   -0.184674390922372d0
    C_ijLM(  3, -1,  1,  0,  2, -1  )=    0.233596680327607d0
    C_ijLM(  3,  0,  1, -1,  2, -1  )=   -0.143048168102669d0
    C_ijLM(  3,  2,  1, -1,  2, -1  )=   -0.184674390922372d0
    C_ijLM(  0,  0,  2, -1,  2, -1  )=    0.282094791773878d0
    C_ijLM(  2, -2,  2,  1,  2, -1  )=   -0.156078347227440d0
    C_ijLM(  2, -1,  2,  0,  2, -1  )=    0.090111875786434d0
    C_ijLM(  2, -1,  2,  2,  2, -1  )=   -0.156078347227440d0
    C_ijLM(  2,  0,  2, -1,  2, -1  )=    0.090111875786434d0
    C_ijLM(  2,  1,  2, -2,  2, -1  )=   -0.156078347227440d0
    C_ijLM(  2,  2,  2, -1,  2, -1  )=   -0.156078347227440d0
    C_ijLM(  4, -3,  2,  2,  2, -1  )=    0.168583882836184d0
    C_ijLM(  4, -2,  2,  1,  2, -1  )=   -0.180223751572869d0
    C_ijLM(  4, -1,  2,  0,  2, -1  )=    0.220728115441823d0
    C_ijLM(  4, -1,  2,  2,  2, -1  )=    0.063718718434028d0
    C_ijLM(  4,  0,  2, -1,  2, -1  )=   -0.161197023870787d0
    C_ijLM(  4,  1,  2, -2,  2, -1  )=    0.063718718434028d0
    C_ijLM(  4,  2,  2, -1,  2, -1  )=   -0.180223751572869d0
    C_ijLM(  4,  3,  2, -2,  2, -1  )=    0.168583882836184d0
    C_ijLM(  2,  0,  0,  0,  2,  0  )=    0.282094791773878d0
    C_ijLM(  1, -1,  1, -1,  2,  0  )=   -0.126156626101008d0
    C_ijLM(  1,  0,  1,  0,  2,  0  )=    0.252313252202016d0
    C_ijLM(  1,  1,  1,  1,  2,  0  )=   -0.126156626101008d0
    C_ijLM(  3, -1,  1, -1,  2,  0  )=    0.202300659403421d0
    C_ijLM(  3,  0,  1,  0,  2,  0  )=    0.247766695083476d0
    C_ijLM(  3,  1,  1,  1,  2,  0  )=    0.202300659403421d0
    C_ijLM(  0,  0,  2,  0,  2,  0  )=    0.282094791773878d0
    C_ijLM(  2, -2,  2, -2,  2,  0  )=   -0.180223751572869d0
    C_ijLM(  2, -1,  2, -1,  2,  0  )=    0.090111875786434d0
    C_ijLM(  2,  0,  2,  0,  2,  0  )=    0.180223751572869d0
    C_ijLM(  2,  1,  2,  1,  2,  0  )=    0.090111875786434d0
    C_ijLM(  2,  2,  2,  2,  2,  0  )=   -0.180223751572869d0
    C_ijLM(  4, -2,  2, -2,  2,  0  )=    0.156078347227440d0
    C_ijLM(  4, -1,  2, -1,  2,  0  )=    0.220728115441823d0
    C_ijLM(  4,  0,  2,  0,  2,  0  )=    0.241795535806181d0
    C_ijLM(  4,  1,  2,  1,  2,  0  )=    0.220728115441823d0
    C_ijLM(  4,  2,  2,  2,  2,  0  )=    0.156078347227440d0
    C_ijLM(  2,  1,  0,  0,  2,  1  )=    0.282094791773878d0
    C_ijLM(  1,  0,  1,  1,  2,  1  )=    0.218509686118416d0
    C_ijLM(  1,  1,  1,  0,  2,  1  )=    0.218509686118416d0
    C_ijLM(  3, -2,  1, -1,  2,  1  )=   -0.184674390922372d0
    C_ijLM(  3,  0,  1,  1,  2,  1  )=   -0.143048168102669d0
    C_ijLM(  3,  1,  1,  0,  2,  1  )=    0.233596680327607d0
    C_ijLM(  3,  2,  1,  1,  2,  1  )=    0.184674390922372d0
    C_ijLM(  0,  0,  2,  1,  2,  1  )=    0.282094791773878d0
    C_ijLM(  2, -2,  2, -1,  2,  1  )=   -0.156078347227440d0
    C_ijLM(  2, -1,  2, -2,  2,  1  )=   -0.156078347227440d0
    C_ijLM(  2,  0,  2,  1,  2,  1  )=    0.090111875786434d0
    C_ijLM(  2,  1,  2,  0,  2,  1  )=    0.090111875786434d0
    C_ijLM(  2,  1,  2,  2,  2,  1  )=    0.156078347227440d0
    C_ijLM(  2,  2,  2,  1,  2,  1  )=    0.156078347227440d0
    C_ijLM(  4, -3,  2, -2,  2,  1  )=   -0.168583882836184d0
    C_ijLM(  4, -2,  2, -1,  2,  1  )=   -0.180223751572869d0
    C_ijLM(  4, -1,  2, -2,  2,  1  )=    0.063718718434028d0
    C_ijLM(  4,  0,  2,  1,  2,  1  )=   -0.161197023870787d0
    C_ijLM(  4,  1,  2,  0,  2,  1  )=    0.220728115441823d0
    C_ijLM(  4,  1,  2,  2,  2,  1  )=   -0.063718718434028d0
    C_ijLM(  4,  2,  2,  1,  2,  1  )=    0.180223751572869d0
    C_ijLM(  4,  3,  2,  2,  2,  1  )=    0.168583882836184d0
    C_ijLM(  2,  2,  0,  0,  2,  2  )=    0.282094791773878d0
    C_ijLM(  1, -1,  1, -1,  2,  2  )=   -0.218509686118416d0
    C_ijLM(  1,  1,  1,  1,  2,  2  )=    0.218509686118416d0
    C_ijLM(  3, -3,  1, -1,  2,  2  )=    0.226179013159540d0
    C_ijLM(  3, -1,  1, -1,  2,  2  )=    0.058399170081902d0
    C_ijLM(  3,  1,  1,  1,  2,  2  )=   -0.058399170081902d0
    C_ijLM(  3,  2,  1,  0,  2,  2  )=    0.184674390922372d0
    C_ijLM(  3,  3,  1,  1,  2,  2  )=    0.226179013159540d0
    C_ijLM(  0,  0,  2,  2,  2,  2  )=    0.282094791773878d0
    C_ijLM(  2, -1,  2, -1,  2,  2  )=   -0.156078347227440d0
    C_ijLM(  2,  0,  2,  2,  2,  2  )=   -0.180223751572869d0
    C_ijLM(  2,  1,  2,  1,  2,  2  )=    0.156078347227440d0
    C_ijLM(  2,  2,  2,  0,  2,  2  )=   -0.180223751572869d0
    C_ijLM(  4, -4,  2, -2,  2,  2  )=    0.238413613504448d0
    C_ijLM(  4, -3,  2, -1,  2,  2  )=    0.168583882836184d0
    C_ijLM(  4, -1,  2, -1,  2,  2  )=    0.063718718434028d0
    C_ijLM(  4,  0,  2,  2,  2,  2  )=    0.040299255967697d0
    C_ijLM(  4,  1,  2,  1,  2,  2  )=   -0.063718718434028d0
    C_ijLM(  4,  2,  2,  0,  2,  2  )=    0.156078347227440d0
    C_ijLM(  4,  3,  2,  1,  2,  2  )=    0.168583882836184d0
    C_ijLM(  4,  4,  2,  2,  2,  2  )=    0.238413613504448d0

    return
  End SUBROUTINE getCijLM

END MODULE ForcePSnonLoc2
