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
    integer :: kk1
    integer :: ierr,M_irad,ir0
    integer,allocatable :: ireq(:),istatus(:,:),irad(:,:)
    real(8),parameter :: ep=1.d-8
    real(8) :: err,err0,maxerr,Rx,Ry,Rz
    real(8) :: a1,a2,a3,c1,c2,c3,d1,d2,d3
    real(8) :: x,y,z,r,kr,pi2,c
    real(8) :: tmp,tmp0,tmp1
    real(8) :: yy1,yy2,yy3
    real(8),allocatable :: work2(:,:)
    real(8),allocatable :: dQY(1:3)
    real(8) :: rf_Q(1:3)
    integer,allocatable :: icheck_tmpQ(:)
    integer :: JJ_tmp(0:3)
#ifdef _DRSDFT_
    real(8) :: ztmp
    real(8),allocatable :: rtmp5(:,:,:),rtmp2(:,:)
#else
    complex(8) :: ztmp
    complex(8),allocatable :: rtmp5(:,:,:),rtmp2(:,:)
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

    allocate( rtmp5(0:2,N_nzqr,MSP_0:MSP_1) )
    allocate( a_rank(MI) )
    allocate( ickeck_tmpQ(0:np_grid-1) )

!$OMP parallel

!$OMP workshare
    rtmp5=zero
    a_rank=.false.
    icheck_tmpQ=0
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


    do ia=1,MI
      ik=ki_atom(ia)
      Rx=aa(1,1)*aa_atom(1,a)+aa(1,2)*aa_atom(2,a)+aa(1,3)*aa_atom(3,a)
      Ry=aa(2,1)*aa_atom(1,a)+aa(2,2)*aa_atom(2,a)+aa(2,3)*aa_atom(3,a)
      Rz=aa(3,1)*aa_atom(1,a)+aa(3,2)*aa_atom(2,a)+aa(3,3)*aa_atom(3,a)
      do ik1=1,N_k1(ik)
        kk1=kk1map(ik1,ia)
        if (kk1==0) cycle
        lma1=nzqr_pair(kk1,1)
        lma2=nzqr_pair(kk1,2)
        l_1=lmap(lma1)
        m_1=mmap(lma1)
        l_2=lmap(lma2)
        m_2=mmap(lma2)
        k2=k1_to_k2(ik1,ik)
        k3=k1_to_k3(ik1,ik)
        iorb1=k1_to_iorb(1,ik1,ik)
        iorb2=k1_to_iorb(2,ik1,ik)
        Rps2=max(Rps(iorb1,ik)**2,Rps(iorb2,ik)**2)
        NRc=max(NRps(iorb1,ik),NRps(iorb2,ik))
        rf_Q=0.d0
        do i=1,M_grid_ion
        ! ppp2????
          d1=map_grid_ion(1,i)
          d2=map_grid_ion(2,i)
          d3=map_grid_ion(3,i)
          x = aa(1,1)*d1+aa(1,2)*d2+aa(1,3)*d3-Rx
          y = aa(2,1)*d1+aa(2,2)*d2+aa(2,3)*d3-Ry
          z = aa(3,1)*d1+aa(3,2)*d2+aa(3,3)*d3-Rz
          r2 = x*x+y*y+z*z
          if (r2>Rps2+1.d-10) cycle
          r = sqrt(x*x+y*y+z*z)
          ir0=irad( int(100.d0*r),ik )
          do ir=ir0,NRc
            if ( r<rad1(ir,ik) ) exit
          end do
          dQY(:)=0.d0

          do ll3=1,nl3v(k2,ik)
            L=l3v(ll3,k2,ik)-1
            cJ=0
            do L1=abs(L-1),L+1
              cJ=cJ+1
              tmp0=0.d0
              err0=0.d0
              if ( abs(x)>1.d-14 .or. abs(y)>1.d-14 &
                .or. abs(z)>1.d-14 .or. L1==0 ) then
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
                  tmp0=tmp0/(r*r*r)
                else
                  write(*,*) "force_ps_Q",ir,NRc
                  stop
                end if
                maxerr=max(maxerr,err0)
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
                dQY(1)=dQY(1)+tmp0*yy1
                dQY(2)=dQY(2)+tmp0*yy2
                dQY(3)=dQY(3)+tmp0*yy3
              end if
            end do ! L1
          enddo ! ll3
          if (dQY(1)/=0.d0.or.dQY(2)/=0.d0.or.dQY(3)/=0.d0) then
            rf_Q(1)=rf_Q(1)+dQY(1)
            rf_Q(2)=rf_Q(2)+dQY(2)
            rf_Q(3)=rf_Q(3)-dQY(3)
          endif
        enddo ! i
        call mpi_allreduce(rf_Q,mpi_in_place,3,mpi_real8,mpi_sum,comm_grid,ierr)
        do s=MSP_0,MSP_1
          rtmp5(1,kk1,s)=rtmp5(1,kk1,s)+rf_Q(1)
        enddo
        !vloc?????????

      end do ! j
    end do ! iqr


!#ifndef _SPLINE_
!$OMP single
    deallocate( irad )
!$OMP end single
!#endif


!if (myrank==0) write(230,'(5a5,a20)') 's','k','n','iqr','j','ztmp'
    do s=MSP_0,MSP_1
      do k=MBZ_0,MBZ_1
!$OMP do schedule(dynamic) private( c,i,d1,d2,d3,kr,ztmp,i1,i2,i3 )
        do n=MB_0,MB_1
!          if ( occ(n,k,s) == 0.d0 ) cycle
          if ( occ(n,k,s) < 1.d-10 ) cycle
!          c=-2.d0*occ(n,k,s)*dV*dV
          do iqr=1,N_nzqr
            do j=1,MJJ_MAP_Q(iqr)
              i1=JJ_MAP_Q(1,j,iqr)
              i2=JJ_MAP_Q(2,j,iqr)
              i3=JJ_MAP_Q(3,j,iqr)
              i = i1-a1b + (i2-a2b)*ab1 + (i3-a3b)*ab1*ab2 + ML_0
              ztmp=Vloc(i,s)
              rtmp5(0,iqr,n,k,s)=rtmp5(0,iqr,n,k,s)+dQY(1,j,iqr)*ztmp
              rtmp5(1,iqr,n,k,s)=rtmp5(1,iqr,n,k,s)+dQY(2,j,iqr)*ztmp
              rtmp5(2,iqr,n,k,s)=rtmp5(2,iqr,n,k,s)-dQY(3,j,iqr)*ztmp
!if (myrank==0) write(230,'(5I5,g20.7)') s,k,n,iqr,j,Vloc(i,s)
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
    do s=MSP_0,MSP_1
      do k=MBZ_0,MBZ_1
        do n=MB_0,MB_1,MB_d
          ib1=n
          ib2=min(ib1+MB_d-1,MB_1)
          if ( occ(n,k,s) < 1.d-10 ) cycle
          call threeWayComm(nrqr_xyz,num_2_rank_Q,sendmap_Q,recvmap_Q,qr_nsend,sbufnl_Q,rbufnl_Q,N_nzqr,ib1,ib2,rtmp5(0,1,ib1,k,s),2)
        enddo
      enddo
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
          if ( a_rank(a) ) then
            if (lma1<lma2) stop 'Nzqr_pair is strange'
            tmp1=-dV
!            tmp1=-rtmp5(0,m,n,k,s)*dV
!            tmp2=-rtmp5(1,m,n,k,s)*dV
!            tmp3=-rtmp5(2,m,n,k,s)*dV
#ifdef _DRSDFT_
            if (lma1==lma2) then
              forceQ(1,a)=forceQ(1,a)+sum(occ(MB_0:MB_1,k,s)*real(uVunk_tmp(lma1,MB_0:MB_1,k,s)*uVunk_tmp(lma2,MB_0:MB_1,k,s),8)*rtmp5(0,m,MB_0:MB_1,k,s))*tmp1
              forceQ(2,a)=forceQ(2,a)+sum(occ(MB_0:MB_1,k,s)*real(uVunk_tmp(lma1,MB_0:MB_1,k,s)*uVunk_tmp(lma2,MB_0:MB_1,k,s),8)*rtmp5(1,m,MB_0:MB_1,k,s))*tmp1
              forceQ(3,a)=forceQ(3,a)+sum(occ(MB_0:MB_1,k,s)*real(uVunk_tmp(lma1,MB_0:MB_1,k,s)*uVunk_tmp(lma2,MB_0:MB_1,k,s),8)*rtmp5(2,m,MB_0:MB_1,k,s))*tmp1
!              forceQ(1,a)=forceQ(1,a)+real(uVunk_tmp(lma1,n,k,s)*uVunk_tmp(lma2,n,k,s),8)*tmp1
!              forceQ(2,a)=forceQ(2,a)+real(uVunk_tmp(lma1,n,k,s)*uVunk_tmp(lma2,n,k,s),8)*tmp2
!              forceQ(3,a)=forceQ(3,a)+real(uVunk_tmp(lma1,n,k,s)*uVunk_tmp(lma2,n,k,s),8)*tmp3
            else
              forceQ(1,a)=forceQ(1,a)+sum(occ(MB_0:MB_1,k,s)*real(uVunk_tmp(lma1,MB_0:MB_1,k,s)*uVunk_tmp(lma2,MB_0:MB_1,k,s),8)*rtmp5(0,m,MB_0:MB_1,k,s))*tmp1
              forceQ(2,a)=forceQ(2,a)+sum(occ(MB_0:MB_1,k,s)*real(uVunk_tmp(lma1,MB_0:MB_1,k,s)*uVunk_tmp(lma2,MB_0:MB_1,k,s),8)*rtmp5(1,m,MB_0:MB_1,k,s))*tmp1
              forceQ(3,a)=forceQ(3,a)+sum(occ(MB_0:MB_1,k,s)*real(uVunk_tmp(lma1,MB_0:MB_1,k,s)*uVunk_tmp(lma2,MB_0:MB_1,k,s),8)*rtmp5(2,m,MB_0:MB_1,k,s))*tmp1
              forceQ(1,a)=forceQ(1,a)+sum(occ(MB_0:MB_1,k,s)*real(uVunk_tmp(lma2,MB_0:MB_1,k,s)*uVunk_tmp(lma1,MB_0:MB_1,k,s),8)*rtmp5(0,m,MB_0:MB_1,k,s))*tmp1
              forceQ(2,a)=forceQ(2,a)+sum(occ(MB_0:MB_1,k,s)*real(uVunk_tmp(lma2,MB_0:MB_1,k,s)*uVunk_tmp(lma1,MB_0:MB_1,k,s),8)*rtmp5(1,m,MB_0:MB_1,k,s))*tmp1
              forceQ(3,a)=forceQ(3,a)+sum(occ(MB_0:MB_1,k,s)*real(uVunk_tmp(lma2,MB_0:MB_1,k,s)*uVunk_tmp(lma1,MB_0:MB_1,k,s),8)*rtmp5(2,m,MB_0:MB_1,k,s))*tmp1
            endif
#else
            if (lma1==lma2) then
              forceQ(1,a)=forceQ(1,a)+sum(occ(MB_0:MB_1,k,s)*real(conjg(uVunk_tmp(lma1,MB_0:MB_1,k,s))*uVunk_tmp(lma2,MB_0:MB_1,k,s),8)*rtmp5(0,m,MB_0:MB_1,k,s))*tmp1
              forceQ(2,a)=forceQ(2,a)+sum(occ(MB_0:MB_1,k,s)*real(conjg(uVunk_tmp(lma1,MB_0:MB_1,k,s))*uVunk_tmp(lma2,MB_0:MB_1,k,s),8)*rtmp5(1,m,MB_0:MB_1,k,s))*tmp1
              forceQ(3,a)=forceQ(3,a)+sum(occ(MB_0:MB_1,k,s)*real(conjg(uVunk_tmp(lma1,MB_0:MB_1,k,s))*uVunk_tmp(lma2,MB_0:MB_1,k,s),8)*rtmp5(2,m,MB_0:MB_1,k,s))*tmp1
!              forceQ(1,a)=forceQ(1,a)+real(conjg(uVunk_tmp(lma1,n,k,s))*uVunk_tmp(lma2,n,k,s),8)*tmp1
!              forceQ(2,a)=forceQ(2,a)+real(conjg(uVunk_tmp(lma1,n,k,s))*uVunk_tmp(lma2,n,k,s),8)*tmp2
!              forceQ(3,a)=forceQ(3,a)+real(conjg(uVunk_tmp(lma1,n,k,s))*uVunk_tmp(lma2,n,k,s),8)*tmp3
            else
              forceQ(1,a)=forceQ(1,a)+sum(occ(MB_0:MB_1,k,s)*real(conjg(uVunk_tmp(lma1,MB_0:MB_1,k,s))*uVunk_tmp(lma2,MB_0:MB_1,k,s),8)*rtmp5(0,m,MB_0:MB_1,k,s))*tmp1
              forceQ(2,a)=forceQ(2,a)+sum(occ(MB_0:MB_1,k,s)*real(conjg(uVunk_tmp(lma1,MB_0:MB_1,k,s))*uVunk_tmp(lma2,MB_0:MB_1,k,s),8)*rtmp5(1,m,MB_0:MB_1,k,s))*tmp1
              forceQ(3,a)=forceQ(3,a)+sum(occ(MB_0:MB_1,k,s)*real(conjg(uVunk_tmp(lma1,MB_0:MB_1,k,s))*uVunk_tmp(lma2,MB_0:MB_1,k,s),8)*rtmp5(2,m,MB_0:MB_1,k,s))*tmp1
              forceQ(1,a)=forceQ(1,a)+sum(occ(MB_0:MB_1,k,s)*real(conjg(uVunk_tmp(lma2,MB_0:MB_1,k,s))*uVunk_tmp(lma1,MB_0:MB_1,k,s),8)*rtmp5(0,m,MB_0:MB_1,k,s))*tmp1
              forceQ(2,a)=forceQ(2,a)+sum(occ(MB_0:MB_1,k,s)*real(conjg(uVunk_tmp(lma2,MB_0:MB_1,k,s))*uVunk_tmp(lma1,MB_0:MB_1,k,s),8)*rtmp5(1,m,MB_0:MB_1,k,s))*tmp1
              forceQ(3,a)=forceQ(3,a)+sum(occ(MB_0:MB_1,k,s)*real(conjg(uVunk_tmp(lma2,MB_0:MB_1,k,s))*uVunk_tmp(lma1,MB_0:MB_1,k,s),8)*rtmp5(2,m,MB_0:MB_1,k,s))*tmp1
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
    call mpi_allreduce(work2,forceQ,3*MI &
         ,mpi_real8,mpi_sum,mpi_comm_world,ierr)
    deallocate( work2 )
    deallocate( a_rank )
    deallocate( rtmp5 )
!$OMP end single

!$OMP end parallel

    return
  END SUBROUTINE calcForceQ
