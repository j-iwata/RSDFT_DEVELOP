MODULE PSQInit
  use VarPSMember
  use VarPSMemberG
  use parallel_module, only: myrank

  use maskf_module

  use Filtering, only: opFiltering

  implicit none
  
  PRIVATE
  PUBLIC :: initKtoKPSQ,ps_Q_init,ps_Q_init_derivative


CONTAINS

!-----------------------------------------------
  SUBROUTINE initKtoKPSQ
    implicit none
    integer :: ik,l1,l2,i1,i2,m1,m2
    integer :: k1,k2,k3,nr1,nr2
    integer :: mm1,mm2,mmin,mmax
    logical :: disp_switch_local
    
    disp_switch_local=(myrank==0)

#ifdef _SHOWALL_INIT_
write(200+myrank,*) ">>>>> initKtoKPSQ"
#endif
#ifdef _SHOWALL_INIT_
write(650+myrank,*) "allocateKtoK"
write(650+myrank,*) "k1max=",k1max
#endif
    call allocateKtoK( k1max,max_k2,Nelement_,max_Rref,max_Lref )
#ifdef _SHOWALL_INIT_
write(650+myrank,*) "l1,i2,nr1,m1,l2,i2,nr2,m2,k1,k2,k3"
#endif
    do ik=1,Nelement_
      k1=0
      nr1=0
      do l1=1,nlf(ik)
        do i1=1,nrf(l1,ik)
          nr1=nr1+1
          do m1=1,2*l1-1
            nr2=0
            do l2=1,l1
              do i2=1,nrf(l2,ik)
                if (.not. ((l1==l2) .and. (i2>i1))) then
                  nr2=nr2+1
                  k2=(nr1*(nr1-1))/2+nr2
                  do m2=1,2*l2-1
                    if (.not. ((l1==l2) .and. (i1==i2) .and. (m2>m1))) then
                      k1=k1+1
                      mm1=(l1-1)**2+m1
                      mm2=(l2-1)**2+m2
                      mmin=min(mm1,mm2)
                      mmax=max(mm1,mm2)
                      k3=(mmax-1)*mmax/2 + mmin
#ifdef _SHOWALL_INIT_
write(650+myrank,'(11I3)') l1,i1,nr1,m1,l2,i2,nr2,m2,k1,k2,k3
#endif

                      k1_to_k2(k1,ik)=k2
                      k1_to_k3(k1,ik)=k3

                      k1_to_iorb(1,k1,ik)=nr1
                      k1_to_iorb(2,k1,ik)=nr2

                      if ((l1==l2) .and. (m1==1) .and. (m2==1)) then
                        qqc(i1,i2,l1,ik) = qqr(i1,i2,l1,ik)
                        qqc(i2,i1,l1,ik) = qqr(i1,i2,l1,ik)

                        if (disp_switch_local) then
                          write(*,*) "#qqc,qqr= ",i1,i2,l1,k1
                          write(*,*) "#qqc,used(=qqr in ps=12)= ",qqc(i1,i2,l1,ik)
                          write(*,*) "#qqr,read= ",qqr(i1,i2,l1,ik)
                        end if
                      end if
                    end if
                  end do ! m2
                end if
              end do ! i2
            end do ! l2
          end do ! m1
        end do ! i1
      end do ! l1
      N_k1(ik)=k1
    end do ! ik

    do ik=1,Nelement_
      do k1=1,N_k1(ik)
        k2=k1_to_k2(k1,ik)
        if (icheck_k2(k2)==1) cycle
        icheck_k2(k2)=1
        
        N_k2(ik)=N_k2(ik)+1
        k2_to_iorb(1,k2,ik)=k1_to_iorb(1,k1,ik)
        k2_to_iorb(2,k2,ik)=k1_to_iorb(2,k1,ik)
      end do
    end do

#ifdef _SHOWALL_KtoK_
    do ik=1,Nelement_
      write(1100+myrank,*) ik,'N_k1(ik)= ',N_k1(ik),'final'
      write(1320+myrank,*) 'k1,k1_to_k2(k1,ik)',ik
      write(1340+myrank,*) 'k1,k1_to_k3(k1,ik)',ik
      write(1360+myrank,*) 'k1,k1_to_iorb(1,k1,ik)',ik
      write(1380+myrank,*) 'k1,k1_to_iorb(2,k1,ik)',ik
      do k1=1,N_k1(ik)
        write(1320+myrank,*) k1,k1_to_k2(k1,ik)
        write(1340+myrank,*) k1,k1_to_k3(k1,ik)
        write(1360+myrank,*) k1,k1_to_iorb(1,k1,ik)
        write(1380+myrank,*) k1,k1_to_iorb(2,k1,ik)
      enddo
      write(1100+myrank,*) ik,'N_k2(ik)= ',N_k2(ik),'final'
      write(1400+myrank,*) 'k2,k2_to_iorb(1,k2,ik)',ik
      write(1420+myrank,*) 'k2,k2_to_iorb(2,k2,ik)',ik
      do k2=1,N_k2(ik)
        write(1400+myrank,*) k2,k2_to_iorb(1,k2,ik)
        write(1420+myrank,*) k2,k2_to_iorb(2,k2,ik)
      enddo
    enddo
#endif
#ifdef _SHOWALL_INIT_
write(650+myrank,*) "<<<<< initKtoKPSQ"
#endif
    return
  END SUBROUTINE initKtoKPSQ

!-----------------------------------------------
#ifdef _USPP_
  SUBROUTINE ps_Q_init(qcut,rcfac,qcfac,etafac)
    implicit none
    real(8),intent(IN) :: qcut
    real(8),intent(IN) :: rcfac,qcfac,etafac

    integer :: NRc
    integer :: MMr
    real(8) :: Rc

    integer :: ik,k2,i,m,m0,m1,m2,ll3,L
    integer :: iorb1,iorb2

    integer :: iloc(1)

    real(8) :: x,y,dy,y0,dy0
    real(8),parameter :: dr=2.d-3,ep=1.d-14
    real(8) :: maxerr

    real(8) :: qc
    real(8) :: Q_rcfac,Q_etafac

    real(8),allocatable :: Q_wm(:,:,:)
    real(8),allocatable :: vrad(:),tmp(:)

#ifdef _SHOWALL_INIT_
write(200+myrank,*) ">>>>> ps_Q_init"
#endif
Q_rcfac=rcfac
Q_etafac=etafac
k2max=max_k2

    qc = qcut*qcfac
    if ( qc<=0.d0 ) qc=qcut

!qc=4.32987324642663
if (myrank==0) write(200,*) 'qc(qr)= ',qc

    call allocateQRps( k2max,Nelement_ )

    do ik=1,Nelement_
      MMr=Mr(ik)
      do k2=1,N_k2(ik)
        iorb1=k2_to_iorb(1,k2,ik)
        iorb2=k2_to_iorb(2,k2,ik)

        Rc=max( Rps0(iorb1,ik),Rps0(iorb2,ik) )*Q_rcfac
        iloc=minloc( abs(rad(1:MMr,ik)-Rc) )
        NRc=iloc(1) ; if (rad(NRc,ik)<Rc) NRc=NRc+1

        if (NRc>MMr) then
          write(*,*) "NRc,MMr= ",NRc,MMr
          stop "rcfac is too large. (Q-part)"
        end if

        Q_NRps(k2,ik)=NRc
        Q_Rps(k2,ik)=rad(NRc,ik)
      end do
    end do

    NRc=maxval( Q_NRps )
    allocate( Q_wm(NRc,k2max,Nelement_) ) ; Q_wm(:,:,:)=0.d0

    do ik=1,Nelement_
        do k2=1,N_k2(ik)
            NRc=Q_NRps(k2,ik)
            Rc=Q_Rps(k2,ik)
!write(520+myrank,*) '------------------------'
!write(52+myrank,'(3I5,2E15.7e2)') ik,k2,NRc,Rc,Q_etafac
            call makemaskf(Q_etafac)

            maxerr=0.d0
            do i=1,NRc
                x=rad(i,ik)/Rc
                if ( x<=dxm ) then
                  y0=1.d0 ; dy0=0.d0
                else
                  m0=int(x/dxm)
                  dy0=1.d10
                  do m=1,20
                      m1=max(m0-m,1) ; m2=min(m0+m,nmsk)
!write(520+myrank,'(3I5,E15.7e2)') m,m1,m2,x
                      call polint(xm(m1),maskr(m1),m2-m1+1,x,y,dy)
!write(520+myrank,'(I5,3E15.7e2)') m,x,y,dy
                      if ( abs(dy)<dy0 ) then
                        y0=y ; dy0=abs(dy)
                      end if
                  end do
                end if
                Q_wm(i,k2,ik)=y0
!write(520+myrank,'(3I5,E15.7e2)') ik,k2,i,Q_wm(i,k2,ik)
                maxerr=max(maxerr,dy0)
            end do
        end do ! k2
    end do ! ik

    do ik=1,Nelement_
        do k2=1,N_k2(ik)
!            Q_NRps(k2,ik)=Q_Rps(k2,ik)/dr
            Q_NRps(k2,ik)=Q_Rps(k2,ik)/dr+1
!            if ( Q_NRps(k2,ik)>max_psgrd ) stop 'Q_NRps too large'
            if ( (Q_NRps(k2,ik)-1)*dr < Q_Rps(k2,ik) ) then
              Q_NRps(k2,ik)=Q_NRps(k2,ik)+1
            endif
        end do 
    end do

    MMr=max( maxval(Mr),maxval(Q_NRps) )
#ifdef _SHOWALL_QINIT_
    if (myrank==0) write(200,'(A12,I5)') 'rad1(Q)=',MMr
#endif
    call allocateRad1(MMr)
    do ik=1,Nelement_
        do i=1,MMr
            rad1(i,ik)=(i-1)*dr
        end do
    end do

    NRc=maxval(NRps0)
    allocate( vrad(NRc),tmp(NRc) )

    do ik=1,Nelement_
        do k2=1,N_k2(ik)
            iorb1=k2_to_iorb(1,k2,ik)
            iorb2=k2_to_iorb(2,k2,ik)
            NRc=max( NRps0(iorb1,ik),NRps0(iorb2,ik) )

            do ll3=1,nl3v(k2,ik)
                L=l3v(ll3,k2,ik)-1
                vrad(1:NRc)=qrL(1:NRc,ll3,k2,ik)*rab(1:NRc,ik)/Q_wm(1:NRc,k2,ik)

#ifdef _SHOWALL_QR_F_
write(520+myrank,*) "ik,k2,ll3= ",ik,k2,ll3
write(520+myrank,*) "    qc L  NRc NRps            rad           rad1            qrL"
write(520+myrank,'(f6.3,I2,2I5,3E15.7e2)') qc,L,NRc,Q_NRps(k2,ik),rad(1,ik),rad1(1,ik),qrL(1,ll3,k2,ik)
#endif
                call opFiltering( qc,L,NRc,Q_NRps(k2,ik),rad(1,ik),rad1(1,ik),vrad,qrL(1,ll3,k2,ik) )

#ifdef _SHOWALL_QR_F_
write(520+myrank,*) "ik,k2,ll3= ",ik,k2,ll3
write(520+myrank,*) "       rad     rad1      vrad      Q_wm       qrL"
do i=1,10
  write(520+myrank,'(I5,5E10.2e2)') i,rad(i,ik),rad1(i,ik),vrad(i),Q_wm(i,k2,ik),qrL(i,ll3,k2,ik)
end do
do i=NRc-10,NRc
  write(520+myrank,'(I5,5E10.2e2)') i,rad(i,ik),rad1(i,ik),vrad(i),Q_wm(i,k2,ik),qrL(i,ll3,k2,ik)
end do
#endif
            end do ! ll3
        end do ! k2
    end do ! ik

    deallocate( vrad,tmp )

    do ik=1,Nelement_
        do k2=1,N_k2(ik)
            NRc=Q_NRps(k2,ik)
            Rc=Q_Rps(k2,ik)
            call makemaskf(Q_etafac)
            maxerr=0.d0

            do i=1,NRc
                x=(i-1)*dr/Rc
                if ( x<dxm ) then
                  y0=1.d0 ; dy0=0.d0
                else
                  m0=int(x/dxm)
                  dy0=1.d10
                  do m=1,20
                      m1=max(m0-m,1) ; m2=min(m0+m,nmsk)
                      call polint(xm(m1),maskr(m1),m2-m1+1,x,y,dy)
                      if ( abs(dy)<dy0 ) then
                        y0=y ; dy0=abs(dy)
                      end if
                  end do
                end if

                if ( maxerr<dy0 ) maxerr=dy0

                do ll3=1,nl3v(k2,ik)
                    qrL(i,ll3,k2,ik)=y0*qrL(i,ll3,k2,ik)
                end do
            end do ! i
        end do ! k2
    end do ! ik

    deallocate( Q_wm )

#ifdef _SHOWALL_INIT_
write(400+myrank,*) "<<<<< ps_Q_init"
#endif
    return
  END SUBROUTINE ps_Q_init
#else
  SUBROUTINE ps_Q_init(qcut,rcfac,qcfac,etafac)
    implicit none
    real(8),intent(IN) :: qcut
    real(8),intent(IN) :: rcfac,qcfac,etafac
    return
  END SUBROUTINE ps_Q_init
#endif

!-----------------------------------------------
#ifdef _USPP_
  SUBROUTINE ps_Q_init_derivative
    implicit none
    integer :: ik,L,NRc,J,i,m,m1,m2,lm,ll3
    integer :: ik1,k2,k3,iorb1,iorb2,cJ
    real(8) :: maxerr,y,dy,y0,dy0
    real(8) :: pi4,const
    real(8),allocatable :: dwork(:,:,:,:,:)
    real(8),parameter :: sqrt_4pi_3=sqrt(4.d0*acos(-1.d0)/3.d0)

    maxcJ=0
    do ik=1,Nelement_
      do ik1=1,N_k1(ik)
        k2=k1_to_k2(ik1,ik)
        do ll3=1,nl3v(k2,ik)
          L=l3v(ll3,k2,ik)-1
          cJ=0
          do J=abs(L-1),L+1
            cJ=cJ+1
          enddo
          maxcJ=max(cJ,maxcJ)
        enddo
      enddo
    enddo
    NRc=max_psgrd
    k2=max_k2
    l=max_Lref
    ik=Nelement_
    allocate( dqrL(NRc,l,k2,ik,maxcJ) ) ; dqrL=0.d0
    allocate( dwork(NRc,l,k2,ik,3) ) ; dwork=0.d0

    do ik=1,Nelement_
      do ik1=1,N_k1(ik)
        k2=k1_to_k2(ik1,ik)
        k3=k1_to_k3(ik1,ik)
        iorb1=k1_to_iorb(1,ik1,ik)
        iorb2=k1_to_iorb(2,ik1,ik)
        NRc=max(NRps(iorb1,ik),NRps(iorb2,ik))
        do ll3=1,nl3v(k2,ik)
          L=l3v(ll3,k2,ik)-1
          maxerr=0.d0
          do i=1,NRc
            dy0=1.d10
            do m=1,20
              m1=max(i-m,1) ; m2=min(i+m,NRc)
              call dpolint( rad1(m1,ik),qrL(m1,ll3,k2,ik),m2-m1+1,rad1(i,ik),y,dy )
              if ( abs(dy)<dy0 ) then
                 y0=y ; dy0=abs(dy)
              end if
            end do ! m
            dqrL(i,ll3,k2,ik,1)=y0
            maxerr=max(maxerr,dy0)
          end do ! i
        enddo ! ll3
      end do ! ik1
    end do !ik

    do ik=1,Nelement_
      do ik1=1,N_k1(ik)
        k2=k1_to_k2(ik1,ik)
        k3=k1_to_k3(ik1,ik)
        iorb1=k1_to_iorb(1,ik1,ik)
        iorb2=k1_to_iorb(2,ik1,ik)
        NRc=max(NRps(iorb1,ik),NRps(iorb2,ik))
        do ll3=1,nl3v(k2,ik)
          L=l3v(ll3,k2,ik)-1
          cJ=0
          do J=abs(L-1),L+1
            cJ=cJ+1
            const=0.5d0*(2.d0+L*(L+1)-J*(J+1))
            do i=1,NRc
              dwork(i,ll3,k2,ik,cJ)=(rad1(i,ik)**3)*dqrL(i,ll3,k2,ik,1)+const*qrL(i,ll3,k2,ik)*(rad1(i,ik)**2)
            end do
          end do ! j
        end do ! ll3
      enddo ! k1
    end do ! ik

    do ik=1,Nelement_
      do ik1=1,N_k1(ik)
        k2=k1_to_k2(ik1,ik)
        k3=k1_to_k3(ik1,ik)
        iorb1=k1_to_iorb(1,ik1,ik)
        iorb2=k1_to_iorb(2,ik1,ik)
        NRc=max(NRps(iorb1,ik),NRps(iorb2,ik))
        do ll3=1,nl3v(k2,ik)
          L=l3v(ll3,k2,ik)-1
          cJ=0
          do J=abs(L-1),L+1
            cJ=cJ+1
            do i=1,NRc
              dqrL(i,ll3,k2,ik,cJ)=sqrt_4pi_3*dwork(i,ll3,k2,ik,cJ)
            end do
          enddo
        end do
      end do
    end do

    deallocate( dwork )
  END SUBROUTINE ps_Q_init_derivative
#else
  SUBROUTINE ps_Q_init_derivative
    implicit none
    return
  END SUBROUTINE ps_Q_init_derivative
#endif

END MODULE PSQInit
