MODULE PSQInit
  
  use VarPSMember
  use VarPSMemberG
  use parallel_module, only: myrank

  implicit none
  
  PRIVATE
  PUBLIC :: initKtoKPSQ


CONTAINS

!-----------------------------------------------
  SUBROUTINE initKtoKPSQ()
    implicit none
    integer :: ik,l1,l2,i1,i2,m1,m2
    integer :: k1,k2,k3,nr1,nr2
    integer :: mm1,mm2,mmin,mmax
    logical :: disp_switch_local
    
    disp_switch_local=(myrank==0)

    call allocateKtoK( k1max,k2max,Nelement_local,Rrefmax,Lrefmax )

    do ik=1,Nelement_local
      k1=0
      nr1=0
      do l1=1,nl(ik)
        do i1=1,nr(l1,ik)
          nr1=nr1+1
          do m1=1,2*l1-1
            nr2=0
            do l2=1,l1
              do i2=1,nr(l2,ik)
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

    do ik=1,Nelement_local
      do k1=1,N_k1(ik)
        k2=k1_to_k2(k1,ik)
        if (icheck_k2(k2)==1) cycle
        icheck_k2(k2)=1
        
        N_k2(ik)=N_k2(ik)+1
        k2_to_iorb(1,k2,ik)=k1_to_iorb(1,k1,ik)
        k2_to_iorb(2,k2,ik)=k1_to_iorb(2,k1,ik)
      end do
    end do

    return
  END SUBROUTINE initKtoKPSQ

!-----------------------------------------------

  SUBROUTINE ps_Q_init(qcut)
    implicit none
    integer,intent(IN) :: qcut

    integer :: NRc
    integer :: MMr
    real(8) :: Rc

    integer :: ik,k2
    integer :: iorb1,iorb2

    integer :: iloc(1)

    integer,allocatable :: Q_NRps(:,:)
    real(8),allocatable :: Q_Rps(:,:)
    real(8),allocatable :: Q_wm(:,:,:)

    qc = qcut*qcfac
    if ( qc<=0.d0 ) qc=qcut

    NRc=maxval( NRps )
    allocate( Q_NRps(k2max,Nelement_local) ) ; Q_NRps(:,:)=0
    allocate( Q_Rps(k2max,Nelement_local)  ) ; Q_Rps(:,:)=0.d0
    allocate( Q_wm(NRc,k2max,nki)          ) ; Q_wm(:,:,:)=0.d0

    do ik=1,Nelement_local
      MMr=Mr(ik)
      do k2=1,N_k2(ik)
        iorb1=k2_to_iorb(1,k2,ik)
        iorb2=k2_to_iorb(2,k2,ik)

        Rc=max( Rps0(iorb1,ik),Rps0(iorb2,ik) )*rcfac
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




  END SUBROUTINE ps_Q_init

END MODULE PSQInit
