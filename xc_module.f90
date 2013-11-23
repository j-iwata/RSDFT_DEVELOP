MODULE xc_module

  use rgrid_module
  use density_module, only: rho
  use ps_pcc_module
  use parallel_module
  use array_bound_module, only: ML_0,ML_1,MSP,MSP_0,MSP_1
  use pw92_gth_module

  implicit none

  PRIVATE
  PUBLIC :: XCtype,calc_xc,read_xc,Vxc,Exc,E_exchange,E_correlation &
           ,read_oldformat_xc

  character(8) :: XCtype
  real(8),allocatable :: Vxc(:,:)
  real(8) :: Exc,E_exchange,E_correlation

CONTAINS


  SUBROUTINE read_xc(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i
    character(6) :: cbuf,ckey
    if ( rank == 0 ) then
       XCtype = "LDAPZ81"
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey == "XCTYPE" ) then
             backspace(unit)
             read(unit,*) cbuf,XCtype
             exit
          end if
       end do
999    continue
       write(*,*) "XCtype= ",XCtype
    end if
    call send_xc(0)
  END SUBROUTINE read_xc


  SUBROUTINE read_oldformat_xc(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    XCtype=""
    if ( rank == 0 ) then
       read(unit,*) XCtype
       write(*,*) "XCtype= ",XCtype
    end if
    call send_xc(0)
  END SUBROUTINE read_oldformat_xc


  SUBROUTINE send_xc(rank)
    implicit none
    integer,intent(IN) :: rank
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(XCtype,8,MPI_CHARACTER,rank,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_xc


  SUBROUTINE calc_xc
    implicit none
    if ( .not.allocated(Vxc) ) then
       allocate( Vxc(ML_0:ML_1,MSP_0:MSP_1) )
    end if
    select case(XCtype)
    case('LDAPZ81')
       call calc_LDAPZ81
    case('LDAPW92')
       if ( flag_pcc_0 ) stop "PCC is not implemented in LDAPW92" 
       call calc_pw92_gth(ML_0,ML_1,MSP,MSP_0,MSP_1,rho,Exc,Vxc,dV,comm_grid)
    case('GGAPBE96')
       call calc_GGAPBE96
    end select
  END SUBROUTINE calc_xc


  SUBROUTINE calc_LDAPZ81
    implicit none
    real(8),parameter :: gam(1:2)=(/-0.1423d0,-0.0843d0/)
    real(8),parameter :: bet1(1:2)=(/1.0529d0,1.3981d0/)
    real(8),parameter :: bet2(1:2)=(/0.3334d0,0.2611d0/)
    real(8),parameter :: A(1:2)=(/0.0311d0,0.01555d0/)
    real(8),parameter :: B(1:2)=(/-0.048d0,-0.0269d0/)
    real(8),parameter :: C(1:2)=(/0.002d0,0.0007d0/)
    real(8),parameter :: D(1:2)=(/-0.0116d0,-0.0048d0/)
    real(8),parameter :: ep=1.d-25
    real(8),allocatable:: Vxc_tmp(:,:)
    real(8) :: onetwo,onethr,onefou,twothr,thrfou,fouthr,sevsix,thrPi
    real(8) :: cns0,cns1,cns2,cns3,cns4,cns5,cns6
    real(8) :: ctime0,ctime1,etime0,etime1
    real(8) :: s0(2),s1(2)
    real(8) :: trho,rs,rssq,rsln,exd,ecd
    real(8) :: Vx,Vxa,Vxb,Vc,Vca,Vcb,Exa,Exb
    real(8) :: decddrs,decddzet
    real(8) :: zet,f,dfdzet,dfdrhoa,dfdrhob,rhoa,rhob
    real(8) :: ecd0(1:2),decddrs0(1:2)
    integer :: i,ierr,ispin

    onetwo=1.d0/2.d0
    onethr=1.d0/3.d0
    onefou=1.d0/4.d0
    twothr=2.d0/3.d0
    thrfou=3.d0/4.d0
    fouthr=4.d0/3.d0
    sevsix=7.d0/6.d0
    thrPi=3.d0/acos(-1.d0)
    cns0=1.d0/(2.d0**fouthr-2.d0)
    cns1=4.d0/9.d0
    cns2=5.d0/54.d0
    cns3=44.d0/135.d0
    cns4=8.d0/9.d0
    cns5=5.d0/27.d0
    cns6=22.d0/45.d0

    allocate( Vxc_tmp(ML_0:ML_1,1:MSP) )

    E_exchange=0.d0
    Exa=0.d0
    Exb=0.d0
    E_correlation=0.d0

    Vxa=0.d0
    Vxb=0.d0

    do i=ML_0,ML_1
!
! --- Setup ---
!
!!!!! total charge = rhoa+rhob (even if nspin=1)
       rhoa=rho(i,1)*dble(MSP)*onetwo
       rhob=rho(i,MSP)*dble(MSP)*onetwo
       if ( flag_pcc_0 ) then
          rhoa=rhoa+rhoc(i)*onetwo
          rhob=rhob+rhoc(i)*onetwo
       end if
       rhoa=max(rhoa,ep*onetwo)
       rhob=max(rhob,ep*onetwo)
!
! --- Exchange ---
!
!!!!! up-spin (or nspin==1)
       trho=2.d0*rhoa
       exd=-thrfou*(thrPi*trho)**onethr
       Exa=Exa+trho*exd*dV
       Vxa=fouthr*exd
!!!!! down-spin
       if ( MSP==2 ) then
          trho=2.d0*rhob
          exd=-thrfou*(thrPi*trho)**onethr
          Exb=Exb+trho*exd*dV
          Vxb=fouthr*exd
       end if
       E_exchange=(Exa+Exb)/dble(MSP)
!
! --- Correlation ---
!
       trho=rhoa+rhob
       zet=(rhoa-rhob)/trho
       if ( abs(zet)>1.d0 ) stop "calc_ldapz81"
!       if ( abs(zet)>1.d-3 ) then
          f=cns0*((1.d0+zet)**fouthr+(1.d0-zet)**fouthr-2.d0)
          dfdzet=cns0*fouthr*((1.d0+zet)**onethr-(1.d0-zet)**onethr)
!       else
!          f=cns0*cns1*zet**2*(1.d0+cns2*zet**2*(1.d0+cns3*zet**2))
!          dfdzet=cns0*cns4*zet*(1.d0+cns5*zet**2*(1.d0+cns6*zet**2))
!       end if
       rs=(onefou*thrPi/trho)**onethr
       if ( rs >= 1.d0 ) then
          rssq=sqrt(rs)
          ecd0(1:2)=gam(1:2)/(1.d0+bet1(1:2)*rssq+bet2(1:2)*rs)
          decddrs0(1:2)=-(onetwo*bet1(1:2)/rssq+bet2(1:2))*ecd0(1:2) &
               & /(1.d0+bet1(1:2)*rssq+bet2(1:2)*rs)
       else
          if ( rs <= 0.d0 ) stop "calc_ldapz81"
          rsln=log(rs)
          ecd0(1:2)=A(1:2)*rsln+B(1:2)+C(1:2)*rs*rsln+D(1:2)*rs
          decddrs0(1:2)=A(1:2)/rs+C(1:2)*(rsln+1.d0)+D(1:2)
       end if
       ecd=ecd0(1)+f*(ecd0(2)-ecd0(1))
       E_correlation=E_correlation+trho*ecd*dV
       decddrs=decddrs0(1)+f*(decddrs0(2)-decddrs0(1))
       decddzet=dfdzet*(ecd0(2)-ecd0(1))
       Vca=ecd-onethr*rs*decddrs+(1.d0-zet)*decddzet
       Vcb=ecd-onethr*rs*decddrs-(1.d0+zet)*decddzet
!
! --- Exchange-Correlation ---
!
       Vxc_tmp(i,MSP)=Vxb+Vcb
       Vxc_tmp(i,1)=Vxa+Vca
    end do

!    if ( iflag>0 ) then
       Vxc(:,MSP_0:MSP_1)=Vxc_tmp(:,MSP_0:MSP_1)
!    end if

    deallocate( Vxc_tmp )

    s0(1)=E_exchange
    s0(2)=E_correlation
    call mpi_allreduce(s0,s1,2,mpi_real8,mpi_sum,comm_grid,ierr)
    E_exchange=s1(1)
    E_correlation=s1(2)
    Exc=E_exchange+E_correlation

    return
  END SUBROUTINE calc_LDAPZ81

!============================================================== GGAPBE96
!--------1---------2---------3---------4---------5---------6---------7--

  SUBROUTINE calc_GGAPBE96
    use aa_module
    use bb_module
    use bc_module
    use kinetic_module, only: Md,SYStype
    use fd_module
    use electron_module, only: Nspin
    implicit none
    real(8),parameter :: mu=0.21951d0,Kp=0.804d0
    real(8),parameter :: zero_density=1.d-10
    real(8),parameter :: A00  =0.031091d0,A01  =0.015545d0,A02  =0.016887d0
    real(8),parameter :: alp10=0.21370d0 ,alp11=0.20548d0 ,alp12=0.11125d0
    real(8),parameter :: bt10 =7.5957d0  ,bt11 =1.41189d1 ,bt12 =1.0357d1
    real(8),parameter :: bt20 =3.5876d0  ,bt21 =6.1977d0  ,bt22 =3.6231d0
    real(8),parameter :: bt30 =1.6382d0  ,bt31 =3.3662d0  ,bt32 =0.88026d0
    real(8),parameter :: bt40 =0.49294d0 ,bt41 =0.62517d0 ,bt42 =0.49671d0
    real(8),parameter :: C1=2.14611945579107d0,C2=0.031091d0
    integer :: i,j,i1,i2,i3,j1,j2,j3,k1,k2,k3,n1,n2,m,ispin,ierr,itmp
    integer :: Mx,My,Mz,ML1,ML2,ML3,ML,a1b,b1b,a2b,b2b,a3b,b3b
    integer :: l1,l2,l3,m1,m2,m3,irank
    integer,allocatable :: LLL2(:,:,:),LL2(:,:)
    real(8),allocatable :: wtmp(:,:,:),wrho(:,:,:),rtmp(:),gx(:),gy(:),gz(:)
    real(8) :: ctime0,ctime1,etime0,etime1,mem,memax
    real(8) :: g1,g2,g3,b(3,3),sbf(3),rbf(3),H1,H2,H3
    real(8) :: trho,s2,kf,ec_lda,ex_lda,vx_lda,T,Fx
    real(8) :: Hs,A,Ai,dH_dT,dA_dn,dec_dn,dH_dA,rs,tmp,tmp1,tmp2
    real(8) :: ec_U,ec_P,deU_dn,deP_dn,alpc,dac_dn,phi,dphi_dz,pi
    real(8) :: dH_dphi,fz,dfz_dz,dec_dz,const1,const2,srho(2),dz_dn(2),aaL(3)
    real(8),allocatable :: rrrr(:,:),rho_tmp(:),zeta(:),nab(:)
    logical :: flag_alloc

    n1=idisp(myrank)+1
    n2=idisp(myrank)+ircnt(myrank)

    pi=acos(-1.d0)

    E_exchange=0.d0
    E_correlation=0.d0
    Exc=0.d0
    Vxc(:,:)=0.d0
    flag_alloc=.false.

    ML  = Ngrid(0)
    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

    Mx=ML1+Md
    My=ML2+Md
    Mz=ML3+Md

    a1b=Igrid(1,1)
    b1b=Igrid(2,1)
    a2b=Igrid(1,2)
    b2b=Igrid(2,2)
    a3b=Igrid(1,3)
    b3b=Igrid(2,3)

    H1 = Hgrid(1)
    H2 = Hgrid(2)
    H3 = Hgrid(3)

    aaL(1) = sqrt( sum(aa(1:3,1)**2) )
    aaL(2) = sqrt( sum(aa(1:3,2)**2) )
    aaL(3) = sqrt( sum(aa(1:3,3)**2) )

    b(:,:)=0.d0
    b(1:3,1)=aaL(1)*bb(1:3,1)/(2.d0*Pi)/H1
    b(1:3,2)=aaL(2)*bb(1:3,2)/(2.d0*Pi)/H2
    b(1:3,3)=aaL(3)*bb(1:3,3)/(2.d0*Pi)/H3

    allocate( nab(-Md:Md) ) ; nab=0.d0
    call get_coef_nabla_fd(Md,nab)

! negative value check >>>>>>>>>>
    sbf(:)=0.d0
    do ispin=1,nspin
       do i=n1,n2
          if ( rho(i,ispin) < 0.d0 ) sbf(ispin)=sbf(ispin)+rho(i,ispin)
       end do
    end do
    if ( flag_pcc_0 ) then
       do i=n1,n2
          if ( rhoc(i) < 0.d0 ) sbf(3)=sbf(3)+rhoc(i)
       end do
    end if
    call mpi_allreduce(sbf,rbf,3,MPI_REAL8,MPI_SUM,comm_grid,ierr)
    if ( disp_switch_parallel ) then
       write(*,'(1x,"negative charge =",3g20.10)') rbf(1:3)*dV
    end if
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    select case(SYStype)
    case default
       allocate( LLL2(0:ML1-1,0:ML2-1,0:ML3-1) ) ; LLL2=0
       i=0
       irank=-1
       do i3=1,node_partition(3)
       do i2=1,node_partition(2)
       do i1=1,node_partition(1)
          irank=irank+1
          l1=pinfo_grid(1,irank) ; m1=pinfo_grid(2,irank)+l1-1
          l2=pinfo_grid(3,irank) ; m2=pinfo_grid(4,irank)+l2-1
          l3=pinfo_grid(5,irank) ; m3=pinfo_grid(6,irank)+l3-1
          do j3=l3,m3
          do j2=l2,m2
          do j1=l1,m1
             i=i+1
             LLL2(j1,j2,j3)=i
          end do
          end do
          end do
       end do
       end do
       end do
    case(1,2)
       stop "stop@GGAPBE96,MOL"
       allocate( LL2(3,ML) ) ; LL2=0
!       call Make_GridMap_1(LL2,1,ML)
       allocate( LLL2(-Mx:Mx,-My:My,-Mz:Mz) ) ; LLL2=0
!       call Make_GridMap_3(LLL2,-Mx,Mx,-My,My,-Mz,Mz)
    end select

    allocate( rrrr(ML,3) ) ; rrrr=0.d0
    allocate( gx(n1:n2),gy(n1:n2),gz(n1:n2) ) ; gx=0.d0 ; gy=0.d0 ; gz=0.d0
    allocate( rtmp(n1:n2) ) ; rtmp=0.d0
    allocate( rho_tmp(n1:n2) ) ; rho_tmp=0.d0
    allocate( zeta(n1:n2) ) ; zeta=0.d0

    do ispin=1,nspin
       rho_tmp(n1:n2)=rho_tmp(n1:n2)+rho(n1:n2,ispin)
    end do
    if ( flag_pcc_0 ) then
       rho_tmp(n1:n2)=rho_tmp(n1:n2)+rhoc(n1:n2)
    end if
    where( rho_tmp < 0.d0 )
       rho_tmp=0.d0
    end where
    zeta(n1:n2)=rho(n1:n2,1)-rho(n1:n2,nspin)
    j=0
    do i=n1,n2
       if ( rho_tmp(i)==0.d0 ) then
          zeta(i)=0.d0
       else
          zeta(i)=zeta(i)/rho_tmp(i)
       end if
       if ( zeta(i)>1.d0 .or. zeta(i)<-1.d0 ) then
          j=j+1
          if ( disp_switch_parallel ) write(*,*) j,zeta(i),rho(i,1:nspin)
       end if
    end do
    where( zeta >  1.d0 )
       zeta= 1.d0
    end where
    where( zeta < -1.d0 )
       zeta=-1.d0
    end where

    www(:,:,:,:)=0.d0
    select case(SYStype)
    case default
       i=n1-1
       do i3=a3b,b3b
       do i2=a2b,b2b
       do i1=a1b,b1b
          i=i+1
          www(i1,i2,i3,1)=rho_tmp(i)
       end do
       end do
       end do
    case(1,2)
       do i=n1,n2
          www(LL2(1,i),LL2(2,i),LL2(3,i),1)=rho_tmp(i)
       end do
    end select

    call bcset(1,1,Md,0)

    select case(SYStype)
    case default
       i=n1-1
       do i3=a3b,b3b
       do i2=a2b,b2b
       do i1=a1b,b1b
          i=i+1
          g1=0.d0 ; g2=0.d0 ; g3=0.d0
          do m=1,Md
             g1=g1-nab(m)*(www(i1-m,i2,i3,1)-www(i1+m,i2,i3,1))
             g2=g2-nab(m)*(www(i1,i2-m,i3,1)-www(i1,i2+m,i3,1))
             g3=g3-nab(m)*(www(i1,i2,i3-m,1)-www(i1,i2,i3+m,1))
          end do
          gx(i)=b(1,1)*g1+b(1,2)*g2+b(1,3)*g3
          gy(i)=b(2,1)*g1+b(2,2)*g2+b(2,3)*g3
          gz(i)=b(3,1)*g1+b(3,2)*g2+b(3,3)*g3
       end do
       end do
       end do
    case(1,2)
       do i=n1,n2
          i1=LL2(1,i) ; i2=LL2(2,i) ; i3=LL2(3,i)
          g1=0.d0 ; g2=0.d0 ; g3=0.d0
          do m=1,Md
             g1=g1-nab(m)*(www(i1-m,i2,i3,1)-www(i1+m,i2,i3,1))
             g2=g2-nab(m)*(www(i1,i2-m,i3,1)-www(i1,i2+m,i3,1))
             g3=g3-nab(m)*(www(i1,i2,i3-m,1)-www(i1,i2,i3+m,1))
          end do
          gx(i)=b(1,1)*g1+b(1,2)*g2+b(1,3)*g3
          gy(i)=b(2,1)*g1+b(2,2)*g2+b(2,3)*g3
          gz(i)=b(3,1)*g1+b(3,2)*g2+b(3,3)*g3
       end do
    end select

!
! --- Exchange ---
!

    do ispin=MSP_0,MSP_1

       do i=n1,n2

          trho=dble(nspin)*rho(i,ispin) ; if (allocated(rhoc)) trho=trho+rhoc(i)

!          if ( trho==0.d0 ) cycle
          if ( trho <= zero_density ) cycle

          kf=(3.d0*Pi*Pi*trho)**(1.d0/3.d0)

          ex_lda=-3.d0/(4.d0*Pi)*kf
          vx_lda=-1.d0/Pi*kf

          g2=gx(i)*gx(i)+gy(i)*gy(i)+gz(i)*gz(i)

          Fx=1.d0+Kp-4.d0*Kp*Kp*(trho*kf)**2/(4.d0*Kp*(trho*kf)**2+mu*g2)

          E_exchange=E_exchange+trho*ex_lda*Fx

          Vxc(i,ispin)=Vxc(i,ispin)+Fx*vx_lda+(24.d0*Pi*Kp*Kp*mu*trho**3*g2)/(4.d0*Kp*(trho*kf)**2+mu*g2)**2

          rtmp(i)=-18.d0*Pi*Kp*Kp*mu*trho**4/(4.d0*Kp*(trho*kf)**2+mu*g2)**2

       end do ! i

       rrrr(n1:n2,1)=rtmp(n1:n2)*gx(n1:n2)
       call mpi_allgatherv(rrrr(n1,1),ir_grid(myrank_g),mpi_real8,rrrr(1,1),ir_grid,id_grid,mpi_real8,comm_grid,ierr)
       rrrr(n1:n2,2)=rtmp(n1:n2)*gy(n1:n2)
       call mpi_allgatherv(rrrr(n1,2),ir_grid(myrank_g),mpi_real8,rrrr(1,2),ir_grid,id_grid,mpi_real8,comm_grid,ierr)
       rrrr(n1:n2,3)=rtmp(n1:n2)*gz(n1:n2)
       call mpi_allgatherv(rrrr(n1,3),ir_grid(myrank_g),mpi_real8,rrrr(1,3),ir_grid,id_grid,mpi_real8,comm_grid,ierr)

       select case(SYStype)
       case default
          do i3=0,ML3-1
          do i2=0,ML2-1
          do i1=0,ML1-1
             i=LLL2(i1,i2,i3)
             do m=-Md,Md
                j1=i1+m
                k1=j1/ML1 ; if ( j1<0 ) k1=(j1+1)/ML1-1
                j1=j1-k1*ML1
                j=LLL2(j1,i2,i3)
                if ( n1<=j .and. j<=n2 ) then
                   Vxc(j,ispin)=Vxc(j,ispin)+nab(m)*sign(1,m)*( rrrr(i,1)*b(1,1)+rrrr(i,2)*b(2,1)+rrrr(i,3)*b(3,1) )
                end if
                j2=i2+m
                k2=j2/ML2 ; if ( j2<0 ) k2=(j2+1)/ML2-1
                j2=j2-k2*ML2
                j=LLL2(i1,j2,i3)
                if ( n1<=j .and. j<=n2 ) then
                   Vxc(j,ispin)=Vxc(j,ispin)+nab(m)*sign(1,m)*( rrrr(i,1)*b(1,2)+rrrr(i,2)*b(2,2)+rrrr(i,3)*b(3,2) )
                end if
                j3=i3+m
                k3=j3/ML3 ; if ( j3<0 ) k3=(j3+1)/ML3-1
                j3=j3-k3*ML3
                j=LLL2(i1,i2,j3)
                if ( n1<=j .and. j<=n2 ) then
                   Vxc(j,ispin)=Vxc(j,ispin)+nab(m)*sign(1,m)*( rrrr(i,1)*b(1,3)+rrrr(i,2)*b(2,3)+rrrr(i,3)*b(3,3) )
                end if
             end do
          end do
          end do
          end do
       case(1,2)
          do i=1,ML
             i1=LL2(1,i) ; i2=LL2(2,i) ; i3=LL2(3,i)
             do m=-Md,Md
                j=LLL2(i1+m,i2,i3)
                if ( n1<=j .and. j<=n2 ) then
                   Vxc(j,ispin)=Vxc(j,ispin)+nab(m)*sign(1,m)*( rrrr(i,1)*b(1,1)+rrrr(i,2)*b(2,1)+rrrr(i,3)*b(3,1) )
                end if
                j=LLL2(i1,i2+m,i3)
                if ( n1<=j .and. j<=n2 ) then
                   Vxc(j,ispin)=Vxc(j,ispin)+nab(m)*sign(1,m)*( rrrr(i,1)*b(1,2)+rrrr(i,2)*b(2,2)+rrrr(i,3)*b(3,2) )
                end if
                j=LLL2(i1,i2,i3+m)
                if ( n1<=j .and. j<=n2 ) then
                   Vxc(j,ispin)=Vxc(j,ispin)+nab(m)*sign(1,m)*( rrrr(i,1)*b(1,3)+rrrr(i,2)*b(2,3)+rrrr(i,3)*b(3,3) )
                end if
             end do
          end do
       end select

    end do ! ispin

!
! --- Correlation ---
!
    const1=2.d0**(4.d0/3.d0)-2.d0
    const2=9.d0*(2.d0**(1.d0/3.d0)-1.d0)/4.d0

    do i=n1,n2

!       trho=rho(i,ispin) ; if ( allocated(rhoc) ) trho=trho+rhoc(i)
       trho=rho_tmp(i)

!       if ( trho==0.d0 ) cycle
       if ( trho <= zero_density ) cycle

       fz=( (1.d0+zeta(i))**(4.d0/3.d0)+(1.d0-zeta(i))**(4.d0/3.d0)-2.d0 )*const1

       kf=(3.d0*Pi*Pi*trho)**(1.d0/3.d0)

!       rs=( 3.d0/(4.d0*Pi*trho) )**(1.d0/3.d0)
       rs=( 3.d0/(4.d0*Pi*abs(trho)) )**(1.d0/3.d0)

       ec_U=-2.d0*A00*(1.d0+alp10*rs)*log( 1.d0+1.d0/(2.d0*A00*(bt10*sqrt(rs)+bt20*rs+bt30*rs**(3.d0/2.d0)+bt40*rs*rs)) )
       ec_P=-2.d0*A01*(1.d0+alp11*rs)*log( 1.d0+1.d0/(2.d0*A01*(bt11*sqrt(rs)+bt21*rs+bt31*rs**(3.d0/2.d0)+bt41*rs*rs)) )
       alpc=-2.d0*A02*(1.d0+alp12*rs)*log( 1.d0+1.d0/(2.d0*A02*(bt12*sqrt(rs)+bt22*rs+bt32*rs**(3.d0/2.d0)+bt42*rs*rs)) )

       ec_lda=ec_U-alpc*fz*const2*(1.d0-zeta(i)**4)+(ec_P-ec_U)*fz*zeta(i)**4

       phi=0.5d0*( (1.d0+zeta(i))**(2.d0/3.d0)+(1.d0-zeta(i))**(2.d0/3.d0) )

       if ( trho==0.d0 ) then
          A=0.d0
          T=0.d0
          Hs=0.d0
       else
          T=(gx(i)*gx(i)+gy(i)*gy(i)+gz(i)*gz(i))*Pi/(16.d0*phi**2*kf*trho**2)
!          A=C1/(exp(-ec_lda/(C2*phi**3))-1.d0)
!          Hs=C2*phi**3*log( 1.d0+C1*(T+A*T*T)/(1.d0+A*T+A*A*T*T) )
          Ai=( exp(-ec_lda/(C2*phi**3)) - 1.d0 )/C1
          Hs=C2*phi**3*log( 1.d0+C1*( Ai*Ai*T + Ai*T*T )/( Ai*Ai + Ai*T + T*T ) )
!          tmp=exp(-ec_lda/(C2*phi**3))-1.d0
!          Hs=C2*phi**3*log( 1.d0+C1*(tmp*tmp/T+C1*tmp)/(tmp*tmp/(T*T)+C1*tmp/T+C1*C1) )
       end if

       E_correlation=E_correlation+trho*(ec_lda+Hs)

       deU_dn=-4.d0*Pi/9.d0*rs**4*alp10*ec_U/(1.d0+alp10*rs) &
            -4.d0*Pi/9.d0*rs*rs*(1.d0+alp10*rs)*(0.5d0*bt10*sqrt(rs)+bt20*rs+1.5d0*bt30*rs*sqrt(rs)+2.d0*bt40*rs*rs) &
            /(bt10+bt20*sqrt(rs)+bt30*rs+bt40*rs*sqrt(rs))**2 * exp(ec_U/(2.d0*A00*(1.d0+alp10*rs)))
       deP_dn=-4.d0*Pi/9.d0*rs**4*alp11*ec_P/(1.d0+alp11*rs) &
            -4.d0*Pi/9.d0*rs*rs*(1.d0+alp11*rs)*(0.5d0*bt11*sqrt(rs)+bt21*rs+1.5d0*bt31*rs*sqrt(rs)+2.d0*bt41*rs*rs) &
            /(bt11+bt21*sqrt(rs)+bt31*rs+bt41*rs*sqrt(rs))**2 * exp(ec_P/(2.d0*A01*(1.d0+alp11*rs)))
       dac_dn=-4.d0*Pi/9.d0*rs**4*alp12*alpc/(1.d0+alp12*rs) &
            -4.d0*Pi/9.d0*rs*rs*(1.d0+alp12*rs)*(0.5d0*bt12*sqrt(rs)+bt22*rs+1.5d0*bt32*rs*sqrt(rs)+2.d0*bt42*rs*rs) &
            /(bt12+bt22*sqrt(rs)+bt32*rs+bt42*rs*sqrt(rs))**2 * exp(alpc/(2.d0*A02*(1.d0+alp12*rs)))

       dfz_dz=4.d0/3.d0*( (1.d0+zeta(i))**(4.d0/3.d0)-(1.d0-zeta(i))**(4.d0/3.d0) )*const1

       dec_dz=-alpc*dfz_dz*const2*(1.d0-zeta(i)**4)+4.d0*alpc*fz*const2*zeta(i)**3 &
            +(ec_P-ec_U)*dfz_dz*zeta(i)**4+(ec_P-ec_U)*fz*4.d0*zeta(i)**3

       if ( zeta(i) == 1.d0 .or. zeta(i) == -1.d0 ) then
          dphi_dz = 0.d0
       else
          dphi_dz=( (1.d0+zeta(i))**(-1.d0/3.d0)-(1.d0-zeta(i))**(-1.d0/3.d0) )/3.d0
       end if

       if ( trho==0.d0 ) then
          dH_dA=0.d0
          dH_dT=0.d0
       else
!          tmp   = 1.d0 + A*T + A*A*T*T
!          dH_dA = -phi**3*C1*C2*A*T**3*(2.d0+A*T)/( tmp*tmp+C1*T*(1.d0+A*T)*tmp )
!          dH_dT =  phi**3*C1*C2*(1.d0+2.d0*A*T)/( tmp*tmp + C1*T*(1.d0+A*T)*tmp )
          tmp   = Ai*Ai + Ai*T + T*T
          dH_dA = -phi**3*C1*C2*Ai*T**3*(2.d0*Ai+T)/( tmp*tmp + C1*T*(Ai*Ai+Ai*T)*tmp ) ! Ai**2 is canceld in dA_dn
          dH_dT =  phi**3*C1*C2*Ai**3*(Ai+2.d0*T)/( tmp*tmp + C1*T*(Ai*Ai+Ai*T)*tmp )
!          tmp1=(exp(-ec_lda/(C2*phi**3))-1.d0)/C1
!          tmp2=tmp1/T
!          dH_dA=-phi**3*C1*C2*tmp1*(2.d0*tmp2+1.d0)/((tmp2*tmp2+tmp2+1.d0)**2+C1*tmp1*(tmp2+1.d0)*(tmp2*tmp2+tmp2+1.d0))
!          dH_dT=phi**3*C1*C2*tmp2**3*(tmp2+2.d0)/((tmp2*tmp2+tmp2+1.d0)**2+tmp1*C1*(tmp2+1.d0)*(tmp2*tmp2+tmp2+1.d0))
       end if

       dH_dphi=3.d0*Hs/phi

       srho(1)    =rho(i,1)     ; if ( flag_pcc_0 ) srho(1)=srho(1)+rhoc(i)/dble(nspin)
       srho(nspin)=rho(i,nspin) ; if ( flag_pcc_0 ) srho(nspin)=srho(nspin)+rhoc(i)/dble(nspin)

       dz_dn(1)    = 2.d0*srho(nspin)/trho
       dz_dn(nspin)=-2.d0*srho(1)/trho

       do ispin=MSP_0,MSP_1

          dec_dn=deU_dn-dac_dn*fz*const2*(1.d0-zeta(i)**4)+(deP_dn-deU_dn)*fz*zeta(i)**4 + dec_dz*dz_dn(ispin)

!          dA_dn=A*(C1+A)/(C1*C2*phi**3)*( dec_dn - 3.d0*ec_lda/phi*dphi_dz*dz_dn(ispin) )
          tmp=exp(-ec_lda/(phi**3*C2))
          dA_dn=tmp/(C1*C2*phi**3)*( dec_dn - 3.d0*ec_lda/phi*dphi_dz*dz_dn(ispin) )

          Vxc(i,ispin)=Vxc(i,ispin) + ec_lda+Hs + trho*dec_dn + trho*dH_dA*dA_dn - 7.d0*T/3.d0*dH_dT &
               + trho*dH_dphi*dphi_dz*dz_dn(ispin)

       end do ! ispin

       rtmp(i)=dH_dT*Pi/(8.d0*kf*trho)

    end do ! i

    rrrr(n1:n2,1)=rtmp(n1:n2)*gx(n1:n2)
    call mpi_allgatherv(rrrr(n1,1),ir_grid(myrank_g),mpi_real8,rrrr(1,1),ir_grid,id_grid,mpi_real8,comm_grid,ierr)
    rrrr(n1:n2,2)=rtmp(n1:n2)*gy(n1:n2)
    call mpi_allgatherv(rrrr(n1,2),ir_grid(myrank_g),mpi_real8,rrrr(1,2),ir_grid,id_grid,mpi_real8,comm_grid,ierr)
    rrrr(n1:n2,3)=rtmp(n1:n2)*gz(n1:n2)
    call mpi_allgatherv(rrrr(n1,3),ir_grid(myrank_g),mpi_real8,rrrr(1,3),ir_grid,id_grid,mpi_real8,comm_grid,ierr)

    select case(SYStype)
    case default
       do i3=0,ML3-1
       do i2=0,ML2-1
       do i1=0,ML1-1
          i=LLL2(i1,i2,i3)
          do m=-Md,Md
             j1=i1+m
             k1=j1/ML1 ; if ( j1<0 ) k1=(j1+1)/ML1-1
             j1=j1-k1*ML1
             j=LLL2(j1,i2,i3)
             if ( n1<=j .and. j<=n2 ) then
                do ispin=MSP_0,MSP_1
                   Vxc(j,ispin)=Vxc(j,ispin)+nab(m)*sign(1,m)*( rrrr(i,1)*b(1,1)+rrrr(i,2)*b(2,1)+rrrr(i,3)*b(3,1) )
                end do
             end if
             j2=i2+m
             k2=j2/ML2 ; if ( j2<0 ) k2=(j2+1)/ML2-1
             j2=j2-k2*ML2
             j=LLL2(i1,j2,i3)
             if ( n1<=j .and. j<=n2 ) then
                do ispin=MSP_0,MSP_1
                   Vxc(j,ispin)=Vxc(j,ispin)+nab(m)*sign(1,m)*( rrrr(i,1)*b(1,2)+rrrr(i,2)*b(2,2)+rrrr(i,3)*b(3,2) )
                end do
             end if
             j3=i3+m
             k3=j3/ML3 ; if ( j3<0 ) k3=(j3+1)/ML3-1
             j3=j3-k3*ML3
             j=LLL2(i1,i2,j3)
             if ( n1<=j .and. j<=n2 ) then
                do ispin=MSP_0,MSP_1
                   Vxc(j,ispin)=Vxc(j,ispin)+nab(m)*sign(1,m)*( rrrr(i,1)*b(1,3)+rrrr(i,2)*b(2,3)+rrrr(i,3)*b(3,3) )
                end do
             end if
          end do
       end do
       end do
       end do
    case(1,2)
       do i=1,ML
          i1=LL2(1,i) ; i2=LL2(2,i) ; i3=LL2(3,i)
          do m=-Md,Md
             j=LLL2(i1+m,i2,i3)
             if ( n1<=j .and. j<=n2 ) then
                do ispin=1,MSP_0,MSP_1
                   Vxc(j,ispin)=Vxc(j,ispin)+nab(m)*sign(1,m)*( rrrr(i,1)*b(1,1)+rrrr(i,2)*b(2,1)+rrrr(i,3)*b(3,1) )
                end do
             end if
             j=LLL2(i1,i2+m,i3)
             if ( n1<=j .and. j<=n2 ) then
                do ispin=MSP_0,MSP_1
                   Vxc(j,ispin)=Vxc(j,ispin)+nab(m)*sign(1,m)*( rrrr(i,1)*b(1,2)+rrrr(i,2)*b(2,2)+rrrr(i,3)*b(3,2) )
                end do
             end if
             j=LLL2(i1,i2,i3+m)
             if ( n1<=j .and. j<=n2 ) then
                do ispin=MSP_0,MSP_1
                   Vxc(j,ispin)=Vxc(j,ispin)+nab(m)*sign(1,m)*( rrrr(i,1)*b(1,3)+rrrr(i,2)*b(2,3)+rrrr(i,3)*b(3,3) )
                end do
             end if
          end do
       end do
    end select

    sbf(1)=E_exchange*dV/dble(nspin)
    sbf(2)=E_correlation*dV
    call mpi_allreduce(sbf,rbf,2,mpi_real8,mpi_sum,comm_grid,ierr)
    E_exchange=rbf(1)
    E_correlation=rbf(2)
    Exc=E_exchange+E_correlation

    deallocate( zeta )
    deallocate( rho_tmp )
    deallocate( rtmp )
    deallocate( gz,gy,gx )
    deallocate( rrrr )
    deallocate( LLL2 )
    deallocate( nab )
    select case(SYStype)
    case(1,2)
       deallocate(LL2)
    end select

    return
  END SUBROUTINE CALC_GGAPBE96

END MODULE xc_module
