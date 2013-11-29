MODULE ps_local_module

  use rgrid_module
  use ggrid_module
  use atom_module
  use strfac_module
  use bb_module
  use density_module
  use electron_module
  use pseudopot_module
  use parallel_module
  use watch_module
  use ps_local_gth_module
  use simc_module

  implicit none

  PRIVATE
  PUBLIC :: Vion,init_ps_local,construct_ps_local,calc_force_ps_local &
       ,construct_ps_local_ffte,calc_force_ps_local_ffte

  real(8),allocatable :: Rcloc(:),vqlg(:,:),vqlgl(:,:),vqls(:,:)
  integer,allocatable :: NRcloc(:)
  real(8),allocatable :: Vion(:)

  INTERFACE
     FUNCTION bberf(x)
       real(8) :: bberf,x
     END FUNCTION bberf
  END INTERFACE

CONTAINS


  SUBROUTINE init_ps_local
    integer :: i,ig,ik,iorb,MMr,NRc,MKI
    real(8) :: Rc,p1,p2,p3,p4,vlong,Pi,const,x,r,sb,sum0,G,G2
    real(8) :: Vcell
    real(8),allocatable :: vshort(:),tmp(:)

    MKI   = Nelement
    Vcell = Ngrid(0)*dV
    Pi    = acos(-1.d0)
    const = 4.d0*Pi/Vcell

    allocate( vqlg(NMGL,MKI)  ) ; vqlg=0.d0

    if ( pselect == 4 .or. pselect == 5 ) then
       call init_ps_local_gth(NMGL,Nelement,GG,vqlg)
       return
    end if

    MMr=maxval(Mr)
    allocate( vqls(MMr,MKI)   ) ; vqls=0.d0
    allocate( vqlgl(NMGL,MKI) ) ; vqlgl=0.d0
    allocate( NRcloc(MKI)     ) ; NRcloc=0
    allocate( Rcloc(MKI)      ) ; Rcloc=0.d0

    allocate( vshort(MMr) )

    do ik=1,MKI

       MMr = Mr(ik)

       Rc=0.d0
       NRc=0
       do iorb=1,norb(ik)
          Rc=max( Rc, Rps(iorb,ik) )
          NRc=max( NRc, NRps(iorb,ik) )
       end do

       if ( Rc<1.d-8 ) Rc=5.d0
       if ( NRc==0 ) then
          do i=1,MMr
             if ( rad(i,ik)>Rc ) then
                NRc=i
                exit
             end if
          end do
       end if

       call simc(rad(1,ik),vql(1,ik),Rc,Zps(ik),parloc(1,ik),MMr)

       p1=parloc(1,ik) ; p2=sqrt(parloc(2,ik))
       p3=parloc(3,ik) ; p4=sqrt(parloc(4,ik))

       do i=1,MMr
          r=rad(i,ik)
          if ( r<1.d-9 ) then
             vlong=-2.d0*Zps(ik)/sqrt(Pi)*(p1*p2+p3*p4)
          else
             vlong=-Zps(ik)/r*( p1*bberf(p2*r)+p3*bberf(p4*r) )
          end if
          vshort(i)=vql(i,ik)-vlong
          vqls(i,ik)=vql(i,ik)-vlong
          if( NRcloc(ik)==0 )then
             if( abs(vshort(i))<1.d-8 ) then
                NRcloc(ik)=i
                Rcloc(ik)=r
             end if
          end if
       end do

       if ( NRcloc(ik)==0 ) then
          Rcloc(ik)=Rc
          NRcloc(ik)=NRc
       end if

       allocate( tmp(MMr) )

       do ig=1,NMGL
          G=sqrt(GG(ig))
          if ( G == 0.d0 ) then
             do i=1,MMr
                tmp(i)=rad(i,ik)*rad(i,ik)*vshort(i)*rab(i,ik)
             end do
          else
             do i=1,MMr
                x=G*rad(i,ik)
                if ( x<1.d-1 ) then
                   sb=-(1.d0/39916800.d0*x**10-1.d0/362880.d0*x**8 &
                       +1.d0/5040.d0*x**6-1.d0/120.d0*x**4+1.d0/6.d0*x**2-1.d0)
                else
                   sb=sin(x)/x
                end if
                tmp(i)=rad(i,ik)*rad(i,ik)*vshort(i)*sb*rab(i,ik)
             end do
          end if
          call simp(tmp,sum0,2)
          vqlg(ig,ik)=sum0*const
       end do

       p1=-Zps(ik)*parloc(1,ik) ; p2=0.25d0/parloc(2,ik)
       p3=-Zps(ik)*parloc(3,ik) ; p4=0.25d0/parloc(4,ik)
       do ig=1,NMGL
          G2=GG(ig)
          if ( G2 == 0.d0 ) then
             vqlgl(ig,ik)=-(p1*p2+p3*p4)*const
             vqlg(ig,ik)=vqlg(ig,ik)+vqlgl(ig,ik)
          else
             vqlgl(ig,ik)=(p1*exp(-G2*p2)+p3*exp(-G2*p4))/G2*const
             vqlg(ig,ik)=vqlg(ig,ik)+vqlgl(ig,ik)
          end if
       end do
       deallocate( tmp )

    end do ! ik

    deallocate( vshort )

  END SUBROUTINE init_ps_local

  SUBROUTINE simp(f,s,m)
    integer,intent(IN)  :: m
    real(8),intent(IN)  :: f(:)
    real(8),intent(OUT) :: s
    real(8),allocatable :: g(:)
    integer :: i,n,nn,nmax
    n=size(f) ; nmax=int(n/m)*m
    do i=0,m
       nmax=nmax+i ; if ( nmax>=n ) exit
    end do
    allocate( g(nmax) ) ; g(1:n)=f ; if ( nmax>n ) g(n+1:)=0.d0
    select case(m)
    case default
       s = 0.5d0*(f(1)+f(n)) + sum(f(2:n-1))
    case(2)
       s=0.d0
       do i=1,nmax-2,2
          s = s + g(i) + 4.d0*g(i+1) + g(i+2)
       end do
       s=s/3.d0
    case(4)
       s=0.d0
       do i=1,nmax-4,4
          s=s+7*g(i)+32*g(i+1)+12*g(i+2)+32*g(i+3)+7*g(i+4)
       end do
       s=s*2.d0/45.d0
    case(6)
       s=0.d0
       do i=1,nmax-6,6
          s=s+41*g(i)+216*g(i+1)+27*g(i+2)+272*g(i+3) &
               +27*g(i+4)+216*g(i+5)+41*g(i+6)
       end do
       s=s/140.d0
    end select
    deallocate( g )
    return
  END SUBROUTINE simp


  SUBROUTINE construct_ps_local

    integer :: a,i,i1,i2,i3,ik,j,MG
    integer :: ML1,ML2,ML3,ML,ML_0,ML_1
    integer :: ifacx(30),ifacy(30),ifacz(30)
    integer,allocatable :: lx1(:),lx2(:),ly1(:),ly2(:),lz1(:),lz2(:)
    complex(8),allocatable :: fftwork(:),zwork(:,:,:),vg(:)
    complex(8),allocatable :: wsavex(:),wsavey(:),wsavez(:)
    real(8) :: ctt(0:3),ett(0:3)

    MG  = NGgrid(0)
    ML  = Ngrid(0)
    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)
    ML_0= Igrid(1,0)
    ML_1= Igrid(2,0)

    ctt(:)=0.d0
    ett(:)=0.d0

    call watch(ctt(0),ett(0))

    if ( .not.allocated(Vion) ) then
       allocate( Vion(ML_0:ML_1) )
    end if

    allocate( zwork(0:ML1-1,0:ML2-1,0:ML3-1) )

    allocate( vg(MG) )

    do i=MG_0,MG_1
       j=MGL(i)
       vg(i)=vqlg(j,1)*SGK(i,1)
    end do
    do ik=2,Nelement
       do i=MG_0,MG_1
          j=MGL(i)
          vg(i)=vg(i)+vqlg(j,ik)*SGK(i,ik)
       end do
    end do
    call allgatherv_Ggrid(vg)

    call construct_Ggrid(2)

    zwork(:,:,:)=(0.d0,0.d0)
    do i=1,NGgrid(0)
       zwork(LLG(1,i),LLG(2,i),LLG(3,i))=vg(i)
    end do

    call destruct_Ggrid

    deallocate( vg )

    allocate( fftwork(ML) )
    allocate( lx1(ML),lx2(ML),ly1(ML),ly2(ML),lz1(ML),lz2(ML) )
    allocate( wsavex(ML1),wsavey(ML2),wsavez(ML3) )

    call prefft(ML1,ML2,ML3,ML,wsavex,wsavey,wsavez &
         ,ifacx,ifacy,ifacz,lx1,lx2,ly1,ly2,lz1,lz2)

    call watch(ctt(1),ett(1))

    call fft3bx(ML1,ML2,ML3,ML,zwork,fftwork,wsavex,wsavey,wsavez &
         ,ifacx,ifacy,ifacz,lx1,lx2,ly1,ly2,lz1,lz2)

    call watch(ctt(2),ett(2))

    i=ML_0-1
    do i3=Igrid(1,3),Igrid(2,3)
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
       i=i+1
       Vion(i)=real( zwork(i1,i2,i3) )
    end do
    end do
    end do

    deallocate( wsavez,wsavey,wsavex )
    deallocate( lz2,lz1,ly2,ly1,lx2,lx1 )
    deallocate( fftwork )
    deallocate( zwork )

    call watch(ctt(3),ett(3))

    if ( disp_switch_parallel ) then
       write(*,*) "time(const_ps_loc_1)",ctt(1)-ctt(0),ett(1)-ett(0)
       write(*,*) "time(const_ps_loc_2)",ctt(2)-ctt(1),ett(2)-ett(1)
       write(*,*) "time(const_ps_loc_3)",ctt(3)-ctt(2),ett(3)-ett(2)
    end if

  END SUBROUTINE construct_ps_local


  SUBROUTINE calc_force_ps_local(MI,force)
    integer,intent(IN) :: MI
    real(8),intent(OUT) :: force(3,MI)
    integer :: a,i,j,ik,ispin,irank,i1,i2,i3,n,l1,l2,l3,m1,m2,m3,j1,j2,j3
    integer :: MI_0,MI_1,ML1,ML2,ML3,ML_0,ML_1,ML,MG,N_MI,ierr
    integer,allocatable :: icnt(:),idis(:)
    real(8) :: a1,a2,a3,pi2,Gr,Gx,Gy,Gz,Vcell
    real(8),allocatable :: work(:)
    integer :: ifacx(30),ifacy(30),ifacz(30)
    integer,allocatable :: lx1(:),lx2(:),ly1(:),ly2(:),lz1(:),lz2(:)
    complex(8),allocatable :: fftwork(:)
    complex(8),allocatable :: zrho3(:,:,:),zrho(:)
    complex(8),allocatable :: wsavex(:),wsavey(:),wsavez(:)
    complex(8),parameter :: z0=(0.d0,0.d0)
    complex(8) :: zsum1,zsum2,zsum3,ztmp

    force(:,:)=0.d0

    MG    = NGgrid(0)
    ML    = Ngrid(0)
    ML1   = Ngrid(1)
    ML2   = Ngrid(2)
    ML3   = Ngrid(3)
    ML_0  = Igrid(1,0)
    ML_1  = Igrid(2,0)
    pi2   = 2.d0*acos(-1.d0)
    Vcell = ML*dV

    allocate( icnt(0:nprocs-1) )
    allocate( idis(0:nprocs-1) )
    N_MI = Natom/nprocs
    icnt(0:nprocs-1) = N_MI
    n = Natom - N_MI*nprocs
    if ( n>0 ) then
       do irank=0,n-1
          icnt(irank)=icnt(irank)+1
       end do
    end if
    do irank=0,nprocs-1
       idis(irank) = sum( icnt(0:irank) ) - icnt(irank)
    end do
    MI_0 = idis(myrank)+1
    MI_1 = idis(myrank)+icnt(myrank)
    deallocate( idis, icnt )

    allocate( zrho3(0:ML1-1,0:ML2-1,0:ML3-1) )

    allocate( work(ML) )
    work(1:ML)=0.d0
    do ispin=1,Nspin
       work(ML_0:ML_1)=work(ML_0:ML_1)+rho(ML_0:ML_1,ispin)
    end do
    call mpi_allgatherv(work(ML_0),ML_1-ML_0+1,mpi_real8,work &
         ,ir_grid,id_grid,mpi_real8,comm_grid,ierr)
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
          zrho3(j1,j2,j3)=work(i)
       end do
       end do
       end do
    end do
    end do
    end do
    deallocate( work )

    allocate( fftwork(ML) )
    allocate( lx1(ML),lx2(ML),ly1(ML),ly2(ML),lz1(ML),lz2(ML) )
    allocate( wsavex(ML1),wsavey(ML2),wsavez(ML3) )

    call prefft(ML1,ML2,ML3,ML,wsavex,wsavey,wsavez &
         ,ifacx,ifacy,ifacz,lx1,lx2,ly1,ly2,lz1,lz2)

    call fft3fx(ML1,ML2,ML3,ML,zrho3,fftwork,wsavex,wsavey,wsavez &
               ,ifacx,ifacy,ifacz,lx1,lx2,ly1,ly2,lz1,lz2)

    call construct_Ggrid(0)

    allocate( zrho(MG) )

    do i=1,NGgrid(0)
       i1=mod(ML1+LLG(1,i),ML1)
       i2=mod(ML2+LLG(2,i),ML2)
       i3=mod(ML3+LLG(3,i),ML3)
       zrho(i) = conjg( zrho3(i1,i2,i3) )
    end do

    do a=MI_0,MI_1

       ik=ki_atom(a)
       a1=pi2*aa_atom(1,a)
       a2=pi2*aa_atom(2,a)
       a3=pi2*aa_atom(3,a)

       zsum1=z0
       zsum2=z0
       zsum3=z0
       do i=1,NGgrid(0)
          Gx=bb(1,1)*LLG(1,i)+bb(1,2)*LLG(2,i)+bb(1,3)*LLG(3,i)
          Gy=bb(2,1)*LLG(1,i)+bb(2,2)*LLG(2,i)+bb(2,3)*LLG(3,i)
          Gz=bb(3,1)*LLG(1,i)+bb(3,2)*LLG(2,i)+bb(3,3)*LLG(3,i)
          Gr=a1*LLG(1,i)+a2*LLG(2,i)+a3*LLG(3,i)
          j=MGL(i)
          ztmp=-vqlg(j,ik)*dcmplx(sin(Gr),cos(Gr))*zrho(i)
          zsum1=zsum1+Gx*ztmp
          zsum2=zsum2+Gy*ztmp
          zsum3=zsum3+Gz*ztmp
       end do
       force(1,a) = -zsum1*Vcell
       force(2,a) = -zsum2*Vcell
       force(3,a) = -zsum3*Vcell

    end do ! a

    call destruct_Ggrid

    deallocate( wsavez,wsavey,wsavex )
    deallocate( lz2,lz1,ly2,ly1,lx2,lx1 )
    deallocate( fftwork )
    deallocate( zrho )
    deallocate( zrho3 )

    allocate( work(3*Natom) )
    n=0
    do a=1,Natom
       do i=1,3
          n=n+1
          work(n)=force(i,a)
       end do
    end do
    call mpi_allreduce(work,force,n,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    deallocate( work )


  END SUBROUTINE calc_force_ps_local


  FUNCTION Vloc(K,z)
    implicit none
    real(8) :: Vloc
    real(8),intent(IN) :: K
    integer,intent(IN) :: z
    real(8) :: pi,Zion,rloc,C1,C2,C3,C4

    pi = acos(-1.d0)

    rloc = 0.d0
    C1   = 0.d0
    C2   = 0.d0
    C3   = 0.d0
    C4   = 0.d0

    select case(z)
    case(1)
       Zion= 1.d0
       rloc= 0.2d0
       C1  =-4.0663326d0
       C2  = 0.6778322d0
    case default
       write(*,*) "z=",z
       write(*,*) "This atomic number is not available"
       stop
    end select

    if ( K == 0.d0 ) then
       Vloc = 2.d0*pi*Zion*rloc**2
    else
       Vloc = -4.d0*pi*Zion*exp(-0.5d0*(K*rloc)**2)/K**2
    end if

    Vloc = Vloc + sqrt((2.d0*pi)**3)*rloc**3*exp(-0.5d0*(rloc*K)**2) &
         *( C1 &
           +C2*(3.d0-(rloc*K)**2) &
           +C3*(15.d0-10.d0*(rloc*K)**2+(rloc*K)**4) &
           +C4*(105.d0-105.d0*(rloc*K)**2+21.d0*(rloc*K)**4-(rloc*K)**6) )
    return

  END FUNCTION Vloc

#ifdef TEST
  SUBROUTINE construct_ps_density_longloc
    use aa_module
    use hartree_module
    integer :: ik,j,ierr,ML_0,ML_1
    real(8) :: a1,a2,a3,c1,c2,c3,d1,d2,d3,pi,q1,q2,q3,q4,c4
    real(8) :: x,y,z,rr,rrmin,chk_conv,const,x2
    real(8),allocatable :: p1(:),p2(:),p3(:),p4(:)
    integer :: a,i,i1,i2,i3,j1,j2,j3,k1,k2,k3,icell,ncell
    integer,allocatable :: grid_from_origin(:,:)

    ML_0= Igrid(1,0)
    ML_1= Igrid(2,0)
    pi  = acos(-1.d0)

    allocate( rho_ps(ML_0:ML_1) )

    c1=1.d0/Ngrid(1)
    c2=1.d0/Ngrid(2)
    c3=1.d0/Ngrid(3)
    const = 1.d0/pi**1.5d0

    allocate( p1(Nelement),p2(Nelement) )
    allocate( p3(Nelement),p4(Nelement) )

    do ik=1,Nelement
       p2(ik)= parloc(2,ik)
       p1(ik)=-Zps(ik)*parloc(1,ik)*const*p2(ik)**1.5d0
       p4(ik)= parloc(4,ik)
       p3(ik)=-Zps(ik)*parloc(3,ik)*const*p4(ik)**1.5d0
    end do

    rewind 15
    do i=0,2000
       rr=(i*0.01d0)**2
       write(15,*) sqrt(rr) &
            ,-4.d0*pi*(p1(1)*exp(-p2(1)*rr)+p3(1)*exp(-p4(1)*rr) )
    end do

    allocate( grid_from_origin(3,ML_0:ML_1) )

    i=ML_0-1
    do i3=Igrid(1,3),Igrid(2,3)
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
       i=i+1
       rrmin=1.d10
       do j3=-1,1
       do j2=-1,1
       do j1=-1,1
          k1=i1+j1*Ngrid(1)
          k2=i2+j2*Ngrid(2)
          k3=i3+j3*Ngrid(3)
          rr=k1*k1+k2*k2+k3*k3
          if ( rr < rrmin ) then
             rrmin=rr
             grid_from_origin(1,i)=k1
             grid_from_origin(2,i)=k2
             grid_from_origin(3,i)=k3
          end if
       end do
       end do
       end do
    end do
    end do
    end do

    rho_ps(:)=0.d0
    chk_conv=0.d0
    do ncell=0,100
       do icell=-ncell,ncell
          if ( abs(icell) < ncell ) cycle
          do a=1,Natom
             ik=ki_atom(a)
             a1=aa_atom(1,a)
             a2=aa_atom(2,a)
             a3=aa_atom(3,a)+icell
             q1=p1(ik)
             q2=p2(ik)
             q3=p3(ik)
             q4=p4(ik)
             do i=ML_0,ML_1
                d1=grid_from_origin(1,i)*c1-a1
                d2=grid_from_origin(2,i)*c2-a2
                d3=grid_from_origin(3,i)*c3-a3
                x =aa(1,1)*d1+aa(1,2)*d2+aa(1,3)*d3
                y =aa(2,1)*d1+aa(2,2)*d2+aa(2,3)*d3
                z =aa(3,1)*d1+aa(3,2)*d2+aa(3,3)*d3
                rr=x*x+y*y+z*z
                rho_ps(i)=rho_ps(i)+q1*exp(-q2*rr)+q3*exp(-q4*rr)
             end do
!             i=ML_0-1
!             do i3=Igrid(1,3),Igrid(2,3)
!             do i2=Igrid(1,2),Igrid(2,2)
!             do i1=Igrid(1,1),Igrid(2,1)
!                d1=i1*c1-a1
!                d2=i2*c2-a2
!                d3=i3*c3-a3
!                x =aa(1,1)*d1+aa(1,2)*d2+aa(1,3)*d3
!                y =aa(2,1)*d1+aa(2,2)*d2+aa(2,3)*d3
!                z =aa(3,1)*d1+aa(3,2)*d2+aa(3,3)*d3
!                rr=x*x+y*y+z*z
!                i=i+1
!                rho_ps(i)=rho_ps(i)+q1*exp(-q2*rr)+q3*exp(-q4*rr)
!             end do
!             end do
!             end do
          end do ! a
       end do ! icell
       x=sum(rho_ps)*dV
       call mpi_allreduce(x,y,1,MPI_REAL8,MPI_SUM,comm_grid,ierr)
       write(*,*) "sum(rho_ps)=",y,ncell,minval(rho_ps),maxval(rho_ps)
       if ( abs(y-chk_conv) < 1.d-12 ) exit
       chk_conv=y
!       rho_ps(:)=rho_ps(:)/abs(x)*sum(Zps)
    end do ! ncell

    allocate( Vh(ML_0:ML_1) )
    Vh=0.d0
    call calc_hartree(ML_0,ML_1,1,rho_ps,0)

    c1 = parloc(1,1)
    c2 = sqrt( parloc(2,1) )
    c3 = parloc(3,1)
    c4 = sqrt( parloc(4,1) )

    rewind 25
    i=ML_0-1
    do i3=Igrid(1,3),Igrid(2,3)
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
       i=i+1
       x=i1*Hgrid(1)
       x2=abs( Hgrid(1)*(i1-Ngrid(1)) )
       if ( i2==0 .and. i3==40 ) then
          if ( i1 == 0 ) then
             y = -2.d0*Zps(1)/sqrt(pi)*(c1*c2+c3*c4)
          else
             y = -Zps(1)/x*(c1*bberf(c2*x)+c3*bberf(c4*x))
          end if
          if ( x2 == 0.d0 ) then
             z = -2.d0*Zps(1)/sqrt(pi)*(c1*c2+c3*c4)
          else
             z = -Zps(1)/x2*(c1*bberf(c2*x2)+c3*bberf(c4*x2))
          end if
          write(25,'(1x,i5,f10.5,4f20.10)') i1,x,rho_ps(i)*(-4.d0*pi),Vh(i),y,z
       end if
    end do
    end do
    end do
    deallocate( Vh )

    deallocate( grid_from_origin )
    deallocate( p4,p3,p2,p1 )

  END SUBROUTINE construct_ps_density_longloc


  SUBROUTINE construct_ps_density_longloc2
    use kinetic_module, only: op_kinetic
    use array_bound_module
    implicit none
    complex(8),allocatable :: work1(:),work2(:)
    integer :: i,i1,i2,i3
    real(8) :: c
    
    allocate( work1(ML_0:ML_1) )
    allocate( work2(ML_0:ML_1) )
    work1(:)=Vion(:)
    work2=0.d0

    call op_kinetic(0,work1,work2,ML_0,ML_1,1,1)

    c = -2.d0/( 4.d0*acos(-1.d0) )
    work2(:) = -c*work2(:)

    rewind 26
    i=ML_0-1
    do i3=Igrid(1,3),Igrid(2,3)
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
       i=i+1
       if ( i2==0 .and. i3==0 ) then
          write(26,'(1x,i5,f10.5,2f20.15)') &
               i1,i1*Hgrid(1),real(work2(i)),aimag(work2(i))
       end if
    end do
    end do
    end do

    write(*,*) "sum(work2)=",sum(work2)*dV
!    rho_ps(:)=work2(:)

    deallocate( work2 )
    deallocate( work1 )

  END SUBROUTINE construct_ps_density_longloc2

  SUBROUTINE construct_ps_density_longloc3
    use array_bound_module
    real(8),allocatable :: work(:)
    integer :: ifacx(30),ifacy(30),ifacz(30)
    integer,allocatable :: lx1(:),lx2(:),ly1(:),ly2(:),lz1(:),lz2(:)
    complex(8),allocatable :: zwork1(:,:,:),zwork0(:,:,:)
    complex(8),allocatable :: wsavex(:),wsavey(:),wsavez(:)
    integer :: n1,n2,ML1,ML2,ML3
    integer :: i1,i2,i3,j1,j2,j3,i,irank,ierr
    real(8) :: pi4,g2

    n1=ML_0
    n2=ML_1
    ML1=Ngrid(1)
    ML2=Ngrid(2)
    ML3=Ngrid(3)

    pi4 = 4.d0*acos(-1.d0)

    allocate( zwork0(0:ML1-1,0:ML2-1,0:ML3-1) )

    allocate( work(ML) ) ; work=0.d0
    work(n1:n2) = Vion(n1:n2)
    call mpi_allgatherv(work(n1),n2-n1+1,mpi_real8 &
         ,work,ir_grid,id_grid,mpi_real8,comm_grid,ierr)

    irank=-1
    i=0
    do i3=1,node_partition(3)
    do i2=1,node_partition(2)
    do i1=1,node_partition(1)
       irank=irank+1
       do j3=pinfo_grid(5,irank),pinfo_grid(5,irank)+pinfo_grid(6,irank)-1
       do j2=pinfo_grid(3,irank),pinfo_grid(3,irank)+pinfo_grid(4,irank)-1
       do j1=pinfo_grid(1,irank),pinfo_grid(1,irank)+pinfo_grid(2,irank)-1
          i=i+1
          zwork0(j1,j2,j3)=work(i)
       end do
       end do
       end do
    end do
    end do
    end do

    deallocate( work )

    allocate( zwork1(0:ML1-1,0:ML2-1,0:ML3-1) )
    allocate( lx1(ML),lx2(ML),ly1(ML),ly2(ML),lz1(ML),lz2(ML) )
    allocate( wsavex(ML1),wsavey(ML2),wsavez(ML3) )

    call prefft(ML1,ML2,ML3,ML,wsavex,wsavey,wsavez &
         ,ifacx,ifacy,ifacz,lx1,lx2,ly1,ly2,lz1,lz2)

    call fft3fx(ML1,ML2,ML3,ML,zwork0,zwork1,wsavex,wsavey,wsavez &
               ,ifacx,ifacy,ifacz,lx1,lx2,ly1,ly2,lz1,lz2)

    call construct_Ggrid(2)

    zwork1(:,:,:)=(0.d0,0.d0)
    do i=1,NGgrid(0)
       g2=GG(MGL(i))
       i1=LLG(1,i)
       i2=LLG(2,i)
       i3=LLG(3,i)
       zwork1(i1,i2,i3)=zwork0(i1,i2,i3)*g2/pi4
    end do

    call destruct_Ggrid

    call fft3bx(ML1,ML2,ML3,ML,zwork1,zwork0,wsavex,wsavey,wsavez &
         ,ifacx,ifacy,ifacz,lx1,lx2,ly1,ly2,lz1,lz2)

    deallocate( wsavez,wsavey,wsavex )
    deallocate( lz2,lz1,ly2,ly1,lx2,lx1 )

    rewind 27
    i=ML_0-1
    do i3=Igrid(1,3),Igrid(2,3)
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
       if ( i3 == 0 .and. i2 == 0 ) then
          write(27,'(1x,i5,f10.5,2f20.15)') &
               i1,i1*Hgrid(1),real(zwork1(i1,i2,i3)),aimag(zwork1(i1,i2,i3))

       end if
       i=i+1
!       rho_ps(i)=zwork1(i1,i2,i3)
    end do
    end do
    end do

    write(*,*) "sum(zwork1)=",sum(zwork1)*dV

    deallocate( zwork1 )

  END SUBROUTINE construct_ps_density_longloc3
#endif

  SUBROUTINE construct_ps_local_ffte

    integer :: a,i,i1,i2,i3,ik,j,MG
    integer :: ML1,ML2,ML3,ML,ML_0,ML_1
    integer :: ifacx(30),ifacy(30),ifacz(30)
    complex(8),allocatable :: zwork(:,:,:),vg(:)
    complex(8),allocatable :: zwork1(:,:,:),zwork2(:,:,:)
    real(8) :: ctt(0:3),ett(0:3)
    integer :: a1b,b1b,a2b,b2b,a3b,b3b,np1,np2,np3
    integer :: n1,n2,MG1,NG1,MG2,NG2,MG3,NG3
    integer :: comm_fftx, npux, ix
    integer :: comm_ffty, npuy, iy
    integer :: comm_fftz, npuz, iz
    integer :: icolor,ierr

    MG  = NGgrid(0)
    ML  = Ngrid(0)
    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)
    ML_0= Igrid(1,0)
    ML_1= Igrid(2,0)
    a1b = Igrid(1,1)
    b1b = Igrid(2,1)
    a2b = Igrid(1,2)
    b2b = Igrid(2,2)
    a3b = Igrid(1,3)
    b3b = Igrid(2,3)
    np1=node_partition(1)
    np2=node_partition(2)
    np3=node_partition(3)
    n1=ML_0
    n2=ML_1
    MG1=-ML1/2     ; MG2=-ML2/2     ; MG3=-ML3/2
    NG1=(ML1-1)/2 ; NG2=(ML2-1)/2 ; NG3=(ML3-1)/2


    ix=a1b/(ML1/np1); iy=a2b/(ML2/np2); iz=a3b/(ML3/np3);
    icolor=iy+iz*np2
    call mpi_comm_split(comm_grid,icolor, 0, comm_fftx, ierr)
    icolor=iz+ix*nprocs
    call mpi_comm_split(comm_grid,icolor, 0, comm_ffty, ierr)
    icolor=iy+ix*nprocs
    call mpi_comm_split(comm_grid,icolor, 0, comm_fftz, ierr)

    call mpi_comm_size(comm_fftx, npux, ierr)
    call mpi_comm_size(comm_ffty, npuy, ierr)
    call mpi_comm_size(comm_fftz, npuz, ierr)

    ctt(:)=0.d0
    ett(:)=0.d0

    call watch(ctt(0),ett(0))

    if ( .not.allocated(Vion) ) then
       allocate( Vion(ML_0:ML_1) )
    end if

    allocate( zwork(0:ML1-1,0:ML2-1,0:ML3-1) )
    allocate( zwork1(0:ML1-1,a2b:b2b,a3b:b3b))
    allocate( zwork2(0:ML1-1,a2b:b2b,a3b:b3b))

    zwork1(:,:,:)=(0.d0,0.d0)
    zwork2(:,:,:)=(0.d0,0.d0)

    allocate( vg(MG) )

    do i=MG_0,MG_1
       j=MGL(i)
       vg(i)=vqlg(j,1)*SGK(i,1)
    end do
    do ik=2,Nelement
       do i=MG_0,MG_1
          j=MGL(i)
          vg(i)=vg(i)+vqlg(j,ik)*SGK(i,ik)
       end do
    end do
    call allgatherv_Ggrid(vg)

    call construct_Ggrid(2)

    zwork(:,:,:)=(0.d0,0.d0)
    do i=1,NGgrid(0)
       zwork(LLG(1,i),LLG(2,i),LLG(3,i))=vg(i)
       i3=LLG(3,i)
       i2=LLG(2,i)
       i1=LLG(1,i)
       if(a2b<=i2 .and. i2<=b2b .and. a3b<=i3 .and. i3<=b3b) then
          zwork1(i1,i2,i3)=zwork(i1,i2,i3)
       endif
    end do

    call destruct_Ggrid

    deallocate( vg )

!    allocate( fftwork(ML) )
!    allocate( lx1(ML),lx2(ML),ly1(ML),ly2(ML),lz1(ML),lz2(ML) )
!    allocate( wsavex(ML1),wsavey(ML2),wsavez(ML3) )

!    call prefft(ML1,ML2,ML3,ML,wsavex,wsavey,wsavez &
!         ,ifacx,ifacy,ifacz,lx1,lx2,ly1,ly2,lz1,lz2)

    call watch(ctt(1),ett(1))

!    call fft3bx(ML1,ML2,ML3,ML,zwork,fftwork,wsavex,wsavey,wsavez &
!         ,ifacx,ifacy,ifacz,lx1,lx2,ly1,ly2,lz1,lz2)

    call pzfft3dv(zwork1,zwork2,ML1,ML2,ML3,comm_ffty,comm_fftz,npuy,npuz,0)
    call pzfft3dv(zwork1,zwork2,ML1,ML2,ML3,comm_ffty,comm_fftz,npuy,npuz,1)
    zwork1=ML*zwork2

    call watch(ctt(2),ett(2))

    i=ML_0-1
    do i3=Igrid(1,3),Igrid(2,3)
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
       i=i+1
       Vion(i)=real( zwork1(i1,i2,i3) )
    end do
    end do
    end do

!    deallocate( wsavez,wsavey,wsavex )
!    deallocate( lz2,lz1,ly2,ly1,lx2,lx1 )
!    deallocate( fftwork )
    deallocate( zwork )
    deallocate(zwork1,zwork2)

    call mpi_comm_free(comm_fftx,ierr)
    call mpi_comm_free(comm_ffty,ierr)
    call mpi_comm_free(comm_fftz,ierr)

    call watch(ctt(3),ett(3))

    if ( disp_switch_parallel ) then
       write(*,*) "time(const_ps_loc_1)",ctt(1)-ctt(0),ett(1)-ett(0)
       write(*,*) "time(const_ps_loc_2)",ctt(2)-ctt(1),ett(2)-ett(1)
       write(*,*) "time(const_ps_loc_3)",ctt(3)-ctt(2),ett(3)-ett(2)
    end if

  END SUBROUTINE construct_ps_local_ffte


  SUBROUTINE calc_force_ps_local_ffte(MI,force)
    integer,intent(IN) :: MI
    real(8),intent(OUT) :: force(3,MI)
    integer :: a,i,j,ik,ispin,irank,i1,i2,i3,n,l1,l2,l3,m1,m2,m3,j1,j2,j3
    integer :: MI_0,MI_1,ML1,ML2,ML3,ML_0,ML_1,ML,MG,N_MI,ierr
    integer,allocatable :: icnt(:),idis(:)
    integer :: m
    real(8) :: a1,a2,a3,pi2,Gr,Gx,Gy,Gz,Vcell,g2
    real(8),allocatable :: work(:)
    integer :: ifacx(30),ifacy(30),ifacz(30)
!    integer,allocatable :: lx1(:),lx2(:),ly1(:),ly2(:),lz1(:),lz2(:)
!    complex(8),allocatable :: fftwork(:)
    complex(8),allocatable :: zrho3(:,:,:),zrho(:)
!    complex(8),allocatable :: wsavex(:),wsavey(:),wsavez(:)
    complex(8),parameter :: z0=(0.d0,0.d0)
    complex(8),allocatable :: zwork1(:,:,:),zwork2(:,:,:),zwork3(:,:,:),w1(:),w2(:)
    complex(8) :: zsum1,zsum2,zsum3,ztmp
    integer :: comm_fftx, npux, ix
    integer :: comm_ffty, npuy, iy
    integer :: comm_fftz, npuz, iz
    integer :: icolor
    integer :: a1b,b1b,a2b,b2b,a3b,b3b,np1,np2,np3
    integer :: MG1,NG1,MG2,NG2,MG3,NG3,n1,n2,MGmax


    force(:,:)=0.d0

    MG    = NGgrid(0)
    ML    = Ngrid(0)
    ML1   = Ngrid(1)
    ML2   = Ngrid(2)
    ML3   = Ngrid(3)
    ML_0  = Igrid(1,0)
    ML_1  = Igrid(2,0)
    pi2   = 2.d0*acos(-1.d0)
    Vcell = ML*dV
    a1b = Igrid(1,1)
    b1b = Igrid(2,1)
    a2b = Igrid(1,2)
    b2b = Igrid(2,2)
    a3b = Igrid(1,3)
    b3b = Igrid(2,3)
    np1=node_partition(1)
    np2=node_partition(2)
    np3=node_partition(3)
    n1=ML_0
    n2=ML_1
    MG1=-ML1/2     ; MG2=-ML2/2     ; MG3=-ML3/2
    NG1=(ML1-1)/2 ; NG2=(ML2-1)/2 ; NG3=(ML3-1)/2


    ix=a1b/(ML1/np1); iy=a2b/(ML2/np2); iz=a3b/(ML3/np3);
    icolor=iy+iz*np2
    call mpi_comm_split(comm_grid,icolor, 0, comm_fftx, ierr)
    icolor=iz+ix*nprocs
    call mpi_comm_split(comm_grid,icolor, 0, comm_ffty, ierr)
    icolor=iy+ix*nprocs
    call mpi_comm_split(comm_grid,icolor, 0, comm_fftz, ierr)

    call mpi_comm_size(comm_fftx, npux, ierr)
    call mpi_comm_size(comm_ffty, npuy, ierr)
    call mpi_comm_size(comm_fftz, npuz, ierr)

    allocate( icnt(0:nprocs-1) )
    allocate( idis(0:nprocs-1) )
    N_MI = Natom/nprocs
    icnt(0:nprocs-1) = N_MI
    n = Natom - N_MI*nprocs
    if ( n>0 ) then
       do irank=0,n-1
          icnt(irank)=icnt(irank)+1
       end do
    end if
    do irank=0,nprocs-1
       idis(irank) = sum( icnt(0:irank) ) - icnt(irank)
    end do
    MI_0 = idis(myrank)+1
    MI_1 = idis(myrank)+icnt(myrank)
    deallocate( idis, icnt )

    allocate( zrho3(0:ML1-1,0:ML2-1,0:ML3-1) )

    allocate( work(ML) )
    work(1:ML)=0.d0
    do ispin=1,Nspin
       work(ML_0:ML_1)=work(ML_0:ML_1)+rho(ML_0:ML_1,ispin)
    end do
    call mpi_allgatherv(work(ML_0),ML_1-ML_0+1,mpi_real8,work &
         ,ir_grid,id_grid,mpi_real8,comm_grid,ierr)

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
          zrho3(j1,j2,j3)=work(i)
       end do
       end do
       end do
    end do
    end do
    end do
    deallocate( work )

    allocate(zwork1(0:ML1-1,a2b:b2b,a3b:b3b))
    allocate(zwork2(0:ML1-1,a2b:b2b,a3b:b3b))
    zwork1(:,:,:)=z0
    zwork2(:,:,:)=z0

    do i3=a3b,b3b
    do i2=a2b,b2b
    do i1=a1b,b1b
      zwork1(i1,i2,i3)=zrho3(i1,i2,i3)
    end do
    end do
    end do

    call mpi_allreduce(zwork1,zwork2,ML1*(b2b-a2b+1)*(b3b-a3b+1),mpi_complex16,mpi_sum,comm_fftx,ierr)

!    allocate( fftwork(ML) )
!    allocate( lx1(ML),lx2(ML),ly1(ML),ly2(ML),lz1(ML),lz2(ML) )
!    allocate( wsavex(ML1),wsavey(ML2),wsavez(ML3) )
!
!    call prefft(ML1,ML2,ML3,ML,wsavex,wsavey,wsavez &
!         ,ifacx,ifacy,ifacz,lx1,lx2,ly1,ly2,lz1,lz2)
!
!    call fft3fx(ML1,ML2,ML3,ML,zrho3,fftwork,wsavex,wsavey,wsavez &
!               ,ifacx,ifacy,ifacz,lx1,lx2,ly1,ly2,lz1,lz2)


!----------------------------------
      call pzfft3dv(zwork2,zwork1,ML1,ML2,ML3,comm_ffty,comm_fftz,npuy,npuz,0)
      call pzfft3dv(zwork2,zwork1,ML1,ML2,ML3,comm_ffty,comm_fftz,npuy,npuz,-1)
!----------------------------------

    zwork1(:,:,:)=zwork1(:,:,:)/DCMPLX(ML1*ML2*ML3)

    call construct_Ggrid(0)

    allocate(w1(ML)) ; w1(:)=z0
    allocate(w2(ML)) ; w2(:)=z0
    allocate(zwork3(a1b:b1b,a2b:b2b,a3b:b3b))

    m=0
    do i3=MG3,NG3
    do i2=MG2,NG2
    do i1=MG1,NG1
       if ( i1==0 .and. i2==0 .and. i3==0 ) cycle
         m=m+1 
    end do
    end do
    end do

    MGmax = m


    zwork3(:,:,:)=z0
    do i3=MG3,NG3
    do i2=MG2,NG2
    do i1=MG1,NG1
         g2=(bb(1,1)*i1+bb(1,2)*i2+bb(1,3)*i3)**2 &
     &          +(bb(2,1)*i1+bb(2,2)*i2+bb(2,3)*i3)**2 &
     &          +(bb(3,1)*i1+bb(3,2)*i2+bb(3,3)*i3)**2
         if ( g2<=0.d0 .or. (MGmax/=MG .and. Ecut-1.d-14<=g2) ) cycle
         j1=mod(i1+ML1,ML1)
         j2=mod(i2+ML2,ML2)
         j3=mod(i3+ML3,ML3)
         if(a1b <= j1 .and. j1 <= b1b .and. a2b <= j2 .and. j2 <= b2b .and. a3b <= j3 .and. j3 <= b3b) then
            zwork3(j1,j2,j3)=zwork1(j1,j2,j3)
         endif
    end do
    end do
    end do

    i=n1-1
    do i3=a3b,b3b
    do i2=a2b,b2b
    do i1=a1b,b1b
      i=i+1
      w1(i)=zwork3(i1,i2,i3)
    enddo
    enddo
    enddo

    call mpi_allreduce(w1,w2,ML,mpi_complex16,mpi_sum,comm_grid,ierr)

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
          zrho3(j1,j2,j3)=w2(i)
       end do
       end do
       end do
    end do
    end do
    end do

    allocate( zrho(MG) )

    do i=1,NGgrid(0)
       i1=mod(ML1+LLG(1,i),ML1)
       i2=mod(ML2+LLG(2,i),ML2)
       i3=mod(ML3+LLG(3,i),ML3)
       zrho(i) = conjg( zrho3(i1,i2,i3) )
    end do

    do a=MI_0,MI_1

       ik=ki_atom(a)
       a1=pi2*aa_atom(1,a)
       a2=pi2*aa_atom(2,a)
       a3=pi2*aa_atom(3,a)

       zsum1=z0
       zsum2=z0
       zsum3=z0
       do i=1,NGgrid(0)
          Gx=bb(1,1)*LLG(1,i)+bb(1,2)*LLG(2,i)+bb(1,3)*LLG(3,i)
          Gy=bb(2,1)*LLG(1,i)+bb(2,2)*LLG(2,i)+bb(2,3)*LLG(3,i)
          Gz=bb(3,1)*LLG(1,i)+bb(3,2)*LLG(2,i)+bb(3,3)*LLG(3,i)
          Gr=a1*LLG(1,i)+a2*LLG(2,i)+a3*LLG(3,i)
          j=MGL(i)
          ztmp=-vqlg(j,ik)*dcmplx(sin(Gr),cos(Gr))*zrho(i)
          zsum1=zsum1+Gx*ztmp
          zsum2=zsum2+Gy*ztmp
          zsum3=zsum3+Gz*ztmp
       end do
       force(1,a) = -zsum1*Vcell
       force(2,a) = -zsum2*Vcell
       force(3,a) = -zsum3*Vcell

    end do ! a

    call destruct_Ggrid

!    deallocate( wsavez,wsavey,wsavex )
!    deallocate( lz2,lz1,ly2,ly1,lx2,lx1 )
!    deallocate( fftwork )
    deallocate( zrho )
    deallocate( zrho3 )
    deallocate( zwork1,zwork2)
    deallocate(w1,w2)
    deallocate(zwork3)

    allocate( work(3*Natom) )
    n=0
    do a=1,Natom
       do i=1,3
          n=n+1
          work(n)=force(i,a)
       end do
    end do
    call mpi_allreduce(work,force,n,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

!--------------------------------------------------------
!    write(*,*) "LOCAL FORCE"
!    do a=1,Natom
!       do i=1,3
!         write(999,*) force(i,a)
!       enddo
!    enddo
!--------------------------------------------------------
    deallocate( work )

    call mpi_comm_free(comm_fftx,ierr)
    call mpi_comm_free(comm_ffty,ierr)
    call mpi_comm_free(comm_fftz,Ierr)


  END SUBROUTINE calc_force_ps_local_ffte


END MODULE ps_local_module
