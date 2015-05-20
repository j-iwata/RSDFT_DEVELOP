MODULE localpot2_ion_module

  use pseudopot_module, only: Zps,rad,rab,NRps,Rps,parloc,norb,Mr,vql
  use aa_module
  use bb_module
  use atom_module
  use simc_module

  use rgrid_module
  use parallel_module
  use localpot2_variables, only: fecut_loc,Ngrid_dense,Igrid_dense,dV_dense &
       ,Ndens_loc
  use array_bound_module
  use electron_module
  use wf_module

  use watch_module

  use ggrid_module
  use strfac_module
  use bberf_module

  implicit none

  PRIVATE
  PUBLIC :: localpot2_ion, localpot2_calc_eion

CONTAINS


  SUBROUTINE localpot2_ion( MKI, ecut_in, vout )
    implicit none
    integer,intent(IN)  :: MKI
    real(8),intent(IN)  :: ecut_in
    real(8),intent(OUT) :: vout(:,:,:)

    real(8) :: pi,rloc,const,C1,C2,C3,C4,G,v,ecut
    integer :: ig,ielm,mm1,mm2,mm3,ierr
    integer :: m1,m2,m3
    real(8),allocatable :: vqlg(:,:)
    complex(8),allocatable :: vg(:)
    integer :: ifacx(30),ifacy(30),ifacz(30)
    integer,allocatable :: lx1(:),lx2(:),ly1(:),ly2(:),lz1(:),lz2(:)
    complex(8),allocatable :: wsavex(:),wsavey(:),wsavez(:)
    complex(8),allocatable :: zw0(:,:,:),zw1(:,:,:)
    integer :: i1,i2,i3,j1,j2,j3,mm,a,ik
    real(8) :: Gx,Gy,Gz,GG,Gr,a1,a2,a3,pi2,ct0,et0,ct1,et1
    integer,allocatable :: ir(:),id(:)
    integer :: mg2,mg2_0,mg2_1,n,i

    call watch(ct0,et0)

    pi = acos(-1.d0)
    const = 1.d0/abs(Va)

    ecut = ecut_in*fecut_loc

! - structure factor

    m1 = Ngrid_dense(1)
    m2 = Ngrid_dense(2)
    m3 = Ngrid_dense(3)

    mm1 = (m1-1)/2
    mm2 = (m2-1)/2
    mm3 = (m3-1)/2

    call get_cutoff_ggrid_2(mm1,mm2,mm3,ecut,mg2)

    allocate( ir(0:nprocs-1) ) ; ir=0
    allocate( id(0:nprocs-1) ) ; id=0
    do i=1,mg2
       n=mod(i-1,nprocs)
       ir(n)=ir(n)+1
    end do
    do n=0,nprocs-1
       id(n) = sum( ir(0:n) ) - ir(n)
    end do
    MG2_0=id(myrank)+1
    MG2_1=id(myrank)+ir(myrank)

    call construct_ggrid_2(mm1,mm2,mm3,MG2,MG2_0,MG2_1,ecut,1)

    call construct_strfac_2(MG2_0,MG2_1)

! - Fourier component

    allocate( vqlg(MG2_0:MG2_1,MKI) )
    vqlg = 0.0d0

    call init_ps_local_2(MG2_0,MG2_1,MKI,vqlg)

    call destruct_ggrid

    allocate( vg(MG2) )
    vg=(0.0d0,0.0d0)

    do i=MG2_0,MG2_1
       vg(i)=vqlg(i,1)*SGK(i,1)
    end do
    do ik=2,Nelement
       do i=MG2_0,MG2_1
          vg(i)=vg(i)+vqlg(i,ik)*SGK(i,ik)
       end do
    end do
    call mpi_allgatherv(vg(MG2_0),ir(myrank),MPI_COMPLEX16 &
         ,vg,ir,id,MPI_COMPLEX16,MPI_COMM_WORLD,ierr)

    call construct_ggrid_2(mm1,mm2,mm3,MG2,MG2_0,MG2_1,ecut,0)

    allocate( zw0(0:m1-1,0:m2-1,0:m3-1) )
    zw0=(0.0d0,0.0d0)
    do i=1,MG2
       i1=mod(LLG(1,i)+m1,m1)
       i2=mod(LLG(2,i)+m2,m2)
       i3=mod(LLG(3,i)+m3,m3)
       zw0(i1,i2,i3)=vg(i)
    end do

    call destruct_ggrid

    deallocate( vg )

    allocate( zw1(0:m1-1,0:m2-1,0:m3-1) )

    mm=m1*m2*m3
    allocate( lx1(mm),lx2(mm),ly1(mm),ly2(mm),lz1(mm),lz2(mm) )
    allocate( wsavex(m1),wsavey(m2),wsavez(m3) )

    call prefft(m1,m2,m3,mm,wsavex,wsavey,wsavez,ifacx,ifacy,ifacz &
               ,lx1,lx2,ly1,ly2,lz1,lz2)
    call fft3bx(m1,m2,m3,mm,zw0,zw1,wsavex,wsavey,wsavez,ifacx,ifacy,ifacz &
               ,lx1,lx2,ly1,ly2,lz1,lz2)

    deallocate( wsavez,wsavey,wsavex )
    deallocate( lz2,lz1,ly2,ly1,lx2,lx1 )

    j1=Igrid_dense(1,1)-1
    j2=Igrid_dense(1,2)-1
    j3=Igrid_dense(1,3)-1
    do i3=Igrid_dense(1,3),Igrid_dense(2,3)
    do i2=Igrid_dense(1,2),Igrid_dense(2,2)
    do i1=Igrid_dense(1,1),Igrid_dense(2,1)
       vout(i1-j1,i2-j2,i3-j3) = real( zw0(i1,i2,i3) )
    end do
    end do
    end do

    deallocate( zw1,zw0 )

    deallocate( vqlg )
    deallocate( id,ir )

    call watch(ct1,et1) ; write(*,*) "localpot2_ion",ct1-ct0,et1-et0

    return

  END SUBROUTINE localpot2_ion

#ifdef TEST
  SUBROUTINE localpot2_calc_eion(m1,m2,m3,vin,eout)
    implicit none
    integer,intent(IN)  :: m1,m2,m3
    real(8),intent(IN)  :: vin(0:m1-1,0:m2-1,0:m3-1)
    real(8),intent(OUT) :: eout
    real(8),allocatable :: vloc_nl_bak(:,:)
    integer :: n1,n2,n,k,s,i,j,ierr
    real(8) :: sum0,eout0

    n1=size( vloc_nl, 1)
    n2=size( vloc_nl, 2)
    allocate( vloc_nl_bak(n1,n2) ) ; vloc_nl_bak=0.0d0

    vloc_nl_bak(:,:) = vloc_nl(:,:)

    call test2_localpot2(m1,m2,m3,vin)

    eout0=0.0d0

    do s=1,MSP
    do k=1,MBZ
    do n=1,MB

       sum0=0.0d0
       do i=1,Ngrid(0)
       do j=1,MLpot
#ifdef _DRSDFT_
          sum0 = sum0 + unk(i,n,k,s)*vloc_nl(j,i)*unk(Lpot(j,i),n,k,s)
#else
          sum0 = sum0 + conjg(unk(i,n,k,s))*vloc_nl(j,i)*unk(Lpot(j,i),n,k,s)
#endif
       end do
       end do

       eout0 = eout0 + occ(n,k,s)*sum0*dV

    end do
    end do
    end do

    call mpi_allreduce(eout0,eout,1,MPI_REAL8,MPI_SUM,comm_grid,ierr)

    vloc_nl(:,:) = vloc_nl_bak(:,:)

    deallocate( vloc_nl_bak )

  END SUBROUTINE localpot2_calc_eion
#else
  SUBROUTINE localpot2_calc_eion( vin, nin, eout )
    implicit none
    real(8),intent(IN)  :: vin(:,:,:)
    real(8),intent(IN)  :: nin(:,:,:)
    real(8),intent(OUT) :: eout
    integer :: ierr
    real(8) :: eout0

    eout0 = sum( vin*nin )*dV_dense

    call mpi_allreduce(eout0,eout,1,MPI_REAL8,MPI_SUM,comm_grid,ierr)

  END SUBROUTINE localpot2_calc_eion
#endif

  SUBROUTINE init_ps_local(mm1,mm2,mm3,MKI,vqlg)
    implicit none
    integer,intent(IN)  :: mm1,mm2,mm3,MKI
    real(8),intent(OUT) :: vqlg(-mm1:mm1,-mm2:mm2,-mm3:mm3,MKI)
    integer :: i,ig,ik,iorb,MMr,NRc,i1,i2,i3
    real(8) :: Rc,p1,p2,p3,p4,vlong,Pi,const,x,r,sb,sum0,G,G2
    real(8) :: Vcell,Gx,Gy,Gz,vqlgl
    real(8),allocatable :: vshort(:),tmp(:)

    vqlg(:,:,:,:)=0.0d0

    Vcell = Ngrid(0)*dV
    Pi    = acos(-1.d0)
    const = 4.d0*Pi/Vcell

    MMr=maxval(Mr)
    allocate( vshort(MMr) ) ; vshort=0.0d0

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
       end do

       allocate( tmp(MMr) )

       do i3=-mm3,mm3
       do i2=-mm2,mm2
       do i1=-mm1,mm1

          Gx=bb(1,1)*i1+bb(1,2)*i2+bb(1,3)*i3
          Gy=bb(2,1)*i1+bb(2,2)*i2+bb(2,3)*i3
          Gz=bb(3,1)*i1+bb(3,2)*i2+bb(3,3)*i3
          G2=Gx*Gx+Gy*Gy+Gz*Gz
          G=sqrt(G2)

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

          vqlg(i1,i2,i3,ik)=sum0*const

       end do
       end do
       end do

       p1=-Zps(ik)*parloc(1,ik) ; p2=0.25d0/parloc(2,ik)
       p3=-Zps(ik)*parloc(3,ik) ; p4=0.25d0/parloc(4,ik)

       do i3=-mm3,mm3
       do i2=-mm2,mm2
       do i1=-mm1,mm1

          Gx=bb(1,1)*i1+bb(1,2)*i2+bb(1,3)*i3
          Gy=bb(2,1)*i1+bb(2,2)*i2+bb(2,3)*i3
          Gz=bb(3,1)*i1+bb(3,2)*i2+bb(3,3)*i3
          G2=Gx*Gx+Gy*Gy+Gz*Gz
          G=sqrt(G2)

          if ( G2 == 0.d0 ) then
             vqlgl=-(p1*p2+p3*p4)*const
             vqlg(i1,i2,i3,ik)=vqlg(i1,i2,i3,ik)+vqlgl
          else
             vqlgl=(p1*exp(-G2*p2)+p3*exp(-G2*p4))/G2*const
             vqlg(i1,i2,i3,ik)=vqlg(i1,i2,i3,ik)+vqlgl
          end if

       end do
       end do
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


  SUBROUTINE init_ps_local_2(MG2_0,MG2_1,MKI,vqlg)
    implicit none
    integer,intent(IN)  :: MG2_0,MG2_1,MKI
    real(8),intent(OUT) :: vqlg(MG2_0:MG2_1,MKI)
    integer :: i,ig,ik,iorb,MMr,NRc,i1,i2,i3
    real(8) :: Rc,p1,p2,p3,p4,vlong,Pi,const,x,r,sb,sum0,G,G2
    real(8) :: Vcell,Gx,Gy,Gz,vqlgl
    real(8),allocatable :: vshort(:),tmp(:)

    vqlg(:,:)=0.0d0

    Vcell = Ngrid(0)*dV
    Pi    = acos(-1.d0)
    const = 4.d0*Pi/Vcell

    MMr=maxval(Mr)
    allocate( vshort(MMr) ) ; vshort=0.0d0

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
       end do

       allocate( tmp(MMr) )

       do ig=MG2_0,MG2_1

          i1=LLG(1,ig)
          i2=LLG(2,ig)
          i3=LLG(3,ig)

          Gx=bb(1,1)*i1+bb(1,2)*i2+bb(1,3)*i3
          Gy=bb(2,1)*i1+bb(2,2)*i2+bb(2,3)*i3
          Gz=bb(3,1)*i1+bb(3,2)*i2+bb(3,3)*i3
          G2=Gx*Gx+Gy*Gy+Gz*Gz
          G=sqrt(G2)

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

       do ig=MG2_0,MG2_1

          i1 = LLG(1,ig)
          i2 = LLG(2,ig)
          i3 = LLG(3,ig)

          Gx=bb(1,1)*i1+bb(1,2)*i2+bb(1,3)*i3
          Gy=bb(2,1)*i1+bb(2,2)*i2+bb(2,3)*i3
          Gz=bb(3,1)*i1+bb(3,2)*i2+bb(3,3)*i3
          G2=Gx*Gx+Gy*Gy+Gz*Gz
          G=sqrt(G2)

          if ( G2 == 0.d0 ) then
             vqlgl=-(p1*p2+p3*p4)*const
             vqlg(ig,ik)=vqlg(ig,ik)+vqlgl
          else
             vqlgl=(p1*exp(-G2*p2)+p3*exp(-G2*p4))/G2*const
             vqlg(ig,ik)=vqlg(ig,ik)+vqlgl
          end if

       end do

       deallocate( tmp )

    end do ! ik

    deallocate( vshort )

  END SUBROUTINE init_ps_local_2


END MODULE localpot2_ion_module
