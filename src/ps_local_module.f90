module ps_local_module

  use fft_module, only: iswitch_fft
  use rgrid_module, only: Igrid, Ngrid, dV
  use ggrid_module, only: NMGL,GG
  use atom_module, only: ki_atom
  use strfac_module, only: SGK
  use bberf_module, only: bberf
  use ps_local_variables, only: vqlg
  use var_ps_member, only: Zps,Rps,NRps,parloc,rad,norb,rab,Mr,ippform,vql
  use simc_module, only: simc
  use ps_local_gth_module, only: init_ps_local_gth
  use ps_local_1dffte_module, only: init_ps_local_1dffte, construct_ps_local_1dffte
  use ps_local_ffte_module, only: init_ps_local_ffte, construct_ps_local_ffte
  use ps_local_fftw_module, only: construct_ps_local_fftw
  use ps_local_fft0_module, only: construct_ps_local_fft0

  implicit none

  private
  public :: init_ps_local
  public :: construct_ps_local
  public :: simp ! for stress

  real(8),public :: const_ps_local
  real(8),public,allocatable :: Vion(:)

  logical :: flag_zero_ave = .false.

contains

  subroutine init_ps_local
    implicit none
    integer :: i,ig,ik,iorb,MMr,NRc,MKI
    real(8) :: Rc,p1,p2,p3,p4,vlong,Pi,const,x,r,sb,sum0,G,G2
    real(8) :: Vcell
    real(8),allocatable :: vshort(:),tmp(:)

    call write_border( 0, " init_ps_local(start)" )

    MKI   = size(Zps)
    Vcell = Ngrid(0)*dV
    Pi    = acos(-1.0d0)
    const = 4.0d0*Pi/Vcell

    if ( allocated(vqlg) ) deallocate(vqlg)
    allocate( vqlg(NMGL,MKI) ); vqlg=0.0d0

    MMr=maxval(Mr)

    allocate( vshort(MMr) ); vshort=0.0d0

    do ik=1,MKI

      if ( ippform(ik) == 4 ) then
        call init_ps_local_gth( Vcell, NMGL, ik, GG, vqlg(1,ik) )
        cycle
      end if

      MMr = Mr(ik)

      Rc=0.0d0
      NRc=0
      do iorb=1,norb(ik)
        Rc=max( Rc, Rps(iorb,ik) )
        NRc=max( NRc, NRps(iorb,ik) )
      end do

      if ( Rc<1.d-8 ) Rc=5.0d0
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
        if ( r < 1.0d-9 ) then
          vlong=-2.d0*Zps(ik)/sqrt(Pi)*(p1*p2+p3*p4)
        else
          vlong=-Zps(ik)/r*( p1*bberf(p2*r)+p3*bberf(p4*r) )
        end if
        vshort(i)=vql(i,ik)-vlong
      end do

      allocate( tmp(MMr) )

      do ig=1,NMGL
        G=sqrt(GG(ig))
        if ( G == 0.0d0 ) then
          do i=1,MMr
            tmp(i)=rad(i,ik)*rad(i,ik)*vshort(i)*rab(i,ik)
          end do
        else
          do i=1,MMr
            x=G*rad(i,ik)
            if ( x < 1.0d-1 ) then
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
        if ( G2 == 0.0d0 ) then
          vqlg(ig,ik)=vqlg(ig,ik)-(p1*p2+p3*p4)*const
        else
          vqlg(ig,ik)=vqlg(ig,ik)+(p1*exp(-G2*p2)+p3*exp(-G2*p4))/G2*const
        end if
      end do
      deallocate( tmp )

    end do ! ik

    deallocate( vshort )

! --- const_ps_local

    const_ps_local=0.0d0

    if ( flag_zero_ave ) then

      do ig=1,NMGL
        if ( GG(ig) == 0.0d0 ) exit
      end do
      do i=1,size(ki_atom)
        ik=ki_atom(i)
        const_ps_local=const_ps_local+vqlg(ig,ik)    
      end do
      vqlg(ig,:)=0.0d0

    end if

    call write_border( 0, " init_ps_local(end)" )

  end subroutine init_ps_local


  subroutine simp(f,s,m)
    implicit none
    integer,intent(in)  :: m
    real(8),intent(in)  :: f(:)
    real(8),intent(out) :: s
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
  end subroutine simp


  subroutine construct_ps_local
    implicit none

    call write_border( 0, " construct_ps_local(start)" )

    if ( .not.allocated(Vion) ) then
      allocate( Vion(Igrid(1,0):Igrid(2,0)) )
      Vion=0.0d0
    end if

    select case( iswitch_fft )
    case( 'FFTE', 'FFTE1' )

      call init_ps_local_ffte( Ngrid, Igrid )
      call construct_ps_local_ffte( vqlg, SGK, Vion )

    case( 'FFTE2' )

      call init_ps_local_1dffte( Ngrid, Igrid )
      call construct_ps_local_1dffte( vqlg, SGK, Vion )

    case( 'FFTW', 'FFTW1' )

      call construct_ps_local_fftw( vqlg, SGK, Vion )

    case default

      call construct_ps_local_fft0( vqlg, SGK, Vion )

    end select

    call write_border( 0, " construct_ps_local(end)" )

  end subroutine construct_ps_local

end module ps_local_module
