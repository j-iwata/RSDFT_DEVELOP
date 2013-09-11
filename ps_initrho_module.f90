MODULE ps_initrho_module

  use rgrid_module, only: dV,Ngrid,Igrid
  use ggrid_module
  use atom_module, only: Nelement
  use strfac_module, only: SGK
  use density_module, only: rho
  use pseudopot_module, only: Mr,rad,rab,cdd
  use electron_module, only: Nspin

  implicit none

  PRIVATE
  PUBLIC :: init_ps_initrho,construct_ps_initrho

  logical :: flag_initrho_0
  logical,allocatable :: flag_initrho(:)
  real(8),allocatable :: cddg(:,:)

CONTAINS


  SUBROUTINE init_ps_initrho
    implicit none
    integer :: i,ig,ik,MKI
    real(8) :: sum0,G,x,sb,const,Vcell
    real(8),allocatable :: tmp(:)

    allocate( flag_initrho(Nelement) )
    flag_initrho(:)=.false.
    flag_initrho_0 =.false.
    if ( allocated(cdd) ) then
       do ik=1,Nelement
          if ( all(cdd(:,ik)==0.d0) ) cycle
          flag_initrho(ik) = .true.
          flag_initrho_0   = .true.
       end do
    end if

    if ( .not.flag_initrho_0 ) return

    MKI   = Nelement
    Vcell = Ngrid(0)*dV
    const = 1.d0/Vcell
    allocate( cddg(NMGL,MKI) ) ; cddg=0.d0
    do ik=1,MKI
       allocate( tmp(Mr(ik)) )
       do ig=1,NMGL
          G=sqrt(GG(ig))
          if ( G == 0.d0 ) then
             do i=1,Mr(ik)
                tmp(i)=cdd(i,ik)*rab(i,ik)
             end do
          else
             do i=1,Mr(ik)
                x=G*rad(i,ik)
                if ( x<1.d-1 ) then
                   sb=-(1.d0/39916800.d0*x**10-1.d0/362880.d0*x**8 &
                   +1.d0/5040.d0*x**6-1.d0/120.d0*x**4+1.d0/6.d0*x**2-1.d0)
                else
                   sb=sin(x)/x
                end if
                tmp(i)=cdd(i,ik)*sb*rab(i,ik)
             end do
          end if
          call simp(tmp,sum0,2)
          cddg(ig,ik)=sum0*const
       end do
       deallocate( tmp )
    end do ! ik
  END SUBROUTINE init_ps_initrho

  SUBROUTINE simp(f,s,m)
    implicit none
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


  SUBROUTINE construct_ps_initrho
    implicit none
    integer :: a,i,i1,i2,i3,ik,j,MG
    integer :: ML1,ML2,ML3,ML,ML_0,ML_1
    integer :: ifacx(30),ifacy(30),ifacz(30)
    integer,allocatable :: lx1(:),lx2(:),ly1(:),ly2(:),lz1(:),lz2(:)
    complex(8),allocatable :: fftwork(:),zwork(:,:,:),vg(:)
    complex(8),allocatable :: wsavex(:),wsavey(:),wsavez(:)
    real(8) :: c

    if ( .not. flag_initrho_0 ) return

    MG   = NGgrid(0)
    ML   = Ngrid(0)
    ML1  = Ngrid(1)
    ML2  = Ngrid(2)
    ML3  = Ngrid(3)
    ML_0 = Igrid(1,0)
    ML_1 = Igrid(2,0)

    if ( .not. allocated(rho) ) then
       allocate( rho(ML_0:ML_1,Nspin) )
       rho=0.d0
    end if

    allocate( zwork(0:ML1-1,0:ML2-1,0:ML3-1) )

    allocate( vg(MG) )

    do i=MG_0,MG_1
       j=MGL(i)
       vg(i)=cddg(j,1)*SGK(i,1)
    end do
    do ik=2,Nelement
       do i=MG_0,MG_1
          j=MGL(i)
          vg(i)=vg(i)+cddg(j,ik)*SGK(i,ik)
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

    call fft3bx(ML1,ML2,ML3,ML,zwork,fftwork,wsavex,wsavey,wsavez &
         ,ifacx,ifacy,ifacz,lx1,lx2,ly1,ly2,lz1,lz2)

    rho(:,:)=0.d0
    i=ML_0-1
    do i3=Igrid(1,3),Igrid(2,3)
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
       i=i+1
       rho(i,1)=rho(i,1)+real( zwork(i1,i2,i3) )
    end do
    end do
    end do
    if ( minval(rho(:,1)) < 0.d0 ) then
       write(*,*) "WARNING: rho is negative at some points",minval(rho(:,1))
    end if
!    rho=abs(rho)
    where( rho < 0.d0 )
       rho=0.d0
    end where
    deallocate( wsavez,wsavey,wsavex )
    deallocate( lz2,lz1,ly2,ly1,lx2,lx1 )
    deallocate( fftwork )
    deallocate( zwork )

    if ( Nspin>1 ) then
       c=1.d0/Nspin
       rho(:,1)=c*rho(:,1)
       do i=2,Nspin
          rho(:,i)=rho(:,1)
       end do
    end if

  END SUBROUTINE construct_ps_initrho

END MODULE ps_initrho_module
