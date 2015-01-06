MODULE localpot2_iongth_module

  use pseudopot_module
  use ps_gth_module
  use aa_module
  use bb_module
  use atom_module

  use rgrid_module
  use parallel_module
  use localpot2_module
  use array_bound_module
  use electron_module
  use wf_module

  implicit none

  PRIVATE
  PUBLIC :: localpot2_ion, localpot2_calc_eion

CONTAINS


  SUBROUTINE localpot2_ion(m1,m2,m3,MKI,ecut_in,vout)
    implicit none
    integer,intent(IN)  :: m1,m2,m3,MKI
    real(8),intent(IN)  :: ecut_in
    real(8),intent(OUT) :: vout(0:m1-1,0:m2-1,0:m3-1)

    real(8) :: pi,rloc,const,C1,C2,C3,C4,G,v,ecut
    integer :: ig,ielm,mm1,mm2,mm3

    real(8),allocatable :: vqlg(:,:,:,:)
    complex(8),allocatable :: sgk(:,:,:,:)
    integer :: ifacx(30),ifacy(30),ifacz(30)
    integer,allocatable :: lx1(:),lx2(:),ly1(:),ly2(:),lz1(:),lz2(:)
    complex(8),allocatable :: wsavex(:),wsavey(:),wsavez(:)
    complex(8),allocatable :: zw0(:,:,:),zw1(:,:,:)
    integer :: i1,i2,i3,j1,j2,j3,mm,a
    real(8) :: Gx,Gy,Gz,GG,Gr,a1,a2,a3,pi2

    pi = acos(-1.d0)
    const = 1.d0/abs(Va)

    ecut = ecut_in*fecut_loc

! - structure factor

    mm1 = (m1-1)/2
    mm2 = (m2-1)/2
    mm3 = (m3-1)/2

    allocate( sgk(-mm1:mm1,-mm2:mm2,-mm3:mm3,MKI) )
    sgk=(0.0d0,0.0d0)

    pi2 = 2.0d0*acos(-1.0d0)

    do a=1,Natom
       ielm=ki_atom(a)
       a1=pi2*aa_atom(1,a)
       a2=pi2*aa_atom(2,a)
       a3=pi2*aa_atom(3,a)
       do i3=-mm3,mm3
       do i2=-mm2,mm2
       do i1=-mm1,mm1
          Gr=i1*a1+i2*a2+i3*a3
          sgk(i1,i2,i3,ielm)=sgk(i1,i2,i3,ielm)+dcmplx(cos(Gr),-sin(Gr))
       end do
       end do
       end do
    end do

! - Fourier component

    allocate( vqlg(-mm1:mm1,-mm2:mm2,-mm3:mm3,MKI) )
    vqlg = 0.0d0

    do ielm=1,MKI

       rloc = Rcloc(ielm)
       C1   = parloc(1,ielm)
       C2   = parloc(2,ielm)
       C3   = parloc(3,ielm)
       C4   = parloc(4,ielm)

       do i3=-mm3,mm3
       do i2=-mm2,mm2
       do i1=-mm1,mm1

          Gx = bb(1,1)*i1 + bb(1,2)*i2 + bb(1,3)*i3
          Gy = bb(2,1)*i1 + bb(2,2)*i2 + bb(2,3)*i3
          Gz = bb(3,1)*i1 + bb(3,2)*i2 + bb(3,3)*i3

          GG = Gx*Gx + Gy*Gy + Gz*Gz

          if ( GG > ecut ) cycle

          G = sqrt(GG)

          if ( G == 0.0d0 ) then

             vqlg(i1,i2,i3,ielm) = const*( 2.d0*pi*Zps(ielm)*rloc**2 &
                  + sqrt((2.d0*pi)**3)*rloc**3*( C1+C2*3.d0+C3*15.d0+C4*105.d0 ) )

          else

             v = -4.d0*pi*Zps(ielm)*exp(-0.5d0*(G*rloc)**2)/G**2 &
               + sqrt((2.d0*pi)**3)*rloc**3*exp(-0.5d0*(rloc*G)**2) &
               *( C1 &
                 +C2*(3.d0-(rloc*G)**2) &
                 +C3*(15.d0-10.d0*(rloc*G)**2+(rloc*G)**4) &
                 +C4*(105.d0-105.d0*(rloc*G)**2+21.d0*(rloc*G)**4-(rloc*G)**6) )

             vqlg(i1,i2,i3,ielm) = const*v

          end if

       end do ! i1
       end do ! i2
       end do ! i3

    end do ! ielm

! -

    allocate( zw0(0:m1-1,0:m2-1,0:m3-1) )
    allocate( zw1(0:m1-1,0:m2-1,0:m3-1) )
    zw0=(0.0d0,0.0d0)
    do ielm=1,MKI
       do i3=-mm3,mm3
       do i2=-mm2,mm2
       do i1=-mm1,mm1
          j1=mod(i1+m1,m1)
          j2=mod(i2+m2,m2)
          j3=mod(i3+m3,m3)
          zw0(j1,j2,j3)=zw0(j1,j2,j3)+vqlg(i1,i2,i3,ielm)*sgk(i1,i2,i3,ielm)
       end do
       end do
       end do
    end do

    mm=m1*m2*m3
    allocate( lx1(mm),lx2(mm),ly1(mm),ly2(mm),lz1(mm),lz2(mm) )
    allocate( wsavex(m1),wsavey(m2),wsavez(m3) )

    call prefft(m1,m2,m3,mm,wsavex,wsavey,wsavez,ifacx,ifacy,ifacz &
               ,lx1,lx2,ly1,ly2,lz1,lz2)
    call fft3bx(m1,m2,m3,mm,zw0,zw1,wsavex,wsavey,wsavez,ifacx,ifacy,ifacz &
               ,lx1,lx2,ly1,ly2,lz1,lz2)

    deallocate( wsavez,wsavey,wsavex )
    deallocate( lz2,lz1,ly2,ly1,lx2,lx1 )

    do i3=0,m3-1
    do i2=0,m2-1
    do i1=0,m1-1
       vout(i1,i2,i3) = real( zw0(i1,i2,i3) )
    end do
    end do
    end do

    deallocate( zw1,zw0 )

    return

  END SUBROUTINE localpot2_ion


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
          sum0 = sum0 + conjg(unk(i,n,k,s))*vloc_nl(j,i)*unk(Lpot(j,i),n,k,s)
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

END MODULE localpot2_iongth_module
