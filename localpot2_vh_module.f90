MODULE localpot2_vh_module

  use rgrid_module
  use bb_module
  use localpot2_variables, only: fecut_loc, Ngrid_dense, Igrid_dense
  use parallel_module
  use watch_module

  implicit none

  PRIVATE
  PUBLIC :: localpot2_vh

CONTAINS

  SUBROUTINE localpot2_vh(m1_0,m2_0,m3_0,ecut_in,n_in,vout,eh)
    implicit none
    integer,intent(IN) :: m1_0,m2_0,m3_0
    real(8),intent(IN) :: ecut_in,n_in(m1_0,m2_0,m3_0)
    real(8),intent(OUT) :: vout(m1_0,m2_0,m3_0),eh
    integer :: ifacx(30),ifacy(30),ifacz(30),mm,m1,m2,m3
    integer,allocatable :: lx1(:),lx2(:),ly1(:),ly2(:),lz1(:),lz2(:)
    complex(8),allocatable :: wsavex(:),wsavey(:),wsavez(:)
    complex(8),allocatable :: zw0(:,:,:),zw1(:,:,:)
    integer :: i1,i2,i3,j1,j2,j3,mm1,mm2,mm3,ierr
    integer :: i10,i20,i30
    real(8) :: Gx,Gy,Gz,GG,pi4,c,ecut,et0,et1,ct0,ct1

    call watch(ct0,et0)

    m1 = Ngrid_dense(1)
    m2 = Ngrid_dense(2)
    m3 = Ngrid_dense(3)
    mm = m1*m2*m3

    eh=0.0d0
    vout(:,:,:)=0.0d0

    ecut = ecut_in*fecut_loc
    pi4 = 4.0d0*acos(-1.0d0)

    i10 = Igrid_dense(1,1)-1
    i20 = Igrid_dense(1,2)-1
    i30 = Igrid_dense(1,3)-1

    allocate( zw0(0:m1-1,0:m2-1,0:m3-1) ) ; zw0=(0.0d0,0.0d0)
    allocate( zw1(0:m1-1,0:m2-1,0:m3-1) ) ; zw1=(0.0d0,0.0d0)

    do i3=Igrid_dense(1,3),Igrid_dense(2,3)
    do i2=Igrid_dense(1,2),Igrid_dense(2,2)
    do i1=Igrid_dense(1,1),Igrid_dense(2,1)
       zw1(i1,i2,i3)=n_in(i1-i10,i2-i20,i3-i30)
    end do
    end do
    end do

    call mpi_allreduce(zw1,zw0,mm,mpi_complex16,mpi_sum,comm_grid,ierr)

    allocate( lx1(mm),lx2(mm),ly1(mm),ly2(mm),lz1(mm),lz2(mm) )
    allocate( wsavex(m1),wsavey(m2),wsavez(m3) )

    call prefft(m1,m2,m3,mm,wsavex,wsavey,wsavez,ifacx,ifacy,ifacz &
               ,lx1,lx2,ly1,ly2,lz1,lz2)

    call fft3fx(m1,m2,m3,mm,zw0,zw1,wsavex,wsavey,wsavez,ifacx,ifacy,ifacz &
               ,lx1,lx2,ly1,ly2,lz1,lz2)

    mm1=(m1-1)/2
    mm2=(m2-1)/2
    mm3=(m3-1)/2

    zw1(:,:,:)=(0.0d0,0.0d0)
    do i3=-mm3,mm3
    do i2=-mm2,mm2
    do i1=-mm1,mm1

       Gx = bb(1,1)*i1 + bb(1,2)*i2 + bb(1,3)*i3
       Gy = bb(2,1)*i1 + bb(2,2)*i2 + bb(2,3)*i3
       Gz = bb(3,1)*i1 + bb(3,2)*i2 + bb(3,3)*i3
       GG = Gx*Gx + Gy*Gy + Gz*Gz

       if ( GG > ecut .or. GG <= 0.0d0 ) cycle

       j1 = mod(i1+m1,m1)
       j2 = mod(i2+m2,m2)
       j3 = mod(i3+m3,m3)

       zw1(j1,j2,j3) = zw0(j1,j2,j3)*pi4/GG

    end do
    end do
    end do

    call fft3bx(m1,m2,m3,mm,zw1,zw0,wsavex,wsavey,wsavez,ifacx,ifacy,ifacz &
               ,lx1,lx2,ly1,ly2,lz1,lz2)

    do i3=Igrid_dense(1,3),Igrid_dense(2,3)
    do i2=Igrid_dense(1,2),Igrid_dense(2,2)
    do i1=Igrid_dense(1,1),Igrid_dense(2,1)
       vout(i1-i10,i2-i20,i3-i30) = real( zw1(i1,i2,i3) ,8 )
    end do
    end do
    end do

    deallocate( wsavez,wsavey,wsavex )
    deallocate( lz2,lz1,ly2,ly1,lx2,lx1 )

    deallocate( zw1,zw0 )

    c=0.5d0*sum( vout*n_in )*dV*Ngrid(0)/mm
    call mpi_allreduce(c,eh,1,mpi_real8,mpi_sum,comm_grid,ierr)

    call watch(ct1,et1) ; write(*,*) "localpot2_vh",ct1-ct0,et1-et0

  END SUBROUTINE localpot2_vh

END MODULE localpot2_vh_module
