MODULE fock_fft_module

  use parallel_module
  use rgrid_module, only: Ngrid,Igrid
  use ggrid_module, only: NGgrid,Ecut
  use bb_module, only: bb
  use xc_hybrid_module
  use watch_module

  implicit none

  PRIVATE
  PUBLIC :: fock_fft,ct_fock_fft,et_focK_fft

  integer,parameter :: TYPE_MAIN=MPI_COMPLEX16

  real(8) :: ct_fock_fft(10),et_fock_fft(10)

CONTAINS

  SUBROUTINE fock_fft(n1,n2,k_fock,q_fock,trho,tVh,tr)
    implicit none
    integer,intent(IN) :: n1,n2,tr
    real(8),intent(IN) :: k_fock(3),q_fock(3)
    complex(8),intent(IN)    :: trho(n1:n2)
    complex(8),intent(INOUT) :: tVh(n1:n2)
    integer :: i,i1,i2,i3,j1,j2,j3,ierr,irank,a1,a2,a3,b1,b2,b3
    integer :: ML1,ML2,ML3,ML
    real(8) :: pi4,mem,memax,ct0,ct1,et0,et1,g2,pi,const1,const2
    real(8) :: ctt0,ett0
    complex(8),allocatable :: work(:)
    complex(8),parameter :: z0=(0.d0,0.d0)
    complex(8),allocatable :: zwork0(:,:,:),zwork1(:,:,:)
    integer :: ifacx(30),ifacy(30),ifacz(30)
    integer,allocatable :: lx1(:),lx2(:),ly1(:),ly2(:),lz1(:),lz2(:)
    complex(8),allocatable :: wsavex(:),wsavey(:),wsavez(:)

    pi  = acos(-1.0d0)
    pi4 = 4.d0*pi
    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)
    ML  = Ngrid(0)

!
! ---
!

    call watch(ct0,et0)
    ctt0=ct0
    ett0=et0

    allocate( work(ML) ) ; work=z0

    call mpi_allgatherv(trho(n1),ir_grid(myrank_g),TYPE_MAIN &
         ,work,ir_grid,id_grid,TYPE_MAIN,comm_grid,ierr)

    allocate( zwork0(0:ML1-1,0:ML2-1,0:ML3-1) ) ; zwork0=z0

    i=0
    irank=-1
    do i3=1,node_partition(3)
    do i2=1,node_partition(2)
    do i1=1,node_partition(1)
       irank=irank+1
       a1=pinfo_grid(1,irank) ; b1=pinfo_grid(2,irank)+a1-1
       a2=pinfo_grid(3,irank) ; b2=pinfo_grid(4,irank)+a2-1
       a3=pinfo_grid(5,irank) ; b3=pinfo_grid(6,irank)+a3-1
       do j3=a3,b3
       do j2=a2,b2
       do j1=a1,b1
          i=i+1
          zwork0(j1,j2,j3)=work(i)
       end do
       end do
       end do
    end do
    end do
    end do

    deallocate( work )

    call watch(ct1,et1)
    ct_fock_fft(1) = ct_fock_fft(1) + ct1-ct0
    et_fock_fft(1) = et_fock_fft(1) + et1-et0

!
! ---
!

    allocate( zwork1(0:ML1-1,0:ML2-1,0:ML3-1) ) ; zwork1=z0

    allocate( lx1(ML),lx2(ML),ly1(ML),ly2(ML),lz1(ML),lz2(ML) )
    allocate( wsavex(ML1),wsavey(ML2),wsavez(ML3) )

    call watch(ct0,et0)
    ct_fock_fft(6) = ct_fock_fft(6) + ct0-ct1
    et_fock_fft(6) = et_fock_fft(6) + et0-et1

    call prefft(ML1,ML2,ML3,ML,wsavex,wsavey,wsavez &
         ,ifacx,ifacy,ifacz,lx1,lx2,ly1,ly2,lz1,lz2)

    call watch(ct1,et1)
    ct_fock_fft(7) = ct_fock_fft(7) + ct1-ct0
    et_fock_fft(7) = et_fock_fft(7) + et1-et0

    call watch(ct0,et0)

    call fft3fx(ML1,ML2,ML3,ML,zwork0,zwork1,wsavex,wsavey,wsavez &
         ,ifacx,ifacy,ifacz,lx1,lx2,ly1,ly2,lz1,lz2)

    call watch(ct1,et1)
    ct_fock_fft(2) = ct_fock_fft(2) + ct1-ct0
    et_fock_fft(2) = et_fock_fft(2) + et1-et0

!    call construct_Ggrid(2)

    const1 = 0.25d0/(omega*omega)
    const2 = pi/(omega*omega)

    zwork1(:,:,:)=z0
!    do i=1,NGgrid(0)
!       i1=LLG(1,i)
!       i2=LLG(2,i)
!       i3=LLG(3,i)
!       g2=GG(MGL(i))
    do i3=-NGgrid(3),NGgrid(3)
    do i2=-NGgrid(2),NGgrid(2)
    do i1=-NGgrid(1),NGgrid(1)
       g2=(bb(1,1)*i1+bb(1,2)*i2+bb(1,3)*i3+k_fock(1)-q_fock(1))**2 &
&        +(bb(2,1)*i1+bb(2,2)*i2+bb(2,3)*i3+k_fock(2)-q_fock(2))**2 &
&        +(bb(3,1)*i1+bb(3,2)*i2+bb(3,3)*i3+k_fock(3)-q_fock(3))**2
       if ( g2 < 0.0d0 .or. g2 > Ecut ) cycle
       j1=mod(i1+ML1,ML1)
       j2=mod(i2+ML2,ML2)
       j3=mod(i3+ML3,ML3)
!       if ( iflag_hf > 0 .or. iflag_pbe0 > 0 ) then
!          if ( g2 <= 1.d-10 ) then
!             zwork1(j1,j2,j3)=zwork0(j1,j2,j3)*2.d0*Pi*R_hf**2.d0
!          else
!             zwork1(j1,j2,j3)=zwork0(j1,j2,j3)*pi4*(1.d0-cos(sqrt(g2)*R_hf))/g2
!          end if
!       end if
!       if ( iflag_hse > 0 ) then
          if ( g2 <= 1.d-10 ) then
!             zwork1(j1,j2,j3)=zwork0(j1,j2,j3)*Pi/(omega*omega)
             zwork1(j1,j2,j3)=zwork0(j1,j2,j3)*const2
          else
!             zwork1(j1,j2,j3)=zwork0(j1,j2,j3) &
!                  *pi4*( 1.0d0 - exp(-0.25d0*g2/(omega*omega)) )/g2
             zwork1(j1,j2,j3)=zwork0(j1,j2,j3)*pi4*( 1.0d0 - exp(-const1*g2) )/g2
          end if
!       end if
!       if ( iflag_lcwpbe > 0 ) then
!          if ( g2 <= 1.d-10 ) then
!             zwork1(j1,j2,j3)=zwork0(j1,j2,j3) &
!                  *(2.d0*Pi*R_hf*R_hf-Pi/(omega*omega))
!          else
!             zwork1(j1,j2,j3)=zwork0(j1,j2,j3) &
!                  *pi4*(exp(-0.25d0*g2/(omega*omega))-cos(sqrt(g2)*R_hf))/g2
!          end if
!       end if
    end do
    end do
    end do

!    call destruct_Ggrid

    call watch(ct0,et0)
    ct_fock_fft(3) = ct_fock_fft(3) + ct0-ct1
    et_fock_fft(3) = et_fock_fft(3) + et0-et1

    call fft3bx(ML1,ML2,ML3,ML,zwork1,zwork0,wsavex,wsavey,wsavez &
         ,ifacx,ifacy,ifacz,lx1,lx2,ly1,ly2,lz1,lz2)

    call watch(ct1,et1)
    ct_fock_fft(4) = ct_fock_fft(4) + ct1-ct0
    et_fock_fft(4) = et_fock_fft(4) + et1-et0

    deallocate( wsavez,wsavey,wsavex )
    deallocate( lz2,lz1,ly2,ly1,lx2,lx1 )
    deallocate( zwork0 )

    i=n1-1
    do i3=Igrid(1,3),Igrid(2,3)
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
       i=i+1
       tVh(i)=zwork1(i1,i2,i3)
    end do
    end do
    end do

    deallocate( zwork1 )

    call watch(ct0,et0)
    ct_fock_fft(5) = ct_fock_fft(5) + ct0-ct1
    et_fock_fft(5) = et_fock_fft(5) + et0-et1

    ct_fock_fft(9) = ct_fock_fft(9) + ct0-ctt0
    et_fock_fft(9) = et_fock_fft(9) + et0-ett0

!    if ( DISP_SWITCH_PARALLEL ) then
!       write(*,*) "TIME(FOCK_FFT)=",ctime1-ctime0,etime1-etime0
!    end if

    return
  END SUBROUTINE fock_fft

END MODULE fock_fft_module
