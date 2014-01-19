MODULE fock_fft_module

!$ use omp_lib
  use parallel_module
  use rgrid_module, only: Ngrid,Igrid
  use ggrid_module, only: NGgrid,Ecut
  use bb_module, only: bb
  use xc_hybrid_module
  use watch_module
  use bz_module
  use omp_variables

  implicit none

  PRIVATE
  PUBLIC :: fock_fft,ct_fock_fft,et_focK_fft

  integer,parameter :: TYPE_MAIN=MPI_COMPLEX16

  real(8) :: ct_fock_fft(10),et_fock_fft(10)

  integer :: comm_fftx,comm_ffty,comm_fftz
  integer :: npux,npuy,npuz
  real(8) :: pi,pi4
  integer :: ML1,ML2,ML3
  integer :: a1b,b1b,a2b,b2b,a3b,b3b
  real(8) :: const1,const2
  complex(8),allocatable :: zwork0(:,:,:),zwork1(:,:,:)
  logical :: first_time=.true.

  integer :: MG_fock
  integer,allocatable :: LG_fock(:,:)
  real(8),allocatable :: Gi_fock(:,:)

CONTAINS

  SUBROUTINE Fock_fft(n1,n2,k_fock,q_fock,trho,tVh,tr)
    implicit none
    integer,intent(IN) :: n1,n2,tr
    real(8),intent(IN) :: k_fock(3),q_fock(3)
    complex(8),intent(IN)    :: trho(n1:n2)
    complex(8),intent(INOUT) :: tVh(n1:n2)
    integer :: i,j,i1,i2,i3,j1,j2,j3,ierr
    integer :: a1,a2,a3,b1,b2,b3,m
    real(8) :: ct0,ct1,et0,et1,g2,gx,gy,gz

    call watch(ct0,et0)

    if ( first_time ) then
       call init_fock_fft
!       call init2_fock_fft
    end if

    call watch(ct1,et1)
    ct_fock_fft(3) = ct_fock_fft(3) + ct1-ct0
    et_fock_fft(3) = et_fock_fft(3) + et1-et0

!$OMP parallel private(i,m)
    m=0
!$  m=omp_get_thread_num()
    i=n1_omp(m)-1
    do i3=a3b_omp(m),b3b_omp(m)
    do i2=a2b,b2b
       do i1=0,a1b-1
          zwork0(i1,i2,i3)=(0.0d0,0.0d0)
       end do
       do i1=a1b,b1b
          i=i+1
          zwork0(i1,i2,i3)=trho(i)
       end do
       do i1=b1b+1,ML1-1
          zwork0(i1,i2,i3)=(0.0d0,0.0d0)
       end do
    end do
    end do
!$OMP end parallel

    call watch(ct0,et0)
    ct_fock_fft(4) = ct_fock_fft(4) + ct0-ct1
    et_fock_fft(4) = et_fock_fft(4) + et0-et1

    call mpi_allreduce(zwork0,zwork1,ML1*(b2b-a2b+1)*(b3b-a3b+1) &
         ,TYPE_MAIN,mpi_sum,comm_fftx,ierr)

    call watch(ct1,et1)
    ct_fock_fft(5) = ct_fock_fft(5) + ct1-ct0
    et_fock_fft(5) = et_fock_fft(5) + et1-et0

    call pzfft3dv(zwork1,zwork0,ML1,ML2,ML3,comm_ffty,comm_fftz,npuy,npuz,-1)

    call watch(ct0,et0)
    ct_fock_fft(6) = ct_fock_fft(6) + ct0-ct1
    et_fock_fft(6) = et_fock_fft(6) + et0-et1

!    do i3=-NGgrid(3),NGgrid(3)
!       j3=mod(i3+ML3,ML3)
!       write(*,*) i3,j3,j3-j3/(NGgrid(3)+1)*ML3
!    end do
!    stop

    goto 1
    zwork1(:,:,:)=(0.0d0,0.0d0)
    do i3=-NGgrid(3),NGgrid(3)
    do i2=-NGgrid(2),NGgrid(2)
    do i1=-NGgrid(1),NGgrid(1)
       g2=(bb(1,1)*i1+bb(1,2)*i2+bb(1,3)*i3+k_fock(1)-q_fock(1))**2 &
&        +(bb(2,1)*i1+bb(2,2)*i2+bb(2,3)*i3+k_fock(2)-q_fock(2))**2 &
&        +(bb(3,1)*i1+bb(3,2)*i2+bb(3,3)*i3+k_fock(3)-q_fock(3))**2
       if ( g2 < 0.d0 .or. g2 > Ecut ) cycle
       j1=mod(i1+ML1,ML1)
       j2=mod(i2+ML2,ML2)
       j3=mod(i3+ML3,ML3)
       if ( a2b <= j2 .and. j2 <= b2b .and. a3b <= j3 .and. j3 <= b3b ) then
          if ( g2 <= 1.d-10 ) then
             zwork1(j1,j2,j3)=zwork0(j1,j2,j3)*const2
          else
             zwork1(j1,j2,j3)=zwork0(j1,j2,j3)*pi4*( 1.d0-exp(-g2*const1) )/g2
          end if
       end if
    end do
    end do
    end do
1 continue
!    goto 2
!$OMP parallel do private(i1,i2,i3,gx,gy,gz,g2)
    do j3=a3b,b3b
    do j2=a2b,b2b
    do j1=0,ML1-1
       i1=j1-j1/(NGgrid(1)+1)*ML1
       i2=j2-j2/(NGgrid(2)+1)*ML2
       i3=j3-j3/(NGgrid(3)+1)*ML3
!       g2=(bb(1,1)*i1+bb(1,2)*i2+bb(1,3)*i3+k_fock(1)-q_fock(1))**2 &
!         +(bb(2,1)*i1+bb(2,2)*i2+bb(2,3)*i3+k_fock(2)-q_fock(2))**2 &
!         +(bb(3,1)*i1+bb(3,2)*i2+bb(3,3)*i3+k_fock(3)-q_fock(3))**2
       gx = bb(1,1)*i1+bb(1,2)*i2+bb(1,3)*i3+k_fock(1)-q_fock(1)
       gy = bb(2,1)*i1+bb(2,2)*i2+bb(2,3)*i3+k_fock(2)-q_fock(2)
       gz = bb(3,1)*i1+bb(3,2)*i2+bb(3,3)*i3+k_fock(3)-q_fock(3)
       g2 = gx*gx + gy*gy + gz*gz
       if ( g2 > Ecut ) then
          zwork1(j1,j2,j3)=(0.0d0,0.0d0)
          cycle
       else if ( g2 <= 1.d-10 ) then
          zwork1(j1,j2,j3)=zwork0(j1,j2,j3)*const2
       else
          zwork1(j1,j2,j3)=zwork0(j1,j2,j3)*pi4*( 1.d0-exp(-g2*const1) )/g2
       end if
    end do
    end do
    end do
!$OMP end parallel do
2 continue
    goto 3
    zwork1(:,:,:)=(0.0d0,0.0d0)
!    if ( a2b == 0 .and. a3b == 0 ) then
!       if ( all(k_fock==q_fock) ) then
!          zwork1(0,0,0) = zwork0(0,0,0)*const2
!       end if
!    end if
    do i=1,MG_fock
       j1=LG_fock(1,i)
       j2=LG_fock(2,i)
       j3=LG_fock(3,i)
!       i1=j1-j1/(NGgrid(1)+1)*ML1
!       i2=j2-j2/(NGgrid(2)+1)*ML2
!       i3=j3-j3/(NGgrid(3)+1)*ML3
!       g2=(bb(1,1)*i1+bb(1,2)*i2+bb(1,3)*i3+k_fock(1)-q_fock(1))**2 &
!         +(bb(2,1)*i1+bb(2,2)*i2+bb(2,3)*i3+k_fock(2)-q_fock(2))**2 &
!         +(bb(3,1)*i1+bb(3,2)*i2+bb(3,3)*i3+k_fock(3)-q_fock(3))**2
!       gx = bb(1,1)*i1+bb(1,2)*i2+bb(1,3)*i3+k_fock(1)-q_fock(1)
!       gy = bb(2,1)*i1+bb(2,2)*i2+bb(2,3)*i3+k_fock(2)-q_fock(2)
!       gz = bb(3,1)*i1+bb(3,2)*i2+bb(3,3)*i3+k_fock(3)-q_fock(3)
       gx = Gi_fock(1,i) + k_fock(1)-q_fock(1)
       gy = Gi_fock(2,i) + k_fock(2)-q_fock(2)
       gz = Gi_fock(3,i) + k_fock(3)-q_fock(3)
       g2 = gx*gx + gy*gy + gz*gz
       if ( g2 <= 1.d-10 ) then
          zwork1(j1,j2,j3)=zwork0(j1,j2,j3)*const2
       else
          zwork1(j1,j2,j3)=zwork0(j1,j2,j3)*pi4*( 1.d0-exp(-g2*const1) )/g2
       end if
    end do
3 continue

    call watch(ct1,et1)
    ct_fock_fft(7) = ct_fock_fft(7) + ct1-ct0
    et_fock_fft(7) = et_fock_fft(7) + et1-et0

    call pzfft3dv(zwork1,zwork0,ML1,ML2,ML3,comm_ffty,comm_fftz,npuy,npuz,1)

    call watch(ct0,et0)
    ct_fock_fft(8) = ct_fock_fft(8) + ct0-ct1
    et_fock_fft(8) = et_fock_fft(8) + et0-et1

!$OMP parallel private(i,m)
    m=0
!$  m=omp_get_thread_num()
    i=n1_omp(m)-1
    do i3=a3b_omp(m),b3b_omp(m)
    do i2=a2b,b2b
    do i1=a1b,b1b
       i=i+1
       tVh(i)=zwork0(i1,i2,i3)
    end do
    end do
    end do
!$OMP end parallel

    call watch(ct1,et1)
    ct_fock_fft(9) = ct_fock_fft(9) + ct1-ct0
    et_fock_fft(9) = et_fock_fft(9) + et1-et0

    return

  END SUBROUTINE Fock_fft


  SUBROUTINE init_fock_fft
    implicit none
    integer :: ML_check(3),ML_switch
    integer :: i1,i2,i3,irank,i,j,ierr
    integer :: icolor1,icolor2,icolor3
    first_time = .false.
    pi     = acos(-1.0d0)
    pi4    = 4*pi
    const1 = 0.25d0/omega**2
    const2 = pi/omega**2
    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)
    a1b = Igrid(1,1)
    b1b = Igrid(2,1)
    a2b = Igrid(1,2)
    b2b = Igrid(2,2)
    a3b = Igrid(1,3)
    b3b = Igrid(2,3)
    ML_check(1)=ML1
    ML_check(2)=ML2
    ML_check(3)=ML3
    do i=1,3
       do j=1,100
          if ( j > 1 .and. ML_check(i) == 1 ) exit
          if ( mod(ML_check(i),2) == 0 ) then
             ML_check(i)=ML_check(i)/2
             cycle
          else if ( mod(ML_check(i),3) == 0 ) then
             ML_check(i)=ML_check(i)/3
             cycle
          else if ( mod(ML_check(i),5) == 0 ) then
             ML_check(i)=ML_check(i)/5
             cycle
          end if
          write(*,*) "If you would like to use the FFT parallel version,"
          write(*,*) "ML must be multiples of 2, 3, 5 without 3^n, 5^n and 3^nx5^n (n:integer)."
          write(*,*) "ML1, ML2, ML3",ML1,ML2,ML3,ML_check
          stop "fock_fft_parallel"
       end do
    end do

    irank=-1
    do i3=1,node_partition(3)
    do i2=1,node_partition(2)
    do i1=1,node_partition(1)
       irank=irank+1
       if ( irank == myrank_g ) then
          icolor1 = i2 + (i3-1)*node_partition(2)
          icolor2 = i1 + (i3-1)*node_partition(1)
          icolor3 = i1 + (i2-1)*node_partition(1)
       end if
    end do
    end do
    end do

    call mpi_comm_split(comm_grid,icolor1,0,comm_fftx,ierr)
    call mpi_comm_split(comm_grid,icolor2,0,comm_ffty,ierr)
    call mpi_comm_split(comm_grid,icolor3,0,comm_fftz,ierr)

    call mpi_comm_size(comm_fftx,npux,ierr)
    call mpi_comm_size(comm_ffty,npuy,ierr)
    call mpi_comm_size(comm_fftz,npuz,ierr)

!---
!    call mpi_comm_rank(comm_fftx,myrank_x,ierr)
!    allocate( ir(0:npux-1) ) ; ir=0
!    allocate( id(0:npux-1) ) ; id=0

    allocate( zwork0(0:ML1-1,a2b:b2b,a3b:b3b) ) ; zwork0=(0.0d0,0.0d0)
    allocate( zwork1(0:ML1-1,a2b:b2b,a3b:b3b) ) ; zwork1=(0.0d0,0.0d0)

    call pzfft3dv(zwork1,zwork0,ML1,ML2,ML3,comm_ffty,comm_fftz,npuy,npuz,0)

  END SUBROUTINE init_fock_fft


  SUBROUTINE init2_fock_fft
    implicit none
    integer :: i1,i2,i3,k,q,j1,j2,j3,i
    real(8) :: g2,k_fock(3),q_fock(3)
    integer,allocatable :: ichk(:,:,:)
    allocate( ichk(0:ML1-1,a2b:b2b,a3b:b3b) ) ; ichk=0
    do k=1,Nbzsm
       k_fock(1)=bb(1,1)*kbb(1,k)+bb(1,2)*kbb(2,k)+bb(1,3)*kbb(3,k)
       k_fock(2)=bb(2,1)*kbb(1,k)+bb(2,2)*kbb(2,k)+bb(2,3)*kbb(3,k)
       k_fock(3)=bb(3,1)*kbb(1,k)+bb(3,2)*kbb(2,k)+bb(3,3)*kbb(3,k)
       do q=1,k
          q_fock(1)=bb(1,1)*kbb(1,q)+bb(1,2)*kbb(2,q)+bb(1,3)*kbb(3,q)
          q_fock(2)=bb(2,1)*kbb(1,q)+bb(2,2)*kbb(2,q)+bb(2,3)*kbb(3,q)
          q_fock(3)=bb(3,1)*kbb(1,q)+bb(3,2)*kbb(2,q)+bb(3,3)*kbb(3,q)
          do j3=a3b,b3b
          do j2=a2b,b2b
          do j1=0,ML1-1
             i1=j1-j1/(NGgrid(1)+1)*ML1
             i2=j2-j2/(NGgrid(2)+1)*ML2
             i3=j3-j3/(NGgrid(3)+1)*ML3
             g2=( bb(1,1)*i1+bb(1,2)*i2+bb(1,3)*i3+k_fock(1)-q_fock(1) )**2 &
               +( bb(2,1)*i1+bb(2,2)*i2+bb(2,3)*i3+k_fock(2)-q_fock(2) )**2 &
               +( bb(3,1)*i1+bb(3,2)*i2+bb(3,3)*i3+k_fock(3)-q_fock(3) )**2
             if ( g2 <= Ecut ) ichk(j1,j2,j3)=ichk(j1,j2,j3)+1
             g2=( bb(1,1)*i1+bb(1,2)*i2+bb(1,3)*i3-k_fock(1)-q_fock(1) )**2 &
               +( bb(2,1)*i1+bb(2,2)*i2+bb(2,3)*i3-k_fock(2)-q_fock(2) )**2 &
               +( bb(3,1)*i1+bb(3,2)*i2+bb(3,3)*i3-k_fock(3)-q_fock(3) )**2
             if ( g2 <= Ecut ) ichk(j1,j2,j3)=ichk(j1,j2,j3)+1
             g2=( bb(1,1)*i1+bb(1,2)*i2+bb(1,3)*i3+k_fock(1)+q_fock(1) )**2 &
               +( bb(2,1)*i1+bb(2,2)*i2+bb(2,3)*i3+k_fock(2)+q_fock(2) )**2 &
               +( bb(3,1)*i1+bb(3,2)*i2+bb(3,3)*i3+k_fock(3)+q_fock(3) )**2
             if ( g2 <= Ecut ) ichk(j1,j2,j3)=ichk(j1,j2,j3)+1
             g2=( bb(1,1)*i1+bb(1,2)*i2+bb(1,3)*i3-k_fock(1)+q_fock(1) )**2 &
               +( bb(2,1)*i1+bb(2,2)*i2+bb(2,3)*i3-k_fock(2)+q_fock(2) )**2 &
               +( bb(3,1)*i1+bb(3,2)*i2+bb(3,3)*i3-k_fock(3)+q_fock(3) )**2
             if ( g2 <= Ecut ) ichk(j1,j2,j3)=ichk(j1,j2,j3)+1
          end do
          end do
          end do
       end do
    end do

!    if ( a2b == 0 .and. a3b == 0 ) ichk(0,0,0)=0

    i=0
    do i3=a3b,b3b
    do i2=a2b,b2b
    do i1=0,ML1-1
       if ( ichk(i1,i2,i3) > 0 ) then
          i=i+1
       end if
    end do
    end do
    end do
    MG_fock=i
    allocate( LG_fock(3,MG_fock) ) ; LG_fock=0
    allocate( Gi_fock(3,MG_fock) ) ; Gi_fock=0.0d0
    i=0
    do i3=a3b,b3b
    do i2=a2b,b2b
    do i1=0,ML1-1
       if ( ichk(i1,i2,i3) > 0 ) then
          i=i+1
          LG_fock(1,i)=i1
          LG_fock(2,i)=i2
          LG_fock(3,i)=i3
          j1=i1-i1/(NGgrid(1)+1)*ML1
          j2=i2-i2/(NGgrid(2)+1)*ML2
          j3=i3-i3/(NGgrid(3)+1)*ML3
          Gi_fock(1,i)=bb(1,1)*j1+bb(1,2)*j2+bb(1,3)*j3
          Gi_fock(2,i)=bb(2,1)*j1+bb(2,2)*j2+bb(2,3)*j3
          Gi_fock(3,i)=bb(3,1)*j1+bb(3,2)*j2+bb(3,3)*j3
       end if
    end do
    end do
    end do
    deallocate( ichk )
  END SUBROUTINE init2_fock_fft


  SUBROUTINE finalize_fock_fft
    implicit none
    integer :: ierr
    first_time = .true.
    if ( allocated(zwork1) ) deallocate( zwork1 )
    if ( allocated(zwork0) ) deallocate( zwork0 )
    call mpi_comm_free(comm_fftx,ierr)
    call mpi_comm_free(comm_ffty,ierr)
    call mpi_comm_free(comm_fftz,ierr)
  END SUBROUTINE finalize_fock_fft

END MODULE fock_fft_module
