MODULE bz_module

  implicit none

  PRIVATE
  PUBLIC :: nk,mmm,Nbzsm,kbb,weight_bz,read_kgrid_bz,generate_bz

  integer :: nk,mmm(3,2)
  integer :: Nbzsm
  real(8),allocatable :: kbb(:,:)
  real(8),allocatable :: weight_bz(:)

  integer :: MMBZ
  real(8),allocatable :: wbz(:)
  integer :: ndata_read_k=0
  real(8) :: kbb0(3)
  data kbb0/0.d0,0.d0,0.d0/

CONTAINS

  SUBROUTINE read_kgrid_bz(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i
    character(6) :: cbuf,ckey
    nk=2
    mmm(1:3,1)=(/ 2,2,2 /)
    mmm(1:3,2)=(/ 2,2,2 /)
    ndata_read_k=0
    kbb0(1:3)=0.d0
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:2) == "NK" ) then
             backspace(unit)
             read(unit,*) cbuf,nk
          else if ( ckey(1:4) == "MMM1" ) then
             backspace(unit)
             read(unit,*) cbuf,mmm(1:3,1)
          else if ( ckey(1:4) == "MMM2" ) then
             backspace(unit)
             read(unit,*) cbuf,mmm(1:3,2)
          end if
       end do
999    continue
       write(*,*) "nk =",nk
       write(*,*) "mmm1 =",mmm(:,1)
       write(*,*) "mmm2 =",mmm(:,2)
!       write(*,*) "ndata_read_k=",ndata_read_k
!       write(*,*) "kbb0=",kbb0(1:3)
    end if
    call send_kgrid_bz(0)
  END SUBROUTINE read_kgrid_bz

!  SUBROUTINE read_kgrid_bz(unit)
!    integer,intent(IN) :: unit
!    read(unit,*) nk,ndata_read_k,kbb0(1:3)
!    read(unit,*) mmm(1:3,1)
!    read(unit,*) mmm(1:3,2)
!    write(*,*) "nk =",nk
!    write(*,*) "mmm1 =",mmm(:,1)
!    write(*,*) "mmm2 =",mmm(:,2)
!    write(*,*) "ndata_read_k=",ndata_read_k
!    write(*,*) "kbb0=",kbb0(1:3)
!  END SUBROUTINE read_kgrid_bz

  SUBROUTINE send_kgrid_bz(rank)
    implicit none
    integer,intent(IN) :: rank
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(nk,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(mmm,6,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(ndata_read_k,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kbb0,3,MPI_REAL8,rank,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_kgrid_bz

  SUBROUTINE generate_bz(disp_switch)
    logical,intent(IN) :: disp_switch
    integer :: i,j,k,k1,iw,m1,m2,m3,mm1,mm2,mm3,i1,i2,i3,p1(3),p2(3)
    integer,allocatable :: mm(:,:),m(:,:),w(:)

    m1 =mmm(1,1) ; m2 =mmm(2,1) ; m3 =mmm(3,1)
    mm1=mmm(1,2) ; mm2=mmm(2,2) ; mm3=mmm(3,2)

    k=(2*m1+1)*(2*m2+1)*(2*m3+1)*2
    allocate( mm(3,k),m(3,k),w(k) )
    mm=0 ; m=0 ; w=0
    k=0 ; k1=0

    do i1=-m1,m1,mm1
    do i2=-m2,m2,mm2
    loop_A : do i3=-m3,m3,mm3

       do iw=1,-1,-2

          p1(1)=i1*iw ; p1(2)=i2*iw ; p1(3)=i3*iw ; p2(1:3)=p1(1:3)
          do i=1,3
             p1(i)=mod(p2(i),nk)
             if ( p1(i)>  nk/2 ) p1(i)=p1(i)-nk
             if ( p1(i)<=-nk/2 ) p1(i)=p1(i)+nk
          end do
          if ( k1>0 ) then
             do i=1,k1
                if ( mm(1,i)==p1(1) .and. &
                     mm(2,i)==p1(2) .and. &
                     mm(3,i)==p1(3)         ) cycle loop_A
             end do
          end if

          if ( iw==1 ) then
             k=k+1
             m(1:3,k)=p1(1:3)
             w(k)=1
          else
             w(k)=2
          end if

          k1=k1+1
          mm(1:3,k1)=p1(1:3)

       end do ! iw

    end do loop_A
    end do
    end do

    Nbzsm = k
    MMBZ  = k1

    allocate( weight_bz(Nbzsm) )
    allocate( kbb(3,Nbzsm) )

    allocate( wbz(Nbzsm) )

    kbb(1:3,1:Nbzsm)=real(m(1:3,1:Nbzsm),8)/nk
    weight_bz(1:Nbzsm)=real(w(1:Nbzsm),8)/MMBZ
    wbz(1:Nbzsm) = weight_bz(1:Nbzsm)

    if ( Nbzsm == ndata_read_k ) then
       kbb(1:3,1) = kbb0(1:3)
       ndata_read_k = 0
    end if

    deallocate( w,m,mm )

    if ( disp_switch ) then
       write(*,*) "Nbzsm, MMBZ =",Nbzsm,MMBZ
       write(*,'(1x,a4,a30,a12)') "","kbb","weight_bz"
       do k=1,Nbzsm
          write(*,'(1x,i4,3f10.5,f12.5)') k,kbb(:,k),weight_bz(k)
       end do
    end if

  END SUBROUTINE generate_bz

END MODULE bz_module
