MODULE bz_module

  implicit none

  PRIVATE
  PUBLIC :: nk,mmm,Nbzsm,kbb,weight_bz,read_kgrid_bz,generate_bz &
           ,read_kgrid_oldformat_bz,generate_bz_sym, MMBZ

  integer :: nk,mmm(3,2)
  integer :: Nbzsm
  real(8),allocatable :: kbb(:,:)
  real(8),allocatable :: weight_bz(:)

  integer :: MMBZ,npbz
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
    npbz=0
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
          else if ( ckey(1:4) == "NPBZ" ) then
             backspace(unit)
             read(unit,*) cbuf,npbz
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

  SUBROUTINE read_kgrid_oldformat_bz(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    if ( rank == 0 ) then
       read(unit,*) nk !,ndata_read_k,kbb0(1:3)
       read(unit,*) mmm(1:3,1)
       read(unit,*) mmm(1:3,2)
       write(*,*) "nk =",nk
       write(*,*) "mmm1 =",mmm(:,1)
       write(*,*) "mmm2 =",mmm(:,2)
!       write(*,*) "ndata_read_k=",ndata_read_k
!       write(*,*) "kbb0=",kbb0(1:3)
    end if
    call send_kgrid_bz(0)
  END SUBROUTINE read_kgrid_oldformat_bz

  SUBROUTINE send_kgrid_bz(rank)
    implicit none
    integer,intent(IN) :: rank
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(nk,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(mmm,6,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(ndata_read_k,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kbb0,3,MPI_REAL8,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(npbz,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_kgrid_bz

  SUBROUTINE generate_bz(disp_switch)
    implicit none
    logical,intent(IN) :: disp_switch
    integer :: i,j,k,k1,iw,m1,m2,m3,mm1,mm2,mm3,i1,i2,i3,p1(3),p2(3)
    integer,allocatable :: mm(:,:),m(:,:),w(:)

    if ( DISP_SWITCH ) then
       write(*,'(a60," generate_bz_sym")') repeat("-",60)
    end if

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
    if ( npbz > k ) Nbzsm=npbz

    MMBZ  = k1

    allocate( weight_bz(Nbzsm) ) ; weight_bz=0.d0
    allocate( kbb(3,Nbzsm)     ) ; kbb=0.d0
    allocate( wbz(Nbzsm)       ) ; wbz=0.d0

    kbb(1:3,1:k)=real(m(1:3,1:k),8)/nk
    weight_bz(1:k)=real(w(1:k),8)/MMBZ
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


  SUBROUTINE generate_bz_sym( nsym,rgb,disp_switch )
    implicit none
    integer,intent(IN) :: nsym
    real(8),intent(IN) :: rgb(3,3,nsym)
    logical,intent(IN) :: disp_switch
    integer,allocatable :: mm(:,:),m(:,:),w(:),w1(:),w2(:)
    integer :: m1,m2,m3,mm1,mm2,mm3,i1,i2,i3,p1(3),p2(3),MBZ_tmp
    integer :: i,j,k,k1,iw,n,nkmax,p3(3),ns,iqrt,ni,is,nni,ig,iqrt2
    real(8) :: c,tmp(3)

    if ( disp_switch ) then
       write(*,'(a60," generate_bz_sym")') repeat("-",60)
    end if

    m1 =mmm(1,1) ; m2 =mmm(2,1) ; m3 =mmm(3,1)
    mm1=mmm(1,2) ; mm2=mmm(2,2) ; mm3=mmm(3,2)

    nkmax=(2*m1+1)*(2*m2+1)*(2*m3+1)*2
    allocate( mm(3,nkmax),m(3,nkmax),w(nkmax) ) ; mm=0 ; m=0 ; w=0
    allocate( w1(nkmax),w2(nkmax) ) ; w1=0 ; w2=0

    ni=1
    is=1

    do i1=-m1,m1,mm1
    do i2=-m2,m2,mm2
       loop_3 : do i3=-m3,m3,mm3

          p1(1)=i1 ; p1(2)=i2 ; p1(3)=i3 ; p2(1:3)=p1(1:3)

          do i=1,3
             p1(i)=mod(p2(i),nk)
             if ( p1(i) >   nk/2 ) p1(i)=p1(i)-nk
             if ( p1(i) <= -nk/2 ) p1(i)=p1(i)+nk
          end do

          do i=1,ni-1
             if ( all(p1(1:3)==mm(1:3,i)) ) cycle loop_3
          end do

          ns =0
          nni=ni

          do iw=1,-1,-2
!          do iw=1,1
             loop_sym : do ig=1,nsym

                if ( ni>nkmax ) stop "generate_bz_sym(1)"

                tmp(:) = matmul( rgb(:,:,ig),p1(:) )*iw
                p3(:) = nint( tmp(:) )

                do i=1,3
                   p2(i)=mod(p3(i),nk)
                   if ( p2(i) >  nk/2 ) p2(i)=p2(i)-nk
                   if ( p2(i) <=-nk/2 ) p2(i)=p2(i)+nk
                end do

                do i=nni,ni-1
                   if ( all(p2(:)==mm(:,i)) ) cycle loop_sym
                end do
                ns=ns+1
                mm(:,ni)=p2(:)
                ni=ni+1

             end do loop_sym
          end do ! iw

          w(is)=ns

          m(1:3,is)=mm(1:3,nni)

          is=is+1

       end do loop_3
    end do ! i2
    end do ! i1

    is=is-1
    ni=ni-1

    do k=1,ni
       w1(k)=1
       do i=1,3
          i1=mod(i,3)+1
          i2=mod(i+1,3)+1
          if ( abs(mm(i,k)) == nk/2 ) then
             do k1=1,ni
                if ( k == k1 ) cycle
                if ( mm(i ,k1) ==-mm(i ,k) .and. &
                     mm(i1,k1) == mm(i1,k) .and. &
                     mm(i2,k1) == mm(i2,k) ) w1(k)=w1(k)*2
             end do
          end if
       end do ! i
    end do ! k

    do k=1,is
       do k1=1,ni
          if ( all(m(1:3,k)==mm(1:3,k1)) ) then
             w2(k)=w1(k1)
             exit
          end if
       end do
    end do

    Nbzsm = is
    if ( npbz > is ) Nbzsm = npbz

    MMBZ = ni

    allocate( weight_bz(Nbzsm) ) ; weight_bz=0.0d0
    allocate( kbb(3,Nbzsm)     ) ; kbb=0.0d0
    allocate( wbz(Nbzsm)       ) ; wbz=0.0d0

    do k=1,Nbzsm
       kbb(1:3,k)   = real( m(1:3,k), 8 )/real( nk, 8 )
       weight_bz(k) = real( w(k), 8 )/real( w2(k), 8 )
    end do

    c=sum( weight_bz(:) )
    weight_bz(:)=weight_bz(:)/c

    wbz(1:Nbzsm) = weight_bz(1:Nbzsm)

    if ( DISP_SWITCH ) then
       write(*,*) "Nbzsm, MMBZ =",Nbzsm,MMBZ
       write(*,*) "KBB"
       do k=1,Nbzsm
          write(*,'(1x,i4,3f10.5,f12.5,2i5)') k,kbb(:,k),wbz(k),w(k),w2(k)
       end do
       write(*,*) "sum(w)  =",sum(w),sum(w2),sum(w1)
       write(*,*) "sum(wbz)=",sum(wbz)
    end if

    deallocate( w2,w1,w,m,mm )

    if ( DISP_SWITCH ) then
       write(*,'(a60," genetrate_bz_sym(END)")') repeat("-",60)
    end if

    return

  END SUBROUTINE generate_bz_sym


END MODULE bz_module
