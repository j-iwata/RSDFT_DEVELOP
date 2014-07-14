MODULE mixing_module

  use mixing_broyden_module

  implicit none

  PRIVATE
  PUBLIC :: sqerr_out, imix, beta &
           ,init_mixing,read_mixing,perform_mixing,read_oldformat_mixing

  include 'mpif.h'

  integer :: imix,mmix
  real(8) :: beta,scf_conv,sqerr_out(4)
  real(8),allocatable :: Xold(:,:,:)
  complex(8),allocatable :: Xin(:,:,:),Xou(:,:,:)
  real(8) :: beta0
  complex(8),parameter :: zero=(0.d0,0.d0)
  integer :: mmix_count=0
  real(8) :: diff(2),dif0(2),dif1(2)

  integer :: ML, ML0, MF0, MSP
  integer :: comm_grid, comm_spin
  real(8) :: dV
  logical :: disp_switch

CONTAINS


  SUBROUTINE read_mixing(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i
    character(7) :: cbuf,ckey
    imix=0
    mmix=4
    beta=1.0d0
    scf_conv=1.d-15
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:4) == "IMIX" ) then
             backspace(unit)
             read(unit,*) cbuf,imix
          else if ( ckey(1:4) == "MMIX" ) then
             backspace(unit)
             read(unit,*) cbuf,mmix
          else if ( ckey(1:4) == "BETA" ) then
             backspace(unit)
             read(unit,*) cbuf,beta
          else if ( ckey(1:7) == "SCFCONV" ) then
             backspace(unit)
             read(unit,*) cbuf,scf_conv
          end if
       end do
999    continue
       write(*,*) "imix    =",imix
       write(*,*) "mmix    =",mmix
       if ( mmix < 1 ) then
          mmix=1
          write(*,*) "mmix is replaced to 1 : mmix=",mmix
       end if
       write(*,*) "beta    =",beta
       write(*,*) "scf_conv=",scf_conv
    end if
    call send_mixing(0)
  END SUBROUTINE read_mixing


  SUBROUTINE read_oldformat_mixing(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    if ( rank == 0 ) then
       read(unit,*) imix, mmix, beta, scf_conv
       write(*,*) "imix, mmix =",imix,mmix
       if ( mmix < 1 ) then
          mmix=1
          write(*,*) "mmix is replaced to 1 : mmix=",mmix
       end if
       write(*,*) "beta =",beta
       write(*,*) "scf_conv=",scf_conv
    end if
    call send_mixing(0)
  END SUBROUTINE read_oldformat_mixing


  SUBROUTINE send_mixing(rank)
    implicit none
    integer,intent(IN) :: rank
    integer :: ierr
    call mpi_bcast(imix,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(mmix,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(beta,1,MPI_REAL8  ,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(scf_conv,1,MPI_REAL8,rank,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_mixing


  SUBROUTINE init_mixing(ML0_in,MSP_in,nf1,nf2,comm_grid_in,comm_spin_in,dV_in,f,g)
    implicit none
    integer,intent(IN) :: ML0_in,MSP_in,nf1,nf2,comm_grid_in,comm_spin_in
    real(8),intent(IN) :: dV_in,f(ML0_in,nf1:nf2), g(ML0_in,nf1:nf2)
    integer :: m,ierr

    beta0      = 1.0d0-beta
    mmix_count = 0

    ML0 = ML0_in
    MSP = MSP_in
    MF0 = MSP*2
    call mpi_allreduce( ML0, ML, 1,MPI_INTEGER,MPI_SUM,comm_grid,ierr )

    comm_grid = comm_grid_in
    comm_spin = comm_spin_in

    dV = dV_in

    if ( allocated(Xou)  ) deallocate( Xou )
    if ( allocated(Xin)  ) deallocate( Xin )
    if ( allocated(Xold) ) deallocate( Xold )
    allocate( Xold(ML0,MF0,mmix) ) ; Xold=0.0d0
    allocate(  Xin(ML0,MSP,mmix) ) ;  Xin=zero
    allocate(  Xou(ML0,MSP,mmix) ) ;  Xou=zero

    m=ML0*(nf2-nf1+1)
    call mpi_allgather( f(1,nf1),m,MPI_REAL8 &
         ,Xold(1,1,1),m,MPI_REAL8,comm_spin,ierr )
    call mpi_allgather( g(1,nf1),m,MPI_REAL8 &
         ,Xold(1,MSP+1,1),m,MPI_REAL8,comm_spin,ierr )

    if ( mod(imix,2) == 0 ) then
       Xin(:,:,mmix) = Xold(:,1:MSP,1)
    else if ( mod(imix,2) == 1 ) then
       Xin(:,:,mmix) = Xold(:,MSP+1:2*MSP,1)
    end if

    dif0(:) = 0.0d0

  END SUBROUTINE init_mixing


  SUBROUTINE perform_mixing( m, n1, n2, f_io, g_io, flag_conv, disp_sw_in )
    implicit none
    integer,intent(IN)    :: m,n1,n2
    real(8),intent(INOUT) :: f_io(m,n1:n2),g_io(m,n1:n2)
    logical,intent(OUT)   :: flag_conv
    logical,optional,intent(IN) :: disp_sw_in
    real(8) :: err0(2),err(2),sum0(2),beta_bak,beta_min,dif_min(2)
    integer :: n,ierr,loop,max_loop,mmix_count_bak,i,i0
    complex(8),allocatable :: Xou_bak(:,:,:),Xin_bak(:,:,:)
    real(8),allocatable :: f(:,:),g(:,:),h(:,:),h_old(:,:),h_bak(:,:)

    disp_switch = .false.
    if ( present(disp_sw_in) ) disp_switch = disp_sw_in

    if ( disp_switch ) write(*,'("----- mixing")')

    allocate( f(ML0,MSP) ) ; f=0.0d0
    allocate( g(ML0,MSP) ) ; g=0.0d0

    n=m*(n2-n1+1)
    call mpi_allgather( f_io(1,n1),n,MPI_REAL8,f,n,MPI_REAL8,comm_spin,ierr)
    call mpi_allgather( g_io(1,n1),n,MPI_REAL8,g,n,MPI_REAL8,comm_spin,ierr)

    call calc_sqerr( ML0, MSP, f, g, flag_conv )
    if ( flag_conv ) then
       deallocate( g,f )
       return
    end if

    allocate( h(ML0,MSP)     ) ; h=0.0d0
    allocate( h_old(ML0,MSP) ) ; h_old=0.0d0

    if ( mod(imix,2) == 0 ) then
       i0=0
       h_old(:,:) = Xold(:,i0+1:i0+MSP,1)
       h(:,:) = f(:,:)
    else if ( mod(imix,2) == 1 ) then
       i0=MSP
       h_old(:,:) = Xold(:,i0+1:i0+MSP,1)
       h(:,:) = g(:,:)
    end if

    allocate( Xin_bak(ML0,MSP,mmix) ) ; Xin_bak=zero
    allocate( Xou_bak(ML0,MSP,mmix) ) ; Xou_bak=zero
    allocate( h_bak(ML0,MSP)        ) ; h_bak=0.0d0

    Xin_bak(:,:,:) = Xin(:,:,:)
    Xou_bak(:,:,:) = Xou(:,:,:)
    h_bak(:,:)     = h(:,:)
    beta_bak       = beta
    mmix_count_bak = mmix_count
    beta_min       = beta
    dif_min(1:MSP) = 1.d10

    max_loop = 1
    if ( beta >= 1.0d0 ) max_loop = 11

    do loop=1,max_loop

       select case(imix)
       case default
          call simple_mixing( ML0, MSP, h_old, h )
       case(10:19)
          call pulay_r_mixing( ML0, MSP, h )
       case(20:29)
          call pulay_r2_mixing( ML0, MSP, h )
!          write(*,*) "imix=",imix
!          stop "this mixing is not available"
       case(30:39)
          call broyden_mixing &
               ( ML0, MSP, comm_grid, mmix_count, mmix, beta, h, Xin, Xou )
       end select

       sum0(:)=0.0d0
       do i=1,MSP
          sum0(i) = sum( (h(:,i)-Xold(:,i0+i,1))**2 )/ML
       end do
       call mpi_allreduce(sum0,dif1,MSP,MPI_REAL8,MPI_SUM,comm_grid,ierr)

       diff(:)=0.0d0
       do i=1,MSP
          if ( dif0(i) /= 0.0d0 ) diff(i) = dif1(i)/dif0(i)
       end do

       if ( present(disp_sw_in) ) then
          if ( disp_switch ) then
             write(*,'(1x,"diff=",i4,2x,2(4g14.6,2x))') &
                  loop,(diff(i),dif1(i),dif0(i),beta,i=1,MSP)
          end if
       end if

       if ( all( dif1(1:MSP) < dif_min(1:MSP) ) ) then
          beta_min = beta
          dif_min(1:MSP) = dif1(1:MSP)
       end if

       if ( loop < max_loop ) then

          if ( any( diff(1:MSP) > 0.5d0 ) ) then

             Xin(:,:,:) = Xin_bak(:,:,:)
             Xou(:,:,:) = Xou_bak(:,:,:)
             h(:,:)     = h_bak(:,:)
             mmix_count = mmix_count_bak
             beta       = beta*0.5d0
             beta0      = 1.0d0 - beta

          else

             exit

          end if

          if ( loop == max_loop-1 ) then
             beta  = beta_min
             beta0 = 1.0d0-beta
          end if

       end if

    end do ! loop

    beta        = beta_bak
    beta0       = 1.0d0 - beta
    dif0(1:MSP) = dif1(1:MSP)

    if ( mod(imix,2) == 0 ) then
       call chk_chrg( ML0,MSP,h )
       f(:,:) = h(:,:)
    else if ( mod(imix,2) == 1 ) then
       g(:,:) = h(:,:)
    end if

    f_io(:,n1:n2) = f(:,n1:n2)
    g_io(:,n1:n2) = g(:,n1:n2)

    Xold(:,1:MSP,1)       = f(:,:)
    Xold(:,MSP+1:2*MSP,1) = g(:,:)

    deallocate( h_bak )
    deallocate( Xou_bak )
    deallocate( Xin_bak )
    deallocate( h_old )
    deallocate( h,g,f )

  END SUBROUTINE perform_mixing


  SUBROUTINE calc_sqerr(m,n,f,g,flag_conv)
    implicit none
    integer,intent(IN)  :: m,n
    real(8),intent(IN)  :: f(m,n),g(m,n)
    logical,intent(OUT) :: flag_conv 
    real(8) :: err0(2*n),err(2*n)
    integer :: i,ierr

    do i=1,n
       err0(i  )=sum( ( f(:,i)-Xold(:,i  ,1) )**2 )/ML
       err0(i+n)=sum( ( g(:,i)-Xold(:,i+n,1) )**2 )/ML
    end do
    call mpi_allreduce(err0,err,2*n,MPI_REAL8,MPI_SUM,comm_grid,ierr)

    if ( all( err(n+1:2*n) <= scf_conv ) ) then
       flag_conv = .true.
    else
       flag_conv = .false.
    end if

    if ( disp_switch ) then
       if ( MSP == 1 ) then
          write(*,'(1x,"RSQERR=",g12.5,3x,"VSQERR=",g12.5)') err(1:2)
          write(40,*) err(1:2)
       else if ( MSP == 2 ) then
          write(*,'(1x,"RSQERR=",2g12.5,3x,"VSQERR=",2g12.5)') err(1:4)
          write(40,'(1x,4g12.5)') err(1:4)
       else
          write(*,*) "MSP is invalid: MSP=",MSP
          stop "stop@parform_mixing"
       end if
    end if

    sqerr_out(1:n)     = err(1:n)
    sqerr_out(n+1:2*n) = err(n+1:2*n)

  END SUBROUTINE calc_sqerr


  SUBROUTINE chk_chrg( m,n,f )
    implicit none
    integer,intent(IN)    :: m,n
    real(8),intent(INOUT) :: f(m,n)
    real(8) :: chrg(2,n)
    integer :: i,ierr

    chrg=0.0d0
    do i=1,n
       chrg(1,i) = sum( f(:,i) )*dV
       chrg(2,i) = sum( f(:,i),mask=(f(:,i)<0.0d0) )*dV
    end do
    call mpi_allreduce(MPI_IN_PLACE,chrg,2*n,MPI_REAL8,MPI_SUM,comm_grid,ierr)

    if ( any( chrg(2,1:n) < 0.0d0 ) ) then
       if ( disp_switch ) then
          write(*,'(1x,"NEGATIVE density exist (tot,neg)=",4g18.9)') &
               ( chrg(1,i),chrg(2,i), i=1,n )
       end if
    end if

!    where( f < 0.0d0 )
!       f = 0.0d0
!    end where
    f(:,:) = abs( f(:,:) )

  END SUBROUTINE chk_chrg


  SUBROUTINE simple_mixing(m,n,fold,f)
    implicit none
    integer,intent(IN)    :: m,n
    real(8),intent(IN)    :: fold(m,n)
    real(8),intent(INOUT) :: f(m,n)
    f(:,:) = beta0*fold(:,:) + beta*f(:,:)
  END SUBROUTINE simple_mixing


  SUBROUTINE pulay_r_mixing(m,n,f)
    implicit none
    integer,intent(IN)    :: m,n
    real(8),intent(INOUT) :: f(m,n)
    integer :: s,mmix0,ierr,i,i0,j0,mm,j
    integer,allocatable :: ipiv(:)
    real(8),allocatable :: rwork(:,:)
    complex(8) :: zc
    complex(8),allocatable :: A0(:,:),A1(:,:),b1(:),X(:),Y(:)

    Xou(:,:,mmix) = f(:,:)

    mmix_count = mmix_count + 1

    mmix0 = min( mmix_count, mmix )

    if ( mmix == 1 .or. mmix0 < mmix ) then
       do i=2,mmix
          Xin(:,:,i-1)=Xin(:,:,i)
          Xou(:,:,i-1)=Xou(:,:,i)
       end do
       Xin(:,:,mmix) = Xin(:,:,mmix) + beta*( Xou(:,:,mmix)-Xin(:,:,mmix) )
       f(:,:) = Xin(:,:,mmix)
       return
    end if

    allocate( ipiv(mmix0) )
    allocate( X(m),Y(m) )
    allocate( b1(mmix0) )
    allocate( A1(mmix0,mmix0) )
    allocate( A0(mmix0,mmix0) )

    b1(:)   = zero
    A1(:,:) = zero
    A0(:,:) = zero
    mm      = size(A1)

    do s=1,n

       do j0=1 ,mmix0
       do i0=j0,mmix0
          i=mmix-mmix0+i0
          j=mmix-mmix0+j0
          A0(i0,j0)=sum( conjg(Xou(:,s,i)-Xin(:,s,i)) &
                             *(Xou(:,s,j)-Xin(:,s,j)) )
          A0(j0,i0)=conjg( A0(i0,j0) )
       end do
       end do

       call mpi_allreduce(A0,A1,mm,mpi_complex16,mpi_sum,comm_grid,ierr)

       b1(1:mmix0) = (1.d0,0.d0)
       A0(:,:)     = A1(:,:)

       call zgesv(mmix0,1,A1,mmix0,ipiv,b1,mmix0,ierr)

       zc=1.d0/sum( b1(1:mmix0) )
       b1(1:mmix0)=zc*b1(1:mmix0)

       X(:)=zero
       Y(:)=zero

       do i0=1,mmix0
          i=mmix-mmix0+i0
          X(:)=X(:)+b1(i0)*Xin(:,s,i)
          Y(:)=Y(:)+b1(i0)*Xou(:,s,i)
       end do

       do i=max(1,mmix-mmix0),mmix-1
          Xin(:,s,i)=Xin(:,s,i+1)
          Xou(:,s,i)=Xou(:,s,i+1)
       end do

       Xin(:,s,mmix) = X(:) + beta*( Y(:)-X(:) )

    end do ! s

    f(:,:) = real( Xin(:,:,mmix) )

    deallocate( A0,A1,b1,Y,X,ipiv )
    return

  END SUBROUTINE pulay_r_mixing


  SUBROUTINE pulay_r2_mixing(m,n,f)
    implicit none
    integer,intent(IN)    :: m,n
    real(8),intent(INOUT) :: f(m,n)
    integer :: s,mmix0,ierr,i,i0,j0,mm,j
    integer,allocatable :: ipiv(:)
    real(8),allocatable :: rwork(:,:)
    complex(8) :: zc
    complex(8),allocatable :: A0(:,:),A1(:,:),b1(:),X(:),Y(:)

    Xou(:,:,mmix) = f(:,:)

    mmix_count = mmix_count + 1

    mmix0 = min( mmix_count, mmix )

    if ( mmix == 1 .or. mmix0 < mmix ) then
       do i=2,mmix
          Xin(:,:,i-1)=Xin(:,:,i)
          Xou(:,:,i-1)=Xou(:,:,i)
       end do
       Xin(:,:,mmix) = Xin(:,:,mmix) + beta*( Xou(:,:,mmix)-Xin(:,:,mmix) )
       f(:,:) = Xin(:,:,mmix)
       return
    end if

    allocate( ipiv(mmix0) )
    allocate( X(m),Y(m) )
    allocate( b1(mmix0) )
    allocate( A1(mmix0,mmix0) )
    allocate( A0(mmix0,mmix0) )

    b1(:)   = zero
    A1(:,:) = zero
    A0(:,:) = zero
    mm      = size(A1)

    do j0=1 ,mmix0
    do i0=j0,mmix0
       i=mmix-mmix0+i0
       j=mmix-mmix0+j0
       A0(i0,j0)=sum( conjg(Xou(:,:,i)-Xin(:,:,i)) &
                          *(Xou(:,:,j)-Xin(:,:,j)) )
       A0(j0,i0)=conjg( A0(i0,j0) )
    end do
    end do

    call mpi_allreduce(A0,A1,mm,mpi_complex16,mpi_sum,comm_grid,ierr)

    b1(1:mmix0) = (1.d0,0.d0)
    A0(:,:)     = A1(:,:)

    call zgesv(mmix0,1,A1,mmix0,ipiv,b1,mmix0,ierr)

    zc=1.d0/sum( b1(1:mmix0) )
    b1(1:mmix0)=zc*b1(1:mmix0)

    do s=1,n

       X(:)=zero
       Y(:)=zero

       do i0=1,mmix0
          i=mmix-mmix0+i0
          X(:)=X(:)+b1(i0)*Xin(:,s,i)
          Y(:)=Y(:)+b1(i0)*Xou(:,s,i)
       end do

       do i=max(1,mmix-mmix0),mmix-1
          Xin(:,s,i)=Xin(:,s,i+1)
          Xou(:,s,i)=Xou(:,s,i+1)
       end do

       Xin(:,s,mmix) = X(:) + beta*( Y(:)-X(:) )

    end do ! s

    f(:,:) = real( Xin(:,:,mmix) )

    deallocate( A0,A1,b1,Y,X,ipiv )
    return

  END SUBROUTINE pulay_r2_mixing


END MODULE mixing_module
