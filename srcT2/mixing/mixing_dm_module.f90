MODULE mixing_dm_module

  use parallel_module, only: RSDFT_MPI_COMPLEX16
  use io_tools_module

  implicit none

  PRIVATE
  PUBLIC :: init_mixing_dm
  PUBLIC :: pulay_mixing_dm

  integer :: comm
  integer :: mmix_count = 0
  integer :: mmix       = 4
  real(8) :: beta       = 0.1d0
  logical :: flag_init  =.false.

  complex(8),allocatable :: Xin(:,:,:)
  complex(8),allocatable :: Xou(:,:,:)

  logical :: disp_switch

CONTAINS


  SUBROUTINE init_mixing_dm( m, n, comm_in, f, mmix_in, beta_in )

    implicit none
    integer,intent(IN) :: m,n,comm_in
    complex(8),optional,intent(IN) :: f(m,n)
    integer,optional,intent(IN) :: mmix_in
    real(8),optional,intent(IN) :: beta_in

    if ( flag_init ) return
    flag_init=.true.

    call write_border( 0, "init_mixing_dm(start)" )

    comm = comm_in

    if ( present(beta_in) ) then
       beta = beta_in
    else
       call IOTools_readReal8Keyword( "BETA", beta )
    end if
    if ( present(mmix_in) ) then
       mmix = mmix_in
    else
       call IOTools_readIntegerKeyword( "MMIX", mmix )
    end if

    call check_disp_switch( disp_switch, 0 )
    if ( disp_switch ) then
       write(*,*) "Max. # of successive data for mixing: mmix=",mmix
       write(*,*) "Ratio of the mixing of the latest data: beta=",beta
       write(*,*) "Only density mixing is available"
    end if

    allocate( Xin(m,n,mmix) ) ; Xin=(0.0d0,0.0d0)
    allocate( Xou(m,n,mmix) ) ; Xou=(0.0d0,0.0d0)

    if ( present(f) ) Xin(:,:,mmix)=f(:,:)

    call write_border( 0, "init_mixing_dm(end)" )

  END SUBROUTINE init_mixing_dm


  SUBROUTINE pulay_mixing_dm( m,n,f )

    implicit none
    integer,intent(IN)       :: m,n
    complex(8),intent(INOUT) :: f(m,n)
    integer :: s,mmix0,ierr,i,i0,j0,mm,j
    integer,allocatable :: ipiv(:)
    real(8),allocatable :: rwork(:,:)
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    complex(8) :: zc
    complex(8),allocatable :: A0(:,:),A1(:,:),b1(:),X(:),Y(:)
    include 'mpif.h'

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
       if ( disp_switch ) then
          write(*,*) "simple mixing: mmix0/mmix=",mmix0,mmix
          do i=1,mmix
             write(*,*) i,sum(abs(Xin(:,:,i))),sum(abs(Xou(:,:,i)))
          end do
       end if
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

    call mpi_allreduce(A0,A1,mm,RSDFT_MPI_COMPLEX16,mpi_sum,comm,ierr)

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

    f(:,:) = Xin(:,:,mmix)

    deallocate( A0,A1,b1,Y,X,ipiv )
    return

  END SUBROUTINE pulay_mixing_dm


END MODULE mixing_dm_module
