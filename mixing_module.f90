MODULE mixing_module

  use parallel_module

  implicit none

  PRIVATE
  PUBLIC :: imix, mmix, beta, sqerr_out &
       , init_mixing, read_mixing, perform_mixing

  integer :: imix,mmix
  real(8) :: beta,scf_conv,sqerr_out(2)
  real(8),allocatable :: Xold(:,:,:)
!  real(8),allocatable :: Xin(:,:,:),Xou(:,:,:)
  complex(8),allocatable :: Xin(:,:,:),Xou(:,:,:)
  real(8) :: beta0
  complex(8),parameter :: zero=(0.d0,0.d0)
  integer :: mmix_count=0
  real(8) :: diff(2),dif0(2),dif1(2)

CONTAINS

  SUBROUTINE read_mixing(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i
    character(7) :: cbuf,ckey
    imix=0
    mmix=4
    beta=1.d0
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
!  SUBROUTINE read_mixing(unit)
!    implicit none
!    integer,intent(IN) :: unit
!    read(unit,*) imix, mmix, beta, scf_conv
!    write(*,*) "imix, mmix =",imix,mmix
!    if ( mmix < 1 ) then
!       mmix=1
!       write(*,*) "mmix is replaced to 1 : mmix=",mmix
!    end if
!    write(*,*) "beta =",beta
!    write(*,*) "scf_conv=",scf_conv
!  END SUBROUTINE read_mixing


  SUBROUTINE send_mixing(rank)
    implicit none
    integer,intent(IN) :: rank
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(imix,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(mmix,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(beta,1,MPI_REAL8  ,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(scf_conv,1,MPI_REAL8,rank,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_mixing


  SUBROUTINE init_mixing(m,n,f)
    implicit none
    integer,intent(IN) :: m,n
    real(8),intent(IN) :: f(m,n)
    beta0      = 1.d0-beta
    mmix_count = 0
    if ( allocated(Xou)  ) deallocate( Xou )
    if ( allocated(Xin)  ) deallocate( Xin )
    if ( allocated(Xold) ) deallocate( Xold )
    allocate( Xold(m,n,mmix) ) ; Xold=0.d0
    allocate( Xin(m,n,mmix)  ) ; Xin=zero
    allocate( Xou(m,n,mmix)  ) ; Xou=zero
    Xold(:,:,1)=f(:,:)
    Xin(:,:,mmix)=f(:,:)
    dif0(:)=0.d0
  END SUBROUTINE init_mixing


  SUBROUTINE perform_mixing(m,n,f,flag_conv,disp_switch)
    implicit none
    integer,intent(IN)    :: m,n
    real(8),intent(INOUT) :: f(m,n)
    logical,intent(OUT)   :: flag_conv
    logical,optional,intent(IN) :: disp_switch
    real(8) :: err0(2),err(2),sum0(2),beta_bak,beta_min,dif_min(2)
    integer :: nn,ierr,loop,max_loop,mmix_count_bak,i
    complex(8),allocatable :: Xou_bak(:,:,:),Xin_bak(:,:,:)
    real(8),allocatable :: f_bak(:,:)

    allocate( Xin_bak(m,n,mmix) )
    allocate( Xou_bak(m,n,mmix) )
    allocate( f_bak(m,n) )

    Xin_bak(:,:,:) = Xin(:,:,:)
    Xou_bak(:,:,:) = Xou(:,:,:)
    f_bak(:,:) = f(:,:)
    beta_bak = beta
    mmix_count_bak=mmix_count

    beta_min     = beta
    dif_min(1:n) = 1.d10
    max_loop     = 1
    if ( beta >= 1.d0 ) max_loop = 11

    do loop=1,max_loop

    select case(imix)
    case default
       call simple_mixing(m,n,f,err0)
    case(10:19)
       call pulay_r_mixing(m,n,f,err0)
    case(20:29)
       write(*,*) "imix=",imix
       stop "this mixing is not available"
    end select

    sum0(1) = sum( (f(:,1)-Xold(:,1,1))**2 )/m
    sum0(n) = sum( (f(:,n)-Xold(:,n,1))**2 )/m
    call mpi_allreduce(sum0,dif1,n,mpi_real8,mpi_sum,comm_grid,ierr)
    if ( all(dif0(1:n)/=0.d0) ) then
       diff(1:n)=dif1(1:n)/dif0(1:n)
       if ( present(disp_switch) ) then
          if ( disp_switch ) then
             write(*,'(1x,"diff=",i4,2x,2(4g14.6,2x))') &
                  loop,(diff(i),dif1(i),dif0(i),beta,i=1,n)
             write(42,*) ( diff(i),i=1,n ),beta
          end if
       end if
    end if

    if ( all(dif1(1:n)<dif_min(1:n)) ) then
       beta_min = beta
       dif_min(1:n)=dif1(1:n)
    end if

    if ( loop < max_loop ) then
       if ( any(diff(1:n)>0.5d0) ) then
          Xin(:,:,:) = Xin_bak(:,:,:)
          Xou(:,:,:) = Xou_bak(:,:,:)
          f(:,:)     = f_bak(:,:)
          mmix_count = mmix_count_bak
          beta       = beta*0.5d0
          beta0      = 1.d0 - beta
       else
          exit
       end if
       if ( loop == max_loop-1 ) then
          beta  = beta_min
          beta0 = 1.d0-beta
       end if
    end if

    end do
    beta = beta_bak
    beta0 = 1.d0 - beta
    dif0(1:n)=dif1(1:n)

    deallocate( f_bak )
    deallocate( Xou_bak )
    deallocate( Xin_bak )

    call mpi_allgather(err0,n,mpi_real8,err,n,mpi_real8,comm_spin,ierr)
    flag_conv = .false.
    nn=n*np_spin
    if ( all(err(1:nn)<=scf_conv) ) flag_conv = .true.
    if ( present(disp_switch) ) then
       if ( disp_switch ) then
          write(*,*) "SQERR=",err(1:nn)
          write(40,*) err(1:nn)
       end if
    end if

    sqerr_out(:) = err(:)

    if ( disp_switch_parallel ) then
    if ( mod(imix,2)==0 .and. minval( f(:,:) ) < 0.d0 ) then
    write(*,*) "NEGATIVECHR",( sum( f(:,nn),mask=(f(:,nn)<0.d0) ) ,nn=1,n )
!    where( f < 0.d0 )
!       f=0.d0
!    end where
    f=abs(f)
    end if
    end if

    Xold(:,:,1) = f(:,:)

  END SUBROUTINE perform_mixing


  SUBROUTINE calc_sqerr(m,n,f,g,err)
    implicit none
    integer,intent(IN) :: m,n
    real(8),intent(IN) :: f(m,n),g(m,n)
    real(8),intent(OUT) :: err(n)
    real(8) :: err0(n)
    integer :: i,ierr
    do i=1,n
       err0(i)=sum( (f(:,i)-g(:,i))**2 )/m
    end do
    call mpi_allreduce(err0,err,n,MPI_REAL8,MPI_SUM,comm_grid,ierr)
  END SUBROUTINE calc_sqerr


  SUBROUTINE simple_mixing(m,n,f,err)
    implicit none
    integer,intent(IN) :: m,n
    real(8),intent(INOUT) :: f(m,n)
    real(8),intent(OUT) :: err(n)
    call calc_sqerr(m,n,f,Xold,err)
    f(:,:) = beta0*Xold(:,:,1) + beta*f(:,:)
  END SUBROUTINE simple_mixing


  SUBROUTINE pulay_r_mixing(m,n,f,err)
    implicit none
    integer,intent(IN) :: m,n
    real(8),intent(INOUT) :: f(m,n)
    real(8),intent(OUT) :: err(n)
    integer :: s,mmix0,ierr,i,i0,j0,mm,j
    integer,allocatable :: ipiv(:)
    real(8),allocatable :: rwork(:,:)
    complex(8) :: zc
    complex(8),allocatable :: A0(:,:),A1(:,:),b1(:),X(:),Y(:)

    Xou(:,:,mmix) = f(:,:)

    mmix_count = mmix_count + 1

    mmix0 = min( mmix_count,mmix )

    allocate( rwork(m,n) )
    rwork(:,:)=Xin(:,:,mmix)
    call calc_sqerr(m,n,f,rwork,err)
    deallocate( rwork )

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

       do j0=1,mmix0
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

    end do

    f(:,:) = real( Xin(:,:,mmix) )

    deallocate( A0,A1,b1,Y,X,ipiv )
    return

  END SUBROUTINE pulay_r_mixing

!--------1---------2---------3---------4---------5---------6---------7--
! Gauss-Jordan
! (from Numerical Recipes)

  SUBROUTINE gaussj(a,n,np,b,m,mp)
    implicit none
    integer,intent(IN) :: n,np,m,mp
    real(8),intent(INOUT) :: a(np,np),b(np,mp)
    integer,parameter :: NMAX=5000
    integer :: indxc(NMAX),indxr(NMAX),ipiv(NMAX)
    integer :: i,j,k,l,irow,icol,lll
    real(8) :: dum,pivinv,big
    do j=1,n
       ipiv(j)=0
    enddo
    do i=1,n
       big=0.d0
       do j=1,n
          if ( ipiv(j) /= 1 ) then
             do k=1,n
                if ( ipiv(k) == 0 ) then
                   if ( abs(a(j,k)) >= big ) then
                      big=abs(a(j,k))
                      irow=j
                      icol=k
                   endif
                else if ( ipiv(k) > 1 ) then
                   write(*,*) 'singular matrix in gaussj'
                   stop
                endif
             enddo
          endif
       enddo
       ipiv(icol)=ipiv(icol)+1
       if ( irow /= icol ) then
          do l=1,n
             dum=a(irow,l)
             a(irow,l)=a(icol,l)
             a(icol,l)=dum
          enddo
          do l=1,m
             dum=b(irow,l)
             b(irow,l)=b(icol,l)
             b(icol,l)=dum
          enddo
       endif
       indxr(i)=irow
       indxc(i)=icol
       if ( abs(a(icol,icol)) == 0.d0 ) then
          write(*,*) 'singular matrix in gaussj'
          stop
       endif
       pivinv=1.d0/a(icol,icol)
       a(icol,icol)=1.d0
       do l=1,n
          a(icol,l)=a(icol,l)*pivinv
       enddo
       do l=1,m
          b(icol,l)=b(icol,l)*pivinv
       enddo
       do lll=1,n
          if ( lll /= icol ) then
             dum=a(lll,icol)
             a(lll,icol)=0.d0
             do l=1,n
                a(lll,l)=a(lll,l)-a(icol,l)*dum
             enddo
             do l=1,m
                b(lll,l)=b(lll,l)-b(icol,l)*dum
             enddo
          endif
       enddo
    enddo
    do l=n,1,-1
       if ( indxr(l) /= indxc(l) ) then
          do k=1,n
             dum=a(k,indxr(l))
             a(k,indxr(l))=a(k,indxc(l))
             a(k,indxc(l))=dum
          enddo
       endif
    enddo
    return
  END SUBROUTINE gaussj

END MODULE mixing_module
