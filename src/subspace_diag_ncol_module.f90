MODULE subspace_diag_ncol_module

  use rgrid_module, only: dV,zdV
  use hamiltonian_module
  use hamiltonian_ncol_module
  use parallel_module, only: comm_grid, id_bzsm, myrank_k
  use watch_module

  implicit none

  PRIVATE
  PUBLIC :: subspace_diag_ncol

CONTAINS


  SUBROUTINE subspace_diag_ncol( k, n1,n2, psi, esp )

    implicit none
    integer,intent(IN) :: k,n1,n2
    complex(8),intent(INOUT) :: psi(:,:,:,:)
    real(8),intent(INOUT) :: esp(:,:,:)
    integer :: ML0,m,n,s,ierr,MB,MS,k0
    complex(8),allocatable :: work(:)
    integer :: WORK1,WORK2
    integer,save :: LWORK=0,LIWORK,LRWORK
    character(6) :: idiag0
    integer,allocatable :: iwork(:)
    real(8),allocatable :: rwork(:)
    complex(8),allocatable :: Hsub(:,:),Htmp(:,:)
    complex(8),allocatable :: psit(:,:,:)
    complex(8),allocatable :: zw1(:,:),zw2(:,:)
    complex(8),parameter :: zero=(0.0d0,0.0d0), one=(1.0d0,0.0d0)
    real(8) :: rtmp(1)
    integer :: itmp(1)
    integer,allocatable :: ir(:),id(:)
    real(8) :: ct(9),et(9)
    type(time) :: t
    include 'mpif.h'

    call write_border( 1, " subspace_diag_ncol(start)" )
    call start_timer( t )

    ML0 = size( psi, 1 )
    MB  = size( psi, 2 )
    MS  = size( psi, 4 )
    k0  = k-id_bzsm(myrank_k)

    ct(:) = 0.0d0
    et(:) = 0.0d0

    idiag0 = "zheevd"

    call watch(ct(1),et(1))

    allocate( Hsub(MB,MB) ) ; Hsub=zero

    allocate( psit(ML0,MB,MS) ) ; psit=zero

    call watch(ct(2),et(2))

    allocate( zw1(ML0,MS) ) ; zw1=zero
    allocate( zw2(ML0,MS) ) ; zw2=zero

    do s=1,MS
       do m=1,MB
#ifndef _DRSDFT_
          call hamiltonian(k,s,psi(:,m,k0,s),psit(:,m,s),n1,n2,m,m)
#endif
       end do
    end do
    do m=1,MB
       zw1(:,:) = psi(:,m,k0,:)
       zw2(:,:) = psit(:,m,:)
!       call hamiltonian_ncol( k, n1,n2, psi(:,m,k,:), psit(:,m,:) )
       call hamiltonian_ncol( k, n1,n2, zw1, zw2 )
       psit(:,m,:) = zw2(:,:)
    end do

    deallocate( zw2 )
    deallocate( zw1 )

    call watch(ct(3),et(3))

    Hsub(:,:)=zero
    do s=1,MS
       do n=1,MB
       do m=1,n
          Hsub(m,n) = Hsub(m,n) + sum( conjg(psi(:,m,k0,s))*psit(:,n,s) )*dV
       end do
       end do
    end do

    call watch(ct(4),et(4))

    deallocate( psit )

! --- allreduce ---

    allocate( Htmp(MB,MB) )
    Htmp=Hsub

    call watch(ct(5),et(5))

    call mpi_allreduce(Htmp,Hsub,MB*MB,MPI_COMPLEX16,MPI_SUM,comm_grid,ierr)

    call watch(ct(6),et(6))

    deallocate( Htmp )

! --- solve eigenvalue problem ---

    call watch(ct(7),et(7))

    WORK1 = 0
    WORK2 = 0

    select case(idiag0)
    case('zheev')

       LWORK=max(LWORK,2*MB-1)
       LRWORK=3*MB-2
       allocate( work(LWORK),rwork(LRWORK) )

       call zheev('V','U',MB,Hsub,MB,esp(1,k,s),work,2*MB,rwork,ierr)
       if ( ierr==0 ) WORK1=nint( real(work(1)) )

       deallocate( work,rwork )

    case('zheevd')

       LWORK=max(LWORK,2*MB+MB*MB)
       LRWORK=1+5*MB+2*MB*MB ; LIWORK=3+5*MB
       allocate( work(LWORK),rwork(LRWORK),iwork(LIWORK) )

       call zheevd('V','U',MB,Hsub,MB,esp(1,k,1) &
                  ,work,LWORK,rwork,LRWORK,iwork,LIWORK,ierr)
       if ( ierr==0 ) WORK1=nint( real(work(1)) )

       deallocate( work,rwork,iwork )

    end select

    esp(:,k,2) = esp(:,k,1)

    call watch(ct(8),et(8))

! --- Rotation ---

    allocate( psit(ML0,MB,1) ) ; psit=zero

    do s=1,MS
       call zgemm('N','N',ML0,MB,MB,one,psi(1,1,k0,s),ML0 &
            ,Hsub(1,1),MB,zero,psit,ML0)
       psi(:,:,k0,s)=psit(:,:,1)
    end do

    deallocate( psit )

    call watch(ct(9),et(9))

    deallocate( Hsub )

    LWORK=max(LWORK,WORK1)
    LIWORK=max(LIWORK,WORK2)

    call result_timer( t, "sd" )
    call write_border( 1, " subspace_diag_ncol(end)" )

    return

  END SUBROUTINE subspace_diag_ncol

END MODULE subspace_diag_ncol_module
