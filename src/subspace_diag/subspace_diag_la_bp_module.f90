module subspace_diag_la_bp_module

  use wf_module, only: unk,esp
  use hamiltonian_module
  use watch_module
  use vector_tools_module, only: vinfo
  use rsdft_mpi_module

  implicit none

  private
  public :: subspace_diag_la_bp

#ifdef _DRSDFT_
  character(6) :: idiag0 = "dsyev"  ! "dsyevd"
! DSYEVD of Intel MKL 11.2 Update 1 or earlier may have failure
! (see https://software.intel.com/en-us/articles/intel-mkl-112-bug-fixes)
#else
  character(6) :: idiag0 = "zheevd"
#endif

contains

  subroutine subspace_diag_la_bp(k,s,v)
    implicit none
    integer,intent(in) :: k,s
    type(vinfo),intent(in) :: v(2)
    integer :: ML0,n1,n2,m1,m,n,ierr,MB,MB0
    integer :: np_band, comm_band, comm_grid
    complex(8),allocatable :: work(:)
    integer :: WORK1,WORK2
    integer,save :: LWORK=0,LIWORK,LRWORK
    integer,allocatable :: iwork(:)
    real(8),allocatable :: rwork(:)
#ifdef _DRSDFT_
    real(8),allocatable :: psit(:,:),psiu(:,:),Hsub(:,:)
    real(8) :: zz
    real(8),parameter :: zero=0.0d0, one=1.0d0
#else
    complex(8),allocatable :: psit(:,:),psiu(:,:),Hsub(:,:)
    complex(8) :: zz
    complex(8),parameter :: zero=(0.0d0,0.0d0), one=(1.0d0,0.0d0)
#endif
    real(8) :: ct(9),et(9)
    real(8) :: rtmp(1),dV
    integer :: itmp(1)
    integer,allocatable :: ir(:),id(:)
    type(time) :: t

    call write_border( 1, " subspace_diag_la_bp(start)" )
    call start_timer( t )

    ML0 = v(1)%pinfo%ir( v(1)%pinfo%me )
    n1  = v(1)%pinfo%id( v(1)%pinfo%me ) + 1
    n2  = n1 + ML0 - 1
    m1  = v(2)%pinfo%id( v(2)%pinfo%me ) + 1
    MB0 = v(2)%pinfo%ir( v(2)%pinfo%me )
    MB  = sum( v(2)%pinfo%ir )
    dV  = v(1)%factor
    zz  = 0.5d0*dV

    ct(:)=0.d0
    et(:)=0.d0

    call watch(ct(1),et(1))

    allocate( Hsub(MB,MB)  ); Hsub=zero
    allocate( psiu(ML0,MB) ); psiu=zero
    allocate( psit(ML0,MB) ); psit=zero

    call watch(ct(2),et(2))

    np_band   = v(2)%pinfo%np
    comm_band = v(2)%pinfo%comm
    allocate( ir(0:np_band-1) ); ir=ML0*v(2)%pinfo%ir
    allocate( id(0:np_band-1) ); id=ML0*v(2)%pinfo%id
    call rsdft_allgatherv( unk(:,:,k,s), psiu, ir, id, comm_band )
    deallocate( id )
    deallocate( ir )

    do m=1,MB
       call hamiltonian(k,s,psiu(1,m),psit(1,m),n1,n2,m,m)
    end do

    call watch(ct(3),et(3))

#ifdef _DRSDFT_
    call dsyr2k('U','T',MB,ML0,zz,psiu,ML0,psit,ML0,zero,Hsub,MB)
!    call dgemm('T','N',MB,MB,ML0, dV,psiu,ML0,psit,ML0,zero,Hsub,MB)
#else
    call zher2k('U','C',MB,ML0,zz,psiu,ML0,psit,ML0,zero,Hsub,MB)
!    call zgemm('C','N',MB,MB,ML0,zdV,psiu,ML0,psit,ML0,zero,Hsub,MB)
#endif

    call watch(ct(4),et(4))

    deallocate( psit )

! --- allreduce ---

    call watch(ct(5),et(5))

    comm_grid = v(1)%pinfo%comm
    call rsdft_allreduce_sum( Hsub, comm_grid )

    call watch(ct(6),et(6))

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

       call zheevd('V','U',MB,Hsub,MB,esp(1,k,s) &
                  ,work,LWORK,rwork,LRWORK,iwork,LIWORK,ierr)
       if ( ierr==0 ) WORK1=nint( real(work(1)) )

       deallocate( work,rwork,iwork )

    case('dsyev')

       if ( LWORK==0 ) LWORK=3*MB-1

       allocate( rwork(LWORK) )

       call DSYEV('V','U',MB,Hsub,MB,esp(1,k,s),rwork,LWORK,ierr)

       WORK1=nint(rwork(1))

       deallocate( rwork )

    case('dsyevd')

       if ( LWORK==0  ) LWORK=1+6*MB+2*MB*MB
       if ( LIWORK==0 ) LIWORK=3+5*MB
       if ( LWORK == 0 .and. LIWORK == 0 ) then
          call DSYEVD('V','U',MB,Hsub,MB,esp(1,k,s),rtmp,-1,itmp,-1,ierr)
          LWORK=nint(rtmp(1))
          LIWORK=itmp(1)
       end if

       allocate( rwork(LWORK),iwork(LIWORK) )

       call DSYEVD('V','U',MB,Hsub,MB,esp(1,k,s) &
                   ,rwork,LWORK,iwork,LIWORK,ierr)

       WORK1=nint(rwork(1))
       WORK2=iwork(1)

       deallocate( iwork,rwork )

    end select

    call watch(ct(8),et(8))

! --- Rotation ---

#ifdef _DRSDFT_
    call dgemm('N','N',ML0,MB0,MB,one,psiu,ML0,Hsub(1,m1),MB,zero,unk(n1,m1,k,s),ML0)
#else
    call zgemm('N','N',ML0,MB0,MB,one,psiu,ML0,Hsub(1,m1),MB,zero,unk(n1,m1,k,s),ML0)
#endif

    call watch(ct(9),et(9))

    deallocate( psiu )
    deallocate( Hsub )

    LWORK=max(LWORK,WORK1)
    LIWORK=max(LIWORK,WORK2)

    call result_timer( t, "sd" )
    call write_border( 1, " subspace_diag_la_bp(end)" )

    return

  end subroutine subspace_diag_la_bp

end module subspace_diag_la_bp_module
