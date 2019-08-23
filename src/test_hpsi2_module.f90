MODULE test_hpsi2_module

  use array_bound_module, only: MB
  use localpot_module, only: Vloc
  use hamiltonian_module
  use parallel_module
  use watch_module
  use rsdft_mpi_module

  implicit none

  PRIVATE
  PUBLIC :: test_hpsi2

CONTAINS

  SUBROUTINE test_hpsi2(nloop)
    implicit none
    integer,intent(IN) :: nloop
#ifdef _DRSDFT_
    real(8),parameter :: zero=0.0d0, one=1.0d0
    real(8), allocatable :: tpsi(:,:)
    real(8), allocatable :: htpsi(:,:)
#else
    real(8),parameter :: zero=(0.0d0,0.0d0), one=(1.0d0,0.0d0)
    complex(8), allocatable :: tpsi(:,:)
    complex(8), allocatable :: htpsi(:,:)
#endif
    integer :: i, j, n1, n2, nrhs, niter, ierr
    integer :: loop,ii
    real(8), allocatable :: sums(:)
    real(8) :: t0, t1, t2, t3, t4, t5, ttmp(2), ttt(2,0:4)
    logical,save :: flag_allocate=.false.
    real(8) :: time_hmlt_min(2,4),time_hmlt_max(2,4)
    real(8) :: time_kine_min(2,11),time_kine_max(2,11)
    real(8) :: time_nlpp_min(2,3),time_nlpp_max(2,3)

    n1  = idisp(myrank)+1
    n2  = idisp(myrank)+ircnt(myrank)

    if ( .not.allocated(Vloc) ) then
       allocate( Vloc(n1:n2,1) ) ; Vloc=0.0d0
       flag_allocate=.true.
    end if

    nrhs=1

    do ii = 0,10

    time_hmlt(:,:)=0.0d0
    time_kine(:,:)=0.0d0
    time_nlpp(:,:)=0.0d0

    if ( ii > 0.and. MB_d == 1 ) then
       exit
    else if ( (ii>0.and.nrhs>=MB) .or. nrhs>MB_d ) then
       exit
    else if ( ii == 0 ) then
       nrhs=1
    else if ( mod(MB,nrhs*2) == 0 ) then
       nrhs=nrhs*2
    else if ( mod(MB,nrhs*3) == 0 ) then
       nrhs=nrhs*3
    else if ( mod(MB,nrhs*5) == 0 ) then
       nrhs=nrhs*5
    else if ( mod(MB,nrhs*7) == 0 ) then
       nrhs=nrhs*7
    else if ( MB_d >= MB ) then
       nrhs=MB
    end if

!    nrhs = 2**ii
!    if ( nrhs > MB_d ) exit

    niter = MB/nrhs

    allocate( tpsi(n1:n2,nrhs)  ) ; tpsi=one
    allocate( htpsi(n1:n2,nrhs) ) ; htpsi=zero
    allocate( sums(nrhs)        ) ; sums=0.0d0

!    call myrand(tpsi, n2-n1+1, nrhs)
  
    if ( DISP_SWITCH_PARALLEL ) then
       write(*,*) '------------ Start test_hpsi2 ------------'
    end if

    !t0 = mpi_wtime()
    call watchb( ttmp, barrier="on" ); ttt=0.0d0 

    do loop=1,nloop
       do i = 1, niter
          call hamiltonian(1, 1, tpsi, htpsi, n1, n2, 1, nrhs) 
       end do
    end do

    !t1 = mpi_wtime()
    call watchb( ttmp, ttt(:,0), barrier="on" )

    time_kine(:,:)=0.0d0
    time_nlpp(:,:)=0.0d0

    !t2 = mpi_wtime()
    call watchb( ttmp, barrier="on" )

!$OMP parallel private( loop, i )
    do loop=1,nloop
       do i=1,niter
          call op_kinetic( tpsi, htpsi )
       end do
    end do
!$OMP end parallel

    !t3 = mpi_wtime()
    call watchb( ttmp, ttt(:,1), barrier="on" )

!$OMP parallel private( loop, i )
    do loop=1,nloop
       do i=1,niter
          call op_localpot( tpsi, htpsi )
       end do
    end do
!$OMP end parallel

    !t4 = mpi_wtime()
    call watchb( ttmp, ttt(:,2), barrier="on" )

!$OMP parallel private( loop, i )
    do loop=1,nloop
       do i=1,niter
          call op_nonlocal( tpsi, htpsi )
       end do
    end do
!$OMP end parallel

    !t5 = mpi_wtime()
    call watchb( ttmp, ttt(:,3), barrier="on" )

    sums(:) = 0d0
    do i = 1, nrhs
       do j = n1, n2
#ifdef _DRSDFT_
          sums(i) = sums(i) + htpsi(j, i)*htpsi(j, i)
#else
          sums(i) = sums(i) + conjg(htpsi(j, i))*htpsi(j, i)
#endif
       end do
    end do

    call rsdft_allreduce_sum( sums(1:nrhs), comm_grid )

    ttt(1,4) = sum( ttt(1,1:3) )
    ttt(2,4) = sum( ttt(2,1:3) )

    call get_time_min(  4, time_hmlt, time_hmlt_min )
    call get_time_max(  4, time_hmlt, time_hmlt_max )
    call get_time_min( 11, time_kine, time_kine_min )
    call get_time_max( 11, time_kine, time_kine_max )
    call get_time_min(  3, time_nlpp, time_nlpp_min )
    call get_time_max(  3, time_nlpp, time_nlpp_max )

    if ( DISP_SWITCH_PARALLEL ) then
       write(*,*) 'nloop =',nloop
       write(*,*) 'nrhs = ', nrhs
       write(*,*) 'niter = ', niter
       write(*,*) 'niter*nrhs =',niter*nrhs,MB
       write(*,*) 'time(tot) = ', ttt(:,0) !t1 - t0
       call write_watchb( time_hmlt, 4, time_hmlt_indx )
       write(*,*) "(min)"
       call write_watchb( time_hmlt_min, 4, time_hmlt_indx )
       write(*,*) "(max)"
       call write_watchb( time_hmlt_max, 4, time_hmlt_indx )
       write(*,*) 'time(kin) = ', ttt(:,1) !t3 - t2
       call write_watchb( time_kine, 11, time_kine_indx )
       write(*,*) "(min)"
       call write_watchb( time_kine_min, 11, time_kine_indx )
       write(*,*) "(max)"
       call write_watchb( time_kine_max, 11, time_kine_indx )
       write(*,*) 'time(loc) = ', ttt(:,2) !t4 - t3
       write(*,*) 'time(nlc) = ', ttt(:,3) !t5 - t4
       write(*,*) 'time(tot) = ', ttt(:,4) !t5 - t2
       call write_watchb( time_nlpp, 3, time_nlpp_indx )
       write(*,*) "(min)"
       call write_watchb( time_nlpp_min, 3, time_nlpp_indx )
       write(*,*) "(max)"
       call write_watchb( time_nlpp_max, 3, time_nlpp_indx )
       !write(*,*) 'check sum'
       !do i = 1, nrhs
       !   write(*,*) sums(i),sum(sums)
       !end do
    end if

    deallocate( sums )
    deallocate( htpsi )
    deallocate( tpsi )

    end do ! ii

    if ( flag_allocate ) then
       deallocate( Vloc )
       flag_allocate=.false.
    end if

  END SUBROUTINE test_hpsi2

  subroutine myrand(V, nrow, ncol)
    implicit none
#ifdef _DRSDFT_
    real(8), intent(out) :: V(nrow,*)
#else
    complex(8), intent(out) :: V(nrow,*)  
#endif
    real(8), allocatable :: tmpV(:,:)
    integer, intent(in) :: nrow, ncol
    integer :: iseed(4)

    allocate(tmpV(nrow,ncol))
    iseed(1) = modulo(myrank_g-2*4096, 4096) ! must be between 0 and 4095
    iseed(2) = modulo(myrank_g-4096, 4096) ! must be between 0 and 4095
    iseed(3) = modulo(myrank_g, 4096) ! must be between 0 and 4095
    iseed(4) = 1 ! must be between 0 and 4095 and odd
  
    call DLARNV(2, iseed, nrow*ncol, tmpV)
    V(:,1:ncol) = tmpV(:,1:ncol)
  end subroutine myrand

  subroutine get_time_min( n, t_in, t_min )
    implicit none
    integer,intent(IN)  :: n
    real(8),intent(IN)  :: t_in(:,:)
    real(8),intent(OUT) :: t_min(:,:)
    integer :: i
    call mpi_allreduce( t_in , t_min, 2*n, mpi_real8, mpi_min, comm_grid, i )
  end subroutine get_time_min

  subroutine get_time_max( n, t_in, t_max )
    implicit none
    integer,intent(IN)  :: n
    real(8),intent(IN)  :: t_in(:,:)
    real(8),intent(OUT) :: t_max(:,:)
    integer :: i
    call mpi_allreduce( t_in , t_max, 2*n, mpi_real8, mpi_max, comm_grid, i )
  end subroutine get_time_max

END MODULE test_hpsi2_module


