module watch_module

  use io_tools_module

  implicit none

  PRIVATE
  PUBLIC :: read_watch
  PUBLIC :: global_watch
  PUBLIC :: watch,watcht
  PUBLIC :: watcha, write_watcha
  PUBLIC :: watchb, write_watchb, watchb_omp
  PUBLIC :: time_cgpc, time_hmlt, time_kine, time_nlpp, time_bcfd
  PUBLIC :: time_cgpc_indx, time_hmlt_indx, time_kine_indx, time_nlpp_indx
  PUBLIC :: time
  PUBLIC :: init_time_watch, calc_time_watch
  public :: start_timer, result_timer, end_timer

  integer,parameter :: DP=kind(0.0d0)

  real(8) :: ct0=0.d0, ctt=0.d0
  real(8) :: ett=0.d0
  integer :: count0=0

  real(8) :: time_cgpc(2,16)
  real(8) :: time_hmlt(2,5)
  real(8) :: time_kine(2,16)
  real(8) :: time_nlpp(2,8)
  real(8) :: time_bcfd(2,8)
  character(5) :: time_hmlt_indx(5)
  character(5) :: time_kine_indx(16)
  character(5) :: time_nlpp_indx(8)
  character(5) :: time_cgpc_indx(16)
  data time_hmlt_indx(1:5)/"kine","loc","nlc","exx","bak"/
  data time_kine_indx(1:11)/"kine1","kine2","bc","kine3","kine4" &
                     ,"recv","send","spack","waita","final","totbc"/
  data time_nlpp_indx(1:7)/"nlc1","nlcom","nlc2","com1","com2","com3","com4"/
  data time_cgpc_indx(8:13)/"recv","send","spack","waita","final","totbc"/

  real(8) :: etime_limit=1.0d100

  real(8) :: global_ctime0,global_etime0
  logical :: flag_count_start=.true.

  integer,parameter :: max_timea_counter=64
  real(8) :: timea(0:max_timea_counter,2)

  type time
     real(8) :: t0, t1
     real(8) :: tmin, tmax
  end type time

  type(time) :: t_save

CONTAINS


  SUBROUTINE read_watch
    implicit none
    call IOTools_readReal8Keyword( "ETLIMIT", etime_limit )
  END SUBROUTINE read_watch


  SUBROUTINE watch( cpu_time, elapsed_time )
    implicit none
    real(DP),optional,intent(OUT) :: cpu_time
    real(DP),optional,intent(OUT) :: elapsed_time
    if ( present(cpu_time)     ) call get_cpu_time( cpu_time )
    if ( present(elapsed_time) ) call get_elapsed_time( elapsed_time )
  END SUBROUTINE watch

  SUBROUTINE get_cpu_time( ct )
    implicit none
    real(DP),intent(OUT) :: ct
    call cpu_time( ct )
  END SUBROUTINE get_cpu_time

  SUBROUTINE get_elapsed_time( et )
    implicit none
    real(DP),intent(OUT) :: et
    integer :: count, count_rate
    real(DP) :: cn,cd
#ifdef _NOMPI_
    call system_clock( count, count_rate )
    cn = count
    cd = count_rate
    et = cn/cd
#else
    include 'mpif.h'
    et = mpi_wtime()
#endif
  END SUBROUTINE get_elapsed_time


  SUBROUTINE watcht(disp_switch,indx,icnt)
    logical,intent(IN) :: disp_switch
    character(*),intent(IN) :: indx 
    integer,intent(IN) :: icnt
    integer :: count,count_rate
    real(8) :: ct
    call cpu_time(ct)
    call system_clock(count,count_rate)
    if ( icnt == 0 ) then
       ctt=0.d0
       ett=0.d0
    else if ( icnt == 1 ) then
       ctt=ct-ct0
       ett=real(count-count0)/real(count_rate)
       if ( disp_switch ) write(*,*) "timet(",indx,")=",ctt,ett
    end if
    ct0=ct
    count0=count
  END SUBROUTINE watcht

  SUBROUTINE global_watch(disp_switch,flag_timelimit,indx)
    implicit none
    logical,intent(IN) :: disp_switch
    logical,optional,intent(OUT) :: flag_timelimit
    character(*),optional,intent(in) :: indx
    integer :: ierr
    real(8) :: ct,et,s(2),r(2)
    include 'mpif.h'
    call watch( ct, et )
    if ( flag_count_start ) then
       global_ctime0=ct
       global_etime0=et
       flag_count_start=.false.
       return
    end if
    s(1)=ct-global_ctime0
    s(2)=et-global_etime0
    call mpi_allreduce(s,r,2,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
    if ( present(flag_timelimit) ) then
       flag_timelimit=.false.
       if ( r(2) > etime_limit ) flag_timelimit=.true.
    end if
    if ( disp_switch ) then
       if ( present(indx) ) then
          write(*,'(1x,"TIME(",a,")",3f12.5)') trim(indx),r(1:2)
       else
          write(*,'(1x,"TIME(END)",3f12.5)') r(1:2)
       end if
    end if
  END SUBROUTINE global_watch


  SUBROUTINE watcha( icounter )
    implicit none
    integer,intent(INOUT) :: icounter
    integer :: n
    include 'mpif.h'
    n = min( max(icounter+1,0), max_timea_counter )
    call cpu_time( timea(n,1) )
    timea(n,2) = mpi_wtime()
    icounter = n
  END SUBROUTINE watcha

  SUBROUTINE write_watcha( n, indx, unit )
    implicit none
    integer,intent(IN) :: n
    character(*),intent(IN) :: indx
    integer,optional,intent(IN) :: unit
    integer :: i,j,u
    u=6
    if ( present(unit) ) u=unit
    do i=1,n
       write(u,'(1x,i3,2x,a,3x,2x,2f15.5)') &
            i, "timea("//indx//")", (timea(i,j)-timea(i-1,j),j=1,2)
    end do
  END SUBROUTINE write_watcha


  SUBROUTINE watchb( t_tmp, t_out, barrier )
    implicit none
    real(8),intent(INOUT) :: t_tmp(2)
    real(8),optional,intent(INOUT) :: t_out(2)
    character(2),optional,intent(in) :: barrier
    real(8) :: t_now(2)
    integer :: ierr
    include 'mpif.h'
    if ( present(barrier) ) then
       if ( barrier == "on" ) call MPI_Barrier( MPI_COMM_WORLD, ierr )
    end if
    call watch( t_now(1), t_now(2) )
    if ( present(t_out) ) t_out(:) = t_out(:) + t_now(:) - t_tmp(:)
    t_tmp(:) = t_now(:)
  END SUBROUTINE watchb


  SUBROUTINE write_watchb( t_results, n, indx, unit )
    implicit none
    integer,intent(IN) :: n
    real(8),intent(IN) :: t_results(2,n)
    character(*),intent(IN) :: indx(n)
    integer,optional,intent(IN) :: unit
    integer :: i,j,u
    u=6
    if ( present(unit) ) u=unit
    do i=1,n
       write(u,'(1x,i3,2x,a,3x,2x,2f15.5)') &
            i, "timeb("//indx(i)//")", (t_results(j,i),j=1,2)
    end do
  END SUBROUTINE write_watchb


  SUBROUTINE watchb_omp( t_tmp, t_out )
    implicit none
    real(8),intent(INOUT) :: t_tmp(2)
    real(8),optional,intent(INOUT) :: t_out(2)
    real(8) :: tnow(2)
    include 'mpif.h'
!$OMP master
    call cpu_time( tnow(1) )
    tnow(2) = mpi_wtime()
    if ( present(t_out) ) t_out(:) = t_out(:) + tnow(:) - t_tmp(:)
    t_tmp(:) = tnow(:)
!$OMP end master
  END SUBROUTINE watchb_omp


  SUBROUTINE init_time_watch( t )
    implicit none
    type(time),intent(INOUT) :: t
    include 'mpif.h'
    t%t0 = mpi_wtime()
  END SUBROUTINE init_time_watch

  SUBROUTINE calc_time_watch( t, comm_in )
    implicit none
    type(time),intent(INOUT) :: t
    integer,optional,intent(IN) :: comm_in
    integer :: i, comm
    include 'mpif.h'
    t%t1 = mpi_wtime()
    t%t0 = t%t1 - t%t0
    comm=MPI_COMM_WORLD ; if ( present(comm_in) ) comm=comm_in
    call mpi_allreduce( t%t0, t%tmin, 1, MPI_REAL8, MPI_MIN, comm, i )
    call mpi_allreduce( t%t0, t%tmax, 1, MPI_REAL8, MPI_MAX, comm, i )
  END SUBROUTINE calc_time_watch


  subroutine start_timer( indx, t_out )
    implicit none
    character(*),optional,intent(in) :: indx
    type(time),optional,intent(out) :: t_out
    type(time) :: t
    include 'mpif.h'
    call cpu_time( t%t0 )
    t%t1 = mpi_wtime()
    if ( present(t_out) ) then
      t_out%t0=t%t0
      t_out%t1=t%t1
    end if
    t_save%t0=t%t0
    t_save%t1=t%t1
  end subroutine start_timer

  subroutine result_timer( indx, t_inout )
    implicit none
    character(*),optional,intent(in) :: indx
    type(time),optional,intent(inout) :: t_inout
    type(time) :: tnow,t
    character(64) :: mesg
    include 'mpif.h'
    call cpu_time( tnow%t0 )
    tnow%t1 = mpi_wtime()
    t=t_save ; if ( present(t_inout) ) t=t_inout
    t%t0 = tnow%t0 - t%t0
    t%t1 = tnow%t1 - t%t1
    if ( present(indx) ) then
       write(mesg,'(1x,"time(",a,")=",2f10.4)') indx, t%t0, t%t1
    else
       write(mesg,'(1x,"time=",2f10.4)') t%t0, t%t1
    end if
    call write_string_log( mesg )
  end subroutine result_timer

  subroutine end_timer( indx, t_inout )
    implicit none
    character(*),optional,intent(in) :: indx
    type(time),optional,intent(inout) :: t_inout
    type(time) :: t_now, t
    character(64) :: mesg
    include 'mpif.h'
    call cpu_time( t_now%t0 )
    t_now%t1 = mpi_wtime()
    if ( present(t_inout) ) then
      t%t0 = t_now%t0 - t_inout%t0
      t%t1 = t_now%t1 - t_inout%t1
    else
      t%t0 = t_now%t0 - t_save%t0
      t%t1 = t_now%t1 - t_save%t1
    end if
    if ( present(indx) ) then
      write(mesg,'(1x,"time(",a,")=",2f10.4)') indx, t%t0, t%t1
    else
      write(mesg,'(1x,"time=",2f10.4)') t%t0, t%t1
    end if
    call write_string_log( mesg )
  end subroutine end_timer


end module watch_module
