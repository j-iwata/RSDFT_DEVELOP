MODULE cpmdio_module

  use cpmd_variables
  use wf_module, only: unk
  use io_tools_module
  use rsdft_mpi_module
  use io2_module

  implicit none

  PRIVATE
  PUBLIC :: write_data_cpmdio
  PUBLIC :: read_data_cpmdio

  integer :: io_ctrl=0
  logical :: flag_init=.true.

CONTAINS


  SUBROUTINE write_data_cpmdio
    implicit none
    if ( flag_init ) then
       call IOTools_readIntegerKeyword( "IOCTRL", io_ctrl )
       flag_init=.false.
    end if
    select case( io_ctrl )
    case default
       call write_data_cpmd_k_seri
    case( 3 )
       call write_data_cpmd_k_para
    end select
  END SUBROUTINE write_data_cpmdio


  SUBROUTINE read_data_cpmdio
    implicit none
    if ( flag_init ) then
       call IOTools_readIntegerKeyword( "IOCTRL", io_ctrl )
       flag_init=.false.
    end if
    select case( io_ctrl )
    case default
       call read_data_cpmd_k_seri
    case( 3 )
       call read_data_cpmd_k_para
    end select
  END SUBROUTINE read_data_cpmdio

!-------------------------------------------------------------

  SUBROUTINE write_data_cpmd_k_seri

    implicit none
    integer :: nn,s,k,n
#ifdef _DRSDFT_
    real(8),allocatable :: utmp(:)
#else
    complex(8),allocatable :: utmp(:)
#endif

    nn=sum(ircnt)

    allocate( utmp(nn) ) ; utmp=0.0d0

    if ( myrank == 0 ) open(1,file='restart_00000.dat',form='unformatted')

    do s=MSP_0,MSP_1
    do k=MBZ_0,MBZ_1
    do n=MB_0_CPMD,MB_1_CPMD
       call rsdft_allgatherv( unk(:,n,k,s),utmp,ir_grid,id_grid,comm_grid )
       if ( myrank == 0 ) write(1) utmp(:)
    end do
    end do
    end do

    do s=MSP_0,MSP_1
    do k=MBZ_0,MBZ_1
    do n=MB_0_CPMD,MB_1_CPMD
#ifdef _DRSDFT_
       call rsdft_allgatherv( psi_v(:,n,k,s),utmp,ir_grid,id_grid,comm_grid )
#endif
       if ( myrank == 0 ) write(1) utmp(:)
    end do
    end do
    end do

    deallocate( utmp )

    if ( myrank == 0 ) close(1)

  END SUBROUTINE write_data_cpmd_k_seri


  SUBROUTINE write_data_cpmd_k_para
    implicit none
    integer :: n1,n2,ML0,n,k,i,ispin
    character(len=64) :: filename

    write(filename,"('restart_',i5.5,'.dat')") myrank

    open(1,file=filename,form="unformatted")
    n1  = idisp(myrank)+1
    n2  = idisp(myrank)+ircnt(myrank)
    ML0 = ircnt(myrank)

!      write(*,*) "MSP_0,MSP_1", MSP_0,MSP_1
!      write(*,*) "MBZ_0,MBZ_1",MBZ_0,MBZ_1

    do ispin=MSP_0,MSP_1
    do k=MBZ_0,MBZ_1
    do n=MB_0_CPMD,MB_1_CPMD
    do i=n1,n2
       write(1) unk(i,n,k,ispin)
    enddo
    enddo
    enddo
    enddo

    do ispin=MSP_0,MSP_1
    do k=MBZ_0,MBZ_1
    do n=MB_0_CPMD,MB_1_CPMD
    do i=n1,n2
       write(1) psi_v(i,n,k,ispin)
    enddo
    enddo
    enddo
    enddo

    close(1)
    return
  END SUBROUTINE write_data_cpmd_k_para

!-------------------------------------------------------------

  SUBROUTINE read_data_cpmd_k_seri

    implicit none
    integer :: n1,n2,nn,s,k,n,ierr
    real(8),allocatable :: utmp(:)

    n1=idisp(myrank)+1
    n2=idisp(myrank)+ircnt(myrank)
    nn=sum(ircnt)

    allocate( utmp(nn) ) ; utmp=0.0d0

    if ( myrank == 0 ) open(1,file='restart_00000.dat',form='unformatted')

    do s=MSP_0,MSP_1
    do k=MBZ_0,MBZ_1
    do n=MB_0_CPMD,MB_1_CPMD
       if ( myrank == 0 ) read(1) utmp(:)
       call MPI_SCATTERV( utmp,ir_grid,id_grid,MPI_REAL8 &
                        , unk(n1,n,k,s),n2-n1+1,MPI_REAL8,0,comm_grid,ierr )
    end do
    end do
    end do

    do s=MSP_0,MSP_1
    do k=MBZ_0,MBZ_1
    do n=MB_0_CPMD,MB_1_CPMD
       if ( myrank == 0 ) read(1) utmp(:)
       call MPI_SCATTERV( utmp,ir_grid,id_grid,MPI_REAL8 &
                        , psi_v(n1,n,k,s),n2-n1+1,MPI_REAL8,0,comm_grid,ierr )
    end do
    end do
    end do

    deallocate( utmp )

    if ( myrank == 0 ) close(1)

  END SUBROUTINE read_data_cpmd_k_seri


  SUBROUTINE read_data_cpmd_k_para
    implicit none
    integer :: n1,n2,ML0,n,k,i,ispin,ierr
    character(len=64) :: filename

    call read_io2( SYStype, ierr )
    if ( ierr == 0 ) then
#ifdef _DRSDFT_
       call read_data2_io2 &
            ( "restart_",".dat", unk, psi_v, MB_0_CPMD, MB_1_CPMD )
#endif
       return
    end if

    write(filename,"('restart_',i5.5,'.dat')") myrank

    open(1,file=filename,form="unformatted")
    n1  = idisp(myrank)+1
    n2  = idisp(myrank)+ircnt(myrank)
    ML0 = ircnt(myrank)

    do ispin=MSP_0,MSP_1
    do k=MBZ_0,MBZ_1
    do n=MB_0_CPMD,MB_1_CPMD
    do i=n1,n2
       read(1) unk(i,n,k,ispin)
    enddo
    enddo
    enddo
    enddo

    do ispin=MSP_0,MSP_1
    do k=MBZ_0,MBZ_1
    do n=MB_0_CPMD,MB_1_CPMD
    do i=n1,n2
       read(1) psi_v(i,n,k,ispin)
    enddo
    enddo
    enddo
    enddo

    close(1)
    return
  END SUBROUTINE read_data_cpmd_k_para


END MODULE cpmdio_module
