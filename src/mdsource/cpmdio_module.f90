MODULE cpmdio_module

  use cpmd_variables
  use wf_module, only: unk
  use io_tools_module
  use rsdft_mpi_module
  use io1_module, only: read_data_io1
  use io2_module
  use io_write_module, only: simple_wf_io_write
  use rgrid_module, only: Igrid,Ngrid

  implicit none

  PRIVATE
  PUBLIC :: write_data_cpmdio
  PUBLIC :: read_data_cpmdio

  integer :: IO_ctrl_r=0
  integer :: IO_ctrl_w=0
  integer :: OC=0
  logical :: flag_init_r=.true.
  logical :: flag_init_w=.true.
  integer :: node_partition_old(3)

CONTAINS


  SUBROUTINE write_data_cpmdio
    implicit none
    logical :: disp_sw
    integer :: i
    disp_sw=(myrank==0)
    if ( flag_init_w ) then
       if ( ctrl_cpmdio < 100 ) then
          IO_ctrl_w = ctrl_cpmdio
       else
          call IOTools_readIntegerKeyword( "IOCTRL", IO_ctrl_w )
       end if
       call IOTools_readIntegerKeyword( "OC", OC )
       flag_init_w=.false.
    end if
    select case( IO_ctrl_w )
    case default
       call write_data_cpmd_k_seri
    case( 1 )
       i=1
       if ( OC == 0 ) OC=3
       call simple_wf_io_write( "restart", IO_ctrl_w, OC, SYStype &
            , i, MBC, disp_sw, psi_v(:,MB_0_CPMD:MB_1_CPMD,:,:) &
            , MBC, MB_0_CPMD, MB_1_CPMD )
    case( 3 )
       call write_data_cpmd_k_para
    end select
  END SUBROUTINE write_data_cpmdio


  SUBROUTINE read_data_cpmdio
    implicit none
    if ( flag_init_r ) then
       if ( ctrl_cpmdio_r < 100 ) then
          IO_ctrl_r = ctrl_cpmdio_r
       else
          call IOTools_readIntegerKeyword( "IOCTRL", IO_ctrl_r )
       end if
       node_partition_old=1
       call IOTools_readIntegerKeyword( "OLD_PROCS", node_partition_old )
       flag_init_r=.false.
    end if
    select case( IO_ctrl_r )
    case default
       call read_data_cpmd_k_seri
    case( 1 )
       call read_data_io1( "restart", SYStype &
            ,wf_out=psi_v(:,MB_0_CPMD:MB_1_CPMD,:,:) &
            ,MB_in=MBC, MB_0_in=MB_0_CPMD, MB_1_in=MB_1_CPMD )
    case( 3 )
       call read_data_cpmd_k_para
    case( 4 )
       call read_data_cpmd_k_para_tmp
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


  SUBROUTINE read_data_cpmd_k_para_tmp
    implicit none
    integer :: n1,n2,n3,ML0_old(3),n,k,i,j,ispin,ierr
    integer :: i1,i2,i3
    integer :: a1_old,a2_old,a3_old,b1_old,b2_old,b3_old
    character(len=64) :: filename
    integer :: n1_now, nprocs_old, irank_old
    integer,allocatable :: Igrid_old(:,:,:)
    real(8) :: tmp0,tmp1
    real(8),allocatable :: w_old(:,:,:)

    n1_now  = idisp(myrank)+1

    nprocs_old = product( node_partition_old )
    ML0_old(1:3) = Ngrid(1:3)/node_partition_old(1:3)

    allocate( Igrid_old(2,3,0:nprocs_old-1) ); Igrid_old=0

    irank_old=-1
    do i3=0,node_partition_old(3)-1
    do i2=0,node_partition_old(2)-1
    do i1=0,node_partition_old(1)-1
       irank_old=irank_old+1
       Igrid_old(1,1,irank_old)=i1*ML0_old(1)
       Igrid_old(2,1,irank_old)=(i1+1)*ML0_old(1)
       Igrid_old(1,2,irank_old)=i2*ML0_old(2)
       Igrid_old(2,2,irank_old)=(i2+1)*ML0_old(2)
       Igrid_old(1,3,irank_old)=i3*ML0_old(3)
       Igrid_old(2,3,irank_old)=(i3+1)*ML0_old(3)
    end do
    end do
    end do

    allocate( w_old(ML0_old(1),ML0_old(2),ML0_old(3)) ); w_old=zero

    do irank_old=0,nprocs_old-1

       a1_old = Igrid_old(1,1,irank_old)
       a2_old = Igrid_old(1,2,irank_old)
       a3_old = Igrid_old(1,3,irank_old)
       b1_old = Igrid_old(2,1,irank_old)
       b2_old = Igrid_old(2,2,irank_old)
       b3_old = Igrid_old(2,3,irank_old)

       write(filename,"('restart_',i5.5,'.dat')") irank_old
       open(1,file=filename,form="unformatted")

       do ispin=MSP_0,MSP_1 
       do k=MBZ_0,MBZ_1
       do n=MB_0_CPMD,MB_1_CPMD
          do i3=1,ML0_old(3)
          do i2=1,ML0_old(2)
          do i1=1,ML0_old(1)
             read(1) w_old(i1,i2,i3)
          end do
          end do
          end do
          i=n1_now-1
          do i3=Igrid(1,3),Igrid(2,3)
          do i2=Igrid(1,2),Igrid(2,2)
          do i1=Igrid(1,1),Igrid(2,1)
             i=i+1
             if ( a1_old <= i1 .and. i1 <= b1_old .and. &
                  a2_old <= i2 .and. i2 <= b2_old .and. &
                  a3_old <= i3 .and. i3 <= b3_old ) then
                unk(i,n,k,ispin)=w_old(i1-a1_old+1,i2-a2_old+1,i3-a3_old+1)
             end if
          end do
          end do
          end do
       end do
       end do
       end do

       do ispin=MSP_0,MSP_1
       do k=MBZ_0,MBZ_1
       do n=MB_0_CPMD,MB_1_CPMD
          do i3=1,ML0_old(3)
          do i2=1,ML0_old(2)
          do i1=1,ML0_old(1)
             read(1) w_old(i1,i2,i3)
          end do
          end do
          end do
          i=n1_now-1
          do i3=Igrid(1,3),Igrid(2,3)
          do i2=Igrid(1,2),Igrid(2,2)
          do i1=Igrid(1,1),Igrid(2,1)
             i=i+1
             if ( a1_old <= i1 .and. i1 <= b1_old .and. &
                  a2_old <= i2 .and. i2 <= b2_old .and. &
                  a3_old <= i3 .and. i3 <= b3_old ) then
                psi_v(i,n,k,ispin)=w_old(i1-a1_old+1,i2-a2_old+1,i3-a3_old+1)
             end if
          end do
          end do
          end do
       end do
       end do
       end do

       close(1)

    end do ! irank_old

    deallocate( w_old )
    deallocate( Igrid_old )

    return
  END SUBROUTINE read_data_cpmd_k_para_tmp


END MODULE cpmdio_module
