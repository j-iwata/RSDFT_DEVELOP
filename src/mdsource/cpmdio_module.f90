module cpmdio_module

  use cpmd_variables, only: psi_v,MB_0_CPMD,MB_1_CPMD,MBZ_0,MBZ_1,ircnt,idisp,MSP_0,MSP_1,myrank &
                           ,ir_grid,id_grid,zero,comm_grid
  use wf_module, only: unk, occ
  use io_tools_module
  use rsdft_mpi_module
  use io1_module, only: read_data_io1
  use io2_module
  use io_write_wf_simple_module, only: write_wf_simple
  use rgrid_module, only: Igrid,Ngrid
  use parallel_module, only: node_partition

  implicit none

  private
  public :: write_data_cpmdio
  public :: read_data_cpmdio
  public :: read_data_cpmdio_0

  integer :: IO_ctrl_r=0
  integer :: IO_ctrl_w=0
  integer :: OC=0
  logical :: flag_init_r=.true.
  logical :: flag_init_w=.true.
  integer :: node_partition_old(7)

contains


  subroutine write_data_cpmdio
    implicit none
    logical :: disp_sw
    integer :: i1,i2,i3

    call IOTools_readIntegerKeyword( "OC", OC )
    if ( OC <= 0 ) return

    if ( myrank == 0 ) then
       open(1,file='restart_.dat')
       write(1,*) size(occ,1),size(occ,2),size(occ,3)
       do i3=1,size(occ,3)
       do i2=1,size(occ,2)
       do i1=1,size(occ,1)
          write(1,*) i1,i2,i3,occ(i1,i2,i3)
       end do
       end do
       end do
       write(1,'(1x,8i5)') node_partition(:)
       close(1)
    end if

    call write_data_cpmd_k_para

  end subroutine write_data_cpmdio

  subroutine read_data_cpmdio_0
    implicit none
    integer :: ierr,n1,n2,n3,i1,i2,i3,j1,j2,j3
    include 'mpif.h'
    if ( myrank == 0 ) then
       occ=0.0d0
       open(1,file='restart_.dat',status='old')
       read(1,*) n1,n2,n3
       do i3=1,n3
       do i2=1,n2
       do i1=1,n1
          read(1,*) j1,j2,j3,occ(i1,i2,i3)
       end do
       end do
       end do
       read(1,*) node_partition_old(:)
       close(1)
    end if
    call MPI_Bcast( occ,size(occ),MPI_REAL8,0,MPI_COMM_WORLD,ierr )
    call MPI_Bcast( node_partition_old,size(node_partition_old),MPI_INTEGER,0,MPI_COMM_WORLD,ierr )
  end subroutine read_data_cpmdio_0

  subroutine read_data_cpmdio
    implicit none
    if ( all(node_partition==node_partition_old) ) then
       call read_data_cpmd_k_para
    else
       call read_data_cpmd_k_para_tmp
    end if
  end subroutine read_data_cpmdio

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


  subroutine write_data_cpmd_k_para
    implicit none
    integer :: n1,n2,ML0,n,k,i,ispin,i1,i2,i3
    character(len=64) :: filename

    write(filename,"('restart_',i5.5,'.dat')") myrank

    open(1,file=filename,form="unformatted")

    n1  = idisp(myrank)+1
    n2  = idisp(myrank)+ircnt(myrank)
    ML0 = ircnt(myrank)

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
  end subroutine write_data_cpmd_k_para

!-------------------------------------------------------------

  SUBROUTINE read_data_cpmd_k_seri

    implicit none
    integer :: n1,n2,nn,s,k,n,ierr
    real(8),allocatable :: utmp(:)
    include 'mpif.h'

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


  subroutine read_data_cpmd_k_para
    implicit none
    integer :: n1,n2,ML0,n,k,i,ispin,ierr
    character(len=64) :: filename

!    call read_io2( SYStype, ierr )
!    if ( ierr == 0 ) then
!#ifdef _DRSDFT_
!       call read_data2_io2 &
!            ( "restart_",".dat", unk, psi_v, MB_0_CPMD, MB_1_CPMD )
!#endif
!       return
!    end if

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
  end subroutine read_data_cpmd_k_para


  subroutine read_data_cpmd_k_para_tmp
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
       Igrid_old(2,1,irank_old)=(i1+1)*ML0_old(1)-1
       Igrid_old(1,2,irank_old)=i2*ML0_old(2)
       Igrid_old(2,2,irank_old)=(i2+1)*ML0_old(2)-1
       Igrid_old(1,3,irank_old)=i3*ML0_old(3)
       Igrid_old(2,3,irank_old)=(i3+1)*ML0_old(3)-1
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
  end subroutine read_data_cpmd_k_para_tmp


end module cpmdio_module
