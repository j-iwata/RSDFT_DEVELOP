module cpmdio_module

  use cpmd_variables, only: psi_v,MB_0_CPMD,MB_1_CPMD,MBZ_0,MBZ_1, &
  ircnt,idisp,MSP_0,MSP_1,myrank,ir_grid,id_grid,comm_grid
  use wf_module, only: unk, occ
  use io_tools_module
  use rsdft_mpi_module
  ! use io1_module, only: read_data_io1
  ! use io2_module
  ! use io_write_wf_simple_module, only: write_wf_simple
  use rgrid_module, only: Igrid,Ngrid
  use parallel_module, only: node_partition
  use watch_module, only: watchb

  implicit none

  private
  public :: write_data_cpmdio
  public :: read_data_cpmdio
  public :: read_data_cpmdio_0

  integer,allocatable :: node_partition_old(:)

contains

  subroutine write_data_cpmdio
    implicit none
    integer :: i1,i2,i3,OC

    call write_border( 0, ' write_data_cpmdio(start)' )

    call IOTools_readIntegerKeyword( "OC", OC )
    if ( OC <= 0 ) then
      if ( myrank == 0 ) write(*,*) "No data is written !"
      call write_border( 0, ' write_data_cpmdio(return)' )
      return
    end if

    if ( myrank == 0 ) then
      open(1,file='restart_.dat')
      write(1,*) size(occ,1),size(occ,2),size(occ,3)
      do i3 = 1, size(occ,3)
      do i2 = 1, size(occ,2)
      do i1 = 1, size(occ,1)
        write(1,*) i1,i2,i3, occ(i1,i2,i3)
      end do
      end do
      end do
      write(1,'(1x,8i5)') node_partition(:)
      close(1)
    end if

    call write_data_cpmd_k_para

    call write_border( 0, ' write_data_cpmdio(end)' )

  end subroutine write_data_cpmdio

  subroutine read_data_cpmdio_0
    implicit none
    integer :: ierr,n1,n2,n3,i1,i2,i3,j1,j2,j3
    real(8),allocatable :: occ_tmp(:,:,:)
    real(8) :: tmp1,tmp2
    include 'mpif.h'
    call write_border( 0, 'read_data_cpmdio_0(start)' )
    allocate( node_partition_old(size(node_partition)) )
    node_partition_old=0
    if ( myrank == 0 ) then
      occ=0.0d0
      open(1,file='restart_.dat',status='old')
      read(1,*) n1,n2,n3
      write(*,*) "size(occ_old),size(occ)",n1,n2,n3,size(occ,1),size(occ,2),size(occ,3)
      allocate( occ_tmp(n1,n2,n3) ); occ_tmp=0.0d0
      do i3=1,n3
      do i2=1,n2
      do i1=1,n1
        read(1,*) j1,j2,j3,occ_tmp(i1,i2,i3)
      end do
      end do
      end do
      read(1,*) node_partition_old(:)
      close(1)
      write(*,*) "sum,min,max(occ_tmp)=",sum(occ_tmp),minval(occ_tmp),maxval(occ_tmp)
      do i3=1,min(size(occ,3),n3)
      do i2=1,min(size(occ,2),n2)
      do i1=1,min(size(occ,1),n1)
        occ(i1,i2,i3) = occ_tmp(i1,i2,i3)
      end do
      end do
      end do
      write(*,*) "sum,min,max(occ)=",sum(occ),minval(occ),maxval(occ)
      tmp1=sum(occ_tmp)
      tmp2=sum(occ)
      if ( abs(tmp1-tmp2) > 1.0d-8 ) occ=-1.0d0
      deallocate( occ_tmp )
    end if
    call MPI_Bcast( occ,size(occ),MPI_REAL8,0,MPI_COMM_WORLD,ierr )
    call MPI_Bcast( node_partition_old,size(node_partition_old),MPI_INTEGER,0,MPI_COMM_WORLD,ierr )
    call write_border( 0, 'read_data_cpmdio_0(end)' )
  end subroutine read_data_cpmdio_0

  subroutine read_data_cpmdio
    implicit none
    logical :: new_format
    integer :: i
    real(8) :: tt(2,0:1)
    call write_border( 0, ' read_data_cpmdio(start)' )
    call watchb( tt(:,0), barrier='on' )
    new_format = check_data_format()
    if ( all(node_partition==node_partition_old) ) then
      call read_data_cpmd_k_para( new_format )
    else
      call read_data_cpmd_k_para_arb( new_format )
    end if
    call watchb( tt(:,1), barrier='on' )
    if ( myrank == 0 ) write(*,'(1x,"time(read):",2f10.3)') (tt(i,1)-tt(i,0),i=1,2)
    call write_border( 0, ' read_data_cpmdio(end)' )
  end subroutine read_data_cpmdio

!-------------------------------------------------------------

!   SUBROUTINE write_data_cpmd_k_seri

!     implicit none
!     integer :: nn,s,k,n
! #ifdef _DRSDFT_
!     real(8),allocatable :: utmp(:)
! #else
!     complex(8),allocatable :: utmp(:)
! #endif

!     call write_border( 0, 'write_data_cpmd_k_seri(start)' )

!     nn=sum(ircnt)

!     allocate( utmp(nn) ) ; utmp=0.0d0

!     if ( myrank == 0 ) open(1,file='restart_00000.dat',form='unformatted')

!     do s=MSP_0,MSP_1
!     do k=MBZ_0,MBZ_1
!     do n=MB_0_CPMD,MB_1_CPMD
!        call rsdft_allgatherv( unk(:,n,k,s),utmp,ir_grid,id_grid,comm_grid )
!        if ( myrank == 0 ) write(1) utmp(:)
!     end do
!     end do
!     end do

!     do s=MSP_0,MSP_1
!     do k=MBZ_0,MBZ_1
!     do n=MB_0_CPMD,MB_1_CPMD
! #ifdef _DRSDFT_
!        call rsdft_allgatherv( psi_v(:,n,k,s),utmp,ir_grid,id_grid,comm_grid )
! #endif
!        if ( myrank == 0 ) write(1) utmp(:)
!     end do
!     end do
!     end do

!     deallocate( utmp )

!     if ( myrank == 0 ) close(1)

!     call write_border( 0, 'write_data_cpmd_k_seri(end)' )

!   END SUBROUTINE write_data_cpmd_k_seri


  subroutine write_data_cpmd_k_para
    implicit none
    integer :: n1,n2,i,n,k,s
    character(len=17) :: filename
    real(8) :: tt(2,0:1)

    call write_border( 0, ' write_data_cpmd_k_para(start)' )
    call watchb( tt(:,0), barrier='on' )

    write(filename,"('restart_',i5.5,'.dat')") myrank

    open(1,file=filename,form="unformatted")

    n1 = idisp(myrank)+1
    n2 = idisp(myrank)+ircnt(myrank)

    do s = MSP_0, MSP_1
    do k = MBZ_0, MBZ_1
    do n = MB_0_CPMD, MB_1_CPMD
      ! do i = n1, n2
      !   write(1) unk(i,n,k,s)
      ! end do
      write(1) unk(:,n,k,s)
    end do
    end do
    end do

    do s = MSP_0, MSP_1
    do k = MBZ_0, MBZ_1
    do n = MB_0_CPMD, MB_1_CPMD
      ! do i = n1, n2
      !   write(1) psi_v(i,n,k,s)
      ! end do
      write(1) psi_v(:,n,k,s)
    end do
    end do
    end do

    close(1)

    call watchb( tt(:,1), barrier='on' )
    if ( myrank == 0 ) write(*,'(1x,"time(write):",2f10.3)') (tt(i,1)-tt(i,0),i=1,2)

    call write_border( 0, ' write_data_cpmd_k_para(end)' )

    return
  end subroutine write_data_cpmd_k_para

!-------------------------------------------------------------

  ! SUBROUTINE read_data_cpmd_k_seri

  !   implicit none
  !   integer :: n1,n2,nn,s,k,n,ierr
  !   real(8),allocatable :: utmp(:)
  !   include 'mpif.h'

  !   call write_border( 0, 'read_data_cpmd_k_seri(start)' )

  !   n1=idisp(myrank)+1
  !   n2=idisp(myrank)+ircnt(myrank)
  !   nn=sum(ircnt)

  !   allocate( utmp(nn) ) ; utmp=0.0d0

  !   if ( myrank == 0 ) open(1,file='restart_00000.dat',form='unformatted')

  !   do s=MSP_0,MSP_1
  !   do k=MBZ_0,MBZ_1
  !   do n=MB_0_CPMD,MB_1_CPMD
  !      if ( myrank == 0 ) read(1) utmp(:)
  !      call MPI_SCATTERV( utmp,ir_grid,id_grid,MPI_REAL8 &
  !                       , unk(n1,n,k,s),n2-n1+1,MPI_REAL8,0,comm_grid,ierr )
  !   end do
  !   end do
  !   end do

  !   do s=MSP_0,MSP_1
  !   do k=MBZ_0,MBZ_1
  !   do n=MB_0_CPMD,MB_1_CPMD
  !      if ( myrank == 0 ) read(1) utmp(:)
  !      call MPI_SCATTERV( utmp,ir_grid,id_grid,MPI_REAL8 &
  !                       , psi_v(n1,n,k,s),n2-n1+1,MPI_REAL8,0,comm_grid,ierr )
  !   end do
  !   end do
  !   end do

  !   deallocate( utmp )

  !   if ( myrank == 0 ) close(1)

  !   call write_border( 0, 'read_data_cpmd_k_seri(end)' )

  ! END SUBROUTINE read_data_cpmd_k_seri


  logical function check_data_format()
    implicit none
    integer :: n1,n2,i,n,k,s,ierr
    integer,parameter :: u=700
    character(len=17) :: filename
    logical :: flag_new_format
    real(8) :: a

    call write_border( 0, ' check_data_format(start)' )

    write(filename,"('restart_',i5.5,'.dat')") 0

    open(u,file=filename,form="unformatted")

    do i = 1, 3
      do s = MSP_0, MSP_1
      do k = MBZ_0, MBZ_1
      do n = MB_0_CPMD, MB_1_CPMD
        read(u,END=9) a
      end do
      end do
      end do
    end do
    9 continue

    close(u)

    flag_new_format=.true.
    if ( i>2 .and. s>MSP_1 .and. k>MBZ_1 .and. n>MB_1_CPMD ) flag_new_format=.false.

    if ( myrank == 0 ) then
      if ( flag_new_format ) write(*,*) "new data format is assumed"
    end if

    checK_data_format = flag_new_format

    call write_border( 0, ' check_data_format(end)' )
  end function check_data_format


  subroutine read_data_cpmd_k_para( new_format )
    implicit none
    logical,intent(in) :: new_format
    integer :: n1,n2,i,n,k,s,ierr
    integer,parameter :: u=700
    character(len=17) :: filename

    call write_border( 0, ' read_data_cpmd_k_para(start)' )

    write(filename,"('restart_',i5.5,'.dat')") myrank

    open(u,file=filename,form="unformatted")

    n1 = idisp(myrank)+1
    n2 = idisp(myrank)+ircnt(myrank)

    rewind u

    do s = MSP_0, MSP_1
    do k = MBZ_0, MBZ_1
    do n = MB_0_CPMD, MB_1_CPMD
      if ( new_format ) then
        read(u) unk(:,n,k,s)
      else
        do i = n1, n2
          read(u) unk(i,n,k,s)
        end do
      end if
    end do
    end do
    end do

    do s = MSP_0, MSP_1
    do k = MBZ_0, MBZ_1
    do n = MB_0_CPMD, MB_1_CPMD
      if ( new_format ) then
        read(u) psi_v(:,n,k,s)
      else
        do i = n1, n2
          read(u) psi_v(i,n,k,s)
        end do
      end if
    end do
    end do
    end do

    close(u)

    call write_border( 0, ' read_data_cpmd_k_para(end)' )

    return
  end subroutine read_data_cpmd_k_para


  subroutine read_data_cpmd_k_para_arb( new_format )
    use parallel_module, only: load_div_parallel,get_range_parallel, &
    ir_bzsm,ir_spin,myrank_b
    use cpmd_variables, only: id_band_cpmd,ir_band_cpmd
    implicit none
    logical,intent(in) :: new_format
    integer :: ML0_old(3),n,k,i,j,s,ierr
    integer :: i0,i1,i2,i3,i4,i5,i6,a1,a2,b1,b2,c1,c2,n1,n2,k1,k2,s1,s2
    integer :: a1_old,a2_old,b1_old,b2_old,c1_old,c2_old,ng_old
    integer :: n1_old,n2_old,k1_old,k2_old,s1_old,s2_old
    character(len=17) :: filename
    integer :: n1_now, nprocs_old, irank_old
    integer :: Nband, Nbzsm, Nspin
    integer,allocatable :: Igrd1_old(:,:),Igrd2_old(:,:),Igrd3_old(:,:)
    integer,allocatable :: Iband_old(:,:),Ibzsm_old(:,:),Ispin_old(:,:)
    integer,allocatable :: ir(:,:),id(:,:),ndata(:)
    integer,parameter :: u=700
    real(8),allocatable :: w_old(:,:,:), u_old(:)

    call write_border( 0, ' read_data_cpmd_k_para_arb(start)' )

    nprocs_old = product( node_partition_old )

    allocate( Igrd1_old(2,0:nprocs_old-1) ); Igrd1_old=0
    allocate( Igrd2_old(2,0:nprocs_old-1) ); Igrd2_old=0
    allocate( Igrd3_old(2,0:nprocs_old-1) ); Igrd3_old=0
    allocate( Iband_old(2,0:nprocs_old-1) ); Iband_old=0
    allocate( Ibzsm_old(2,0:nprocs_old-1) ); Ibzsm_old=0
    allocate( Ispin_old(2,0:nprocs_old-1) ); Ispin_old=0

    n = maxval( node_partition_old )
    allocate( ir(0:n-1,size(node_partition_old)) ); ir=0
    allocate( id(0:n-1,size(node_partition_old)) ); id=0
    allocate( ndata(size(node_partition_old))    ); ndata=0

    Nband = sum( ir_band_cpmd )
    Nbzsm = sum( ir_bzsm )
    Nspin = sum( ir_spin )

    ndata(1:6) = (/ Ngrid(1:3), Nband, Nbzsm, Nspin /)

    do i = 1, size(node_partition_old)
      n = node_partition_old(i)
      call load_div_parallel( ir(0:n-1,i), id(0:n-1,i), ndata(i) )
    end do

    deallocate( ndata )

    irank_old=-1
    do i6 = 0, node_partition_old(6)-1
    do i5 = 0, node_partition_old(5)-1
    do i4 = 0, node_partition_old(4)-1
      do i3 = 0, node_partition_old(3)-1
      do i2 = 0, node_partition_old(2)-1
      do i1 = 0, node_partition_old(1)-1
        irank_old = irank_old + 1
        Igrd1_old(1,irank_old) = id(i1,1)
        Igrd1_old(2,irank_old) = id(i1,1) + ir(i1,1) - 1
        Igrd2_old(1,irank_old) = id(i2,2)
        Igrd2_old(2,irank_old) = id(i2,2) + ir(i2,2) - 1
        Igrd3_old(1,irank_old) = id(i3,3)
        Igrd3_old(2,irank_old) = id(i3,3) + ir(i3,3) - 1
        Iband_old(1,irank_old) = id(i4,4) + 1
        Iband_old(2,irank_old) = id(i4,4) + ir(i4,4)
        Ibzsm_old(1,irank_old) = id(i5,5) + 1
        Ibzsm_old(2,irank_old) = id(i5,5) + ir(i5,5)
        Ispin_old(1,irank_old) = id(i6,6) + 1
        Ispin_old(2,irank_old) = id(i6,6) + ir(i6,6)
      end do
      end do
      end do
    end do
    end do
    end do

    ML0_old = (/ maxval(ir(:,1)), maxval(ir(:,2)), maxval(ir(:,3)) /)

    deallocate( id, ir )

    allocate( w_old(ML0_old(1),ML0_old(2),ML0_old(3)) ); w_old=0.0d0

    if ( new_format ) then
      allocate( u_old(size(w_old)) ); u_old=0.0d0
    end if

    a1 = Igrid(1,1) ; a2 = Igrid(2,1)
    b1 = Igrid(1,2) ; b2 = Igrid(2,2)
    c1 = Igrid(1,3) ; c2 = Igrid(2,3)
    ! call get_range_parallel( n1, n2, 'b' )
    n1 = id_band_cpmd(myrank_b) + 1
    n2 = n1 + ir_band_cpmd(myrank_b) - 1
    call get_range_parallel( k1, k2, 'k' )
    call get_range_parallel( s1, s2, 's' )

    do irank_old = 0, nprocs_old-1

      a1_old = Igrd1_old(1,irank_old)
      a2_old = Igrd1_old(2,irank_old)
      b1_old = Igrd2_old(1,irank_old)
      b2_old = Igrd2_old(2,irank_old)
      c1_old = Igrd3_old(1,irank_old)
      c2_old = Igrd3_old(2,irank_old)
      ng_old = (a2_old-a1_old+1)*(b2_old-b1_old+1)*(c2_old-c1_old+1)
      n1_old = Iband_old(1,irank_old)
      n2_old = Iband_old(2,irank_old)
      k1_old = Ibzsm_old(1,irank_old)
      k2_old = Ibzsm_old(2,irank_old)
      s1_old = Ispin_old(1,irank_old)
      s2_old = Ispin_old(2,irank_old)

      if ( a2_old < a1 .or. a2 < a1_old .or. &
           b2_old < b1 .or. b2 < b1_old .or. &
           c2_old < c1 .or. c2 < c1_old .or. &
           n2_old < n1 .or. n2 < n1_old .or. &
           k2_old < k1 .or. k2 < k1_old .or. &
           s2_old < s1 .or. s2 < s1_old ) cycle

      write(filename,"('restart_',i5.5,'.dat')") irank_old
      open(u,file=filename,form="unformatted")

      do s = s1_old, s2_old
      do k = k1_old, k2_old
      do n = n1_old, n2_old

        if ( new_format ) then
          read(u) u_old(1:ng_old)
          i=0
          do i3 = 1, c2_old-c1_old+1
          do i2 = 1, b2_old-b1_old+1
          do i1 = 1, a2_old-a1_old+1
            i=i+1
            w_old(i1,i2,i3) = u_old(i)
          end do
          end do
          end do
        else
          do i3 = 1, c2_old-c1_old+1
          do i2 = 1, b2_old-b1_old+1
          do i1 = 1, a2_old-a1_old+1
            read(u) w_old(i1,i2,i3)
          end do
          end do
          end do
        end if

        if ( MSP_0<=s .and. s<=MSP_1 .and. &
             MBZ_0<=k .and. k<=MBZ_1 .and. MB_0_CPMD<=n .and. n<=MB_1_CPMD ) then
          i=Igrid(1,0)-1
          do i3 = Igrid(1,3), Igrid(2,3)
          do i2 = Igrid(1,2), Igrid(2,2)
          do i1 = Igrid(1,1), Igrid(2,1)
            i=i+1
            if ( a1_old <= i1 .and. i1 <= a2_old .and. &
                 b1_old <= i2 .and. i2 <= b2_old .and. &
                 c1_old <= i3 .and. i3 <= c2_old ) then
              unk(i,n,k,s) = w_old(i1-a1_old+1,i2-b1_old+1,i3-c1_old+1)
            end if
          end do
          end do
          end do
        end if

      end do !n
      end do !k
      end do !s

      do s = s1_old, s2_old
      do k = k1_old, k2_old
      do n = n1_old, n2_old

        if ( new_format ) then
          read(u) u_old(1:ng_old)
          i=0
          do i3 = 1, c2_old-c1_old+1
          do i2 = 1, b2_old-b1_old+1
          do i1 = 1, a2_old-a1_old+1
            i=i+1
            w_old(i1,i2,i3) = u_old(i)
          end do
          end do
          end do
        else
          do i3 = 1, c2_old-c1_old+1
          do i2 = 1, b2_old-b1_old+1
          do i1 = 1, a2_old-a1_old+1
            read(u) w_old(i1,i2,i3)
          end do
          end do
          end do
        end if

        if ( MSP_0<=s .and. s<=MSP_1 .and. &
             MBZ_0<=k .and. k<=MBZ_1 .and. MB_0_CPMD<=n .and. n<=MB_1_CPMD ) then
          i=Igrid(1,0)-1
          do i3 = Igrid(1,3), Igrid(2,3)
          do i2 = Igrid(1,2), Igrid(2,2)
          do i1 = Igrid(1,1), Igrid(2,1)
            i=i+1
            if ( a1_old <= i1 .and. i1 <= a2_old .and. &
                 b1_old <= i2 .and. i2 <= b2_old .and. &
                 c1_old <= i3 .and. i3 <= c2_old ) then
              psi_v(i,n,k,s) = w_old(i1-a1_old+1,i2-b1_old+1,i3-c1_old+1)
            end if
          end do
          end do
          end do
        end if

      end do !n
      end do !k
      end do !s

      close(u)

    end do ! irank_old

    if ( allocated(u_old) ) deallocate( u_old )
    deallocate( w_old )
    deallocate( Ispin_old, Ibzsm_old, Iband_old )
    deallocate( Igrd3_old, Igrd2_old, Igrd1_old )

    call write_border( 0, ' read_data_cpmd_k_para_arb(end)' )

    return
  end subroutine read_data_cpmd_k_para_arb

end module cpmdio_module
