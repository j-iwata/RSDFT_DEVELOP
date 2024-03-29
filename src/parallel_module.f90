module parallel_module

  implicit none

  include 'mpif.h'

!  PRIVATE
  public :: init_parallel
  public :: start_mpi_parallel
  public :: end_mpi_parallel

  public :: para, construct_para
  public :: get_np_parallel
  public :: node_partition &
           ,comm_grid, comm_band, comm_bzsm, comm_spin &
           ,myrank,myrank_g,myrank_b,myrank_k,myrank_s &
           ,nprocs,nprocs_g,nprocs_k,nprocs_s &
           ,np_band,np_spin,np_bzsm,np_grid &
           ,id_class,ircnt,idisp,ir_spin,id_spin &
           ,ir_band,id_band,ir_grid,id_grid,ir_bzsm,id_bzsm &
           ,pinfo_grid
  public :: MB_d, MB_d_nl
  public :: disp_switch_parallel
  public :: comm_fkmb, myrank_f, np_fkmb, ir_fkmb, id_fkmb
  public :: get_range_parallel
  public :: construct_id_ir_parallel
  public :: get_ParaInfo
  public :: load_div_parallel
  public :: load_modify_parallel

  private :: dr, construct_dr
  private :: read_parallel
  private :: set_np_parallel
  private :: construct_ParaInfo

#ifdef _NO_MPI_COMPLEX16_
  integer,parameter,public :: RSDFT_MPI_COMPLEX16 = MPI_DOUBLE_COMPLEX
#else
  integer,parameter,public :: RSDFT_MPI_COMPLEX16 = MPI_COMPLEX16
#endif
  integer,parameter,public :: RSDFT_MPI_REAL8 = MPI_REAL8

  integer,parameter :: max_parallel=7
  integer :: node_partition(max_parallel)
  integer :: nprocs, myrank
  integer :: comm_grid, myrank_g, nprocs_g, np_grid
  integer :: comm_spin, myrank_s, nprocs_s, np_spin
  integer :: comm_band, myrank_b, nprocs_b, np_band
  integer :: comm_bzsm, myrank_k, nprocs_k, np_bzsm
  integer :: comm_fkmb, myrank_f, nprocs_f, np_fkmb
  integer :: comm_bks
  integer :: comm_gb
  integer :: MB_d, MB_d_nl
  integer,allocatable :: id_class(:,:),ircnt(:),idisp(:)
  integer,allocatable :: ir_grid(:),id_grid(:)
  integer,allocatable :: ir_band(:),id_band(:)
  integer,allocatable :: ir_bzsm(:),id_bzsm(:)
  integer,allocatable :: ir_spin(:),id_spin(:)
  integer,allocatable :: ir_fkmb(:),id_fkmb(:)
  integer,allocatable :: pinfo_grid(:,:)
  logical :: disp_switch_parallel=.false.

  type dr
    integer :: np
    integer,allocatable :: id(:)
    integer,allocatable :: ir(:)
  end type dr

  type para
    integer  :: np(0:7)
    type(dr) :: grid
    type(dr) :: band
    type(dr) :: bzsm
    type(dr) :: spin
  end type para

  type pinfo
    integer :: comm
    integer :: myrank
    integer :: nprocs
    integer,allocatable :: ir(:)
    integer,allocatable :: id(:)
  end type pinfo

  type(pinfo),private :: ParaInfo(6)
  integer,private :: nprocs_bks, nprocs_gb

contains


  subroutine get_np_parallel( np )
    implicit none
    integer,intent(out) :: np(0:)
    integer :: n
    np(0)=nprocs
    n=size(np)-1
    np(1:n)=node_partition(1:n)
  end subroutine get_np_parallel


  SUBROUTINE construct_para( ng, nb, nk, ns, p )
    implicit none
    integer,intent(IN) :: ng, nb, nk, ns
    type(para),intent(INOUT) :: p
    call construct_dr( ng, p%np(1)*p%np(2)*p%np(3), p%grid )
    call construct_dr( nb, p%np(4)                , p%band )
    call construct_dr( nk, p%np(5)                , p%bzsm )
    call construct_dr( ns, p%np(6)                , p%spin )
  END SUBROUTINE construct_para


  SUBROUTINE construct_dr( nn, np, p )
    implicit none
    integer,intent(IN) :: nn, np
    type(dr),intent(INOUT) :: p
    integer :: i,j
    p%np=np
    allocate( p%id(0:np-1) ) ; p%id=0
    allocate( p%ir(0:np-1) ) ; p%ir=0
    do i=0,nn-1
       j=mod( i, np )
       p%ir(j)=p%ir(j)+1
    end do
    do j=0,np-1
       p%id(j) = sum( p%ir(0:j) ) - p%ir(j)
    end do
  END SUBROUTINE construct_dr


  subroutine start_mpi_parallel
    integer :: ierr,iprovided
    call MPI_Init(ierr)
    !call MPI_Init_thread( MPI_THREAD_MULTIPLE, iprovided, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD,nprocs,ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)
  end subroutine start_mpi_parallel


  subroutine end_mpi_parallel
    integer :: ierr
    call MPI_Finalize(ierr)
  end subroutine end_mpi_parallel


  subroutine read_parallel( flag_exist )
    use io_tools_module, only: IOTools_readIntegerKeyword, IOTools_findKeyword
    implicit none
    logical,intent(out) :: flag_exist
    integer :: itmp(2)
    call write_border( 0, ' read_parallel(start)' )
    node_partition(:)=1
    MB_d=0
    MB_d_nl=0
    itmp(:)=0
    call IOTools_findKeyword( 'PROCS', flag_exist, flag_bcast=.true. )
    call IOTools_readIntegerKeyword( 'PROCS', node_partition )
    call IOTools_readIntegerKeyword( 'MBD', itmp )
    MB_d = itmp(1)
    MB_d_nl = itmp(2)
    if ( MB_d_nl == 0 ) MB_d_nl = MB_d
    if ( myrank == 0 ) then
      write(*,'(1x,"node_partition(1:7)=",7i4)') node_partition(1:7)
      write(*,*) "MB_d   =",MB_d
      write(*,*) "MB_d_nl=",MB_d_nl
    end if
    call write_border( 0, ' read_parallel(end)' )
  end subroutine read_parallel


  subroutine set_np_parallel(Ngrid,Nband,Nbzsm,Nspin)
    use gcd_lcm_pf_module, only: gcd,lcm,prime_factorization
    implicit none
    integer,intent(in) :: Ngrid(0:3),Nspin,Nbzsm
    integer,intent(inout) :: Nband
    integer :: m,n,k,ng(3),i,j,mg(3)
    integer,allocatable :: ifact(:), plst(:)
    logical :: disp_on

    call write_border( 0, " set_np_parallel(start)" )
    call check_disp_switch( disp_on, 0 )

    if ( all(node_partition==1) ) then
      n = nprocs
      m = gcd( Nspin, n ); node_partition(6)=m
      n = n/m
      m = gcd( Nbzsm, n ); node_partition(5)=m

      n = n/m
      allocate( ifact(n) ); ifact=0
      call prime_factorization( n, ifact )

      k = sum( ifact )
      allocate( plst(k) ); plst=0

      k=0
      do j = n, 1, -1
        do i = 1, ifact(j)
          k=k+1
          plst(k)=j
        end do
      end do

      ng(1:3) = Ngrid(1:3)

      do
        mg = ng
        do j = 1, 3
          i = maxloc( ng, 1 )
          do k = 1, size(plst)-1
            if ( plst(k) == 0 ) cycle
            if ( mod(ng(i),plst(k)) == 0 ) then
              ng(i) = ng(i)/plst(k)
              plst(k) = 0
              exit
            end if
          end do !k
        end do !j
        if ( all(mg==ng) ) exit
      end do

      node_partition(1:3) = Ngrid(1:3)/ng(1:3)

      deallocate( plst )
      deallocate( ifact )

      m = product( node_partition(1:3) )

      n = n/m
      if ( mod(Nband,n) == 0 ) then
        node_partition(4) = n
      else
        k = (Nband+n-1)/n * n
        if ( disp_on ) then
          write(*,*) "Nband =",Nband
          write(*,*) "Nband is replaced to ",k
        end if
        Nband = k
        node_partition(4) = n
      end if
    end if
    if ( disp_on ) then
      write(*,'(1x,"node_partition= ",8i6)') node_partition(:)
    end if
    call write_border( 0, " set_np_parallel(end)" )
  end subroutine set_np_parallel


  subroutine init_parallel( Ngrid, Nband, Nbzsm, Nspin )
    implicit none
    integer,intent(in) :: Ngrid(0:3),Nspin,Nbzsm
    integer,intent(inout) :: Nband
    integer :: m,n,ierr,irank,jrank,itags,itagr,nreq
    integer :: istatus(MPI_STATUS_SIZE,123)
    integer :: ip,fp,nc,ns,a1,a2,a3,a4,a5,a6,a7,b1,b2,b3,b4,b5,b6
    integer :: i1,i2,i3,i4,ib,j0,j1,j2,j3,m0,m1,m2,m3,i,j,k,id,ie,ir
    integer :: ML1,ML2,ML3,np1,np2,np3,MLI(2,3),icolor
    integer,allocatable :: mtmp(:),ntmp(:,:),ireq(:)
    integer,allocatable :: map_grid_2_pinfo(:,:,:,:)
    logical :: exist_procs_parameter

    call write_border( 0, " init_parallel(start)" )

    call read_parallel( exist_procs_parameter )
    if ( .not.exist_procs_parameter ) call set_np_parallel( Ngrid,Nband,Nbzsm,Nspin )

    n = product( node_partition )

    if ( n /= nprocs ) then
      write(*,'(1x,"myrank,n,nprocs,node_partition=",10i5)') myrank,n,nprocs,node_partition
      call stop_program('stop@init_parallel')
    end if

    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

    np1 = node_partition(1)
    np2 = node_partition(2)
    np3 = node_partition(3)

    if ( disp_switch_parallel ) then
      write(*,*) "np1,np2,np3=",np1,np2,np3
      write(*,*) "np_grid    =",np1*np2*np3
      write(*,*) "np_band    =",node_partition(4)
      write(*,*) "np_bz      =",node_partition(5)
      write(*,*) "np_spin    =",node_partition(6)
      write(*,*) "np_fkmb    =",node_partition(7)
      write(*,*) "nprocs     =",nprocs
    end if
    if ( np1>ML1 .or. np2>ML2 .or. np3>ML3 .or. &
         node_partition(4)>Nband .or. node_partition(6)>Nspin ) then
      write(*,'(1x,12i6)') node_partition(1:6),ML1,ML2,ML3,Nband,nspin
      call stop_program("stop@init_parallel")
    end if
    if ( node_partition(5)>Nbzsm ) then
      write(*,'(1x,i6)') Nbzsm
      call stop_program("stop@init_parallel")
    end if

! --- class ---

    allocate( id_class(0:nprocs-1,0:max_parallel) )
    id_class=-1

    j=-1
    do a7=0,node_partition(7)-1
    do a6=0,node_partition(6)-1
    do a5=0,node_partition(5)-1
    do a4=0,node_partition(4)-1
      i=-1
      do a3=0,node_partition(3)-1
      do a2=0,node_partition(2)-1
      do a1=0,node_partition(1)-1
        i=i+1
        j=j+1
        id_class(j,0)=i
        id_class(j,1)=a1
        id_class(j,2)=a2
        id_class(j,3)=a3
        id_class(j,4)=a4
        id_class(j,5)=a5
        id_class(j,6)=a6
        id_class(j,7)=a7
      end do ! a1
      end do ! a2
      end do ! a3
    end do ! a4
    end do ! a5
    end do ! a6
    end do ! a7

    ! if ( disp_switch_parallel ) then
    !   write(*,'(1x,5a5)') "irank","grid","band","bz","spin"
    !   do i=0,nprocs-1
    !     write(*,'(1x,6(2x,i3))') i,id_class(i,0),(id_class(i,j),j=4,7)
    !   end do
    ! end if

! --- fock-band parallel ---

    if ( disp_switch_parallel ) write(*,'("--- fock-band-parallel ---")')

    np_fkmb = node_partition(7)

    if ( np_fkmb > Nband ) then
      write(*,*) "np_fkmb is too large!"
      write(*,*) ",np_fkmb, Nband=",np_fkmb,Nband
      call stop_program( "stop@init_parallel" )
    end if

    allocate( id_fkmb(0:np_fkmb-1) ) ; id_fkmb=0
    allocate( ir_fkmb(0:np_fkmb-1) ) ; ir_fkmb=0

    do i=0,Nband-1
      j=mod(i,np_fkmb)
      ir_fkmb(j)=ir_fkmb(j)+1
    end do

    do i=0,np_fkmb-1
      id_fkmb(i)=sum( ir_fkmb(0:i) )-ir_fkmb(i)
    end do

    if ( disp_switch_parallel ) then
      write(*,*) "Nband  =",Nband
      write(*,*) "np_fkmb=",np_fkmb
      ! write(*,'(1x,3a8)') "i","id_fkmb","ir_fkmb"
      ! do i=0,np_fkmb-1
      !   write(*,'(1x,3i8)') i,id_fkmb(i),ir_fkmb(i)
      ! end do
    end if

    nprocs_f = count( id_class(:,7)==id_class(myrank,7) )

    if ( disp_switch_parallel ) write(*,*) "nprocs_f =",nprocs_f

! --- spin parallel ---

    if ( disp_switch_parallel ) write(*,'("--- spin-parallel ---")')

    np_spin = node_partition(6)

    if ( np_spin>nspin ) then
      write(*,*) "np_spin is too large!"
      write(*,*) ",np_spin,nspin=",np_spin,nspin
      call stop_program('np_spin,parallel')
    end if

    allocate( id_spin(0:np_spin-1) ) ; id_spin=0
    allocate( ir_spin(0:np_spin-1) ) ; ir_spin=0

    do i=0,nspin-1
      j=mod(i,np_spin)
      ir_spin(j)=ir_spin(j)+1
    end do

    do i=0,np_spin-1
      id_spin(i)=sum( ir_spin(0:i) )-ir_spin(i)
    end do

    if ( disp_switch_parallel ) then
      write(*,*) "nspin  =",nspin
      write(*,*) "np_spin=",np_spin
      ! write(*,'(1x,3a8)') "i","id_spin","ir_spin"
      ! do i=0,np_spin-1
      !   write(*,'(1x,3i8)') i,id_spin(i),ir_spin(i)
      ! end do
    end if

    nprocs_s = count( id_class(:,6)==id_class(myrank,6) .and. &
                      id_class(:,7)==id_class(myrank,7) )

    if ( disp_switch_parallel ) write(*,*) "nprocs_s =",nprocs_s

! --- k parallel ---

    if ( disp_switch_parallel ) write(*,'("--- k-parallel ---")')

    np_bzsm = node_partition(5)

    allocate( id_bzsm(0:np_bzsm-1) ) ; id_bzsm=0
    allocate( ir_bzsm(0:np_bzsm-1) ) ; ir_bzsm=0

    do i=0,max(Nbzsm,np_bzsm)-1
      j=mod(i,np_bzsm)
      ir_bzsm(j)=ir_bzsm(j)+1
    end do

    do i=0,np_bzsm-1
      id_bzsm(i)=sum( ir_bzsm(0:i) )-ir_bzsm(i)
    end do

    if ( disp_switch_parallel ) then
      write(*,*) "Nbzsm    =",Nbzsm
      write(*,*) "np_bzsm=",np_bzsm
      ! write(*,'(1x,3a8)') "i","id_bzsm","ir_bzsm"
      ! do i=0,np_bzsm-1
      !   write(*,'(1x,3i8)') i,id_bzsm(i),ir_bzsm(i)
      ! end do
    end if

    nprocs_k = count( id_class(:,5)==id_class(myrank,5) .and. &
                      id_class(:,6)==id_class(myrank,6) .and. &
                      id_class(:,7)==id_class(myrank,7) )

    if ( disp_switch_parallel ) write(*,*) "nprocs_k =",nprocs_k

! --- band parallel ---

    if ( disp_switch_parallel ) write(*,'("--- band-parallel ---")')

    np_band = node_partition(4)

    allocate( id_band(0:np_band-1) ) ; id_band=0
    allocate( ir_band(0:np_band-1) ) ; ir_band=0

    do i=0,Nband-1
      j=mod(i,np_band)
      ir_band(j)=ir_band(j)+1
    end do

    do i=0,np_band-1
      id_band(i)=sum( ir_band(0:i) )-ir_band(i)
    end do

    if ( disp_switch_parallel ) then
      write(*,*) "Nband  =",Nband
      write(*,*) "np_band=",np_band
      ! write(*,'(1x,3a8)') "i","id_band","ir_band"
      ! do i=0,np_band-1
      !   write(*,'(1x,3i8)') i,id_band(i),ir_band(i)
      ! end do
    end if

    nprocs_b = count( id_class(:,4)==id_class(myrank,4) .and.  &
                      id_class(:,5)==id_class(myrank,5) .and.  &
                      id_class(:,6)==id_class(myrank,6) .and.  &
                      id_class(:,7)==id_class(myrank,7)       )

    if ( disp_switch_parallel ) write(*,*) "nprocs_b =",nprocs_b

! --- grid parallel ---

    if ( disp_switch_parallel ) write(*,'("--- grid-parallel ---")')

    np_grid = node_partition(1)*node_partition(2)*node_partition(3)
    if ( np_spin*np_bzsm*np_band*np_grid*np_fkmb /= nprocs ) call stop_program("@init_parallel")

    allocate( id_grid(0:np_grid-1)      ) ; id_grid=0
    allocate( ir_grid(0:np_grid-1)      ) ; ir_grid=0
    allocate( pinfo_grid(8,0:np_grid-1) ) ; pinfo_grid=0

    nprocs_g = count( id_class(:,0)==id_class(myrank,0) .and. &
                      id_class(:,4)==id_class(myrank,4) .and. &
                      id_class(:,5)==id_class(myrank,5) .and. &
                      id_class(:,6)==id_class(myrank,6) .and. &
                      id_class(:,7)==id_class(myrank,7)       )

! ---

    if ( disp_switch_parallel ) then
      write(*,*) "nprocs_g =",nprocs_g
      write(*,*) "nprocs_b =",nprocs_b
      write(*,*) "nprocs_k =",nprocs_k
      write(*,*) "nprocs_s =",nprocs_s
      write(*,*) "nprocs_f =",nprocs_f
      write(*,*) "nprocs   =",nprocs
    end if

! --- communicators ---

    m=0
    i=0
    do i4=0,np_fkmb-1
    do i3=0,np_spin-1
    do i2=0,np_bzsm-1
    do i1=0,np_band-1
      i=i+1
      if ( id_class(myrank,4)==i1 .and. &
           id_class(myrank,5)==i2 .and. &
           id_class(myrank,6)==i3 .and. &
           id_class(myrank,7)==i4       ) m=i
    end do
    end do
    end do
    end do
    call mpi_comm_split(mpi_comm_world,m,myrank,comm_grid,ierr)
    call mpi_comm_rank(comm_grid,myrank_g,ierr)
    call mpi_comm_size(comm_grid,nprocs_g,ierr)

    m=0
    i=0
    do i4=0,np_fkmb-1
    do i3=0,np_spin-1
    do i2=0,np_bzsm-1
    do i1=0,np_grid-1
      i=i+1
      if ( id_class(myrank,0)==i1 .and. &
           id_class(myrank,5)==i2 .and. &
           id_class(myrank,6)==i3 .and. &
           id_class(myrank,7)==i4       ) m=i
    end do
    end do
    end do
    end do
    call mpi_comm_split(mpi_comm_world,m,myrank,comm_band,ierr)
    call mpi_comm_rank(comm_band,myrank_b,ierr)
    call mpi_comm_size(comm_band,nprocs_b,ierr)

    m=0
    i=0
    do i4=0,np_fkmb-1
    do i3=0,np_spin-1
    do i2=0,np_band-1
    do i1=0,np_grid-1
      i=i+1
      if ( id_class(myrank,0)==i1 .and. &
           id_class(myrank,4)==i2 .and. &
           id_class(myrank,6)==i3 .and. &
           id_class(myrank,7)==i4 ) m=i
    end do
    end do
    end do
    end do
    call mpi_comm_split(mpi_comm_world,m,myrank,comm_bzsm,ierr)
    call mpi_comm_rank(comm_bzsm,myrank_k,ierr)
    call mpi_comm_size(comm_bzsm,nprocs_k,ierr)

    m=0
    i=0
    do i4=0,np_fkmb-1
    do i3=0,np_bzsm-1
    do i2=0,np_band-1
    do i1=0,np_grid-1
      i=i+1
      if ( id_class(myrank,0)==i1 .and. &
           id_class(myrank,4)==i2 .and. &
           id_class(myrank,5)==i3 .and. &
           id_class(myrank,7)==i4 ) m=i
    end do
    end do
    end do
    end do
    call mpi_comm_split(mpi_comm_world,m,myrank,comm_spin,ierr)
    call mpi_comm_rank(comm_spin,myrank_s,ierr)
    call mpi_comm_size(comm_spin,nprocs_s,ierr)

    m=0
    i=0
    do i4=0,np_spin-1
    do i3=0,np_bzsm-1
    do i2=0,np_band-1
    do i1=0,np_grid-1
      i=i+1
      if ( id_class(myrank,0)==i1 .and. &
           id_class(myrank,4)==i2 .and. &
           id_class(myrank,5)==i3 .and. &
           id_class(myrank,6)==i4 ) m=i
    end do
    end do
    end do
    end do
    call mpi_comm_split(mpi_comm_world,m,myrank,comm_fkmb,ierr)
    call mpi_comm_rank(comm_fkmb,myrank_f,ierr)
    call mpi_comm_size(comm_fkmb,nprocs_f,ierr)

! --- Communicator of band + k + spin

    icolor = myrank_g
    call MPI_Comm_split( MPI_COMM_WORLD, icolor, myrank, comm_bks, ierr )
    call construct_ParaInfo( ParaInfo(6), comm_bks )
    call MPI_Comm_size( comm_bks, nprocs_bks, ierr )

! --- Communicator of grid + band

    icolor = myrank_k + myrank_s*nprocs_k
    call MPI_Comm_split( MPI_COMM_WORLD, icolor, myrank, comm_gb, ierr )
    call MPI_Comm_size( comm_gb, nprocs_gb, ierr )

    if ( disp_switch_parallel ) then
      write(*,*) "comm_world, nprocs     =",MPI_COMM_WORLD,nprocs
      write(*,*) "comm_grid,  nprocs_g   =",comm_grid, nprocs_g
      write(*,*) "comm_band,  nprocs_b   =",comm_band, nprocs_b
      write(*,*) "comm_bzsm,  nprocs_k   =",comm_bzsm, nprocs_k
      write(*,*) "comm_spin,  nprocs_s   =",comm_spin, nprocs_s
      write(*,*) "comm_fkmb,  nprocs_f   =",comm_fkmb, nprocs_f
      write(*,*) "comm_bks ,  nprocs_bks =",comm_bks , nprocs_bks
      write(*,*) "comm_gb  ,  nprocs_gb  =",comm_gb  , nprocs_gb
    end if

    if ( myrank_g /= id_class(myrank,0) ) then
      write(*,*) "myrank,myrank_g=",myrank,myrank_g
      call stop_program('@parallel')
    end if

! --- band bundle ---

    n = min( MB_d,ir_band(id_class(myrank,4)) )
    MB_d = max( n,1 )

    n = min( MB_d_nl,ir_band(id_class(myrank,4)) )
    MB_d_nl = max( n,1 )

! ---

    allocate( idisp(0:nprocs-1) ) ; idisp=-1
    allocate( ircnt(0:nprocs-1) ) ; ircnt=0

    call write_border( 0, " init_parallel(end)" )

  end subroutine init_parallel


  subroutine get_range_parallel( n1, n2, indx )
    implicit none
    integer,intent(out) :: n1,n2
    character(*),intent(in) :: indx
    select case( indx )
    case( 'g','grid' )
       n1=id_grid(myrank_g)+1
       n2=n1+ir_grid(myrank_g)-1
    case( 'b','band' )
       n1=id_band(myrank_b)+1
       n2=n1+ir_band(myrank_b)-1
    case( 'k','bzsm' )
       n1=id_bzsm(myrank_k)+1
       n2=n1+ir_bzsm(myrank_k)-1
    case( 's','spin' )
       n1=id_spin(myrank_s)+1
       n2=n1+ir_spin(myrank_s)-1
    end select
  end subroutine get_range_parallel


  subroutine construct_id_ir_parallel( id, ir, nn, comm_in, n0, n1 )
    implicit none
    integer,allocatable,intent(inout) :: id(:), ir(:)
    integer,intent(in) :: nn
    integer,optional,intent(in) :: comm_in
    integer,optional,intent(out) :: n0, n1
    integer :: np, ierr, i, j, comm, irank
    comm=MPI_COMM_WORLD
    if ( present(comm_in) ) comm=comm_in
    call MPI_Comm_size( comm, np, ierr )
    allocate( id(0:np-1) ); id=0
    allocate( ir(0:np-1) ); ir=0
    do i=0,nn-1
       j=mod( i, np )
       ir(j)=ir(j)+1
    end do
    do j=0,np-1
       id(j) = sum( ir(0:j) ) - ir(j)
    end do
    if ( present(n0) .or. present(n1) ) then
       call MPI_Comm_rank( comm, irank, ierr )
       if ( present(n0) ) n0=id(irank)+1
       if ( present(n1) ) n1=id(irank)+ir(irank)
    end if
  end subroutine construct_id_ir_parallel


  subroutine construct_ParaInfo( pi, comm )
    implicit none
    type(pinfo),intent(inout) :: pi
    integer,intent(in) :: comm
    integer :: ierr
    call MPI_Comm_rank( comm, pi%myrank, ierr )
    call MPI_Comm_size( comm, pi%nprocs, ierr )
    allocate( pi%ir(0:nprocs-1) ); pi%ir=0
    allocate( pi%id(0:nprocs-1) ); pi%id=0
  end subroutine construct_ParaInfo


  subroutine get_ParaInfo( pi, indx )
    implicit none
    type(pinfo),intent(inout) :: pi
    character(*) :: indx
    integer :: np
    select case( indx )
    case( 'bks','BKS' )
      np = ParaInfo(6)%nprocs
      if ( .not.allocated(pi%ir) ) then
        allocate( pi%ir(0:np-1) ); pi%ir=0
      end if
      if ( .not.allocated(pi%id) ) then
        allocate( pi%id(0:np-1) ); pi%id=0
      end if
      pi = ParaInfo(6)
    case default
      write(*,*) "indx=",indx," is not implemented yet"
      call stop_program('@get_ParaInfo')
    end select
  end subroutine get_ParaInfo


  subroutine load_div_parallel( ir, id, ndat )
    implicit none
    integer,intent(out) :: ir(0:),id(0:)
    integer,intent(in) :: ndat
    integer :: np, nres, i, j
    np=size(ir)
    ir=ndat/np
    nres = ndat - sum(ir)
    do i = 1, nres
      j = mod(i-1,np)
      ir(j) = ir(j) + 1
    end do
    id=0
    do i=0,np-1
      id(i) = sum( ir(0:i) ) - ir(i)
    end do
  end subroutine load_div_parallel


  subroutine load_modify_parallel( ng, nb, nk, ns )
    implicit none
    integer,optional,intent(in) :: ng,nb,nk,ns
    if ( present(ng) ) call load_div_parallel( ir_grid, id_grid, ng )
    if ( present(nb) ) call load_div_parallel( ir_band, id_band, nb )
    if ( present(nk) ) call load_div_parallel( ir_bzsm, id_bzsm, nk )
    if ( present(ns) ) call load_div_parallel( ir_spin, id_spin, ns )
  end subroutine load_modify_parallel

end module parallel_module
