module io_module

  use parallel_module
  use rgrid_module, only: Ngrid,Igrid
  use array_bound_module, only: ML,ML_0,ML_1,MSP,MSP_0,MSP_1
  use rgrid_mol_module, only: LL
  use kinetic_module, only: SYStype

  use io1_module, only: read_data_io1
  use io2_module
  use io_read_wf_simple_module, only: read_wf_simple
  use io_write_wf_simple_module, only: write_wf_simple
  use rsdft_mpi_module

  use density_module, only: rho
  use xc_module, only: Vxc
  use localpot_module, only: Vloc
  use hartree_variables, only: Vh

  use aa_module, only: aa
  use atom_module, only: Natom, Nelement, ki_atom, zn_atom, aa_atom
  use grid_module, only: construct_map_1d_to_3d_grid
  use fermi_module, only: efermi
  use io_ctrl_parameters, only: IC, OC, OC2, IO_ctrl, MBwr1, MBwr2, lat_new, lat_old, &
  icount, flag_overwrite_io, read_io_ctrl_parameters, file_wf0,file_wf1,file_wf2, &
  file_vrho0, file_vrho1, file_vrho2
  use watch_module
  use rsdft_mpi_module

  implicit none

  private
  public :: read_io
  public :: write_data
  public :: read_data
  public :: GetParam_IO
  public :: Init_IO
  public :: read_vrho_data

  ! character(30) :: file_wf0   ="wf.dat1"
  ! character(30) :: file_vrho0 ="vrho.dat1"
  ! character(30) :: file_wf1   ="wf.dat1"
  ! character(30) :: file_vrho1 ="vrho.dat1"
  ! character(30) :: file_wf2   ="wf.dat1"
  ! character(30) :: file_vrho2 ="vrho.dat1"

  character(64),parameter :: version="version3.0, comment_length=64"
  character(64) :: comment

contains


  subroutine read_io
    implicit none
    call read_io_ctrl_parameters
  end subroutine read_io


  subroutine write_data( disp_switch, flag, ctrl_wr, suffix )
    implicit none
    logical,intent(in) :: flag, disp_switch
    character(*),optional,intent(in) :: ctrl_wr
    character(*),optional,intent(in) :: suffix
    integer :: i,j,k,n,n1,n2,ierr,i1,i2,i3,j1,j2,j3,isym,ir,ie,id,i0
    integer :: a1,a2,a3,b1,b2,b3,ML0,irank,s,lt(3),jd
    integer :: istatus(MPI_STATUS_SIZE,123)
    integer,allocatable :: irc(:),ids(:),irtmp(:),itmp3(:,:,:)
    integer,allocatable :: idtmp(:),LL2(:,:)
    real(8) :: c,fs,ct0,ct1,et0,et1,mem,memax
    real(8),allocatable :: rtmp(:)
    integer :: ML1,ML2,ML3
    type(time) :: tt
    logical :: flag_wr_pot, flag_wr_wf

    if ( OC2 < 1 .or. OC < 1 .or. OC > 15 ) return

    icount=icount+1
    if ( .not.(flag .or. icount==OC2) ) return

    call write_border( 0, " write_data(start)" )
    call start_timer( t_out=tt )

    n1  = idisp(myrank)+1
    n2  = idisp(myrank)+ircnt(myrank)
    ML0 = ircnt(myrank)
    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

    call construct_gridmap_sub( LL2 )

! ---

    flag_wr_pot = .true.
    flag_wr_wf  = .true.
    if ( present(ctrl_wr) ) then
       flag_wr_pot = ( ctrl_wr=="pot" .or. ctrl_wr=="pot_and_wf" )
       flag_wr_wf  = ( ctrl_wr=="wf"  .or. ctrl_wr=="pot_and_wf" )
    end if

! ---

    if ( present(suffix) .and. .not.flag_overwrite_io ) call Init_IO( suffix )

! ---

!
! --- density and potentials ---
!

    if ( ( OC ==  2 .or. OC ==  3 .or. OC ==  5 .or. &
           OC == 12 .or. OC == 13 ) .and. flag_wr_pot ) then

       allocate( rtmp(ML) ) ; rtmp=0.0d0

       if ( myrank == 0 ) then
          open(2,file=file_vrho1,form='unformatted')
          write(2) version
          write(2) ML,ML1,ML2,ML3
          write(2) LL2(:,:)
          comment="lattice aa(3,3)"
          write(2) comment
          write(2) aa
          comment="atom: Nelement Natom / aa_atom / ki_atom /zn_atom"
          write(2) comment
          write(2) Nelement, Natom
          write(2) aa_atom
          write(2) ki_atom
          write(2) zn_atom
          comment="Fermi energy"
          write(2) comment
          write(2) efermi
       end if

       do s=1,MSP

          do irank=0,nprocs-1
             i1=id_band(id_class(irank,4))+1
             j1=id_band(id_class(irank,4))+ir_band(id_class(irank,4))
             i2=id_bzsm(id_class(irank,5))+1
             j2=id_bzsm(id_class(irank,5))+ir_bzsm(id_class(irank,5))
             i3=id_spin(id_class(irank,6))+1
             j3=id_spin(id_class(irank,6))+ir_spin(id_class(irank,6))
             if ( id_grid(id_class(irank,0))==0 .and. &
                  i1<=1 .and. 1<=j1 .and. &
                  i2<=1 .and. 1<=j2 .and. &
                  i3<=s .and. s<=j3 ) exit
          end do
          if ( irank>=nprocs ) then
             write(*,*) "ERROR(write_data)",myrank
             stop
          end if

          if ( MSP_0<=s .and. s<=MSP_1 ) then
             call mpi_gatherv(rho(n1,s),ML0,mpi_real8,rtmp, &
                  ir_grid,id_grid,mpi_real8,0,comm_grid,ierr)
          end if

          call mpi_barrier(mpi_comm_world,ierr)

          if ( irank /= 0 .and. myrank_f == 0 ) then
             if ( irank==myrank ) then
                call mpi_send(rtmp,ML,mpi_real8,0,0,mpi_comm_world,ierr)
             end if
             if ( myrank==0 ) then
                call mpi_recv(rtmp,ML,mpi_real8,irank,0, &
                     mpi_comm_world,istatus,ierr)
             end if
          end if
          if ( myrank==0 ) then
             write(2) rtmp(:)
          end if

          if ( MSP_0<=s .and. s<=MSP_1 ) then
             call mpi_gatherv(Vloc(n1,s),ML0,mpi_real8,rtmp, &
                  ir_grid,id_grid,mpi_real8,0,comm_grid,ierr)
          end if

          call mpi_barrier(mpi_comm_world,ierr)

          if ( irank /= 0 .and. myrank_f == 0 ) then
             if ( irank==myrank ) then
                call mpi_send(rtmp,ML,mpi_real8,0,0,mpi_comm_world,ierr)
             end if
             if ( myrank==0 ) then
                call mpi_recv(rtmp,ML,mpi_real8,irank,0, &
                     & mpi_comm_world,istatus,ierr)
             end if
          end if
          if ( myrank==0 ) then
             write(2) rtmp(:)
          end if

          if ( s==1 ) then

             call mpi_gatherv(Vh(n1),ML0,mpi_real8,rtmp,ir_grid, &
                  id_grid,mpi_real8,0,comm_grid,ierr)

             if ( myrank==0 ) then
                write(2) rtmp(:)
             end if

          end if

          if ( MSP_0<=s .and. s<=MSP_1 ) then
             call mpi_gatherv(Vxc(n1,s),ML0,mpi_real8,rtmp, &
                  & ir_grid,id_grid,mpi_real8,0,comm_grid,ierr)
          end if

          call mpi_barrier(mpi_comm_world,ierr)

          if ( irank /= 0 .and. myrank_f == 0 ) then
             if ( irank==myrank ) then
                call mpi_send(rtmp,ML,mpi_real8,0,0,mpi_comm_world, &
                     & ierr)
             end if
             if ( myrank==0 ) then
                call mpi_recv(rtmp,ML,mpi_real8,irank,0, &
                     & mpi_comm_world,istatus,ierr)
             end if
          end if
          if ( myrank==0 ) then
             write(2) rtmp(:)
          end if

       end do ! s

       if ( myrank == 0 ) then
          close(2)
       end if

       if (DISP_SWITCH) then
          write(*,*) "write to ",file_vrho1
       end if

       deallocate( rtmp )

    end if

!
! --- Wave function ---
!
    

    if ( ( OC ==  1 .or. OC ==  3 .or. OC ==  4 .or. OC ==  5 .or. &
           OC == 11 .or. OC == 13 ) .and. flag_wr_wf ) then

       call write_wf_simple &
            ( file_wf1, IO_ctrl, OC, SYStype, MBwr1, MBwr2, disp_switch )

    end if

    deallocate( LL2 )

    icount=0

    call Init_IO( "" )  ! remove suffix from the I/O file names

    call result_timer( "write_data", tt )
    call write_border( 0, " write_data(end)" )

    return
  end subroutine write_data

!--------1---------2---------3---------4---------5---------6---------7--

  subroutine read_data(disp_switch)
    use var_sys_parameter, only: use_real8_wf
    use wf_module, only: unk, occ
    implicit none
    logical,intent(in) :: disp_switch
    integer :: k,n,i,j,ML_tmp,n1,n2,ML0,irank
    integer :: ML1_tmp,ML2_tmp,ML3_tmp,ierr,i1,i2,i3,j1,j2,j3,s,i0
    integer :: ML1,ML2,ML3,mx,my,mz,itmp(7)
    real(8) :: fs,mem,memax,ct0,et0,ct1,et1
    real(8),allocatable :: rtmp(:),rtmp3(:,:,:)
    character(len=64) :: cbuf
    logical :: flag_versioned
    type(time) :: tt

    if ( IC <= 0 ) return

    call write_border( 0, " read_data(start)" )
    call start_timer( t_out=tt )

    n1    = idisp(myrank)+1
    n2    = idisp(myrank)+ircnt(myrank)
    ML0   = ircnt(myrank)

    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

    allocate( lat_old(3,ML) ) ; lat_old=0

    call construct_gridmap_sub( lat_new )

!
! --- Read VRHO ---
!
    if ( IC==2 .or. IC==3 .or. IC==4 .or. IC==5 .or. IC==12 .or. IC==13 ) then

       if ( myrank==0 ) then
          open(80,file=file_vrho2,form='unformatted')
          read(80) cbuf(1:7)
          rewind(80)
          if ( cbuf(1:7) == "version" ) then
             read(80) cbuf
             write(*,*) "file format version: "//cbuf
             flag_versioned=.true.
          else
             flag_versioned=.false.
          end if
          read(80) ML_tmp,ML1_tmp,ML2_tmp,ML3_tmp
          itmp=0
          itmp(1)=ML_tmp
          itmp(2)=ML1_tmp
          itmp(3)=ML2_tmp
          itmp(4)=ML3_tmp
       end if
       call mpi_bcast(itmp,4,mpi_integer,0,mpi_comm_world,ierr)
       ML_tmp=itmp(1)
       ML1_tmp=itmp(2)
       ML2_tmp=itmp(3)
       ML3_tmp=itmp(4)

       if ( ML_tmp/=ML ) then
          write(*,*) "ML,ML_tmp =",ML,ML_tmp
          call stop_program("stop@read_data(1)")
       end if
       if ( ML1_tmp/=ML1 .or. ML2_tmp/=ML2 .or. ML3_tmp/=ML3 ) then
          if (DISP_SWITCH) then
             write(*,*) "ML1_tmp,ML2_tmp,ML3_tmp =", &
                  & ML1_tmp,ML2_tmp,ML3_tmp
             write(*,*) "ML1,ML2,ML3 =",ML1,ML2,ML3
          end if
          call stop_program("stop@read_data(2)")
       end if

       if ( myrank==0 ) then
          read(80) lat_old(:,:)
          if ( flag_versioned ) then
             read(80) cbuf
             read(80)
             read(80) cbuf
             read(80)
             read(80)
             read(80)
             read(80)
             read(80) cbuf
             read(80) efermi
          end if
       end if
       call mpi_bcast(lat_old,3*ML,mpi_integer,0,mpi_comm_world,ierr)
       call mpi_bcast(efermi,1,mpi_real8,0,mpi_comm_world,ierr)

       i=sum(abs(lat_old(:,:)-lat_new(:,:)))
       if ( i/=0 ) then
          if (DISP_SWITCH) then
             write(*,*) "LL and lat_old is different"
          end if
       end if

       if ( SYStype == 0 ) then

          allocate( rtmp3(0:ML1-1,0:ML2-1,0:ML3-1) )

       else if ( SYStype == 1 ) then

          mx = (Ngrid(1)-1)/2
          my = (Ngrid(2)-1)/2
          mz = (Ngrid(3)-1)/2
          allocate( rtmp3(-mx:mx,-my:my,-mz:mz) )

       end if

       allocate( rtmp(ML) )

       rtmp(:)=0.d0
       rtmp3(:,:,:)=0.d0

       do s=1,MSP

          if ( myrank==0 ) then
             read(80) rtmp(:)
          end if
          call mpi_bcast(rtmp,ML,mpi_real8,0,mpi_comm_world,ierr)
          do i=1,ML
             i1=lat_old(1,i) ; i2=lat_old(2,i) ; i3=lat_old(3,i)
             rtmp3(i1,i2,i3)=rtmp(i)
          end do

          do i=n1,n2
             i1=lat_new(1,i) ; i2=lat_new(2,i) ; i3=lat_new(3,i)
             rho(i,s)=rtmp3(i1,i2,i3)
          end do

          if ( myrank==0 ) then
             read(80) rtmp(:)
          end if
          call mpi_bcast(rtmp,ML,mpi_real8,0,mpi_comm_world,ierr)
          do i=1,ML
             i1=lat_old(1,i) ; i2=lat_old(2,i) ; i3=lat_old(3,i)
             rtmp3(i1,i2,i3)=rtmp(i)
          end do
          if ( MSP_0<=s .and. s<=MSP_1 ) then
             do i=n1,n2
                i1=lat_new(1,i) ; i2=lat_new(2,i) ; i3=lat_new(3,i)
                Vloc(i,s)=rtmp3(i1,i2,i3)
             end do
          end if

          if ( s==1 ) then
             if ( myrank==0 ) then
                read(80) rtmp(:)
             end if
             call mpi_bcast(rtmp,ML,mpi_real8,0,mpi_comm_world,ierr)
             do i=1,ML
                i1=lat_old(1,i) ; i2=lat_old(2,i) ; i3=lat_old(3,i)
                rtmp3(i1,i2,i3)=rtmp(i)
             end do
             do i=n1,n2
                i1=lat_new(1,i) ; i2=lat_new(2,i) ; i3=lat_new(3,i)
                Vh(i)=rtmp3(i1,i2,i3)
             end do
          end if

          if ( myrank==0 ) then
             read(80) rtmp(:)
          end if
          call mpi_bcast(rtmp,ML,mpi_real8,0,mpi_comm_world,ierr)
          do i=1,ML
             i1=lat_old(1,i) ; i2=lat_old(2,i) ; i3=lat_old(3,i)
             rtmp3(i1,i2,i3)=rtmp(i)
          end do
          if ( MSP_0<=s .and. s<=MSP_1 ) then
             do i=n1,n2
                i1=lat_new(1,i) ; i2=lat_new(2,i) ; i3=lat_new(3,i)
                Vxc(i,s)=rtmp3(i1,i2,i3)
             end do
          end if

       end do ! s

       deallocate( rtmp,rtmp3 )

       if ( myrank==0 ) then
          close(80)
       end if

       if (DISP_SWITCH) then
          write(*,*) "read from ",file_vrho2
       end if

    end if ! IC = 2,3, 4,5, 12,13

!
! --- Read WF ---
!

    if ( IC == 2 .or. IC == 4 .or. IC == 12 ) return

    if ( IO_ctrl == 2 .or. IO_ctrl == 3 ) then
       IO_ctrl=3
       call read_io2( SYStype, ierr )
       if ( ierr == 0 ) IO_ctrl=2
    end if

    select case( IO_ctrl )
    case( 0, 3 )
      if ( use_real8_wf() ) then
#ifdef _DRSDFT_
        call read_wf_simple( d_wf_out=unk, occ_out=occ )
#endif
      else
#ifndef _DRSDFT_
        call read_wf_simple( z_wf_out=unk, occ_out=occ )
#endif
      end if
    case( 1 )
       call read_data_io1( file_wf2, SYStype )
    case( 2 )
       call read_data_io2
    end select

    call result_timer( "read_data", tt )
    call write_border( 0, " read_data(end)" )

    return

  end subroutine read_data


  FUNCTION GetParam_IO(i)
    implicit none
    integer :: GetParam_IO
    integer,intent(IN) :: i
    select case(i)
    case( 1 )
       GetParam_IO = IC
    case( 2 )
       GetParam_IO = OC
    case( 3 )
       GetParam_IO = OC2
    case default
       GetParam_IO = -1
    end select
  END FUNCTION GetParam_IO


  SUBROUTINE Init_IO( index )
    implicit none
    character(*),intent(IN) :: index
    if ( index /= "" ) then
       file_wf1   = trim(file_wf0)//"."//index
       file_vrho1 = trim(file_vrho0)//"."//index
    else
       file_wf1   = file_wf0
       file_vrho1 = file_vrho0
    end if
  END SUBROUTINE Init_IO


  subroutine construct_gridmap_sub( LL2 )
    implicit none
    integer,allocatable,intent(inout) :: LL2(:,:)
    integer :: ierr
    include 'mpif.h'

    if ( SYStype == 0 ) then

       call construct_map_1d_to_3d_grid( LL2 )

    else if ( SYStype == 1 ) then

       allocate( LL2(3,Ngrid(0)) ) ; LL2=0

       LL2(1:3,Igrid(1,0):Igrid(2,0)) = LL(1:3,Igrid(1,0):Igrid(2,0))

       call rsdft_allreduce_sum( LL2, comm_grid )

    end if

  end subroutine construct_gridmap_sub


  subroutine read_vrho_data( file_in, Vh_out, Ves_out )
    implicit none
    character(30),intent(in) :: file_in
    real(8),optional,intent(out) :: Vh_out(:)
    real(8),optional,intent(out) :: Ves_out(:)
    integer :: k,n,i,j,ML_tmp,n1,n2,ML0,irank
    integer :: ML1_tmp,ML2_tmp,ML3_tmp,ierr,i1,i2,i3,j1,j2,j3,s,i0
    integer :: ML1,ML2,ML3,mx,my,mz,itmp(7)
    real(8) :: fs,mem,memax,ct0,et0,ct1,et1
    real(8),allocatable :: rtmp(:),rtmp3(:,:,:)
    character(len=64) :: cbuf
    logical :: flag_versioned
    type(time) :: tt
    logical :: disp_switch

    call write_border( 0, " read_vrho_data(start)" )
    call start_timer( t_out=tt )
    call check_disp_switch( disp_switch, 0 )

    n1  = idisp(myrank)+1
    n2  = idisp(myrank)+ircnt(myrank)
    ML0 = ircnt(myrank)

    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

    if ( allocated(lat_old) ) deallocate( lat_old )
    allocate( lat_old(3,ML) ); lat_old=0

    call construct_gridmap_sub( lat_new )

!
! --- Read VRHO ---
!
    if ( myrank == 0 ) then
       open(80,file=file_in,form='unformatted')
       read(80) cbuf(1:7)
       rewind(80)
       if ( cbuf(1:7) == "version" ) then
          read(80) cbuf
          write(*,*) "file format version: "//cbuf
          flag_versioned=.true.
       else
          flag_versioned=.false.
       end if
       read(80) ML_tmp,ML1_tmp,ML2_tmp,ML3_tmp
       itmp=0
       itmp(1)=ML_tmp
       itmp(2)=ML1_tmp
       itmp(3)=ML2_tmp
       itmp(4)=ML3_tmp
    end if
    call mpi_bcast(itmp,4,mpi_integer,0,mpi_comm_world,ierr)
    ML_tmp=itmp(1)
    ML1_tmp=itmp(2)
    ML2_tmp=itmp(3)
    ML3_tmp=itmp(4)

    if ( ML_tmp/=ML ) then
       write(*,*) "ML,ML_tmp =",ML,ML_tmp
       stop
    end if
    if ( ML1_tmp/=ML1 .or. ML2_tmp/=ML2 .or. ML3_tmp/=ML3 ) then
       if ( disp_switch ) then
          write(*,*) "ML1_tmp,ML2_tmp,ML3_tmp =",ML1_tmp,ML2_tmp,ML3_tmp
          write(*,*) "ML1,ML2,ML3 =",ML1,ML2,ML3
       end if
       stop
    end if

    if ( myrank == 0 ) then
       read(80) lat_old(:,:)
       if ( flag_versioned ) then
          read(80) cbuf
          read(80)
          read(80) cbuf
          read(80)
          read(80)
          read(80)
          read(80)
          read(80) cbuf
          read(80) efermi
       end if
    end if
    call mpi_bcast(lat_old,3*ML,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(efermi,1,mpi_real8,0,mpi_comm_world,ierr)

    i=sum(abs(lat_old(:,:)-lat_new(:,:)))
    if ( i/=0 ) then
       if ( disp_switch ) then
          write(*,*) "LL and lat_old is different"
       end if
    end if

    if ( SYStype == 0 ) then

       allocate( rtmp3(0:ML1-1,0:ML2-1,0:ML3-1) )

    else if ( SYStype == 1 ) then

       mx = (Ngrid(1)-1)/2
       my = (Ngrid(2)-1)/2
       mz = (Ngrid(3)-1)/2
       allocate( rtmp3(-mx:mx,-my:my,-mz:mz) )

    end if

    allocate( rtmp(ML) )

    rtmp(:)=0.d0
    rtmp3(:,:,:)=0.d0

    do s=1,1 !MSP

       if ( myrank==0 ) then
          read(80) rtmp(:)
       end if
       !call mpi_bcast(rtmp,ML,mpi_real8,0,mpi_comm_world,ierr)
       !do i=1,ML
       !   i1=lat_old(1,i) ; i2=lat_old(2,i) ; i3=lat_old(3,i)
       !   rtmp3(i1,i2,i3)=rtmp(i)
       !end do
       !do i=n1,n2
       !   i1=lat_new(1,i) ; i2=lat_new(2,i) ; i3=lat_new(3,i)
       !   rho(i,s)=rtmp3(i1,i2,i3)
       !end do

       if ( myrank == 0 ) then
          read(80) rtmp(:)
       end if
       call mpi_bcast(rtmp,ML,mpi_real8,0,mpi_comm_world,ierr)
       do i=1,ML
          i1=lat_old(1,i) ; i2=lat_old(2,i) ; i3=lat_old(3,i)
          rtmp3(i1,i2,i3)=rtmp(i)
       end do
       if ( present(Ves_out) ) then
          do i=n1,n2
             i1=lat_new(1,i) ; i2=lat_new(2,i) ; i3=lat_new(3,i)
             Ves_out(i-n1+1)=rtmp3(i1,i2,i3)
          end do
       end if

       if ( s == 1 ) then
          if ( myrank == 0 ) then
             read(80) rtmp(:)
          end if
          call mpi_bcast(rtmp,ML,mpi_real8,0,mpi_comm_world,ierr)
          do i=1,ML
             i1=lat_old(1,i) ; i2=lat_old(2,i) ; i3=lat_old(3,i)
             rtmp3(i1,i2,i3)=rtmp(i)
          end do
          if ( present(Vh_out) ) then
             do i=n1,n2
                i1=lat_new(1,i) ; i2=lat_new(2,i) ; i3=lat_new(3,i)
                Vh_out(i-n1+1)=rtmp3(i1,i2,i3)
             end do
          end if
       end if

       if ( myrank == 0 ) then
          read(80) rtmp(:)
       end if
       call mpi_bcast(rtmp,ML,mpi_real8,0,mpi_comm_world,ierr)
       do i=1,ML
          i1=lat_old(1,i) ; i2=lat_old(2,i) ; i3=lat_old(3,i)
          rtmp3(i1,i2,i3)=rtmp(i)
       end do
       if ( present(Ves_out) ) then
          do i=n1,n2
             i1=lat_new(1,i) ; i2=lat_new(2,i) ; i3=lat_new(3,i)
             Ves_out(i-n1+1)=Ves_out(i-n1+1)-rtmp3(i1,i2,i3)
          end do
       end if

    end do ! s

    deallocate( rtmp,rtmp3 )

    if ( myrank==0 ) then
       close(80)
    end if

    call write_border( 0, "read_vrho_data(end)" )

  end subroutine read_vrho_data


end module io_module
