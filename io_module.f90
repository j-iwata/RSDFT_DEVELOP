MODULE io_module

  use rgrid_module, only: Ngrid,Igrid
  use wf_module, only: unk,occ
  use density_module, only: rho
  use xc_module, only: Vxc
  use localpot_module, only: Vloc
  use hartree_variables, only: Vh
  use parallel_module
  use array_bound_module, only: ML,ML_0,ML_1,MB,MB_0,MB_1 &
                               ,MBZ,MBZ_0,MBZ_1,MSP,MSP_0,MSP_1

  use rgrid_mol_module, only: LL
  use kinetic_module, only: SYStype

  use io2_module
  use io_read_module
  use io_write_module

  implicit none

  PRIVATE
  PUBLIC :: read_io, write_data, read_data, read_oldformat_io &
           ,GetParam_IO, Init_IO

  integer :: IO_ctrl=0
  integer :: IC,OC,OC2
  integer :: MBwr1=0
  integer :: MBwr2=0
  character(30) :: file_wf0   ="wf.dat1"
  character(30) :: file_vrho0 ="vrho.dat1"
  character(30) :: file_wf1   ="wf.dat1"
  character(30) :: file_vrho1 ="vrho.dat1"
  character(30) :: file_wf2   ="wf.dat1"
  character(30) :: file_vrho2 ="vrho.dat1"

#ifdef _DRSDFT_
  integer,parameter :: TYPE_MAIN=MPI_REAL8
  real(8),allocatable :: utmp(:)
  real(8),allocatable :: utmp3(:,:,:)
  real(4),allocatable :: utmpSP(:)
  real(8),parameter :: zero=0.d0
  real(4),parameter :: zeroSP=0.0
#else
  integer,parameter :: TYPE_MAIN=MPI_COMPLEX16
  complex(8),allocatable :: utmp(:)
  complex(8),allocatable :: utmp3(:,:,:)
  complex(4),allocatable :: utmpSP(:)
  complex(8),parameter :: zero=(0.d0,0.d0)
  complex(4),parameter :: zeroSP=(0.0,0.0)
#endif

  integer,save :: icount=0

CONTAINS

  SUBROUTINE read_io(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i
    character(6) :: cbuf,ckey
    IC  = 0
    OC  = 0
    OC2 = 100
    IO_ctrl = 0
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:2) == "IC" ) then
             backspace(unit)
             read(unit,*) cbuf,IC
          else if ( ckey(1:3) == "OC2" ) then
             backspace(unit)
             read(unit,*) cbuf,OC2
          else if ( ckey(1:3) == "OC" ) then
             backspace(unit)
             read(unit,*) cbuf,OC
          else if ( ckey(1:6) == "IOCTRL" ) then
             backspace(unit)
             read(unit,*) cbuf,IO_ctrl
          else if ( ckey(1:4) == "MBWR" ) then
             backspace(unit)
             read(unit,*) cbuf,MBwr1,MBwr2
          end if
       end do
999    continue
       write(*,*) "IC =",IC
       write(*,*) "OC =",OC
       write(*,*) "OC2=",OC2
       write(*,*) "IO_ctrl=",IO_ctrl
       write(*,*) "MBwr1,MBwr2=",MBwr1,MBwr2
    end if
    call send_io(0)
  END SUBROUTINE read_io

  SUBROUTINE read_oldformat_io(rank,unit)
    integer,intent(IN) :: rank,unit
    if ( rank == 0 ) then
       read(unit,*) IC,OC2,OC
       write(*,*) "IC,OC,OC2=",IC,OC,OC2
    end if
    call send_io(0)
  END SUBROUTINE read_oldformat_io

  SUBROUTINE send_io(rank)
    integer,intent(IN) :: rank
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(IC ,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(OC ,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(OC2,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(IO_ctrl,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(MBwr1,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(MBwr2,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_io

  SUBROUTINE write_data(disp_switch,flag)
    logical,intent(IN) :: flag,disp_switch
    integer :: i,j,k,n,n1,n2,ierr,i1,i2,i3,j1,j2,j3,isym,ir,ie,id,i0
    integer :: a1,a2,a3,b1,b2,b3,ML0,irank,s,lt(3),jd
    integer :: istatus(MPI_STATUS_SIZE,123)
    integer,allocatable :: irc(:),ids(:),irtmp(:),itmp3(:,:,:)
    integer,allocatable :: idtmp(:),LL2(:,:)
    real(8) :: c,fs,ct0,ct1,et0,et1,mem,memax
    real(8),allocatable :: rtmp(:)
    logical :: flag_related
    integer :: ML1,ML2,ML3
    character(len=5) :: cmyrank
    character(len=32) :: file_wf_split

    if ( OC2 < 1 .or. OC < 1 .or. OC > 15 ) return

    icount=icount+1
    if ( .not.(flag .or. icount==OC2) ) return

    if ( DISP_SWITCH ) then
       write(*,*)
       write(*,'(a60," write_data")') repeat("-",60)
    end if

    n1  = idisp(myrank)+1
    n2  = idisp(myrank)+ircnt(myrank)
    ML0 = ircnt(myrank)

    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

    if ( IO_ctrl==3 ) then
       write(cmyrank,'(i5.5)') myrank
       file_wf_split = trim(file_wf1)//"."//trim(adjustl(cmyrank))
    end if

    if ( MBwr1<1 .or. MBwr2<MBwr1 .or. MB<MBwr1 ) MBwr1=1
    if ( MBwr2<1 .or. MBwr2<MBwr1 .or. MB<MBwr2 ) MBwr2=MB

    allocate( LL2(3,ML) ) ; LL2=0

    if ( SYStype == 0 ) then

       i=n1-1
       do i3=Igrid(1,3),Igrid(2,3)
       do i2=Igrid(1,2),Igrid(2,2)
       do i1=Igrid(1,1),Igrid(2,1)
          i=i+1
          LL2(1,i)=i1
          LL2(2,i)=i2
          LL2(3,i)=i3
       end do
       end do
       end do

    else if ( SYStype == 1 ) then

       LL2(1:3,n1:n2) = LL(1:3,n1:n2)

    end if

    allocate( irc(0:np_grid-1),ids(0:np_grid-1) )
    irc(0:np_grid-1)=3*ir_grid(0:np_grid-1)
    ids(0:np_grid-1)=3*id_grid(0:np_grid-1)

    call mpi_allgatherv(LL2(1,n1),irc(myrank_g),mpi_integer, &
         LL2,irc,ids,mpi_integer,comm_grid,ierr)

    deallocate( ids, irc )

!
! --- density and potentials ---
!

    if ( OC ==  2 .or. OC ==  3 .or. OC ==  5 .or. &
         OC == 12 .or. OC == 13 .or. OC == 15 ) then

       allocate( rtmp(ML) )

       if ( myrank==0 ) then
          open(2,file=file_vrho1,form='unformatted')
          write(2) ML,ML1,ML2,ML3
          write(2) LL2(:,:)
       end if

       do s=1,MSP

          do irank=0,nprocs-1
             i1=id_band(id_class(irank,4))+1
             j1=id_band(id_class(irank,4))+ir_band(id_class(irank,4))
             i2=id_bzsm(id_class(irank,5))+1
             j2=id_bzsm(id_class(irank,5))+ir_bzsm(id_class(irank,5))
             i3=id_spin(id_class(irank,6))+1
             j3=id_spin(id_class(irank,6))+ir_spin(id_class(irank,6))
             if ( id_grid(id_class(irank,0))==0 .and. i1<=1 .and. &
                  1<=j1 .and. i2<=1 .and. 1<=j2 .and. i3<=s &
                  .and. s<=j3 ) exit
          end do
          if ( irank>=nprocs ) then
             write(*,*) "ERROR(read_data)",myrank
             stop
          end if

          if ( MSP_0<=s .and. s<=MSP_1 ) then
             call mpi_gatherv(rho(n1,s),ML0,mpi_real8,rtmp, &
                  ir_grid,id_grid,mpi_real8,0,comm_grid,ierr)
          end if

          call mpi_barrier(mpi_comm_world,ierr)

          if ( irank/=0 ) then
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

          if ( irank/=0 ) then
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

          if ( irank/=0 ) then
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

       if ( myrank==0 ) then
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
    

    if ( OC ==  1 .or. OC ==  3 .or. OC ==  4 .or. OC ==  5 .or. &
         OC == 11 .or. OC == 13 .or. OC == 14 .or. OC == 15 ) then

       if ( OC == 4 ) OC=14
       if ( OC == 5 ) OC=15

       call simple_wf_io_write &
            ( file_wf1, IO_ctrl, OC, SYStype, MBwr1, MBwr2, disp_switch )

goto 100
       if (DISP_SWITCH) then
          select case(OC)
#ifdef _DRSDFT_
          case default
             print*,'WF: real(8) -> real(8)'
          case(4,5)
             print*,'WF: real(8) -> real(4)'
          end select
#else
          case default
             print*,'WF: complex(8) -> complex(8)'
          case(4,5)
             print*,'WF: complex(8) -> complex(4)'
          end select
#endif
       end if

       allocate( utmp(ML) ) ; utmp=zero
       if ( OC==4 .or. OC==5 ) then
          allocate( utmpSP(ML) )
          utmpSP=zero
       end if

       if ( myrank==0 ) then
          open(1,file=file_wf1,form="unformatted")
          write(1) ML,ML1,ML2,ML3
          if ( OC > 10 ) then
             write(1) MB,MBwr1,MBwr2,MBZ,MSP
          else
             write(1) MB,MBwr1,MBwr2
          end if
          write(1) LL2(:,:)
          write(1) occ(:,:,:)
       end if
       if ( IO_ctrl==3 ) then
          if ( myrank==0 ) then
             close(1)
          end if
          open(1,file=file_wf_split,form="unformatted")
       end if

       do s=1,MSP
       do k=1,MBZ
       do n=MBwr1,MBwr2

          do irank=0,nprocs-1
             i1=id_band(id_class(irank,4))+1
             j1=id_band(id_class(irank,4))+ir_band(id_class(irank,4))
             i2=id_bzsm(id_class(irank,5))+1
             j2=id_bzsm(id_class(irank,5))+ir_bzsm(id_class(irank,5))
             i3=id_spin(id_class(irank,6))+1
             j3=id_spin(id_class(irank,6))+ir_spin(id_class(irank,6))
             if ( id_grid(id_class(irank,0))==0 .and. i1<=n .and. &
                  n<=j1 .and. i2<=k .and. k<=j2 .and. i3<=s .and. &
                  s<=j3 ) exit
          end do

          if ( irank>=nprocs ) then
             write(*,*) "ERROR(write_data)",myrank
             stop
          end if

          flag_related = .false.
          if ( MBZ_0<=k .and. k<=MBZ_1 .and. MB_0 <=n .and. n<=MB_1 &
               & .and. MSP_0<=s .and. s<=MSP_1 ) then
             flag_related=.true.
          end if

          if ( IO_ctrl==0 ) then

             if ( flag_related ) then
                call mpi_gatherv(unk(n1,n,k,s),ML0,TYPE_MAIN, &
                     utmp,ir_grid,id_grid,TYPE_MAIN,0, &
                     comm_grid,ierr)
             end if

             call mpi_barrier(mpi_comm_world,ierr)

             if ( irank/=0 ) then
                if ( irank==myrank ) then
                   call mpi_send(utmp,ML,TYPE_MAIN,0,0, &
                        & mpi_comm_world,ierr)
                end if
                if ( myrank==0 ) then
                   call mpi_recv(utmp,ML,TYPE_MAIN,irank,0, &
                        & mpi_comm_world,istatus,ierr)
                end if
             end if
             if ( myrank==0 ) then
#ifdef _DRSDFT_
                select case(OC)
                case default
                   write(1) utmp(:)
                case(4,5)
                   utmpSP(:)=utmp(:)
                   write(1) utmpSP(:)
                end select
#else
                select case(OC)
                case default
                   write(1) utmp(:)
                case(4,5)
                   utmpSP(:)=utmp(:)
                   write(1) utmpSP(:)
                end select
#endif
             end if

          else !( IO_ctrl /= 0 )

             if ( flag_related ) then
#ifdef _DRSDFT_
                select case(OC)
                case default
                   write(1) unk(n1:n2,n,k,s)
                case(4,5)
                   utmpSP(n1:n2)=unk(n1:n2,n,k,s)
                   write(1) utmpSP(n1:n2)
                end select
#else
                select case(OC)
                case default
                   write(1) unk(n1:n2,n,k,s)
                case(4,5)
                   utmpSP(n1:n2)=unk(n1:n2,n,k,s)
                   write(1) utmpSP(n1:n2)
                end select
#endif
             end if

          end if

       end do ! n
       end do ! k
       end do ! s

       if ( IO_ctrl==0 ) then

          if ( myrank==0 ) then
             close(1)
          end if

       else

          close(1)

       end if

       if (DISP_SWITCH) then
          write(*,*) "write to ",file_wf1
          write(*,*) "MBwr1,MBwr2 =",MBwr1,MBwr2
       end if

       if ( allocated(utmpSP) ) deallocate( utmpSP )
       deallocate( utmp )
100 continue

    end if

    deallocate( LL2 )

    icount=0

    return
  END SUBROUTINE write_data

!--------1---------2---------3---------4---------5---------6---------7--
!1:wf, 2:vrho, 3:wf&vrho, 4:wf(Real->Comp), 5:wf(Real->Comp)&vrho

  SUBROUTINE read_data(disp_switch)
    logical,intent(IN) :: disp_switch
    integer :: k,n,i,j,ML_tmp,MB_tmp,MB1_tmp,MB2_tmp,n1,n2,ML0,irank
    integer :: ML1_tmp,ML2_tmp,ML3_tmp,ierr,i1,i2,i3,j1,j2,j3,s,i0
    integer :: itmp(7),a1,a2,a3,b1,b2,b3,istatus(MPI_STATUS_SIZE,123)
    integer :: lt(3),isym,ML1,ML2,ML3,mx,my,mz
    integer,allocatable :: LL_tmp(:,:),ir(:),id(:),itmp3(:,:,:)
    integer,allocatable :: idtmp(:),irtmp(:),LL2(:,:)
    real(8) :: fs,mem,memax,ct0,et0,ct1,et1
    real(8),allocatable :: rtmp(:),rtmp3(:,:,:)
    logical :: flag_related
    integer(kind=4) :: int4
    logical :: flag_SP,flag_R2C
    real(kind=4),allocatable :: rtmpSP(:)
    real(kind=8),allocatable :: rtmpDP(:)
    character(len=5) :: cmyrank
    character(len=32) :: file_wf_split

    if ( IC <= 0 ) return

    n1    = idisp(myrank)+1
    n2    = idisp(myrank)+ircnt(myrank)
    ML0   = ircnt(myrank)

    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

    if ( IO_ctrl==3 ) then
       write(cmyrank,'(i5.5)') myrank
       file_wf_split = trim(file_wf2)//"."//trim(adjustl(cmyrank))
    end if

    if (DISP_SWITCH) then
       write(*,'(a60," read_data")') repeat("-",60)
    end if

    allocate( LL2(3,ML)    ) ; LL2=0
    allocate( LL_tmp(3,ML) ) ; LL_tmp=0

    if ( TYPE_MAIN==mpi_complex16 .and. (IC==4 .or. IC==5) ) then
       flag_R2C=.true.
    else
       flag_R2C=.false.
    end if

    if ( SYStype == 0 ) then
       i=n1-1
       do i3=Igrid(1,3),Igrid(2,3)
       do i2=Igrid(1,2),Igrid(2,2)
       do i1=Igrid(1,1),Igrid(2,1)
          i=i+1
          LL2(1,i)=i1
          LL2(2,i)=i2
          LL2(3,i)=i3
       end do
       end do
       end do

    else if ( SYStype == 1 ) then

       LL2(1:3,n1:n2) = LL(1:3,n1:n2)

    end if

    allocate( ir(0:np_grid-1),id(0:np_grid-1) )

    ir(0:np_grid-1)=3*ir_grid(0:np_grid-1)
    id(0:np_grid-1)=3*id_grid(0:np_grid-1)
    call mpi_allgatherv(LL2(1,n1),ir(myrank_g),mpi_integer &
         ,LL2,ir,id,mpi_integer,comm_grid,ierr)

    deallocate( id,ir )

!
! --- Read VRHO ---
!
    if ( IC==2 .or. IC==3 .or. IC==5 .or. IC==6 .or. IC==13 ) then

       if ( myrank==0 ) then
          open(80,file=file_vrho2,form='unformatted')
          read(80) ML_tmp,ML1_tmp,ML2_tmp,ML3_tmp
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
          if (DISP_SWITCH) then
             write(*,*) "ML1_tmp,ML2_tmp,ML3_tmp =", &
                  & ML1_tmp,ML2_tmp,ML3_tmp
             write(*,*) "ML1,ML2,ML3 =",ML1,ML2,ML3
          end if
          stop
       end if

       if ( myrank==0 ) then
          read(80) LL_tmp(:,:)
       end if
       call mpi_bcast(LL_tmp,3*ML,mpi_integer,0,mpi_comm_world,ierr)

       i=sum(abs(LL_tmp(:,:)-LL2(:,:)))
       if ( i/=0 ) then
          if (DISP_SWITCH) then
             write(*,*) "LL and LL_tmp is different"
          end if
!          if ( IO_ctrl /= 0 ) stop
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
             i1=LL_tmp(1,i) ; i2=LL_tmp(2,i) ; i3=LL_tmp(3,i)
             rtmp3(i1,i2,i3)=rtmp(i)
          end do

          do i=n1,n2
             i1=LL2(1,i) ; i2=LL2(2,i) ; i3=LL2(3,i)
             rho(i,s)=rtmp3(i1,i2,i3)
          end do

          if ( myrank==0 ) then
             read(80) rtmp(:)
          end if
          call mpi_bcast(rtmp,ML,mpi_real8,0,mpi_comm_world,ierr)
          do i=1,ML
             i1=LL_tmp(1,i) ; i2=LL_tmp(2,i) ; i3=LL_tmp(3,i)
             rtmp3(i1,i2,i3)=rtmp(i)
          end do
          if ( MSP_0<=s .and. s<=MSP_1 ) then
             do i=n1,n2
                i1=LL2(1,i) ; i2=LL2(2,i) ; i3=LL2(3,i)
                Vloc(i,s)=rtmp3(i1,i2,i3)
             end do
          end if

          if ( s==1 ) then
             if ( myrank==0 ) then
                read(80) rtmp(:)
             end if
             call mpi_bcast(rtmp,ML,mpi_real8,0,mpi_comm_world,ierr)
             do i=1,ML
                i1=LL_tmp(1,i) ; i2=LL_tmp(2,i) ; i3=LL_tmp(3,i)
                rtmp3(i1,i2,i3)=rtmp(i)
             end do
             do i=n1,n2
                i1=LL2(1,i) ; i2=LL2(2,i) ; i3=LL2(3,i)
                Vh(i)=rtmp3(i1,i2,i3)
             end do
          end if

          if ( myrank==0 ) then
             read(80) rtmp(:)
          end if
          call mpi_bcast(rtmp,ML,mpi_real8,0,mpi_comm_world,ierr)
          do i=1,ML
             i1=LL_tmp(1,i) ; i2=LL_tmp(2,i) ; i3=LL_tmp(3,i)
             rtmp3(i1,i2,i3)=rtmp(i)
          end do
          if ( MSP_0<=s .and. s<=MSP_1 ) then
             do i=n1,n2
                i1=LL2(1,i) ; i2=LL2(2,i) ; i3=LL2(3,i)
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

    end if

!
! --- Read WF ---
!

    if ( IC == 2 ) return

    call simple_wf_io_read( file_wf2, SYStype, IO_ctrl, disp_switch )

goto 100

    if ( IC == 6 ) then
       if ( allocated(LL_tmp) ) deallocate(LL_tmp)
       if ( allocated(LL2) ) deallocate(LL2)
       call read_data_io2( disp_switch )
       return
    end if

    if ( IC==1 .or. IC==3 .or. IC==4 .or. IC==5 ) then

!check format (single/double)
       if ( myrank==0 ) then
          inquire(iolength=i) int4
          if ( IO_ctrl==0 ) then
             open(3,file=file_wf2,form='unformatted', &
                  & access='direct',recl=i)
             read(3,rec=(2+4)+(2+3)+(2+3*ML)+(2+MB*MBZ*MSP*2)+1) i
             close(3)
             i=i/ML
          else if ( IO_ctrl==3 ) then
             open(3,file=file_wf_split,form='unformatted', &
                  & access='direct',recl=i)
             read(3,rec=1) i
             close(3)
             i=i/ML0
          end if
          select case(i)
          case default
             print*,'wrong WF file:',i
             stop
#ifdef _DRSDFT_
          case(4)
             print*,'WF: real(4) -> real(8)'
             flag_SP=.true.
          case(8)
             print*,'WF: real(8) -> real(8)'
             flag_SP=.false.
#else
          case(4)
             if ( flag_R2C ) then
                print*,'WF: real(4) -> complex(8)'
                flag_SP=.true.
             else
                print*,'WF: real(4)?'
                stop
             end if
          case(8)
             if ( flag_R2C ) then
                print*,'WF: real(8) -> complex(8)'
                flag_SP=.false.
             else
                print*,'WF: complex(4) -> complex(8)'
                flag_SP=.true.
             end if
          case(16)
             if ( flag_R2C ) then
                print*,'WF: complex(8)?'
                stop
             else
                print*,'WF: complex(8) -> complex(8)'
                flag_SP=.false.
             end if
#endif
          end select
       end if

       call mpi_bcast(flag_SP,1,mpi_logical,0,mpi_comm_world,ierr)

       if ( myrank==0 ) then
          open(3,file=file_wf2,form='unformatted')
          read(3) ML_tmp,ML1_tmp,ML2_tmp,ML3_tmp
          read(3) MB_tmp,MB1_tmp,MB2_tmp
          itmp(1)=ML_tmp
          itmp(2)=ML1_tmp
          itmp(3)=ML2_tmp
          itmp(4)=ML3_tmp
          itmp(5)=MB_tmp
          itmp(6)=MB1_tmp
          itmp(7)=MB2_tmp
       end if
       call mpi_bcast(itmp,7,mpi_integer,0,mpi_comm_world,ierr)
       ML_tmp=itmp(1)
       ML1_tmp=itmp(2)
       ML2_tmp=itmp(3)
       ML3_tmp=itmp(4)
       MB_tmp=itmp(5)
       MB1_tmp=itmp(6)
       MB2_tmp=itmp(7)

       if (DISP_SWITCH) then
          write(*,*) "ML,ML_tmp  =",ML,ML_tmp
          write(*,*) "MB,MB_tmp=",MB,MB_tmp
       end if
       if ( ML_tmp/=ML ) then
          write(*,*) "ML,ML_tmp  =",ML,ML_tmp
          stop
       end if
       if ( ML1_tmp/=ML1 .or. ML2_tmp/=ML2 .or. ML3_tmp/=ML3 ) then
          if ( DISP_SWITCH ) then
             write(*,*) "ML1_tmp,ML2_tmp,ML3_tmp =", &
                  & ML1_tmp,ML2_tmp,ML3_tmp
             write(*,*) "ML1,ML2,ML3 =",ML1,ML2,ML3
          end if
          stop
       end if

       if ( 1/=MB1_tmp .or. MB2_tmp/=MB ) then
          if (DISP_SWITCH) then
             write(*,*) "******** WARNING! ********"
             write(*,*) "MB1_tmp,MB2_tmp=",MB1_tmp,MB2_tmp
             write(*,*) "stop!"
          end if
          stop
       end if

       if ( myrank==0 ) then
          read(3) LL_tmp(:,:)
       end if
       call mpi_bcast(LL_tmp,3*ML,mpi_integer,0,mpi_comm_world,ierr)

       if(DISP_SWITCH)then
          write(*,*) "minval(LL_tmp),maxval(LL_tmp)"
          write(*,*) minval(LL_tmp(1,1:ML) ),maxval( LL_tmp(1,1:ML))
          write(*,*) minval(LL_tmp(2,1:ML) ),maxval( LL_tmp(2,1:ML))
          write(*,*) minval(LL_tmp(3,1:ML) ),maxval( LL_tmp(3,1:ML))
       endif

       i=sum(abs(LL_tmp(:,:)-LL2(:,:)))
       if ( i/=0 ) then
          if (DISP_SWITCH) then
             write(*,*) "LL and LL_tmp is different"
          end if
          if ( IO_ctrl /= 0 ) stop
       end if

       if ( myrank==0 ) then
          read(3) occ(:,:,:)
       end if
       call mpi_bcast(occ,MB*MBZ*MSP,mpi_real8,0,mpi_comm_world,ierr)

       if ( myrank==0 ) then
          write(*,*) "sum(occ)=",sum(occ)
       end if

       if ( SYStype == 0 ) then

          allocate( utmp3(0:ML1-1,0:ML2-1,0:ML3-1) ) ; utmp3=zero

       else if ( SYStype == 1 ) then

          mx=(Ngrid(1)-1)/2
          my=(Ngrid(1)-1)/2
          mz=(Ngrid(1)-1)/2
          allocate( utmp3(-mx:mx,-my:my,-mz:mz) ) ; utmp3=0.d0

       end if

       if ( IO_ctrl==3 ) then
          if ( myrank==0 ) then
             close(3)
          end if
          open(3,file=file_wf_split,form="unformatted")
       end if

       allocate( utmp(ML) ) ; utmp=zero

       if ( IO_ctrl==0 ) then
          if ( flag_R2C .and. flag_SP ) then
             allocate( rtmpSP(ML) ) ; rtmpSP=0.e0
          elseif ( flag_R2C ) then
             allocate( rtmpDP(ML) ) ; rtmpDP=0.d0
          elseif ( flag_SP ) then
             allocate( utmpSP(ML) ) ; utmpSP=zeroSP
          end if
       elseif ( IO_ctrl==3 ) then
          if ( flag_R2C .and. flag_SP ) then
             allocate( rtmpSP(ML0) ) ; rtmpSP=0.e0
          elseif ( flag_R2C ) then
             allocate( rtmpDP(ML0) ) ; rtmpDP=0.d0
          elseif ( flag_SP ) then
             allocate( utmpSP(ML0) ) ; utmpSP=zeroSP
          end if
       end if

       do s=1,MSP
       do k=1,MBZ
       do n=MB1_tmp,min(MB2_tmp,MB)

          do irank=0,nprocs-1
             i1=id_band(id_class(irank,4))+1
             j1=id_band(id_class(irank,4))+ir_band(id_class(irank,4))
             i2=id_bzsm(id_class(irank,5))+1
             j2=id_bzsm(id_class(irank,5))+ir_bzsm(id_class(irank,5))
             i3=id_spin(id_class(irank,6))+1
             j3=id_spin(id_class(irank,6))+ir_spin(id_class(irank,6))
             if ( id_grid(id_class(irank,0))==0 .and. &
                  i1<=n .and. n<=j1 .and. &
                  i2<=k .and. k<=j2 .and. &
                  i3<=s .and. s<=j3 ) exit
          end do
          if ( irank>=nprocs ) then
             write(*,*) "ERROR(read_data)",myrank
             stop
          end if
            
          flag_related = .false.
          if ( MBZ_0<=k .and. k<=MBZ_1 .and. MB_0 <=n .and. n<=MB_1 &
               & .and. MSP_0<=s .and. s<=MSP_1 ) then
             flag_related=.true.
          end if

          if ( IO_ctrl==0 ) then

             if ( myrank==0 ) then
#ifdef _DRSDFT_
                if ( flag_R2C .and. flag_SP ) then
                   read(3) rtmpSP
                   utmp=rtmpSP
                elseif ( flag_R2C ) then
                   read(3) rtmpDP
                   utmp=rtmpDP
                elseif ( flag_SP ) then
                   read(3) utmpSP
                   utmp=utmpSP
                else
                   read(3) utmp
                endif
#else
                if ( flag_R2C .and. flag_SP ) then
                   read(3) rtmpSP
                   utmp=rtmpSP
                elseif ( flag_R2C ) then
                   read(3) rtmpDP
                   utmp=rtmpDP
                elseif ( flag_SP ) then
                   read(3) utmpSP
                   utmp=utmpSP
                else
                   read(3) utmp
                endif
#endif
                do i=1,ML
                   i1=LL_tmp(1,i) ; i2=LL_tmp(2,i) ; i3=LL_tmp(3,i)
                   utmp3(i1,i2,i3)=utmp(i)
                end do
                do i=1,ML
                   i1=LL2(1,i) ; i2=LL2(2,i) ; i3=LL2(3,i)
                   utmp(i)=utmp3(i1,i2,i3)
                end do
             end if

             call mpi_barrier(mpi_comm_world,ierr)

             if ( irank/=0 ) then
                if ( myrank==0 ) then
                   call mpi_send(utmp,ML,TYPE_MAIN,irank,0, &
                        & mpi_comm_world,ierr)
                end if
                if ( myrank==irank ) then
                   call mpi_recv(utmp,ML,TYPE_MAIN,0,0, &
                        & mpi_comm_world,istatus,ierr)
                end if
             end if

             call mpi_barrier(mpi_comm_world,ierr)

             if ( flag_related ) then
                call mpi_scatterv(utmp,ir_grid,id_grid,TYPE_MAIN &
                     & ,unk(n1,n,k,s),ML0,TYPE_MAIN,0,comm_grid,ierr)
             end if

          elseif ( IO_ctrl==3 ) then

             if ( flag_related ) then
#ifdef _DRSDFT_
                if ( flag_R2C .and. flag_SP ) then
                   read(3) rtmpSP(1:ML0)
                   unk(n1:n2,n,k,s)=rtmpSP(1:ML0)
                elseif ( flag_R2C ) then
                   read(3) rtmpDP(1:ML0)
                   unk(n1:n2,n,k,s)=rtmpDP(1:ML0)
                elseif ( flag_SP ) then
                   read(3) utmpSP(1:ML0)
                   unk(n1:n2,n,k,s)=utmpSP(1:ML0)
                else
                   read(3) unk(n1:n2,n,k,s)
                end if
#else
                if ( flag_R2C .and. flag_SP ) then
                   read(3) rtmpSP(1:ML0)
                   unk(n1:n2,n,k,s)=rtmpSP(1:ML0)
                elseif ( flag_R2C ) then
                   read(3) rtmpDP(1:ML0)
                   unk(n1:n2,n,k,s)=rtmpDP(1:ML0)
                elseif ( flag_SP ) then
                   read(3) utmpSP(1:ML0)
                   unk(n1:n2,n,k,s)=utmpSP(1:ML0)
                else
                   read(3) unk(n1:n2,n,k,s)
                end if
#endif
             end if

          end if

       end do ! n
       end do ! k
       end do ! s

       if ( IO_ctrl==0 ) then
          if ( myrank==0 ) then
             close(3)
          end if
       else
          close(3)
       end if

       deallocate( utmp,utmp3 )
       if ( flag_R2C .and. flag_SP ) then
          deallocate( rtmpSP )
       elseif ( flag_R2C ) then
          deallocate( rtmpDP )
       elseif ( flag_SP ) then
          deallocate( utmpSP )
       end if

       if (DISP_SWITCH) then
          write(*,*) "read from ",file_wf2
       end if

    end if

100 continue

    deallocate( LL_tmp )
    deallocate( LL2 )

    return

  END SUBROUTINE read_data


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


END MODULE io_module
