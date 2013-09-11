MODULE parameters_module

  use global_variables

  implicit none

  PRIVATE
  PUBLIC :: read_parameters

  integer,parameter :: unit=1

CONTAINS

  SUBROUTINE read_parameters
    implicit none
    integer :: i
    character(7) :: label,cbuf,ckey

!    if ( myrank == 0 ) then
!       do i=1,1000
!          read(unit,'(a7)') label
!          if ( label=='# RSDFT' ) exit
!       end do
!       if ( i > 1000 ) then
!          write(*,*) "'# RSDFT' is not found"
!          call mpi_abort(i)
!          stop
!       end if
!    end if

    call read_xc(myrank,unit)

    call read_aa(myrank,unit)

    call read_electron(myrank,unit)

    call read_kgrid_bz(myrank,unit)

    call read_cg(myrank,unit)

    call read_rgrid(myrank,unit)

    call read_kinetic(myrank,unit)

    call read_ppname_pseudopot(myrank,unit)

    call read_mixing(myrank,unit)

    call read_param_pseudopot(myrank,unit)

    call read_ps_nloc1(myrank,unit)
    call read_ps_nloc2_init(myrank,unit)

    call read_parallel(myrank,unit)

    call read_gram_schmidt_t(myrank,unit)

    call read_cgpc(myrank,unit)

    call read_fermi(myrank,unit)

    call read_io(myrank,unit)

    call read_watch(myrank,unit)

    call read_atom(myrank,970)

    Diter  = 100
    Nsweep = 0
    if ( myrank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:5) == "DITER" ) then
             backspace(unit)
             read(unit,*) cbuf,Diter
          else if ( ckey(1:6) == "NSWEEP" ) then
             backspace(unit)
             read(unit,*) cbuf,Nsweep
          end if
       end do
999    continue
       write(*,*) "Diter =",Diter
       write(*,*) "Nsweep=",Nsweep
    end if

    iswitch_scf  = 1
    iswitch_opt  = 0
    iswitch_band = 0
    if ( myrank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=990) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:5) == "SWSCF" ) then
             backspace(unit)
             read(unit,*) cbuf,iswitch_scf
          else if ( ckey(1:5) == "SWOPT" ) then
             backspace(unit)
             read(unit,*) cbuf,iswitch_opt
          else if ( ckey(1:6) == "SWBAND" ) then
             backspace(unit)
             read(unit,*) cbuf,iswitch_band
          end if
       end do
990    continue
       write(*,*) "iswitch_scf =",iswitch_scf
       write(*,*) "iswitch_opt =",iswitch_opt
       write(*,*) "iswitch_band=",iswitch_band
    end if

    call send_parameters(0)

    call read_atomopt(myrank,unit)

  END SUBROUTINE read_parameters


  SUBROUTINE send_parameters(rank)
    implicit none
    integer,intent(IN) :: rank
    integer :: ierr

!    call send_aa(0)
!    call send_electron(0)
!    call send_rgrid(0)
!    call send_kinetic(0)
!    call send_xc(0)

!    call send_cg(0)
!    call send_mixing(0)

!    call send_kgrid_bz(0)
!    call send_parallel(0)

!    call send_param_pseudopot(0)
!    call send_ps_nloc1(0)
!    call send_ps_nloc2(0)

!    call send_gram_schmidt_t(0)

!    call send_cgpc(0)

!    call send_fermi(0)

! myrank/=0 must allocate some arrays in send_atom
!    call send_atom(myrank)

!    call send_io(0)

    call mpi_bcast(Diter ,1,mpi_integer,rank,mpi_comm_world,ierr)
    call mpi_bcast(Nsweep,1,mpi_integer,rank,mpi_comm_world,ierr)
    call mpi_bcast(iswitch_scf ,1,mpi_integer,rank,mpi_comm_world,ierr)
    call mpi_bcast(iswitch_opt ,1,mpi_integer,rank,mpi_comm_world,ierr)
    call mpi_bcast(iswitch_band,1,mpi_integer,rank,mpi_comm_world,ierr)

!    call send_watch(0)

  END SUBROUTINE send_parameters

END MODULE parameters_module
