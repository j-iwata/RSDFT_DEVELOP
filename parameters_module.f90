MODULE parameters_module

  use global_variables
  use atomopt_module

  implicit none

  PRIVATE
  PUBLIC :: read_parameters,send_parameters,read_parameters_2

  integer,parameter :: unit=1

CONTAINS

  SUBROUTINE read_parameters
    implicit none
    integer :: i
    character(7) :: label

    do i=1,1000
       read(unit,'(a7)') label
       if ( label=='# RSDFT' ) exit
    end do
    if ( i>1000 ) stop "'# RSDFT' is not found"

    call read_xc(unit)

    call read_aa(unit)

    call read_atom(970)

    call read_electron(unit)

    call read_kgrid_bz(unit)

    call read_cg(unit)

    call read_rgrid(unit)

    call read_kinetic(unit)

    call read_ppname_pseudopot(Nelement,unit)

    call read_mixing(unit)

    call read_param_pseudopot(unit)
    call read_ps_nloc1(unit)
    call read_ps_nloc2(unit)

    call read_parallel(unit)

    call read_gram_schmidt_t(unit)

    call read_cgpc(unit)

    call read_fermi(unit)

    read(unit,*) Diter, Nsweep
    write(*,*) "Diter =",Diter
    write(*,*) "Nsweep=",Nsweep

    call read_io(unit)

    read(unit,*) iswitch_scf,iswitch_opt,iswitch_band
    write(*,*) "iswitch_scf, iswitch_opt, iswitch_band =" &
         ,iswitch_scf,iswitch_opt,iswitch_band

    call read_watch(unit)

  END SUBROUTINE read_parameters


  SUBROUTINE send_parameters
    implicit none
    integer :: ierr

    call send_aa(0)
    call send_electron(0)
    call send_rgrid(0)
    call send_kinetic(0)
    call send_xc(0)

    call send_cg(0)
    call send_mixing(0)

    call send_kgrid_bz(0)
    call send_parallel(0)

    call send_param_pseudopot(0)
    call send_ps_nloc1(0)
    call send_ps_nloc2(0)

    call send_gram_schmidt_t(0)

    call send_cgpc(0)

    call send_fermi(0)

! myrank/=0 must allocate some arrays in send_atom
    call send_atom(myrank)

    call send_io(0)

    call mpi_bcast(Diter ,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(Nsweep,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(iswitch_scf ,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(iswitch_opt ,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(iswitch_band,1,mpi_integer,0,mpi_comm_world,ierr)

    call send_watch(0)

  END SUBROUTINE send_parameters


  SUBROUTINE read_parameters_2
    implicit none
    call read_atomopt(myrank,unit)
  END SUBROUTINE read_parameters_2


END MODULE parameters_module
