module parameters_module

  use global_variables, only: iswitch_tddft, iswitch_scf, iswitch_opt, iswitch_band, iswitch_dos, iswitch_latopt, iswitch_test
  use io_tools_module
  use band_unfold_module, only: read_band_unfold
  use xc_hybrid_module, only: read_xc_hybrid
  use vdw_grimme_module, only: read_vdw_grimme
  use bz_module
  use ps_nloc2_init_module
  use gram_schmidt_t_module
  use cg_module
  use cgpc_module
  use mixing_module
  use io_module
  use watch_module
  use atomopt_module
  use xc_module
  use sweep_module
  use scf_module

  implicit none

  private
  public :: read_parameters

  integer,parameter :: unit=1

contains

  subroutine read_parameters
    implicit none
    integer :: ierr,myrank
    include 'mpif.h'

    call write_border( 0, " read_parameters(start)" )

    call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr )

    call read_bz

    call read_ps_nloc2_init(myrank,unit)

    call read_gram_schmidt_t(myrank,unit)

    call read_cg

    call read_cgpc(myrank,unit)

    call read_mixing

    call read_io

    call read_watch

    iswitch_scf   = 0 
    iswitch_opt   = 0
    iswitch_latopt= 0
    iswitch_band  = 0
    iswitch_test  = 0
    iswitch_tddft = 0
    iswitch_dos   = 0
    call IOTools_readIntegerKeyword( "SWSCF"  ,iswitch_scf   )
    call IOTools_readIntegerKeyword( "SWOPT"  ,iswitch_opt   )
    call IOTools_readIntegerKeyword( "LATOPT" ,iswitch_latopt)
    call IOTools_readIntegerKeyword( "SWBAND" ,iswitch_band  )
    call IOTools_readIntegerKeyword( "SWTEST" ,iswitch_test  )
    call IOTools_readIntegerKeyword( "SWTDDFT",iswitch_tddft )
    call IOTools_readIntegerKeyword( "SWDOS"  ,iswitch_dos   )

    call read_atomopt(myrank,unit)

    call read_xc
    call read_xc_hybrid
    call read_vdw_grimme

    call read_sweep
    call read_scf

    call read_band_unfold( myrank, unit )

    call write_border( 0," read_parameters(end)" )

  end subroutine read_parameters

end module parameters_module
