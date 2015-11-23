MODULE parameters_module

  use global_variables
  use info_module
  use io_tools_module

  implicit none

  PRIVATE
  PUBLIC :: read_parameters

  integer,parameter :: unit=1

CONTAINS

  SUBROUTINE read_parameters
    implicit none

    if ( disp_switch_parallel ) then
       write(*,'(a50," read_parameters(START)")') repeat("-",50)
    end if

    call read_electron(myrank,unit)

    call Read_RgridSol(myrank,unit)

    call read_kinetic(myrank,unit)

    call read_ps_nloc1(myrank,unit)
    call read_ps_nloc2_init(myrank,unit)

    call read_parallel(myrank,unit)

    call read_gram_schmidt_t(myrank,unit)

    call read_cgpc(myrank,unit)

    call read_fermi(myrank,unit)

    call read_io(myrank,unit)

    call read_watch(myrank,unit)

    iswitch_scf   = 1
    iswitch_opt   = 0
    iswitch_band  = 0
    iswitch_test  = 0
    iswitch_tddft = 0
    call IOTools_readIntegerKeyword( "SWSCF"  ,iswitch_scf   )
    call IOTools_readIntegerKeyword( "SWOPT"  ,iswitch_opt   )
    call IOTools_readIntegerKeyword( "SWBAND" ,iswitch_band  )
    call IOTools_readIntegerKeyword( "SWTEST" ,iswitch_test  )
    call IOTools_readIntegerKeyword( "SWTDDFT",iswitch_tddft )

    call read_atomopt(myrank,unit)

    call read_symmetry( myrank, unit )

    call read_sweep( myrank, unit )

    select case( iswitch_scf )
    case default
       call read_scf( myrank, unit )
    case( 2 )
       call read_scf_chefsi( myrank, unit )
    end select

    if ( disp_switch_parallel ) then
       write(*,'(a50," read_parameters(END)")') repeat("-",50)
    end if

  END SUBROUTINE read_parameters

END MODULE parameters_module
