MODULE parameters_module

  use global_variables

  implicit none

  PRIVATE
  PUBLIC :: read_parameters, read_oldformat_parameters

  integer,parameter :: unit=1, unit_atom=970

  integer :: atom_format

CONTAINS

  SUBROUTINE read_parameters
    implicit none
    integer :: i
    character(7) :: label,cbuf,ckey
    real(8) :: Ratom(3)

    call read_atom(myrank,unit_atom)

    call read_xc(myrank,unit)

    call read_aa(myrank,unit)

    call read_electron(myrank,unit)

    call read_kgrid_bz(myrank,unit)

    call read_cg(myrank,unit)

    call read_rgrid(myrank,unit)

    call read_kinetic(myrank,unit)

    call read_ppname_pseudopot(Nelement,myrank,unit)

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
          else if ( ckey(1:6) == "NDIAG" ) then
             backspace(unit)
             read(unit,*) cbuf,Ndiag
          end if
       end do
999    continue
       write(*,*) "Diter =",Diter
       write(*,*) "Nsweep=",Nsweep
    end if

    iswitch_scf  = 1
    iswitch_opt  = 0
    iswitch_test = 0
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
          else if ( ckey(1:6) == "SWTEST" ) then
             backspace(unit)
             read(unit,*) cbuf,iswitch_test
          end if
          end if
       end do
990    continue
       write(*,*) "iswitch_scf =",iswitch_scf
       write(*,*) "iswitch_opt =",iswitch_opt
       write(*,*) "iswitch_band=",iswitch_band
       write(*,*) "iswitch_test=",iswitch_test
    end if

    atom_format = 0
    if ( myrank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=980) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:3) == "XYZ" ) then
             atom_format = 1
             write(*,*) "atom_format=",atom_format,ckey(1:3)
             exit
          end if
       end do
980    continue
    end if

    call send_parameters(0)

    call read_atomopt(myrank,unit)

    if ( SYStype == 1 ) then
       if ( atom_format == 0 ) then
          aa=ax*aa
          do i=1,Natom
             Ratom(1:3) = matmul( aa, aa_atom(:,i) )
             aa_atom(1:3,i) = Ratom(1:3)
             if ( myrank == 0 ) write(*,'(1x,i6,3f20.15)') i,aa_atom(1:3,i)
          end do
       end if
       ax=1.0d0
       aa=0.0d0
       aa(1,1)=1.0d0
       aa(2,2)=1.0d0
       aa(3,3)=1.0d0
       if ( myrank == 0 ) then
          write(*,*) "Molecule"
          write(*,*) "aa & aa_atom are modified"
       end if
    else if ( SYStype == 0 ) then
       if ( atom_format == 1 ) then
          call construct_bb(aa)
          bb(:,:)=transpose(bb(:,:))/(ax*2.0d0*acos(-1.0d0))
          do i=1,Natom
             !aa_atom(1,i) = aa_atom(1,i)/( ax*aa(1,1) )
             !aa_atom(2,i) = aa_atom(2,i)/( ax*aa(2,2) )
             !aa_atom(3,i) = aa_atom(3,i)/( ax*aa(3,3) )
             aa_atom(1:3,i) = matmul( bb(1:3,1:3), aa_atom(1:3,i) )
          end do
       end if
    end if

  END SUBROUTINE read_parameters


  SUBROUTINE read_oldformat_parameters
    implicit none
    integer :: i
    character(7) :: label,cbuf,ckey

    call read_atom(myrank,unit_atom)

    if ( myrank == 0 ) then
       do i=1,1000
          read(unit,'(a7)') label
          if ( label=='# RSDFT' ) exit
       end do
       if ( i > 1000 ) then
          write(*,*) "'# RSDFT' is not found"
          call mpi_abort(i)
          stop
       end if
    end if

    call read_oldformat_xc(myrank,unit)

    call read_oldformat_aa(myrank,unit)

    call read_oldformat_electron(myrank,unit)

    call read_kgrid_oldformat_bz(myrank,unit)

    call read_oldformat_cg(myrank,unit)

    call read_oldformat_rgrid(myrank,unit)

    call read_oldformat_kinetic(myrank,unit)

    call read_ppname_oldformat_pseudopot(Nelement,myrank,unit)

    call read_oldformat_mixing(myrank,unit)

    call read_param_oldformat_pseudopot(myrank,unit)
    call read_oldformat_ps_nloc1(myrank,unit)
    call read_oldformat_ps_nloc2_init(myrank,unit)

    call read_oldformat_parallel(myrank,unit)

    call read_oldformat_gram_schmidt_t(myrank,unit)

    call read_oldformat_cgpc(myrank,unit)

    call read_oldformat_fermi(myrank,unit)

    if ( myrank == 0 ) then
       read(unit,*) Diter, Nsweep, Ndiag
       write(*,*) "Diter =",Diter
       write(*,*) "Nsweep=",Nsweep
       write(*,*) "Nsweep=",Ndiag
    end if

    call read_oldformat_io(myrank,unit)

    if ( myrank == 0 ) then
       read(unit,*) iswitch_scf,iswitch_opt,iswitch_band,iswitch_test
       write(*,*) "iswitch_scf =",iswitch_scf
       write(*,*) "iswitch_opt =",iswitch_opt
       write(*,*) "iswitch_band=",iswitch_band
       write(*,*) "iswitch_test=",iswitch_test
    end if

    call read_oldformat_watch(myrank,unit)

    call read_oldformat_atomopt(myrank,unit)

    call send_parameters(0)

  END SUBROUTINE read_oldformat_parameters


  SUBROUTINE send_parameters(rank)
    implicit none
    integer,intent(IN) :: rank
    integer :: ierr
    call mpi_bcast(Diter ,1,mpi_integer,rank,mpi_comm_world,ierr)
    call mpi_bcast(Nsweep,1,mpi_integer,rank,mpi_comm_world,ierr)
    call mpi_bcast(Ndiag ,1,mpi_integer,rank,mpi_comm_world,ierr)
    call mpi_bcast(iswitch_scf ,1,mpi_integer,rank,mpi_comm_world,ierr)
    call mpi_bcast(iswitch_opt ,1,mpi_integer,rank,mpi_comm_world,ierr)
    call mpi_bcast(iswitch_band,1,mpi_integer,rank,mpi_comm_world,ierr)
    call mpi_bcast(iswitch_test,1,mpi_integer,rank,mpi_comm_world,ierr)
    call mpi_bcast(atom_format,1,mpi_integer,rank,mpi_comm_world,ierr)
  END SUBROUTINE send_parameters


END MODULE parameters_module
