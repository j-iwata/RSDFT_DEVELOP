MODULE rgrid_module

  use rgrid_variables
  use rgrid_sol_module
  use rgrid_mol_module
 !use esm_rgrid_module
  use parallel_module, only: node_partition, np_grid, pinfo_grid, myrank, myrank_g &
       , id_grid, comm_grid, ircnt, idisp, ir_grid
  use lattice_module, only: lattice
  use grid_module, only: grid_info

  implicit none

  PRIVATE
  PUBLIC :: Igrid, Ngrid, Hgrid, dV, zdV
  PUBLIC :: Init_Rgrid, InitParallel_Rgrid
  PUBLIC :: set_grid_info_rgrid

  integer :: SYStype=0

CONTAINS


  SUBROUTINE Init_Rgrid( SYStype_in, aa )
    implicit none
    integer,intent(IN) :: SYStype_in
    type(lattice),intent(INOUT) :: aa
    logical :: disp

    call write_border( 0," init_rgrid(start)")

    SYStype = SYStype_in

    select case( SYStype )

    case(0) ! ----- SOL Sol sol

       call Read_RgridSol

       call Init_RgridSol( aa%LatticeVector )

    case(1) ! ----- MOL Mol mol

       call Read_RgridMol

       call GetNumGrids_RgridMol(Ngrid)

       call GetSimBox_RgridMol( aa%LatticeVector )

       aa%LatticeConstant = 1.0d0

       call GetGridSize_RgridMol(Hgrid)

       dV  = Hgrid(1)*Hgrid(2)*Hgrid(3)
       zdV = dV

!    case(2) ! ----- ESM Esm esm
!
!       call Init_RgridSol(aa)
!
!       call Read_RgridESM(myrank,unit)
!
!       call Init_RgridESM(aa,Ngrid,Md)

    end select

    call check_disp_switch( disp, 0 )
    if ( disp ) then
       write(*,*) "SYStype=",SYStype
       write(*,'(1x,"Ngrid(1:3)=",3i5)') Ngrid(1:3)
       write(*,'(1x,"Ngrid(0)  =",i8 )') Ngrid(0)
       write(*,'(1x,"Hgrid(1:3)=",3f15.8)') Hgrid(1:3)
    end if

    call write_border( 0," init_rgrid(end)")

  END SUBROUTINE Init_Rgrid


  SUBROUTINE InitParallel_Rgrid
    implicit none
    integer :: ierr, Nshift(3)
    include 'mpif.h'

    call write_border( 0, " InitParallel_Rgrid(start)" )

    Nshift(:) = 0

    select case(SYStype)

    case(0) ! ----- SOL sol

       call InitParallel_RgridSol( node_partition, np_grid, pinfo_grid )

    case(1) ! ----- MOL mol

       call InitParallel_RgridMol(node_partition,np_grid,pinfo_grid,myrank==0)

       Nshift(1:3) = -1

!    case(2)
!
!       call InitParallel_RgridSol( node_partition, np_grid, pinfo_grid )
!       call InitParallel_RgridESM( np_grid, pinfo_grid, Nshift )

    end select

    Igrid(1,0) = pinfo_grid(7,myrank_g) + 1
    Igrid(2,0) = pinfo_grid(7,myrank_g) + pinfo_grid(8,myrank_g)
    Igrid(1,1) = pinfo_grid(1,myrank_g)
    Igrid(2,1) = pinfo_grid(1,myrank_g) + pinfo_grid(2,myrank_g) - 1
    Igrid(1,2) = pinfo_grid(3,myrank_g)
    Igrid(2,2) = pinfo_grid(3,myrank_g) + pinfo_grid(4,myrank_g) - 1
    Igrid(1,3) = pinfo_grid(5,myrank_g)
    Igrid(2,3) = pinfo_grid(5,myrank_g) + pinfo_grid(6,myrank_g) - 1

    Igrid(1,1:3) = Igrid(1,1:3) - Nshift(1:3)
    Igrid(2,1:3) = Igrid(2,1:3) - Nshift(1:3)

    id_grid(:) = pinfo_grid(7,:)
    ir_grid(:) = pinfo_grid(8,:)

    idisp(myrank) = id_grid(myrank_g)
    ircnt(myrank) = ir_grid(myrank_g)
    call mpi_allgather(id_grid(myrank_g),1,mpi_integer &
         ,idisp,1,mpi_integer,mpi_comm_world,ierr)
    call mpi_allgather(ir_grid(myrank_g),1,mpi_integer &
         ,ircnt,1,mpi_integer,mpi_comm_world,ierr)

    call write_border( 0, " InitParallel_Rgrid(end)" )

  END SUBROUTINE InitParallel_Rgrid


  SUBROUTINE set_grid_info_rgrid( rgrid, lattice_vector, Md )
    implicit none
    type(grid_info),intent(INOUT) :: rgrid
    real(8),intent(IN) :: lattice_vector(:,:)
    integer,intent(IN) :: Md
    rgrid%head(0:3) = Igrid(1,0:3)
    rgrid%tail(0:3) = Igrid(2,0:3)
    rgrid%size(0:3) = Igrid(2,0:3)-Igrid(1,0:3)+1
    rgrid%spacing(1:3) = Hgrid(1:3)
    rgrid%dV = dV
    rgrid%lattice(1:3,1:3) = lattice_vector(1:3,1:3)
    if ( SYStype == 1 ) then
      call GetGridSize_RgridMol( Rsize_out=rgrid%Rsize, Zsize_out=rgrid%Zsize )
    else
      rgrid%Rsize=0.0d0
      rgrid%Zsize=0.0d0
    end if
    rgrid%comm = comm_grid
    rgrid%Norder_FD = Md
  END SUBROUTINE set_grid_info_rgrid


END MODULE rgrid_module
