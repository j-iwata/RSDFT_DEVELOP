MODULE rtddft_mol_module

  use wf_module
  use rgrid_mol_module, only: Hsize, LL
  use grid_module, only: grid, get_range_rgrid
  use hamiltonian_module
  use parallel_module, only: comm_grid
  use density_module
  use total_energy_module
  use hartree_module, only: calc_hartree
  use xc_module, only: calc_xc
  use localpot_module, only: update_localpot
  use watch_module
  use io_module

  implicit none

  PRIVATE
  PUBLIC :: rtddft_mol, init_rtddft_mol

  type td
     real(8) :: dt
     real(8) :: tmax
     integer :: nt
     integer :: ialg
     integer :: nalg
     real(8) :: field(3)
     real(8) :: strength
  END type td

  type(td) :: tddft
  real(8),allocatable :: dipole(:,:)
  integer,parameter :: unit=91
  character(14),parameter :: file_tddft = 'rtddft_dat_mol'

  real(8),parameter :: Tau2fs = 2.418884326505d-2

  integer :: nt_wr=1000

CONTAINS


  SUBROUTINE init_rtddft_mol( unit, rank )
     implicit none
     integer,intent(IN) :: unit, rank
     integer :: i
     real(8) :: sbuf(7)
     character(5) :: cbuf,ckey
     include 'mpif.h'
     call write_border( 0, " init_rtddft_mol(start)" )
     tddft%dt  =0.0d0
     tddft%nt  =0
     tddft%ialg=1
     tddft%nalg=4
     tddft%field=0.0d0
     tddft%strength=0.0d0
     if ( rank == 0 ) then
        rewind unit
        do i=1,100000
           read(unit,*,END=900) cbuf
           call convert_capital(cbuf,ckey)
           if ( ckey == "TDDFT" ) then
              backspace(unit)
              read(unit,*) cbuf, tddft%dt, tddft%nt, tddft%ialg, tddft%nalg
           else if ( ckey == "FIELD" ) then
              backspace(unit)
              read(unit,*) cbuf, tddft%field(1:3)
           end if
        end do ! i
900     continue
     end if
     sbuf(1)=tddft%dt
     sbuf(2)=tddft%nt
     sbuf(3)=tddft%ialg
     sbuf(4)=tddft%nalg
     sbuf(5:7)=tddft%field(1:3)
     call MPI_BCAST(sbuf,7,MPI_REAL8,0,MPI_COMM_WORLD,i)
     tddft%dt         = sbuf(1)
     tddft%nt         = nint( sbuf(2) )
     tddft%ialg       = nint( sbuf(3) )
     tddft%nalg       = nint( sbuf(4) )
     tddft%field(1:3) = sbuf(5:7)

     tddft%tmax = tddft%dt * tddft%nt
     tddft%strength = sqrt(sum(tddft%field(:)**2))

     if ( rank == 0 ) then
        write(*,*) "dt=",tddft%dt
        write(*,*) "nt=",tddft%nt
        write(*,*) "tmax=",tddft%tmax
        write(*,*) "ialgorithm=",tddft%ialg
        write(*,*) "nalgorithm=",tddft%nalg
        write(*,'("field=",3f15.8)') tddft%field(1:3)
        write(*,'("strength=",3f15.8)') tddft%strength
     end if

     call write_border( 0, "init_rtddft_mol(end)" )

  END SUBROUTINE init_rtddft_mol


  SUBROUTINE rtddft_mol( job_ctrl )
    implicit none
    integer,intent(IN) :: job_ctrl
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    complex(8),parameter :: zi=(0.0d0,1.0d0)
    complex(8),allocatable :: tpsi(:,:),hpsi(:,:),zcoef(:)
    integer :: itaylor,i,n,k,s,it,ierr,t_1,t_0,myrank
    real(8) :: c,t,ct(0:9),et(0:9)
    type(grid) :: rgrid
    logical :: disp_sw,flag_end,flag_wr
    real(8) :: Etot
    include 'mpif.h'

#ifdef _DRSDFT_

    write(*,*) "rtddft_mol is available only for COMPLEX16 calculations"
    return

#else

    call write_border( 0, " rtddft_mol(start)" )
    call check_disp_switch( disp_sw, 0 )
!    call check_disp_switch( .false., 1 )

    call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )

    ct(:)=0.0d0
    et(:)=0.0d0

    call get_range_rgrid( rgrid )

! ---

    call Init_IO( "tddft" )

! ---

    if ( job_ctrl == 1 ) then

       if ( myrank == 0 ) open(unit,file=file_tddft,status="replace",position="rewind")

       call initial_condition

       call calc_total_energy( .true., Etot )

       t_0 = 1
       t_1 = tddft%nt
       allocate( dipole(0:3,0:tddft%nt) ) ; dipole=0.0d0

       call calc_dipole( dipole(0,0), rgrid%VolumeElement )

       if ( disp_sw ) then
          write(*,'(1x,i6,1x,6f16.10,1x,2f10.5)') 0,0.0,dipole(1:3,0),Etot,dipole(0,0)
       end if
       if ( myrank == 0 ) then
          write(unit,'(1x,i6,1x,6f22.16)') 0,0.0,dipole(1:3,0),Etot,dipole(0,0)
       end if

    else if ( job_ctrl == 2 ) then

       if ( myrank == 0 ) then
          open(unit,file=file_tddft,status="old",position="append")
          backspace(unit)
          read(unit,*) t_0
       end if
       call MPI_BCAST( t_0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )

       t_0 = t_0 + 1
       t_1 = t_0 + tddft%nt - 1
       allocate( dipole(0:3,t_0:t_1) ) ; dipole=0.0d0

    end if

    allocate( tpsi(ML_0_WF:ML_1_WF,1) ) ; tpsi=zero
    allocate( hpsi(ML_0_WF:ML_1_WF,1) ) ; hpsi=zero
    allocate( zcoef(tddft%nalg)       ) ; zcoef=zero

    do itaylor=1,tddft%nalg
       c=1.0d0
       do i=1,itaylor
          c=c*i
       end do
       zcoef(itaylor) = (-zi*tddft%dt)**itaylor/c
    end do

! ---

    do it=t_0,t_1

       call watch(ct(0),et(0))

       t = it*tddft%dt

       do s=MS_0_WF,MS_1_WF
       do k=MK_0_WF,MK_1_WF
       do n=MB_0_WF,MB_1_WF

          tpsi(:,1) = unk(:,n,k,s)
          do itaylor=1,tddft%nalg
             call hamiltonian(k,s,tpsi,hpsi,ML_0_WF,ML_1_WF,1,1)
             unk(:,n,k,s) = unk(:,n,k,s) + zcoef(itaylor)*hpsi(:,1)
             tpsi(:,:) = hpsi(:,:)
          end do

       end do ! n
       end do ! k
       end do ! s

       call calc_density
       call calc_hartree( rho )
       call calc_xc
       call update_localpot

       call calc_total_energy( .true., Etot )

       call calc_dipole( dipole(0,it), rgrid%VolumeElement )

       call watch(ct(1),et(1))
       call global_watch( .false., flag_end )

       if ( disp_sw ) then
          write(*,'(1x,i6,1x,6f16.10,1x,2f10.5)') &
               it,t,dipole(1:3,it),Etot,dipole(0,it),ct(1)-ct(0),et(1)-et(0)
       end if
       if ( myrank == 0 ) then
          write(unit,'(1x,i6,1x,6f22.16)') it,t,dipole(1:3,it),Etot,dipole(0,it)
       end if

       flag_wr = ( mod(it-t_0+1,nt_wr)==0 .or. it==t_1 .or. flag_end )
       call write_data( disp_sw, flag_wr, suffix="tddft" )

       if ( flag_end ) exit

    end do ! it

    deallocate( zcoef  )
    deallocate( hpsi   )
    deallocate( tpsi   )
    deallocate( dipole )

! ---

    if ( myrank == 0 ) close(unit)

! ---

    call write_border( 0, " rtddft_mol(end)" )

#endif

  END SUBROUTINE rtddft_mol


  SUBROUTINE initial_condition
    implicit none
    integer :: i,n,k,s
    real(8) :: x,y,z,kr,kx,ky,kz

    kx = tddft%field(1)
    ky = tddft%field(2)
    kz = tddft%field(3)

    do s=MS_0_WF,MS_1_WF
    do k=MK_0_WF,MK_1_WF
    do n=MB_0_WF,MB_1_WF

       do i=ML_0_WF,ML_1_WF
          x=LL(1,i)*Hsize
          y=LL(2,i)*Hsize
          z=LL(3,i)*Hsize
          kr=x*kx+y*ky+z*kz
          unk(i,n,k,s) = dcmplx(cos(kr),sin(kr))*unk(i,n,k,s)
       end do ! i

    end do ! n
    end do ! k
    end do ! s

  END SUBROUTINE initial_condition


  SUBROUTINE calc_dipole( d, dV )
    implicit none
    real(8),intent(OUT) :: d(0:3)
    real(8),intent(IN) :: dV
    integer :: i
    real(8) :: x,y,z,c,d0(0:3),trho
    include 'mpif.h'

    c=0.5d0 ; if ( MS_WF == 2 ) c=1.0d0

    d0(:)=0.0d0
    do i=ML_0_WF,ML_1_WF
       x=LL(1,i)*Hsize
       y=LL(2,i)*Hsize
       z=LL(3,i)*Hsize
       trho=c*( rho(i,1) + rho(i,MS_WF) )
       d0(0) = d0(0) + trho
       d0(1) = d0(1) + x*trho
       d0(2) = d0(2) + y*trho
       d0(3) = d0(3) + z*trho
    end do ! i
    d0(:)=d0(:)*dV

    call MPI_ALLREDUCE(d0,d,4,MPI_REAL8,MPI_SUM,comm_grid,i)

  END SUBROUTINE calc_dipole


END MODULE rtddft_mol_module
