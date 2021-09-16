MODULE rtddft_sol_module

  use io_tools_module
  use lattice_module, only: get_inverse_lattice
  use aa_module, only: aa, Va
  use bb_module, only: bb
  use bz_module, only: kbb
  use nonlocal_module, only: update_k_dependence_nonlocal
  use array_bound_module, only: MBZ_0,MBZ_1
  use current_module
  use kinetic_module
  use density_module
  use localpot_module, only: update_localpot
  use hartree_module, only: calc_hartree
  use xc_module, only: calc_xc
  use wf_module
  use ps_local_module
  use total_energy_module
  use hamiltonian_module
  use io_module
  use io_ctrl_parameters, only: OC2
  use watch_module

  implicit none

  PRIVATE
  PUBLIC :: rtddft_sol

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
  logical :: init_flag = .true.
  logical,parameter :: flag_momentum = .true.

  complex(8),allocatable :: zc(:)
  real(8) :: bb_inv(3,3)
  real(8),allocatable :: kbb_bak(:,:)

  real(8),allocatable :: Amac(:,:) ! macroscopic vector potential
  real(8) :: jav(3)  ! macroscopic (cell-averaged) current density

  real(8) :: dV
  complex(8),allocatable :: tpsi(:,:)
  complex(8),allocatable :: hpsi(:,:)

  integer,parameter :: unit_rt1=91
  integer,parameter :: unit_rt2=92
  integer,parameter :: unit_rt3=93

  integer :: myrank
  integer :: job_ctrl=0
  integer :: t_0
  integer :: nt_wr=1000

CONTAINS


  SUBROUTINE init_rtddft_sol

    implicit none
    logical :: flag
    integer :: unit,rank,ierr,n,i
    character(8) :: cbuf,posi,stat
    real(8) :: sbuf(4),Atmp(3,0:1),tx
    include 'mpif.h'
    complex(8),parameter :: zi=(0.0d0,1.0d0)

    call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )
    call IOTools_findKeyword( "TDDFT", flag, unit, .true. )
    if ( flag .and. rank == 0 ) then
       backspace(unit)
       read(unit,*) cbuf, tddft%dt, tddft%nt, tddft%ialg, tddft%nalg
       sbuf(1)=tddft%dt
       sbuf(2)=tddft%nt
       sbuf(3)=tddft%ialg
       sbuf(4)=tddft%nalg
    end if
    call MPI_BCAST( sbuf, 4, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
    tddft%dt=sbuf(1)
    tddft%nt=nint(sbuf(2))
    tddft%ialg=nint(sbuf(3))
    tddft%nalg=nint(sbuf(4))
    call IOTools_readReal8Keyword( "FIELD", tddft%field )
    call IOTools_readIntegerKeyword( "OC2TD", nt_wr )

    tddft%tmax = tddft%dt * tddft%nt
    tddft%strength = sqrt(sum(tddft%field(:)**2))

    if ( rank == 0 ) then
       write(*,*) "dt        =",tddft%dt
       write(*,*) "nt        =",tddft%nt
       write(*,*) "tmax      =",tddft%tmax
       write(*,*) "ialgorithm=",tddft%ialg
       write(*,*) "nalgorithm=",tddft%nalg
       write(*,'(1x,"field     =",3f15.8)') tddft%field(1:3)
       write(*,'(1x,"strength  =",3f15.8)') tddft%strength
       write(*,'(1x,"Output interval =",i6)') nt_wr
    end if

    OC2 = nt_wr

! ---

    allocate( zc(tddft%nalg) ) ; zc=(0.0d0,0.0d0)

    do i=1,tddft%nalg
       zc(i)=(-zi*tddft%dt)**i
       do n=2,i
          zc(i)=zc(i)/n
       end do
    end do

    call get_inverse_lattice( bb, bb_inv )

! ---

    dV = Va/ML_WF

    n=size(unk,1)
    allocate( tpsi(n,1) ) ; tpsi=(0.0d0,0.0d0)
    allocate( hpsi(n,1) ) ; hpsi=(0.0d0,0.0d0)

    myrank=rank

    if ( rank == 0 ) then
       stat="replace" ; if ( job_ctrl == 2 ) stat="old"
       posi="rewind"  ; if ( job_ctrl == 2 ) posi="append"
       open(unit_rt1,file="rtddft_dat_amac" ,status=stat,position=posi)
       open(unit_rt2,file="rtddft_dat_jemac",status=stat,position=posi)
       open(unit_rt3,file="rtddft_dat_etot" ,status=stat,position=posi)
    end if

    if ( job_ctrl == 2 ) then

       if ( rank == 0 ) then
          backspace(unit_rt1)
          read(unit_rt1,*) t_0, tx, Atmp(:,1), Atmp(:,0)
       end if
       call MPI_BCAST( t_0 ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr )
       call MPI_BCAST( Atmp,6,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr )

       allocate( Amac(3,t_0-1:t_0+tddft%nt) ) ; Amac=0.0d0
       Amac(:,t_0-1)=Atmp(:,0)
       Amac(:,t_0  )=Atmp(:,1)

    else

       t_0 = 1
       allocate( Amac(3,t_0-1:t_0+tddft%nt) ) ; Amac=0.0d0

    end if

    init_flag = .false.

  END SUBROUTINE init_rtddft_sol


  SUBROUTINE rtddft_sol( job_ctrl_in )

    implicit none
    integer,intent(IN) :: job_ctrl_in
    real(8) :: dkbb(3),const1,const2,const3,tx,tx1
    real(8) :: E_field, Emac(3), Etot, pi, Nele(2)
    integer :: k, nt_ini, t, s, t_1
    logical :: disp_sw, disp_off=.false.
    logical :: flag_end, flag_wr

    call write_border( 0, " rtddft_sol(start)" )
    call check_disp_switch( disp_sw , 0 )
    call check_disp_switch( disp_off, 1 )

    job_ctrl = job_ctrl_in

    if ( init_flag ) call init_rtddft_sol

    call Init_IO( "tddft" )

! ---

    allocate( kbb_bak(size(kbb,1),size(kbb,2)) )
    kbb_bak=kbb

! ---

    if ( job_ctrl == 1 ) then

       t_0 = 1

       jav(:)    = 0.0d0
       Amac(:,:) = 0.0d0

       Amac(1:3,0) = tddft%field(1:3)
       Amac(1:3,1) = Amac(1:3,0)

       dkbb(1:3) = matmul( bb_inv, Amac(:,0) )
       do k=1,size(kbb,2)
          kbb(:,k) = kbb_bak(:,k) + dkbb(:)
       end do

       call init_kinetic( aa,bb,size(kbb,2),kbb )
       call update_k_dependence_nonlocal( MBZ_0, MBZ_1, kbb, flag_momentum )

       call calc_macro_current( jav )

       if ( myrank == 0 ) then
          write(unit_rt1,'(1x,i6,f12.6,6f22.16)') 1,tddft%dt,Amac(:,1),Amac(:,0)
       end if

    end if

! ---

    t_1 = t_0 + tddft%nt - 1

    pi = acos(-1.0d0)
    const1 = 4.0d0*pi*(tddft%dt)**2
    const2 = 1.0d0/(2.0d0*tddft%dt)
    const3 = Va/(8.0d0*pi)

    do t=t_0,t_1

       dkbb(1:3) = matmul( bb_inv, Amac(:,t) )
       do k=1,size(kbb,2)
          kbb(:,k) = kbb_bak(:,k) + dkbb(:)
       end do

       call init_kinetic( aa,bb,size(kbb,2),kbb )
       call update_k_dependence_nonlocal( MBZ_0, MBZ_1, kbb, flag_momentum )

#ifndef _DRSDFT_
       call time_evolution( unk )
#endif

       call calc_density( Nele(1:MS_WF) )
       call calc_hartree( rho )
       call calc_xc
       call update_localpot

       call calc_macro_current( jav )
       where( Amac(:,t) == 0.0d0 ) jav=0.0d0

       Amac(:,t+1) = 2.0d0*Amac(:,t) - Amac(:,t-1) - const1*jav(:)

       Emac(:) = -( Amac(:,t+1)-Amac(:,t-1) )*const2

       call calc_total_energy( .true., Etot )
       E_field = const3*sum(Emac*Emac)

       if ( disp_sw ) then
          write(*,'(1x,i5,f12.6,3g18.9,1x,5f20.12)') &
               t,t*tddft%dt,jav(:),Etot+E_field,Etot,E_field,Nele(1:MS_WF)
       end if

       tx =t*tddft%dt
       tx1=(t+1)*tddft%dt
       if ( myrank == 0 ) then
          write(unit_rt1,'(1x,i6,f12.6,6f22.16)') t+1,tx1,Amac(:,t+1),Amac(:,t)
          write(unit_rt2,'(1x,i6,f12.6,6f22.16)') t  ,tx ,jav(:),Emac(:)
          write(unit_rt3,'(1x,i6,f12.6,5f22.16)') t  ,tx ,Etot+E_field,Etot,E_field,Nele(1:MS_WF)
       end if

       call global_watch( .false., flag_end )
       flag_wr = ( mod(t-t_0+1,nt_wr)==0 .or. t==t_1 .or. flag_end )
       call write_data( disp_sw, flag_wr, suffix="tddft" )

       if ( flag_end ) exit

    end do ! t

! ---

    if ( myrank == 0 ) then
       close(unit_rt3)
       close(unit_rt2)
       close(unit_rt1)
    end if

! ---

    call write_border( 0, " rtddft_sol(end)" )

  END SUBROUTINE rtddft_sol


  SUBROUTINE time_evolution( psi )
    implicit none
    complex(8),intent(INOUT) :: psi(:,:,:,:)
    integer :: n,k,s,m,i,s0,k0,n0
    complex(8),parameter :: z0=(0.0d0,0.0d0)
#ifndef _DRSDFT_
    do s=MS_0_WF,MS_1_WF
       s0=s-MS_0_WF+1
    do k=MK_0_WF,MK_1_WF
       k0=k-MK_0_WF+1
    do n=MB_0_WF,MB_1_WF
       n0=n-MB_0_WF+1
       if ( abs(occ(n,k,s)) < 1.d-10 ) cycle
       tpsi(:,1)=psi(:,n0,k0,s0)
       do i=1,tddft%nalg
          call hamiltonian( tpsi, hpsi, n,k,s )
          psi(:,n0,k0,s0)=psi(:,n0,k0,s0)+zc(i)*hpsi(:,1)
          tpsi=hpsi
       end do
    end do ! n
    end do ! k
    end do ! s
#endif
  END SUBROUTINE time_evolution

END MODULE rtddft_sol_module
