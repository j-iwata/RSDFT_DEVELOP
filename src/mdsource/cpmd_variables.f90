!--------1---------2---------3---------4---------5---------6---------7---------8

MODULE cpmd_variables

  use array_bound_module
  use parallel_module
  use watch_module

!_switch
  integer :: iswitch_md
!_constants
  real(8),parameter :: TOANGS=0.52917724924d+00   ! (a.u./Ang)
  real(8),parameter :: TOMETR=TOANGS*1.0D-10      ! (Ang/m)
  real(8),parameter :: TOKCAL=627.5131018d+00     ! (kcal/mol/a.u.)
  real(8),parameter :: TOJOUL=4.35975D-18         ! (J/a.u.)
  real(8),parameter :: AVOGAD=6.022045d+26        ! (N/mole)
  real(8),parameter :: PLANCK=6.626176D-34        ! h
  real(8),parameter :: CLIGHT=2.99792458d+10      ! (m/s)
  real(8),parameter :: FEMTO=1.0D-15              ! (fs)
  real(8),parameter :: KB=1.380658d-23            ! (J/K)
  real(8),parameter :: FS_TO_AU=41.3413737896d0   ! (a.u./fs)
  real(8),parameter :: AMU=1822.8880d0            ! (atomic mass unit)
  real(8),parameter :: Pi2= 6.283185307179586d0   ! Pi*2
!_for Born-Oppenheimer molecular dynamics
  integer :: nstep
  real(8) :: dt
  real(8) :: deltat
  real(8) :: temp,trange,dsettemp
  real(8) :: omegan
  real(8) :: emass
!  integer,allocatable :: iatom(:)          ! ordering   of nuclei
  real(8),allocatable :: Velocity(:,:)     ! velocity   of nuclei
!_for Car-Parrinello molecualr dynamics
!#ifdef _DRSDFT_
  real(8),allocatable :: psi_v(:,:,:,:) ! velocity   of wavefunction
  real(8),allocatable :: psi_n(:,:,:,:) ! updated    of wavefunction 
!#else
!  complex(8),allocatable :: psi_v(:,:,:,:) ! velocity   of wavefunction
!  complex(8),allocatable :: psi_n(:,:,:,:) ! updated    of wavefunction 
!#endif
  real(8),allocatable :: zam(:,:)
  real(8),allocatable :: tau(:,:)
  real(8),allocatable :: sig(:,:)
  real(8),allocatable :: gam(:,:)
  real(8),allocatable :: gamn(:,:)
  real(8),allocatable :: scr(:,:)
  real(8),allocatable :: wrk(:,:)
  integer,allocatable :: mstocck(:,:)
  integer :: MBC,MBT
!--- band-parallel
  integer :: MB_0_CPMD, MB_1_CPMD
  integer :: MB_0_SCF, MB_1_SCF
  integer,allocatable :: ir_band_cpmd(:),id_band_cpmd(:)
! massive Nose-Hoover chain
  real(8) :: wt2(7),wt4(7),wt8(7)
  integer, parameter :: nch=4
  integer, parameter :: nche=4
  integer, parameter :: ndc=2
  integer, parameter :: nys=7
  real(8),allocatable :: xmns(:)
  real(8),allocatable :: vmns(:)
  real(8),allocatable :: gmns(:)
  real(8),allocatable :: qmns(:)
  real(8) :: fkt
  real(8) :: fnkt
! ewald
  integer,parameter :: maxkap=10
  real(8),allocatable :: Clm(:),Slm(:)
  real(8),allocatable :: elc(:,:),els(:,:)
  real(8),allocatable :: emc(:,:),ems(:,:)
  real(8),allocatable :: enc(:,:),ens(:,:)
  real(8),allocatable :: Ake(:),Cks(:),Ckc(:)
!_common
  real(8) :: vw(4),pmass(36)
  character*2  :: batm(36)
  data pmass/1.00794d0,4.003d0,6.941d0,9.012d0,10.81d0,&
       12.00d0,14.01d0,16.00d0,19.00d0,20.18d0,&
       22.99d0,24.31d0,26.98d0,28.0855d0,30.97d0,&
       32.07d0,35.45d0,39.95d0,39.10d0,40.08d0,&
       44.96d0,47.87d0,50.94d0,52.00d0,54.94d0,&
       55.84d0,58.93d0,58.69d0,63.55d0,65.39d0,&
       69.72d0,72.61d0,74.92d0,78.96d0,79.90d0,&
       83.80d0/
!_atom
  data batm/'H' ,'He','Li','Be','B' ,&
       'C' ,'N' ,'O' ,'F' ,'Ne',&
       'Na','Mg','Al','Si','P' ,&
       'S' ,'Cl','Ar','K' ,'Ca',&
       'Sc','Ti','V' ,'Cr','Mn',&
       'Fe','Co','Ni','Cu','Zn',&
       'Ga','Ge','As','Se','Br',&
       'Kr'/
!_control
  logical :: inivel,lcpmd,lbath,lmeta
!_end
  real(8) :: dtsuz(9)
  real(8) :: gkt(1,1)
  real(8) :: qnospc(nch,1)
  real(8) :: etap1(nch,1)
  real(8) :: fetapv(nch)
  real(8) :: etap1dot(nch,1)
  logical :: lbere,lbathnew,lscale,lscaleele,lquench,lforce_fast,lbathnewe
  logical :: lrestart_p,lrestart_e,linitnose,linitnosee
  real(8) :: ekin1, ekin2, ekinw
!_electron
  real(8) :: qnosee(nche)
  real(8) :: etadot(nche,1)
  real(8) :: etap(nche,1)
  integer :: nedof
  real(8) :: feta(nche+1)
  real(8) :: wnosee,wnose0

  integer :: Ndof
  integer :: MI
!  integer :: myrank
  integer,allocatable :: mkd(:)
  character(30),allocatable :: file_ps(:)
  real(8),allocatable :: asi(:,:),Rion(:,:)
  real(8),allocatable :: aa(:,:)
  real(8),allocatable :: occ(:,:,:)
!  integer :: MSP_0,MSP_1,MBZ_0,MBZ_1,MB_0,MB_1
!  integer,allocatable :: ircnt(:),idisp(:)
!  integer :: nprocs
  real(8),allocatable :: Force(:,:)
!  integer :: MB
  real(8) :: pi
!  real(8),allocatable :: unk(:,:,:,:)
  real(8) :: Etot
  logical :: DISP_SWITCH
  real(8),allocatable :: aaL(:)
  real(8),allocatable :: esp(:,:,:)
  real(8),allocatable :: Vloc(:,:)
  real(8),allocatable :: bb(:,:)
  real(8),allocatable :: Vh(:)
  real(8),allocatable :: Vxc(:,:)
  real(8),allocatable :: Vion(:)
  integer :: Diter,Diter1
  integer :: TYPE_MAIN
  integer :: MBSIZE,NBSIZE
  integer :: NBLK2
  real(8),allocatable :: vtmp2(:,:)
  integer :: NBLK1
  integer :: NBLK
!  integer :: MBZ
  real(8) :: zero
  real(8),allocatable :: utmp2(:,:)
  integer :: TRANSA,TRANSB
  real(8) :: dV
!  integer :: comm_grid
!  integer,allocatable :: id_class(:,:)
!  integer :: np_band
!  integer :: myrank_b
!  integer :: comm_band
!  integer,allocatable :: id_band(:)
!  integer,allocatable :: ir_band(:)
!  integer :: np_grid
!  integer,allocatable :: id_grid(:),ir_grid(:)
!  integer,allocatable :: id_spin(:),ir_spin(:)
!  integer,allocatable :: id_bzsm(:),ir_bzsm(:)
  integer,allocatable :: irlabel(:,:),irdim(:)
  integer,allocatable :: irlabel2n(:,:,:)
!  integer :: myrank_g
  real(8),allocatable :: utmp3(:,:,:),utmp(:)
  integer,allocatable :: LL2(:,:),LL(:,:)
!  integer :: ML
  integer :: ML1,ML2,ML3
  integer :: ML_irreducible
  integer,allocatable :: rga(:,:,:)
  integer :: isymmetry,nsym
  real(8),allocatable :: Rir(:,:,:,:)
  integer :: MBwr1,MBwr2
!  integer :: MB_d
  integer :: nspin
  integer :: SYStype
  integer :: IO_ctrl
  character(30) :: FILE_WF2,FILE_WF
!  include 'mpif.h'
!--------------for blue-moon------------------------
!  real(8),allocatable :: ylagr(:),xlagr(:)
!  integer,allocatable :: ipvt(:)
!  real(8),allocatable :: DT2BYM(:),DTB2MI(:)
!  real(8),allocatable :: pm_dim(:)
!  real(8),allocatable :: Rion0(:,:)
!  real(8),allocatable :: fc(:),fv(:)
!  real(8),allocatable :: anorm(:,:)
  logical :: lblue
!  integer :: mcnstr,nodim,ityp
!  integer :: ia(10),ib(10),ic(10)
!  real(8) :: cnpar(2,10)
!  real(8) :: cval(10)
!  integer :: index(10)

  integer :: itime
  integer :: trjstep=1
  integer :: wrtstep=0
  integer :: all_traj=0

END MODULE cpmd_variables

