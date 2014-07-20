MODULE VarPSMember
use parallel_module, only: myrank
  implicit none

  integer :: Nelement_PP
  integer :: Nelement_
  integer,allocatable :: ippform(:)
  character(30),allocatable :: file_ps(:)
  real(8),allocatable :: rad(:,:),rab(:,:),rad1(:,:)
  real(8),allocatable :: rabr2(:,:)
  real(8),allocatable :: vql(:,:),viod(:,:,:)   ! dvql(:,:),dviod(:,:,:)
  real(8),allocatable :: cdc(:,:),cdd(:,:)
  real(8),allocatable :: anorm(:,:)

  real(8),allocatable :: Rps(:,:)   ! Rps0(:,:),Rps1(:,:)
  integer,allocatable :: lo(:,:),inorm(:,:),NRps(:,:)   ! NRps0(:,:),NRps1(:,:)

  integer,allocatable :: Mr(:),norb(:)
  real(8),allocatable :: parloc(:,:)
  real(8),allocatable :: Zps(:),Zelement(:)
  real(8),allocatable :: cdd_coef(:,:,:)

  integer,allocatable :: NRps0(:,:)
  real(8),allocatable :: Rps0(:,:)

  integer :: max_psgrd=0,max_psorb=0,max_ngauss=0

  integer,allocatable :: nlf(:)
  integer,allocatable :: nrf(:,:)
  integer,allocatable :: no(:,:)
  


CONTAINS
!------------------------------------------
  SUBROUTINE allocateRps
    implicit none
    integer :: m
if (myrank==0) write(400+myrank,*) ">>>> allocateRps"
    if (allocated( NRps0 )) then
      call deallocateRps
    end if
if (myrank==0) write(400+myrank,*) "inside allocateRps 1"
    m=max_psorb
if (myrank==0) write(400+myrank,*) "inside allocateRps 2"
    allocate( NRps0(m,Nelement_) ) ; NRps0=0
if (myrank==0) write(400+myrank,*) "inside allocateRps 3"
    allocate( Rps0(m,Nelement_)  ) ; Rps0=0.d0
    NRps0(:,:)=NRps(:,:)
    Rps0(:,:)=Rps(:,:)

if (myrank==0) write(400+myrank,*) "<<<< allocateRps"
    return
  END SUBROUTINE allocateRps

!------------------------------------------
  SUBROUTINE deallocateRps
    implicit none

    deallocate( NRps0 )
    deallocate( Rps0  )

    return
  END SUBROUTINE deallocateRps

!------------------------------------------
  SUBROUTINE send_pseudopot(myrank)
    implicit none
    integer,intent(IN) :: myrank
    integer :: m,n,ierr
    include 'mpif.h'
    m=max_psgrd
    n=max_psorb
    call mpi_bcast(m,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Nelement_PP,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if ( myrank /= 0 ) then
       call ps_allocate(m,n)
    end if
    call mpi_bcast(Mr    ,Nelement_PP,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(norb  ,Nelement_PP,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Zps   ,Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(parloc,4*Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(anorm ,n*Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(inorm ,n*Nelement_PP,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Rps   ,n*Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(NRps  ,n*Nelement_PP,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(lo    ,n*Nelement_PP,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(vql   ,m*Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(cdd   ,m*Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(cdc   ,m*Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(rad   ,m*Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(rab   ,m*Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
#ifdef _USPP_
    call mpi_bcast(nlf   ,Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nrf   ,n*Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(no    ,n*Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(rabr2 ,m*Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
#endif
    call mpi_bcast(viod  ,m*n*Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
!
    call mpi_bcast(max_ngauss,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if ( max_ngauss /= 0 ) then
       if ( myrank /= 0 ) then
          allocate( cdd_coef(3,max_ngauss,Nelement_PP) ) ; cdd_coef(:,:,:)=0.0d0
       end if
       call mpi_bcast(cdd_coef,3*max_ngauss*Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    end if
  END SUBROUTINE send_pseudopot

!-------------------------------------------------------
  SUBROUTINE ps_allocate(n_grd,n_orb)
    implicit none
    integer,intent(IN) :: n_grd,n_orb
    integer :: mg,mo
    real(8),allocatable :: vql_tmp(:,:),cdd_tmp(:,:),rad_tmp(:,:)
    real(8),allocatable :: cdc_tmp(:,:),viod_tmp(:,:,:)
    real(8),allocatable :: anorm_tmp(:,:),Rps_tmp(:,:),rab_tmp(:,:)
    real(8),allocatable :: rabr2_tmp(:,:)
    integer,allocatable :: inorm_tmp(:,:),lo_tmp(:,:),NRps_tmp(:,:)
    integer,allocatable :: nrf_tmp(:,:),no_tmp(:,:)
    if ( max_psgrd==0 .or. max_psorb==0 ) then
       allocate( Mr(Nelement_PP)   ) ; Mr(:)=0
       allocate( norb(Nelement_PP) ) ; norb(:)=0
       allocate( nlf(Nelement_PP)   ) ; nlf(:)=0
       allocate( nrf(n_orb,Nelement_PP) ) ; nrf(:,:)=0
       allocate( no(n_orb,Nelement_PP) ) ; no(:,:)=0
       allocate( Zps(Nelement_PP)      ) ; Zps(:)=0.d0
       allocate( Zelement(Nelement_PP) ) ; Zelement(:)=0.d0
       allocate( parloc(4,Nelement_PP)    ) ; parloc(:,:)=0.d0
       allocate( anorm(n_orb,Nelement_PP) ) ; anorm(:,:)=0.d0
       allocate( inorm(n_orb,Nelement_PP) ) ; inorm(:,:)=0
       allocate( Rps(n_orb,Nelement_PP)   ) ; Rps(:,:)=0.d0
       allocate( NRps(n_orb,Nelement_PP)  ) ; NRps(:,:)=0
       allocate( lo(n_orb,Nelement_PP)    ) ; lo(:,:)=0
       allocate( vql(n_grd,Nelement_PP)   ) ; vql(:,:)=0.d0
       allocate( cdd(n_grd,Nelement_PP)   ) ; cdd(:,:)=0.d0
       allocate( cdc(n_grd,Nelement_PP)   ) ; cdc(:,:)=0.d0
       allocate( rad(n_grd,Nelement_PP)   ) ; rad(:,:)=0.d0
       allocate( rad1(n_grd,Nelement_PP)  ) ; rad1(:,:)=0.d0
       allocate( rab(n_grd,Nelement_PP)   ) ; rab(:,:)=0.d0
       allocate( rabr2(n_grd,Nelement_PP) ) ; rabr2(:,:)=0.d0
       allocate( viod(n_grd,n_orb,Nelement_PP) ) ; viod(:,:,:)=0.d0
       max_psgrd=n_grd
       max_psorb=n_orb
if (myrank==0) write(400+myrank,*) "max_psgrd,psorb=",max_psgrd,max_psorb
       return
    end if
    mg = max( max_psgrd, n_grd )
    mo = max( max_psorb, n_orb )
    if ( max_psgrd < mg ) then
       allocate( vql_tmp(mg,Nelement_PP) ) ; vql_tmp(:,:)=0.d0
       allocate( cdd_tmp(mg,Nelement_PP) ) ; cdd_tmp(:,:)=0.d0
       allocate( rad_tmp(mg,Nelement_PP) ) ; rad_tmp(:,:)=0.d0
       allocate( rab_tmp(mg,Nelement_PP) ) ; rab_tmp(:,:)=0.d0
       allocate( rabr2_tmp(mg,Nelement_PP) ) ; rabr2_tmp(:,:)=0.d0
       allocate( cdc_tmp(mg,Nelement_PP) ) ; cdc_tmp(:,:)=0.d0
       vql_tmp(1:max_psgrd,1:Nelement_PP) = vql(1:max_psgrd,1:Nelement_PP)
       cdd_tmp(1:max_psgrd,1:Nelement_PP) = cdd(1:max_psgrd,1:Nelement_PP)
       rad_tmp(1:max_psgrd,1:Nelement_PP) = rad(1:max_psgrd,1:Nelement_PP)
       rab_tmp(1:max_psgrd,1:Nelement_PP) = rab(1:max_psgrd,1:Nelement_PP)
       rabr2_tmp(1:max_psgrd,1:Nelement_PP) = rabr2_tmp(1:max_psgrd,1:Nelement_PP)
       cdc_tmp(1:max_psgrd,1:Nelement_PP) = cdc(1:max_psgrd,1:Nelement_PP)
       deallocate( cdc )
       deallocate( rab )
       deallocate( rabr2 )
       deallocate( rad )
       deallocate( rad1 )
       deallocate( cdd )
       deallocate( vql )
       allocate( vql(mg,Nelement_PP) ) ; vql(:,:)=0.d0
       allocate( cdd(mg,Nelement_PP) ) ; cdd(:,:)=0.d0
       allocate( rad(mg,Nelement_PP) ) ; rad(:,:)=0.d0
       allocate( rad1(mg,Nelement_PP) ) ; rad1(:,:)=0.d0
       allocate( rab(mg,Nelement_PP) ) ; rab(:,:)=0.d0
       allocate( rabr2(mg,Nelement_PP) ) ; rabr2(:,:)=0.d0
       allocate( cdc(mg,Nelement_PP) ) ; cdc(:,:)=0.d0
       vql(:,:)=vql_tmp(:,:)
       cdd(:,:)=cdd_tmp(:,:)
       rad(:,:)=rad_tmp(:,:)
       rab(:,:)=rab_tmp(:,:)
       rabr2(:,:)=rabr2_tmp(:,:)
       cdc(:,:)=cdc_tmp(:,:)
       deallocate( cdc_tmp )
       deallocate( rab_tmp )
       deallocate( rabr2_tmp )
       deallocate( rad_tmp )
       deallocate( cdd_tmp )
       deallocate( vql_tmp )
       allocate( viod_tmp(mg,mo,Nelement_PP) ) ; viod_tmp(:,:,:)=0.d0
       viod_tmp(1:max_psgrd,1:max_psorb,1:Nelement_PP) &
            = viod(1:max_psgrd,1:max_psorb,1:Nelement_PP)
       deallocate( viod )
       allocate( viod(mg,mo,Nelement_PP) ) ; viod(:,:,:)=0.d0
       viod(:,:,:)=viod_tmp(:,:,:)
       deallocate( viod_tmp )
    end if
    if ( max_psorb < mo ) then
       allocate( anorm_tmp(mo,Nelement_PP) ) ; anorm_tmp(:,:)=0.d0
       allocate( inorm_tmp(mo,Nelement_PP) ) ; inorm_tmp(:,:)=0
       allocate( lo_tmp(mo,Nelement_PP) ) ; lo_tmp(:,:)=0
       allocate( Rps_tmp(mo,Nelement_PP) ) ; Rps_tmp(:,:)=0.d0
       allocate( NRps_tmp(mo,Nelement_PP) ) ; NRps_tmp(:,:)=0
       allocate( nrf_tmp(mo,Nelement_PP) ) ; nrf_tmp(:,:)=0
       allocate( no_tmp(mo,Nelement_PP) ) ; no_tmp(:,:)=0
       anorm_tmp(1:max_psorb,1:Nelement_PP) = anorm(1:max_psorb,1:Nelement_PP)
       inorm_tmp(1:max_psorb,1:Nelement_PP) = inorm(1:max_psorb,1:Nelement_PP)
       lo_tmp(1:max_psorb,1:Nelement_PP) = lo(1:max_psorb,1:Nelement_PP)
       Rps_tmp(1:max_psorb,1:Nelement_PP) = Rps(1:max_psorb,1:Nelement_PP)
       NRps_tmp(1:max_psorb,1:Nelement_PP) = NRps(1:max_psorb,1:Nelement_PP)
       nrf_tmp(1:max_psorb,1:Nelement_PP) = nrf(1:max_psorb,1:Nelement_PP)
       no_tmp(1:max_psorb,1:Nelement_PP) = no(1:max_psorb,1:Nelement_PP)
       deallocate( NRps )
       deallocate( Rps )
       deallocate( lo )
       deallocate( inorm )
       deallocate( anorm )
       deallocate( nrf )
       deallocate( no )
       allocate( anorm(mo,Nelement_PP) ) ; anorm(:,:)=0.d0
       allocate( inorm(mo,Nelement_PP) ) ; inorm(:,:)=0
       allocate( lo(mo,Nelement_PP)    ) ; lo(:,:)=0
       allocate( Rps(mo,Nelement_PP)   ) ; Rps(:,:)=0.d0
       allocate( NRps(mo,Nelement_PP)  ) ; NRps(:,:)=0
       allocate( nrf(mo,Nelement_PP)    ) ; nrf(:,:)=0
       allocate( no(mo,Nelement_PP)    ) ; no(:,:)=0
       anorm(:,:) = anorm_tmp(:,:)
       inorm(:,:) = inorm_tmp(:,:)
       lo(:,:) = lo_tmp(:,:)
       Rps(:,:) = Rps_tmp(:,:)
       NRps(:,:) = NRps_tmp(:,:)
       nrf(:,:) = nrf_tmp(:,:)
       no(:,:) = no_tmp(:,:)
       deallocate( NRps_tmp )
       deallocate( Rps_tmp )
       deallocate( lo_tmp )
       deallocate( inorm_tmp )
       deallocate( anorm_tmp )
       deallocate( nrf_tmp )
       deallocate( no_tmp )
       if ( max_psgrd >= mg ) then
          allocate( viod_tmp(mg,mo,Nelement_PP) ) ; viod_tmp(:,:,:)=0.d0
          viod_tmp(1:max_psgrd,1:max_psorb,1:Nelement_PP) &
               = viod(1:max_psgrd,1:max_psorb,1:Nelement_PP)
          deallocate( viod )
          allocate( viod(mg,mo,Nelement_PP) ) ; viod(:,:,:)=0.d0
          viod(:,:,:)=viod_tmp(:,:,:)
          deallocate( viod_tmp )
       end if
    end if
    max_psgrd = mg
    max_psorb = mo
  END SUBROUTINE ps_allocate

END MODULE VarPSMember
