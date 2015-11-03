MODULE pseudopot_module

  use parallel_module, only:myrank
  use VarPSMember
  use PSreadPSV
!#ifdef _USPP_
  use VarPSMemberG, only: sendPSG
!#endif
  use ps_read_TM_module
  use ps_read_YB_module
  use ps_read_UPF_module
  use ps_gth_module
  use io_tools_module

  implicit none

  PRIVATE
  PUBLIC :: ippform,file_ps,inorm,NRps,norb,Npseudopot &
           ,Mr,lo,no,vql,cdd,cdc,rad,anorm,viod,Rps,Zps,parloc,rab &
           ,cdd_coef,ps_type,Rcloc,hnml,knml,hnl,knl &
           ,read_ppname_pseudopot,read_pseudopot &
           ,read_param_pseudopot


  integer,PUBLIC :: pselect = 2

  integer :: Npseudopot
  integer :: unit_ps,ielm

CONTAINS


!-------------------------------------------------------
  SUBROUTINE read_ppname_pseudopot(MKI,rank,unit)
    implicit none
    integer,intent(IN) :: MKI,rank,unit
    integer :: i,j,ierr
    character(2) :: cbuf,ckey
    Nelement_PP=0
    Nelement_=MKI
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=990) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:2) == "PP" ) Nelement_PP=Nelement_PP+1
       end do
990    continue
    end if
    call send_ppname_1(0)
    allocate( ippform(Nelement_PP),file_ps(Nelement_PP) )
    ippform(:)=0
    file_ps(:)=""
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:2) == "PP" ) then
             backspace(unit)
             do j=1,Nelement_PP
                read(unit,*) cbuf,ippform(j),file_ps(j)
             end do ! j
             exit
          end if
       end do
999    continue
       write(*,*) "Nelement_PP, Nelement=",Nelement_PP,MKI
       do i=1,Nelement_PP
          write(*,'(1x,"ippform, file_ps = ",i3,2x,a30,3x,3f10.5)') &
               ippform(i),file_ps(i)
       end do
       if ( Nelement_PP > MKI ) then
          write(*,'("WARNING(Nelement_PP>Nelement):First ",i1," potentials are used")') MKI
       end if
    end if
    call send_ppname_2(0)
    ierr=0
    if ( Nelement_PP < MKI ) ierr=-1
    do i=1,Nelement_PP
       if ( ippform(i) <= 0 .or. file_ps(i)=="" ) ierr=-1
       if ( ippform(1) /= ippform(i) ) ierr=-1
    end do
    if ( ierr == -1 ) stop "stop@read_ppname_pseudopot"
    Nelement_PP=MKI
    Npseudopot =MKI
  END SUBROUTINE read_ppname_pseudopot

!-------------------------------------------------------

  SUBROUTINE send_ppname_1(rank)
    implicit none
    integer,intent(IN) :: rank
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(Nelement_PP,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)    
  END SUBROUTINE send_ppname_1

!-------------------------------------------------------

  SUBROUTINE send_ppname_2(rank)
    implicit none
    integer,intent(IN) :: rank
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(file_ps,30*Nelement_PP,MPI_CHARACTER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(ippform,Nelement_PP,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_ppname_2

!-------------------------------------------------------

  SUBROUTINE read_param_pseudopot
    implicit none
    call IOTools_readIntegerKeyword( "PSELECT", pselect )
    if ( .not.( pselect==2 .or. pselect==3 .or. pselect==102 ) ) then
       stop "invalid pselect(stop@read_param_pseudopot)"
    end if
  END SUBROUTINE read_param_pseudopot

!-------------------------------------------------------

  SUBROUTINE read_pseudopot(rank)

    implicit none
    integer,intent(IN) :: rank
    real(8),allocatable :: psi_(:,:,:),phi_(:,:,:),bet_(:,:,:)
    real(8),allocatable :: ddi_(:,:,:),qqr_(:,:,:)
    integer :: i,j,io,jo,lj

#ifdef _SHOWALL_
    if ( rank == 0 ) write(200+rank,*) '>>>>>>>> read_pseudopot'
#endif

    allocate( ps(Nelement_PP) )

    if ( rank == 0 ) then

       write(*,'(a60," read_pseudopot")') repeat("-",60) 

       max_psgrd=0
       max_psorb=0

       do ielm=1,Nelement_PP

          unit_ps=33+ielm
          open(unit_ps,FILE=file_ps(ielm),STATUS='old')

          select case(ippform(ielm))
          case(1)

             close(unit_ps)
             open(unit_ps,FILE=file_ps(ielm),form='unformatted',STATUS='old')

             call ps_read_TM(unit_ps)
             call ps_allocate(ps_tm%nrr,ps_tm%norb)
             Mr(ielm)                 = ps_tm%nrr
             norb(ielm)               = ps_tm%norb
             Zps(ielm)                = ps_tm%znuc
             anorm(1:norb(ielm),ielm) = ps_tm%anorm(1:norb(ielm))
             inorm(1:norb(ielm),ielm) = ps_tm%inorm(1:norb(ielm))
             Rps(1:norb(ielm),ielm)   = ps_tm%Rc(1:norb(ielm))
             NRps(1:norb(ielm),ielm)  = ps_tm%NRc(1:norb(ielm))
             lo(1:norb(ielm),ielm)    = ps_tm%lo(1:norb(ielm))
             vql(1:Mr(ielm),ielm)     = ps_tm%vql(1:Mr(ielm))
             cdd(1:Mr(ielm),ielm)     = ps_tm%cdd(1:Mr(ielm))
             cdc(1:Mr(ielm),ielm)     = ps_tm%cdc(1:Mr(ielm))
             rad(1:Mr(ielm),ielm)     = ps_tm%rr(1:Mr(ielm))
             rab(1:Mr(ielm),ielm)     = ps_tm%rx(1:Mr(ielm))
             viod(1:Mr(ielm),1:norb(ielm),ielm) &
                                      = ps_tm%vps(1:Mr(ielm),1:norb(ielm))

          case( 2, 102 )

             call read_PSV( unit_ps,ielm,ddi_,qqr_,psi_,phi_,bet_,ps(ielm) )

          case(3)

             call ps_read_YB(unit_ps)
             call ps_allocate(ps_yb%nrr,ps_yb%norb)
             Mr(ielm)                 = ps_yb%nrr
             norb(ielm)               = ps_yb%norb
             Zps(ielm)                = ps_yb%znuc
             anorm(1:norb(ielm),ielm) = ps_yb%anorm(1:norb(ielm))
             inorm(1:norb(ielm),ielm) = ps_yb%inorm(1:norb(ielm))
             Rps(1:norb(ielm),ielm)   = ps_yb%Rc(1:norb(ielm))
             NRps(1:norb(ielm),ielm)  = ps_yb%NRc(1:norb(ielm))
             lo(1:norb(ielm),ielm)    = ps_yb%lo(1:norb(ielm))
             vql(1:Mr(ielm),ielm)     = ps_yb%vql(1:Mr(ielm))
             cdd(1:Mr(ielm),ielm)     = ps_yb%cdd(1:Mr(ielm))
             cdc(1:Mr(ielm),ielm)     = ps_yb%cdc(1:Mr(ielm))
             rad(1:Mr(ielm),ielm)     = ps_yb%rr(1:Mr(ielm))
             rab(1:Mr(ielm),ielm)     = ps_yb%rx(1:Mr(ielm))
             viod(1:Mr(ielm),1:norb(ielm),ielm) &
                                      = ps_yb%vps(1:Mr(ielm),1:norb(ielm))

          case( 4 )

             call read_ps_gth( unit_ps, ippform(ielm) )
             call ps_allocate(1,ps_gth%norb)
             norb(ielm)             = ps_gth%norb
             Zps(ielm)              = ps_gth%znuc
             Rps(1:norb(ielm),ielm) = ps_gth%Rc(1:norb(ielm))
             lo(1:norb(ielm),ielm)  = ps_gth%lo(1:norb(ielm))
             no(1:norb(ielm),ielm)  = ps_gth%no(1:norb(ielm))
             parloc(1:4,ielm)       = ps_gth%parloc(1:4)
             Rcloc(ielm)            = ps_gth%Rcloc
             hnl(:,:,ielm)          = ps_gth%hnl(:,:)
             knl(:,:,ielm)          = ps_gth%knl(:,:)
             hnml(:,:,:,ielm)       = ps_gth%hnml(:,:,:)
             knml(:,:,:,ielm)       = ps_gth%knml(:,:,:)
             inorm(:,ielm)          = ps_gth%inorm(:)

             if ( any( hnml /= 0.0d0 ) ) ps_type=1

          case( 5 )

             call ps_read_UPF(unit_ps)
             call ps_allocate(ps_upf%nrr,ps_upf%norb)
             Mr(ielm)                 = ps_upf%nrr
             norb(ielm)               = ps_upf%norb
             Zps(ielm)                = ps_upf%znuc
             anorm(1:norb(ielm),ielm) = ps_upf%anorm(1:norb(ielm))
             inorm(1:norb(ielm),ielm) = ps_upf%inorm(1:norb(ielm))
             Rps(1:norb(ielm),ielm)   = ps_upf%Rc(1:norb(ielm))
             NRps(1:norb(ielm),ielm)  = ps_upf%NRc(1:norb(ielm))
             lo(1:norb(ielm),ielm)    = ps_upf%lo(1:norb(ielm))
             no(1:norb(ielm),ielm)    = ps_upf%no(1:norb(ielm))
             vql(1:Mr(ielm),ielm)     = ps_upf%vql(1:Mr(ielm))
             cdd(1:Mr(ielm),ielm)     = ps_upf%cdd(1:Mr(ielm))
             cdc(1:Mr(ielm),ielm)     = ps_upf%cdc(1:Mr(ielm))
             rad(1:Mr(ielm),ielm)     = ps_upf%rr(1:Mr(ielm))
             rab(1:Mr(ielm),ielm)     = ps_upf%rx(1:Mr(ielm))
             viod(1:Mr(ielm),1:norb(ielm),ielm) &
                                      = ps_upf%vps(1:Mr(ielm),1:norb(ielm))
             if ( any( ps_upf%Dij /= 0.0d0 ) ) then
                ps_type = 1
                do j=1,norb(ielm)
                   jo=no(j,ielm)
                   lj=lo(j,ielm)
                do i=1,norb(ielm)
                   io=no(i,ielm)
                   if ( lo(i,ielm) /= lj ) cycle
                   hnml(io,jo,lj,ielm) = ps_upf%Dij(i,j)
                end do
                end do
             end if

          case default

             stop "ippform error"

          end select ! ippform

          close(unit_ps)

       end do ! ielm

       write(*,*) "ps_type = ",ps_type

    end if ! [ rank == 0 ]

! --- bcast pseudopotential data

    call send_pseudopot(rank)
!#ifdef _USPP_
    if ( all(ippform == 102) ) call sendPSG(rank,Nelement_PP)
!#endif
    do ielm=1,Nelement_PP
       call ps_send_ps1d( ps(ielm) )
    end do

! ---

!    call chk_pot(1,rank)

#ifdef _SHOWALL_
    if ( rank == 0 ) write(200+rank,*) '<<<<<<<< read_pseudopot'
#endif

  END SUBROUTINE read_pseudopot


  SUBROUTINE chk_pot(iflag,rank)
    implicit none
    integer,intent(IN) :: iflag,rank
    integer :: u,ielm,i,j
    if ( rank == 0 ) then
       do ielm=1,Nelement_PP
          u=9+ielm
          rewind u
          do i=1,Mr(ielm)
             if ( iflag == 2 ) then
                write(u,'(1x,5f20.10)') rad(i,ielm),cdd(i,ielm),cdc(i,ielm)
             else
                write(u,'(1x,5f20.10)') &
                     rad(i,ielm),vql(i,ielm),(viod(i,j,ielm),j=1,norb(ielm))
             end if
          end do
          write(*,'(1x,"chk_pot(",i1,"): fort.",i2)') iflag,u
       end do
    end if
    stop "stop@chk_pot"
  END SUBROUTINE chk_pot


  SUBROUTINE send_pseudopot_1(myrank)
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
    call mpi_bcast(no    ,n*Nelement_PP,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(vql   ,m*Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(cdd   ,m*Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(cdc   ,m*Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(rad   ,m*Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(rab   ,m*Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(viod  ,m*n*Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Rcloc ,Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
!
    call mpi_bcast(max_ngauss,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if ( max_ngauss /= 0 ) then
       if ( myrank /= 0 ) then
          allocate( cdd_coef(3,max_ngauss,Nelement_PP) ) ; cdd_coef=0.0d0
       end if
       call mpi_bcast(cdd_coef,size(cdd_coef),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    end if
!
    call mpi_bcast(hnl ,size(hnl) ,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(knl ,size(knl) ,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(hnml,size(hnml),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(knml,size(knml),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(ps_type,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  END SUBROUTINE send_pseudopot_1


  SUBROUTINE ps_allocate_1(n_grd,n_orb)
    implicit none
    integer,intent(IN) :: n_grd,n_orb
    integer :: mg,mo
    real(8),allocatable :: vql_tmp(:,:),cdd_tmp(:,:),rad_tmp(:,:)
    real(8),allocatable :: cdc_tmp(:,:),viod_tmp(:,:,:)
    real(8),allocatable :: anorm_tmp(:,:),Rps_tmp(:,:),rab_tmp(:,:)
    integer,allocatable :: inorm_tmp(:,:),lo_tmp(:,:),NRps_tmp(:,:)
    integer,allocatable :: no_tmp(:,:)
    if ( .not.allocated(hnl) ) then
       allocate( hnl(3,0:2,Nelement_PP) ) ; hnl=0.0d0
       allocate( knl(3,1:2,Nelement_PP) ) ; knl=0.0d0
    end if
    if ( .not.allocated(hnml) ) then
       allocate( hnml(3,3,0:2,Nelement_PP) ) ; hnml=0.0d0
       allocate( knml(3,3,1:2,Nelement_PP) ) ; knml=0.0d0
    end if
    if ( max_psgrd==0 .or. max_psorb==0 ) then
       allocate( Mr(Nelement_PP)   ) ; Mr=0
       allocate( norb(Nelement_PP) ) ; norb=0
       allocate( Zps(Nelement_PP)  ) ; Zps=0.d0
       allocate( parloc(4,Nelement_PP)    ) ; parloc=0.d0
       allocate( anorm(n_orb,Nelement_PP) ) ; anorm=0.d0
       allocate( inorm(n_orb,Nelement_PP) ) ; inorm=0
       allocate( Rps(n_orb,Nelement_PP)   ) ; Rps=0.d0
       allocate( NRps(n_orb,Nelement_PP)  ) ; NRps=0
       allocate( lo(n_orb,Nelement_PP)    ) ; lo=0
       allocate( no(n_orb,Nelement_PP)    ) ; no=0
       allocate( vql(n_grd,Nelement_PP)   ) ; vql=0.d0
       allocate( cdd(n_grd,Nelement_PP)   ) ; cdd=0.d0
       allocate( cdc(n_grd,Nelement_PP)   ) ; cdc=0.d0
       allocate( rad(n_grd,Nelement_PP)   ) ; rad=0.d0
       allocate( rab(n_grd,Nelement_PP)   ) ; rab=0.d0
       allocate( viod(n_grd,n_orb,Nelement_PP) ) ; viod=0.d0
       allocate( Rcloc(Nelement_PP) ) ; Rcloc=0.0d0
       max_psgrd=n_grd
       max_psorb=n_orb
       return
    end if
    mg = max( max_psgrd, n_grd )
    mo = max( max_psorb, n_orb )
    if ( max_psgrd < mg ) then
       allocate( vql_tmp(mg,Nelement_PP) ) ; vql_tmp=0.0d0
       allocate( cdd_tmp(mg,Nelement_PP) ) ; cdd_tmp=0.0d0
       allocate( rad_tmp(mg,Nelement_PP) ) ; rad_tmp=0.0d0
       allocate( rab_tmp(mg,Nelement_PP) ) ; rab_tmp=0.0d0
       allocate( cdc_tmp(mg,Nelement_PP) ) ; cdc_tmp=0.0d0
       vql_tmp(1:max_psgrd,1:Nelement_PP) = vql(1:max_psgrd,1:Nelement_PP)
       cdd_tmp(1:max_psgrd,1:Nelement_PP) = cdd(1:max_psgrd,1:Nelement_PP)
       rad_tmp(1:max_psgrd,1:Nelement_PP) = rad(1:max_psgrd,1:Nelement_PP)
       rab_tmp(1:max_psgrd,1:Nelement_PP) = rab(1:max_psgrd,1:Nelement_PP)
       cdc_tmp(1:max_psgrd,1:Nelement_PP) = cdc(1:max_psgrd,1:Nelement_PP)
       deallocate( cdc )
       deallocate( rab )
       deallocate( rad )
       deallocate( cdd )
       deallocate( vql )
       allocate( vql(mg,Nelement_PP) ) ; vql=0.d0
       allocate( cdd(mg,Nelement_PP) ) ; cdd=0.d0
       allocate( rad(mg,Nelement_PP) ) ; rad=0.d0
       allocate( rab(mg,Nelement_PP) ) ; rab=0.d0
       allocate( cdc(mg,Nelement_PP) ) ; cdc=0.d0
       vql(:,:)=vql_tmp(:,:)
       cdd(:,:)=cdd_tmp(:,:)
       rad(:,:)=rad_tmp(:,:)
       rab(:,:)=rab_tmp(:,:)
       cdc(:,:)=cdc_tmp(:,:)
       deallocate( cdc_tmp )
       deallocate( rab_tmp )
       deallocate( rad_tmp )
       deallocate( cdd_tmp )
       deallocate( vql_tmp )
       allocate( viod_tmp(mg,mo,Nelement_PP) )
       viod_tmp(1:max_psgrd,1:max_psorb,1:Nelement_PP) &
            = viod(1:max_psgrd,1:max_psorb,1:Nelement_PP)
       deallocate( viod )
       allocate( viod(mg,mo,Nelement_PP) ) ; viod=0.d0
       viod(:,:,:)=viod_tmp(:,:,:)
       deallocate( viod_tmp )
    end if
    if ( max_psorb < mo ) then
       allocate( anorm_tmp(mo,Nelement_PP) ) ; anorm_tmp=0.0d0
       allocate( inorm_tmp(mo,Nelement_PP) ) ; inorm_tmp=0
       allocate( lo_tmp(mo,Nelement_PP) ) ; lo_tmp=0
       allocate( no_tmp(mo,Nelement_PP) ) ; no_tmp=0
       allocate( Rps_tmp(mo,Nelement_PP) ) ; Rps_tmp=0.0d0
       allocate( NRps_tmp(mo,Nelement_PP) ) ; NRps_tmp=0
       anorm_tmp(1:max_psorb,1:Nelement_PP) = anorm(1:max_psorb,1:Nelement_PP)
       inorm_tmp(1:max_psorb,1:Nelement_PP) = inorm(1:max_psorb,1:Nelement_PP)
       lo_tmp(1:max_psorb,1:Nelement_PP) = lo(1:max_psorb,1:Nelement_PP)
       no_tmp(1:max_psorb,1:Nelement_PP) = no(1:max_psorb,1:Nelement_PP)
       Rps_tmp(1:max_psorb,1:Nelement_PP) = Rps(1:max_psorb,1:Nelement_PP)
       NRps_tmp(1:max_psorb,1:Nelement_PP) = NRps(1:max_psorb,1:Nelement_PP)
       deallocate( NRps )
       deallocate( Rps )
       deallocate( no )
       deallocate( lo )
       deallocate( inorm )
       deallocate( anorm )
       allocate( anorm(mo,Nelement_PP) ) ; anorm=0.d0
       allocate( inorm(mo,Nelement_PP) ) ; inorm=0
       allocate( lo(mo,Nelement_PP)    ) ; lo=0
       allocate( no(mo,Nelement_PP)    ) ; no=0
       allocate( Rps(mo,Nelement_PP)   ) ; Rps=0.d0
       allocate( NRps(mo,Nelement_PP)  ) ; NRps=0
       anorm(:,:) = anorm_tmp(:,:)
       inorm(:,:) = inorm_tmp(:,:)
       lo(:,:) = lo_tmp(:,:)
       no(:,:) = no_tmp(:,:)
       Rps(:,:) = Rps_tmp(:,:)
       NRps(:,:) = NRps_tmp(:,:)
       deallocate( NRps_tmp )
       deallocate( Rps_tmp )
       deallocate( lo_tmp )
       deallocate( inorm_tmp )
       deallocate( anorm_tmp )
       if ( max_psgrd >= mg ) then
          allocate( viod_tmp(mg,mo,Nelement_PP) )
          viod_tmp(1:max_psgrd,1:max_psorb,1:Nelement_PP) &
               = viod(1:max_psgrd,1:max_psorb,1:Nelement_PP)
          deallocate( viod )
          allocate( viod(mg,mo,Nelement_PP) ) ; viod=0.d0
          viod(:,:,:)=viod_tmp(:,:,:)
          deallocate( viod_tmp )
       end if
    end if
    max_psgrd = mg
    max_psorb = mo
  END SUBROUTINE ps_allocate_1


END MODULE pseudopot_module
