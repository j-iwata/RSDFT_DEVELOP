MODULE pseudopot_module

  use ps_read_YB_module
  use VarPSMember
  use PSreadPSV
#ifdef _USPP_
  use VarPSMemberG
  use PSReadPSVG
#endif

  implicit none

  PRIVATE
  PUBLIC :: pselect,ippform,file_ps,inorm,NRps,norb,Npseudopot &
           ,Mr,lo,vql,cdd,cdc,rad,anorm,viod,Rps,Zps,parloc,rab &
           ,read_ppname_pseudopot,read_pseudopot &
           ,cdd_coef,read_param_pseudopot,read_param_oldformat_pseudopot &
           ,read_ppname_oldformat_pseudopot

  integer :: pselect,Npseudopot
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
    end do
    if ( ierr == -1 ) stop "stop@read_ppname_pseudopot"
    Nelement_PP=MKI
    Npseudopot =MKI
  END SUBROUTINE read_ppname_pseudopot

!-------------------------------------------------------
  SUBROUTINE read_ppname_oldformat_pseudopot(MKI,rank,unit)
    implicit none
    integer,intent(IN) :: MKI,rank,unit
    integer :: i
    Nelement_PP=MKI
    Npseudopot =MKI
    allocate( ippform(Nelement_PP),file_ps(Nelement_PP) )
    if ( rank == 0 ) then
       do i=1,Nelement_PP
          read(unit,*) ippform(i),file_ps(i)
       end do
       do i=1,Nelement_PP
          write(*,'(1x,"ippform, file_ps = ",i3,2x,a30,3x,3f10.5)') &
               ippform(i),file_ps(i)
       end do
    end if
    call send_ppname_2(0)
  END SUBROUTINE read_ppname_oldformat_pseudopot

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
  SUBROUTINE read_param_pseudopot(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    character(7) :: cbuf,ckey
    integer :: i
    pselect=2
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:7) == "PSELECT" ) then
             backspace(unit)
             read(unit,*) cbuf,pselect
          end if
       end do
999    continue
       write(*,*) "pslect=",pselect
    end if
    call send_param_pseudopot(0)
  END SUBROUTINE read_param_pseudopot

!-------------------------------------------------------
  SUBROUTINE read_param_oldformat_pseudopot(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    if ( rank == 0 ) then
       read(unit,*) pselect
       write(*,*) "pselect=",pselect
    end if
    call send_param_pseudopot(0)
  END SUBROUTINE read_param_oldformat_pseudopot

  SUBROUTINE send_param_pseudopot(rank)
    implicit none
    integer,intent(IN) :: rank
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(pselect,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_param_pseudopot

!-------------------------------------------------------
  SUBROUTINE read_pseudopot(rank)
    implicit none
    integer,intent(IN) :: rank
    real(8),allocatable :: psi_(:,:,:),phi_(:,:,:),bet_(:,:,:)
    real(8),allocatable :: ddi_(:,:,:),qqr_(:,:,:)

    if ( rank == 0 ) then
       max_psgrd=0
       max_psorb=0
       do ielm=1,Nelement_PP
          unit_ps=33+ielm
          select case(ippform(ielm))
          case(1)
             call read_TM
          case(2)
             open(unit_ps,FILE=file_ps(ielm),STATUS='old')
             call read_PSV( unit_ps,ielm,ddi_,qqr_,psi_,phi_,bet_ )
             close(unit_ps)
          case(3)
             open(unit_ps,FILE=file_ps(ielm),STATUS='old')
             call ps_read_YB(unit_ps)
             close(unit_ps)
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
#ifdef _USPP_
          case(102)
            open(unit_ps,FILE=file_ps(ielm),STATUS='old')
            call read_PSV( unit_ps,ielm,ddi_,qqr_,psi_,phi_,bet_ )
            write(*,*) 'normal PSV finished'
            call readPSVG( unit_ps,ielm,ddi_,qqr_,psi_,phi_,bet_ )
            write(*,*) 'new PSV finished'
            close(unit_ps)

#endif
          case default
             stop "ippform error"
          end select
       end do
    end if

    if (pselect==2) then
      call send_pseudopot(rank)
#ifdef _USPP_
    elseif (pselect==102) then
      call send_pseudopot(rank)
      call sendPSG(rank,Nelement_PP)
#endif
    else
      stop 'pselect must = 2(NCPP),102(USPP)'
    end if
  END SUBROUTINE read_pseudopot


!-------------------------------------------------------
  SUBROUTINE read_TM
    implicit none
    stop
  END SUBROUTINE read_TM

END MODULE pseudopot_module
