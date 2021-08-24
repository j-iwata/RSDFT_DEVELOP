module hamiltonian_module

  use kinetic_module, only: op_kinetic
  use localpot_module, only: op_localpot, Vloc
  use nonlocal_module, only: op_nonlocal
  use fock_module, only: op_fock
  use watch_module, only: watchb, watchb_omp,time_hmlt

  implicit none

  private
  public :: hamiltonian
  public :: op_kinetic
  public :: op_localpot
  public :: op_nonlocal
!  public :: op_hamiltonian
!  public :: init_op_hamiltonian
!  public :: send_status_to_hamiltonian
!  public :: reset_op_hamiltonian
!  public :: backup_op_hamiltonian
!  public :: prep_work_op_hamiltonian
!  public :: release_work_op_hamiltonian
!  public :: rotv_op_hamiltonian

#ifdef _DRSDFT_
  real(8), public, allocatable :: hpsi_backup(:,:,:,:)
#else
  complex(8), public, allocatable :: hpsi_backup(:,:,:,:)
#endif

  interface op_hamiltonian
    module procedure d_op_hamiltonian, z_op_hamiltonian
  end interface

  interface backup_op_hamiltonian
    module procedure d_backup_op_hamiltonian, z_backup_op_hamiltonian
  end interface

  interface rotv_op_hamiltonian
    module procedure d_rotv_op_hamiltonian, z_rotv_op_hamiltonian
  end interface

  interface hamiltonian
    module procedure d_hamiltonian1, z_hamiltonian1
  end interface

  logical :: USE_BACKUP_AND_RESTORE=.false.
  logical :: BACKUP_AVAILABLE=.false.
  logical :: VLOC_IS_MODIFIED=.false.
  real(8),allocatable :: vloc_backup(:,:)
  real(8),allocatable :: dvloc(:,:)

#ifdef _DRSDFT_
  real(8), allocatable :: hpsi_work(:,:)
#else
  complex(8), allocatable :: hpsi_work(:,:)
#endif

contains

  subroutine hamiltonian0(n,k,s,tpsi,htpsi,n1,n2,ib1,ib2)
    implicit none
    integer,intent(in) :: n,k,s,n1,n2,ib1,ib2
#ifdef _DRSDFT_
    real(8),intent(in)  :: tpsi(n1:,ib1:)
    real(8),intent(out) :: htpsi(n1:,ib1:)
    real(8),parameter :: zero=0.0d0
#else
    complex(8),intent(in)  :: tpsi(n1:,ib1:)
    complex(8),intent(out) :: htpsi(n1:,ib1:)
    complex(8),parameter :: zero=(0.0d0,0.0d0)
#endif
    real(8) :: ttmp(2)

!$omp parallel

!$omp workshare
    htpsi=zero
!$omp end workshare

    !call watchb_omp( ttmp )

! --- Kinetic energy ---

    call op_kinetic( tpsi, htpsi, k, Vloc(:,s) )

!$omp barrier

    !call watchb_omp( ttmp, time_hmlt(1,1) )

! --- local potential ---

    !call op_localpot( tpsi, htpsi, s )

!$omp barrier

    !call watchb_omp( ttmp, time_hmlt(1,2) )

! --- nonlocal potential ---

    call op_nonlocal( tpsi, htpsi, n,k,s )

!$omp barrier

    !call watchb_omp( ttmp, time_hmlt(1,3) )

!$omp end parallel

    !call watchb( ttmp )

    call op_fock( tpsi, htpsi, n,k,s )

    !call watchb( ttmp, time_hmlt(1,4) )

  end subroutine hamiltonian0

  subroutine d_hamiltonian1( tpsi, htpsi, n,k,s )
    implicit none
    integer,intent(in) :: n,k,s
    real(8),intent(in) :: tpsi(:,:)
    real(8),intent(inout) :: htpsi(:,:)
!$omp parallel
    !call watchb_omp( ttmp )
    call op_kinetic( tpsi, htpsi, k, Vloc(:,s) )
!$omp barrier
    !call watchb_omp( ttmp, time_hmlt(1,1) )
    !call op_localpot( tpsi, htpsi, s )
!!$omp barrier
    !call watchb_omp( ttmp, time_hmlt(1,2) )
    call op_nonlocal( tpsi, htpsi, n,k,s )
!$omp barrier
    !call watchb_omp( ttmp, time_hmlt(1,3) )
!$omp end parallel
    !call watchb( ttmp )
    call op_fock( tpsi, htpsi, n,k,s )
    !call watchb( ttmp, time_hmlt(1,4) )
  end subroutine d_hamiltonian1

  subroutine z_hamiltonian1( tpsi, htpsi, n,k,s )
    implicit none
    integer,intent(in) :: n,k,s
    complex(8),intent(in) :: tpsi(:,:)
    complex(8),intent(inout) :: htpsi(:,:)
    !$omp parallel
    !call watchb_omp( ttmp )
    call op_kinetic( tpsi, htpsi, k, Vloc(:,s) )
    !$omp barrier
    !call watchb_omp( ttmp, time_hmlt(1,1) )
    !call op_localpot( tpsi, htpsi, s )
    !!$omp barrier
    !call watchb_omp( ttmp, time_hmlt(1,2) )
    call op_nonlocal( tpsi, htpsi, n,k,s )
    !$omp barrier
    !call watchb_omp( ttmp, time_hmlt(1,3) )
    !$omp end parallel
    !call watchb( ttmp )
    call op_fock( tpsi, htpsi, n,k,s )
    !call watchb( ttmp, time_hmlt(1,4) )
  end subroutine z_hamiltonian1


  subroutine init_op_hamiltonian( lswitch, ng,nb,nk,ns )
    implicit none
    logical, intent(in) :: lswitch
    integer, intent(in) :: ng, nb, nk(2), ns(2)
    call write_border( 0, "init_op_hamiltonian(start)" )
    USE_BACKUP_AND_RESTORE = lswitch
    if( USE_BACKUP_AND_RESTORE )then
       allocate( hpsi_backup(ng,nb,nk(1):nk(2),ns(1):ns(2)) )
       hpsi_backup=(0.0d0,0.0d0)
       allocate( vloc_backup(ng,ns(1):ns(2)) )
       vloc_backup=0.0d0
    end if
    call write_border( 0, "init_op_hamiltonian(end)" )
  end subroutine init_op_hamiltonian


  subroutine d_op_hamiltonian( n,k,s, tpsi, hpsi, backup, restore )
    implicit none
    integer,intent(in)  :: n,k,s
    real(8),intent(in)  :: tpsi(:,:)
    real(8),intent(out) :: hpsi(:,:)
    logical,optional,intent(in) :: backup
    logical,optional,intent(in) :: restore
    real(8) :: ttmp(2)
    integer :: n1,n2,ib,jb
#ifdef _DRSDFT_
    if( USE_BACKUP_AND_RESTORE .and. present(restore) )then
       if( restore .and. BACKUP_AVAILABLE )then
          n1=n
          n2=n+size(hpsi,2)-1
          if( VLOC_IS_MODIFIED )then
             do ib = n1, n2
                jb = ib - n + 1 
                !$omp parallel workshare
                hpsi(:,jb)=hpsi_backup(:,ib,k,s)+dvloc(:,s)*tpsi(:,jb)
                !$omp end parallel workshare
             end do
          else
             !$omp parallel workshare
             hpsi(:,:) = hpsi_backup(:,n1:n2,k,s)
             !$omp end parallel workshare
          end if
          return
       end if
    end if

!$omp parallel

!$omp workshare
    !hpsi=0.0d0
!$omp end workshare

    !call watchb_omp( ttmp )

! --- Kinetic energy ---

    call op_kinetic( tpsi, hpsi, k, Vloc(:,s) )

!$omp barrier

    !call watchb_omp( ttmp, time_hmlt(1,1) )

! --- local potential ---

    !call op_localpot( tpsi, hpsi, s )

!$omp barrier

    !call watchb_omp( ttmp, time_hmlt(1,2) )

! --- nonlocal potential ---

    call op_nonlocal( tpsi, hpsi, k, s )

!$omp barrier

    !call watchb_omp( ttmp, time_hmlt(1,3) )

!$omp end parallel

    !call watchb( ttmp )

    call op_fock( tpsi, hpsi, n,k,s )

    !call watchb( ttmp, time_hmlt(1,4) )

    if( USE_BACKUP_AND_RESTORE .and. present(backup) )then
       if ( backup )then
          n1=n
          n2=n+size(hpsi,2)-1
!$omp parallel workshare
          hpsi_backup(:,n1:n2,k,s) = hpsi(:,:)
!$omp end parallel workshare
          BACKUP_AVAILABLE = .true.
          if( VLOC_IS_MODIFIED .eqv. .true. )then
             vloc_backup(:,:) = Vloc(:,:)
             VLOC_IS_MODIFIED = .false.
          end if
       end if
    end if

    !call watchb( ttmp, time_hmlt(1,5) )

#else
    hpsi=0.0d0
#endif

  end subroutine d_op_hamiltonian

  subroutine z_op_hamiltonian( n,k,s, tpsi, hpsi, backup, restore )
    implicit none
    integer,intent(in) :: n,k,s
    complex(8),intent(in)  :: tpsi(:,:)
    complex(8),intent(out) :: hpsi(:,:)
    logical,optional,intent(in) :: backup
    logical,optional,intent(in) :: restore
    real(8) :: ttmp(2)
    integer :: n1,n2,ib,jb
#ifndef _DRSDFT_
    if( USE_BACKUP_AND_RESTORE .and. present(restore) )then
       if ( restore .and. BACKUP_AVAILABLE )then
          n1=n
          n2=n+size(hpsi,2)-1
          if( VLOC_IS_MODIFIED )then
             do ib = n1, n2
                jb = ib - n + 1 
                !$omp parallel workshare
                hpsi(:,jb)=hpsi_backup(:,ib,k,s)+dvloc(:,s)*tpsi(:,jb)
                !$omp end parallel workshare
             end do
          else
             !$omp parallel workshare
             hpsi(:,:) = hpsi_backup(:,n1:n2,k,s)
             !$omp end parallel workshare
          end if
          return
       end if
    end if

!$omp parallel

!$omp workshare
    !hpsi=(0.0d0,0.0d0)
!$omp end workshare

    !call watchb_omp( ttmp )

! --- Kinetic energy ---

    call op_kinetic( tpsi, hpsi, k, Vloc(:,s) )

!$omp barrier

    !call watchb_omp( ttmp, time_hmlt(1,1) )

! --- local potential ---

    !call op_localpot( tpsi, hpsi, s )

!$omp barrier

    !call watchb_omp( ttmp, time_hmlt(1,2) )

! --- nonlocal potential ---

    call op_nonlocal( tpsi, hpsi, k, s )

!$omp barrier

    !call watchb_omp( ttmp, time_hmlt(1,3) )

!$omp end parallel

    !call watchb( ttmp )

    call op_fock( tpsi, hpsi, n,k,s )

    !call watchb( ttmp, time_hmlt(1,4) )

    if( USE_BACKUP_AND_RESTORE .and. present(backup) )then
       if( backup )then
          n1=n
          n2=n+size(hpsi,2)-1
          !$omp parallel workshare
          hpsi_backup(:,n1:n2,k,s) = hpsi(:,:)
          !$omp end parallel workshare
          BACKUP_AVAILABLE = .true.
          if( VLOC_IS_MODIFIED .eqv. .true. )then
             vloc_backup(:,:) = Vloc(:,:)
             VLOC_IS_MODIFIED = .false.
          end if
       end if
    end if

    !call watchb( ttmp, time_hmlt(1,5) )

#else
    hpsi=(0.0d0,0.0d0)
#endif
  end subroutine z_op_hamiltonian

  subroutine send_status_to_hamiltonian( msg )
    implicit none
    character(*),intent(in) :: msg
    integer :: a,b,c
    if ( .not.USE_BACKUP_AND_RESTORE ) return
    select case( msg )
    case( "vloc is modified" )
       if( .not.allocated(dvloc) )then
          a=size(Vloc,1)
          b=lbound(Vloc,2)
          c=ubound(Vloc,2)
          allocate( dvloc(a,b:c) ); dvloc=0.0d0
       end if
       VLOC_IS_MODIFIED=.true.
       dvloc = Vloc - vloc_backup
       vloc_backup = Vloc
    case( "wf is modified" )
       BACKUP_AVAILABLE=.false.
    case default
       write(*,*) "msg= ",msg
       call stop_program("stop@send_status_to_hamiltonian")
    end select
  end subroutine send_status_to_hamiltonian

  subroutine reset_op_hamiltonian
    implicit none
    BACKUP_AVAILABLE=.false.
    VLOC_IS_MODIFIED=.false.
    if ( allocated(hpsi_backup) ) hpsi_backup=(0.0d0,0.0d0)
    if ( allocated(vloc_backup) ) vloc_backup=0.0d0
    if ( allocated(dvloc) ) dvloc=0.0d0
  end subroutine reset_op_hamiltonian

  subroutine d_backup_op_hamiltonian( n,k,s,hpsi )
    implicit none
    integer, intent(in) :: n, k, s
    real(8), intent(in) :: hpsi(:)
    if ( .not.USE_BACKUP_AND_RESTORE ) return
    !$omp parallel workshare
    hpsi_backup(:,n,k,s) = hpsi(:)
    !$omp end parallel workshare
    BACKUP_AVAILABLE = .true.
  end subroutine d_backup_op_hamiltonian

  subroutine z_backup_op_hamiltonian( n,k,s,hpsi )
    implicit none
    integer, intent(in) :: n, k, s
    complex(8), intent(in) :: hpsi(:)
    if ( .not.USE_BACKUP_AND_RESTORE ) return
    !$omp parallel workshare
    hpsi_backup(:,n,k,s) = hpsi(:)
    !$omp end parallel workshare
    BACKUP_AVAILABLE = .true.
  end subroutine z_backup_op_hamiltonian


  subroutine prep_work_op_hamiltonian( ml0, mb0, mb1 )
    implicit none
    integer,intent(in) :: ml0, mb0, mb1
    if ( .not.(USE_BACKUP_AND_RESTORE.and.BACKUP_AVAILABLE) ) return
    allocate( hpsi_work(ml0,mb0:mb1) ); hpsi_work=(0.0d0,0.0d0)
  end subroutine prep_work_op_hamiltonian

  subroutine release_work_op_hamiltonian( k, s )
    implicit none
    integer,intent(in) :: k, s
    integer :: mb0,mb1
    if ( .not.USE_BACKUP_AND_RESTORE ) return
    if ( allocated(hpsi_work) ) then
       mb0=lbound(hpsi_work,2)
       mb1=ubound(hpsi_work,2)
       hpsi_backup(:,mb0:mb1,k,s)=hpsi_work
       BACKUP_AVAILABLE=.true.
       deallocate(hpsi_work)
    end if
  end subroutine release_work_op_hamiltonian

  subroutine d_rotv_op_hamiltonian( k,s,i1,ii,ms,ns,RotMat )
    implicit none
    integer,intent(in) :: k,s,i1,ii,ms,ns
    real(8),intent(in) :: RotMat(:,:)
    real(8),parameter :: one=1.0d0
    integer :: mm,nn,ll
    if ( .not.(USE_BACKUP_AND_RESTORE.and.BACKUP_AVAILABLE) ) return
    call write_border( 0, "d_rotv_op_hamiltonian" )
    mm=size(RotMat,1)
    nn=size(RotMat,2)
    ll=size(hpsi_backup,1)
    call dgemm('N','N',ii,nn,mm,one,hpsi_backup(i1,ms,k,s),ll &
               ,RotMat,mm,one,hpsi_work(i1,ns),ll)
  end subroutine d_rotv_op_hamiltonian

  subroutine z_rotv_op_hamiltonian( k,s,i1,ii,ms,ns,RotMat )
    implicit none
    integer,intent(in) :: k,s,i1,ii,ms,ns
    complex(8),intent(in) :: RotMat(:,:)
    complex(8),parameter :: one=(1.0d0,0.0d0)
    integer :: mm,nn,ll
    if ( .not.(USE_BACKUP_AND_RESTORE.and.BACKUP_AVAILABLE) ) return
    mm=size(RotMat,1)
    nn=size(RotMat,2)
    ll=size(hpsi_backup,1)
    call zgemm('N','N',ii,nn,mm,one,hpsi_backup(i1,ms,k,s),ll &
               ,RotMat,mm,one,hpsi_work(i1,ns),ll)
  end subroutine z_rotv_op_hamiltonian

end module hamiltonian_module
