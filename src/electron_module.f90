module electron_module

  use io_tools_module

  implicit none

  private
  integer,public :: Nband
  real(8),public :: Nelectron
  integer,public :: Nspin
  real(8),public :: Nelectron_spin(2)
  real(8),public :: Ndspin
  integer,public :: Nfixed

  public :: read_electron
  public :: count_electron
  public :: check_Nband_electron
  public :: get_atomic_site_spin_diff_electron

  integer,parameter :: unit_sc=980
  real(8),allocatable :: dspin(:,:)
  real(8) :: Next_electron

contains

  subroutine read_electron
    implicit none
    integer :: i
    character(6) :: cbuf,ckey
    logical :: disp
    call write_border( 0, " read_electron(start)" )
    call check_disp_switch( disp, 0 )
    Nband=0
    Nspin=1
    Nfixed=0
    Ndspin=0.0d0
    Next_electron=0.0d0
    call IOTools_readIntegerKeyword('NBAND' ,Nband)
    call IOTools_readIntegerKeyword('NSPIN' ,Nspin)
    call IOTools_readIntegerKeyword('NFIXED',Nfixed)
    call IOTools_readReal8Keyword('NDSPIN',Ndspin)
    call IOTools_readReal8Keyword('NEXTE' ,Next_electron)
    if( disp )then
      write(*,*) "# of Bnads: Nband=",Nband
      write(*,*) "# of spin degeree of freedom: Nspin=",Nspin
      write(*,*) "Absolute # of spin difference: Ndspin=",Ndspin
      write(*,*) "(Spin configuration is read from fort.980 if Ndspin is negative)"
      write(*,*) "# of iteration steps keeping Ndspin: Nfixed=",Nfixed
      write(*,*) "Extra # of electrons: Next_electron=",Next_electron
    end if
    if ( Nspin == 1 .and. Ndspin /= 0.0d0 ) then
      Ndspin=0.0d0
      if ( disp ) write(*,*) "Ndspin is replaced to 0.0: Ndspin=",Ndspin
    end if
    if ( Ndspin < 0.0d0 ) call Read_SpinConf( dspin )
    call write_border( 0, " read_electron(end)" )
  end subroutine read_electron


  subroutine Read_SpinConf( dspin )
    implicit none
    real(8),allocatable,intent(inout) :: dspin(:,:)
    integer :: ndat,i1,i2,i,ierr,myrank
    real(8) :: diff_ele
    include 'mpif.h'
    call write_border( 0, " Read_SpinConf(start)" )
    call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr )
    if ( myrank == 0 ) then
      rewind unit_sc
      ndat=0
      do
        read(unit_sc,*,END=90) i1,i2,diff_ele
        ndat = ndat + 1
      end do
90    continue
    end if
    call MPI_BCAST(ndat,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)      
    allocate( dspin(3,ndat) ); dspin=0.0d0
    if ( myrank == 0 ) then
      rewind unit_sc
      do i = 1, ndat
        read(unit_sc,*) dspin(1:3,i)
      end do
    end if
    call MPI_BCAST(dspin,size(dspin),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    Ndspin=sum( dspin(3,:) )
    if ( myrank == 0 ) write(*,*) "Ndspin is replaced: Ndspin=",Ndspin
    call write_border( 0, " Read_SpinConf(end)" )
  end subroutine Read_SpinConf


  subroutine get_atomic_site_spin_diff_electron( d, iatom )
    implicit none
    real(8),intent(out) :: d
    integer,intent(in)  :: iatom
    integer :: i, i1, i2
    d=0.0d0
    if ( Nspin == 1 .or. Ndspin == 0.0d0 .or. .not.allocated(dspin) ) return
    do i = 1, size(dspin,2)
      i1 = nint( dspin(1,i) )
      i2 = nint( dspin(2,i) )
      if ( i1 <= iatom .and. iatom <= i2 ) then
        d = dspin(3,i)
        return
      end if
    end do
  end subroutine get_atomic_site_spin_diff_electron


  subroutine count_electron( Zps, element_id )
    implicit none
    real(8),intent(in) :: Zps(:)
    integer,intent(in) :: element_id(:)
    integer :: iatom,natom
    call write_border( 0, " count_electron(start)" )
    natom=size(element_id)
    Nelectron = 0.0d0
    do iatom = 1, natom
       Nelectron = Nelectron + Zps( element_id(iatom) )
    end do
    Nelectron = Nelectron + Next_electron
    Nelectron_spin(1) = 0.5d0*Nelectron + 0.5d0*Ndspin
    Nelectron_spin(2) = 0.5d0*Nelectron - 0.5d0*Ndspin
    call write_border( 0, " count_electron(end)" )
  end subroutine count_electron


  subroutine check_Nband_electron
    implicit none
    real(8),parameter :: factor=1.5d0
    character(30) :: mesg
    call write_border( 0, " check_Nband_electron(start)" )
    if ( dble(Nband) < 0.5d0*Nelectron ) then
       Nband = nint( 0.5d0*Nelectron * factor )
       Nband = max( Nband, 8 )
       write(mesg,'(1x,"Nband is replaced to ",i8)') Nband
       call write_string( mesg )
    end if
    call write_border( 0, " check_Nband_electron(end)" )
  end subroutine check_Nband_electron


end module electron_module
