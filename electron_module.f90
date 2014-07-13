MODULE electron_module

  use atom_module, only: Natom,Nelement,ki_atom
  use pseudopot_module, only: Zps
  use bz_module, only: Nbzsm,weight_bz

  implicit none

  PRIVATE
  PUBLIC :: Nband,Nspin,Nelectron,Next_electron &
           ,Ndspin,Nfixed &
           ,read_electron,count_electron &
           ,read_oldformat_electron &
           ,NSelectron, dspin

  integer :: Nband, Nspin, Nfixed
  real(8) :: Nelectron,Next_electron,Ndspin
  integer :: MB_0,MB_1,MSP_0,MSP_1

  integer,parameter :: unit_sc=980
  real(8) :: NSelectron(2)
  real(8),allocatable :: dspin(:)

CONTAINS

  SUBROUTINE read_electron(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i
    character(6) :: cbuf,ckey
    Nband=0
    Nspin=1
    Nfixed=0
    Next_electron=0.d0
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:5) == "NBAND" ) then
             backspace(unit)
             read(unit,*) cbuf,Nband
          else if ( ckey(1:5) == "NSPIN" ) then
             backspace(unit)
             read(unit,*) cbuf,Nspin
          else if ( ckey(1:6) == "NDSPIN" ) then
             backspace(unit)
             read(unit,*) cbuf,Ndspin
          else if ( ckey(1:6) == "NFIXED" ) then
             backspace(unit)
             read(unit,*) cbuf,Nfixed
          else if ( ckey(1:6) == "NEXTE" ) then
             backspace(unit)
             read(unit,*) cbuf,Next_electron
          end if
       end do
999    continue
       write(*,*) "Nband=",Nband
       write(*,*) "Nspin=",Nspin
       write(*,*) "Ndspin=",Ndspin
       write(*,*) "Nfixed=",Nfixed
       write(*,*) "Next_electron=",Next_electron
       if ( Nspin == 1 .and. Ndspin /= 0.d0 ) then
          Ndspin=0.d0
          write(*,*) "Ndspin is replaced to 0.0: Ndspin=",Ndspin
       end if
    end if
    call send_electron(0)
    if ( Ndspin < 0.0d0 ) then
       allocate( dspin(Natom) ) ; dspin=0.0d0
       call Read_SpinConf(rank)
    end if
  END SUBROUTINE read_electron


  SUBROUTINE Read_SpinConf(rank)
    implicit none
    integer,intent(IN) :: rank
    integer :: i1,i2,i,ierr
    real(8) :: diff_ele
    include 'mpif.h'
    if ( rank == 0 ) then
       rewind unit_sc
       do i=1,Natom
          read(unit_sc,*,END=90) i1,i2,diff_ele
          dspin(i1:i2)=diff_ele
       end do
90     continue
    end if
    call MPI_BCAST(dspin,Natom,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    Ndspin=sum(dspin)
    if ( rank == 0 ) write(*,*) "Ndspin is replaced: Ndspin=",Ndspin
  END SUBROUTINE Read_SpinConf


  SUBROUTINE read_oldformat_electron(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    if ( rank == 0 ) then
       Nband=0
       Nspin=0
       read(unit,*) Nband, Nspin
       read(unit,*) Next_electron, Ndspin, Nfixed
       write(*,*) "Nband=",Nband
       write(*,*) "Nspin=",Nspin
       write(*,*) "Next_electron=",Next_electron
       write(*,*) "Ndspin=",Ndspin
       write(*,*) "Nfixed=",Nfixed
       if ( Nspin == 1 .and. Ndspin /= 0.d0 ) then
          Ndspin=0.d0
          write(*,*) "Ndspin is replaced to 0.0: Ndspin=",Ndspin
       end if
    end if
    call send_electron(0)
  END SUBROUTINE read_oldformat_electron

  SUBROUTINE send_electron(rank)
    integer,intent(IN) :: rank
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(Nband,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Nspin,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Next_electron,1,MPI_REAL8,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Ndspin,1,MPI_REAL8,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Nfixed,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_electron


  SUBROUTINE count_electron
    integer :: i
    Nelectron=0.d0
    do i=1,Natom
       Nelectron=Nelectron+Zps(ki_atom(i))
    end do
    Nelectron=Nelectron+Next_electron
    if ( nint(Nelectron)>2*Nband ) then
       write(*,*) "Nband is too small!!!"
       stop 'count_eletron'
    end if
    NSelectron(1)=0.5d0*Nelectron + 0.5d0*Ndspin
    NSelectron(2)=0.5d0*Nelectron - 0.5d0*Ndspin
  END SUBROUTINE count_electron


#ifdef TEST
  SUBROUTINE init_occupation
    integer :: n,k,s
    real(8) :: sum0
    allocate( occ(Nband,Nbzsm,Nspin) )
    occ=0.d0
    sum0=0.d0
    do n=1,Nband
       if ( sum0+2.d0 > Nelectron ) exit
       sum0=sum0+2.d0
       occ(n,1,s)=2.d0
    end do
    if ( sum0 < Nelectron ) then
       if ( n > Nband ) then
          stop "Nband is small (init_occupation)"
       else
          occ(n,1,1)=Nelectron-sum0
       end if
    end if
    do k=2,Nbzsm
       do n=1,Nband
          occ(n,k,1)=occ(n,1,1)*weight_bz(k)
       end do
    end do
    do n=1,Nband
       occ(n,1,1)=occ(n,1,1)*weight_bz(1)
    end do
    if ( Nspin > 1 ) then
       occ(:,:,1)=occ(:,:,1)/dble(Nspin)
       do s=2,Nspin
          occ(:,:,s)=occ(:,:,1)
       end do
    end if
  END SUBROUTINE init_occupation

  SUBROUTINE init_occ_electron
    integer :: n,k,s
    real(8) :: sum0,d,Nel
    allocate( occ(Nband,Nbzsm,Nspin) )
    occ=0.d0
    d=2.d0/Nspin
    do s=1,Nspin
       Nel = 0.5d0*d*Nelectron + (3-2*s)*0.5d0*Ndspin
       if ( Nel < 0.d0 ) stop "Ndspin is too large !!!"
       sum0=0.d0
       do n=1,Nband
          if ( sum0+d > Nel ) exit
          sum0=sum0+d
          occ(n,1,s)=d
       end do
       if ( sum0 < Nel ) then
          if ( n > Nband ) then
             stop "Nband is small (init_occupation)"
          else
             occ(n,1,s) = Nel-sum0
             write(*,*) Nel,sum0
          end if
       end if
       do k=2,Nbzsm
          do n=1,Nband
             occ(n,k,s)=occ(n,1,s)*weight_bz(k)
          end do
       end do
       do n=1,Nband
          occ(n,1,s)=occ(n,1,s)*weight_bz(1)
       end do
    end do
    sum0=sum(occ)
    if ( abs(sum0-Nelectron)>1.d-10 ) then
       write(*,'(1x,"sum(occ), Nelectron =",2g30.20)') sum(occ),Nelectron
       stop "sum(occ) /= Nelectron !!!"
    end if
  END SUBROUTINE init_occ_electron
#endif


END MODULE electron_module
