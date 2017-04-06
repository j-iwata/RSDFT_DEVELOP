MODULE fermi_ncol_module

  use bberf_module

  implicit none

  PRIVATE
  PUBLIC :: calc_fermi_ncol
  PUBLIC :: get_parameters_fermi_ncol

  real(8) :: ekbt, efermi, Eentropy
  integer :: mb1,mb2,kinteg

  integer :: nsetocc,isetocc(2)
  real(8) :: setocc(20)

  logical :: first_time = .true.
  real(8),allocatable :: factor(:)

CONTAINS


  SUBROUTINE read_fermi(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i,n1,n2
    character(6) :: cbuf,ckey
    ekbt=1.d-5
    kinteg=5
    nsetocc=0
    isetocc=0
    setocc=0.0d0
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:4) == "EKBT" ) then
             backspace(unit)
             read(unit,*) cbuf,ekbt
          else if ( ckey(1:6) == "KINTEG" ) then
             backspace(unit)
             read(unit,*) cbuf,kinteg
          else if ( ckey(1:6) == "SETOCC" ) then
             backspace(unit)
             read(unit,*) cbuf,n1,n2,setocc(nsetocc+1:nsetocc+n2-n1+1)
             nsetocc=nsetocc+n2-n1+1
             if ( all( isetocc(:) == 0 ) ) then
                isetocc(1)=n1
                isetocc(2)=n2
             else
                if ( isetocc(1) /= n1 .or. isetocc(2) /= n2 ) then
                   write(*,*) "isetocc(:) are different !!!"
                   write(*,*) "isetocc(:)=",isetocc(:)
                   write(*,*) " /= n1,n2 =",n1,n2
                   stop "stop@read_fermi"
                end if
             end if
          end if
       end do
999    continue
       write(*,*) "ekbt   =",ekbt
       write(*,*) "kinteg =",kinteg
       write(*,*) "nsetocc=",nsetocc,setocc(1:nsetocc)
       write(*,*) "isetocc=",isetocc(1:2)
    end if
    call send_fermi(0)
  END SUBROUTINE read_fermi


  SUBROUTINE send_fermi(rank)
    implicit none
    integer,intent(IN) :: rank
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(ekbt,1,MPI_REAL8,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kinteg,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nsetocc,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(isetocc,2,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(setocc,nsetocc,MPI_REAL8,rank,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_fermi


  SUBROUTINE calc_fermi_ncol(iter,nfixed,MB,MBZ,MSP,znel,dspn,esp,wbz,occ)

    implicit none
    integer,intent(IN)  :: iter,nfixed,MB,MBZ,MSP
    real(8),intent(IN)  :: esp(MB,MBZ,MSP),wbz(MBZ),znel,dspn
    real(8),intent(OUT) :: occ(MB,MBZ,MSP)
    real(8) :: ef1,ef2,ef,ef0
    integer :: id,n,k,s,efconv,nnk,myrank,ierr
    real(8) :: zne,octmp,ff,xx
    real(8),parameter :: eps=0.d0
    integer,parameter :: mxcycl=1000
    logical :: disp_switch
    include 'mpif.h'

    call write_border( 1, " calc_fermi_ncol(start)" )

    call check_disp_switch( disp_switch, 0 )

    if ( first_time ) then
       call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )
       call read_fermi( myrank, 1 )
       first_time = .false.
       if ( kinteg > 0 ) then
          allocate( factor(kinteg) )
          factor(1)=-1.d0/(4.d0*sqrt(acos(-1.d0)))
          do n=2,kinteg
             factor(n)=-factor(n-1)/(4.d0*n)
          end do
       end if
    end if

    mb1 = 1
    mb2 = MB

! Set upper & lower boundarires of Fermi energy

    ef1 = minval( esp(mb1:mb2,1:MBZ,1) )
    ef2 = maxval( esp(mb1:mb2,1:MBZ,1) )
    if ( ef1 == ef2 ) then
       ef1 = ef1 - 0.01d0
       ef2 = ef2 + 0.01d0
    end if

!C Safety margin for highly degenerate systems & artificial fault
!C
    ef2 = ef2 + min( ekbt*1.d2, 0.1d0 )

    if ( .true. ) then

       zne = znel - 1.0d0*(mb1-1)

       ef0 = 1.d10

       do id=1,mxcycl

          ef = 0.5d0*( ef1 + ef2 )
          if ( ef == ef0 ) goto 100
          octmp = 0.0d0

          do s=1,1
          do k=1,MBZ
          do n=mb1,mb2

             xx=(esp(n,k,s)-ef)/ekbt
             ff=ff0(kinteg,xx)

             octmp = octmp + ff*wbz(k)
             occ(n,k,s)=ff

          end do
          end do
          end do

          if ( octmp-zne > eps ) then
             ef2=ef
          else if ( octmp-zne < -eps ) then
             ef1=ef
          else
             goto 100
          end if

          ef0 = ef

       end do ! id

    end if

    if ( abs(octmp-zne) > 1.d-10 ) then
       if ( disp_switch ) then
          write(6,*)' EF IS NOT CONVERGED'
          write(6,*)' Check the # of electron, mb1, and mb2'
          write(6,*)' EF1 & EF2=',ef1,ef2
          write(6,*)' octmp,zne=',octmp,zne,octmp-zne
          do s=1,1
             do k=1,MBZ
                write(*,*) "s,k =",s,k
                do n=mb1,mb2
                   write(*,*) n,occ(n,k,s),esp(n,k,s)
                end do
             end do
          end do
       end if
       stop'FERMI'
    end if

100 continue

    efermi   = ef
    Eentropy = 0.0d0

    do s=1,1
    do k=1,MBZ
       do n=1,mb1-1
          occ(n,k,s)=1.0d0*wbz(k)
       end do
       do n=mb1,mb2
          occ(n,k,s)=1.0d0*occ(n,k,s)*wbz(k)
       end do
       do n=mb2+1,MB
          occ(n,k,s)=0.0d0
       end do
    end do
    end do
    occ(:,:,MSP)=occ(:,:,1)

    call write_border( 1, " calc_fermi_mp(end)" )

  END SUBROUTINE calc_fermi_ncol

  FUNCTION ff0(n,x)
    implicit none
    integer :: n,i
    real(8) :: x,ff,ff0,hp0,hp1,hp2,hp3
    ff0 = 0.5d0*(1.d0-bberf(x))
    if ( n <= 0 ) return
    hp0 = 1.d0
    hp1 = 2.d0*x
    ff  = factor(1)*hp1
    do i=2,n
       hp2 = 2.d0*x*hp1 - 2.d0*(2*i-3)*hp0
       hp3 = 2.d0*x*hp2 - 2.d0*(2*i-2)*hp1
       ff  = ff + factor(i)*hp3
       hp0 = hp2
       hp1 = hp3
    end do
    ff0 = ff0 + ff*exp(-x*x)
    return
  END FUNCTION ff0


  SUBROUTINE get_parameters_fermi_ncol &
       ( ekbt_out, kinteg_out, efermi_out, Eentropy_out )
    implicit none
    real(8),optional,intent(OUT) :: ekbt_out, efermi_out, Eentropy_out
    integer,optional,intent(OUT) :: kinteg_out
    if ( present(ekbt_out)     ) ekbt_out=ekbt
    if ( present(efermi_out)   ) efermi_out=efermi
    if ( present(Eentropy_out) ) Eentropy_out=Eentropy
    if ( present(kinteg_out)   ) kinteg_out=kinteg
  END SUBROUTINE get_parameters_fermi_ncol


END MODULE fermi_ncol_module
