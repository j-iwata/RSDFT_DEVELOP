MODULE psv_initrho_module

  use pseudopot_module, only: cdd_coef, ippform, file_ps

  implicit none

  PRIVATE
  PUBLIC :: read_psv_initrho

  integer :: max_ngauss
  integer :: iswitch_initrho = 0

  include 'mpif.h'

CONTAINS


  SUBROUTINE read_psv_initrho( Nelement, rank, unit, info )
    implicit none
    integer,intent(IN) :: Nelement, rank, unit
    integer,intent(OUT) :: info
    integer :: i, ierr
    character(7) :: cbuf, ckey
    iswitch_initrho = 0
    if ( rank == 0 ) then
       write(*,'(a50," psv_initrho")') repeat("-",50)
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey == "INITRHO" ) then
             backspace(unit)
             read(unit,*) cbuf, iswitch_initrho
             exit
          end if
       end do
999    continue
       write(*,*) "iswitch_initrho=",iswitch_initrho
    end if
    call mpi_bcast(iswitch_initrho,1,mpi_integer,0,mpi_comm_world,ierr)

    if ( iswitch_initrho == 2 ) then
       call read_coef_psv_initrho( Nelement, rank )
    end if

    info = iswitch_initrho

  END SUBROUTINE read_psv_initrho


  SUBROUTINE read_coef_psv_initrho( Nelement, rank )
    implicit none
    integer,intent(IN) :: Nelement, rank
    integer,parameter :: max_loop=1000000
    integer :: unit_ps,ierr
    integer :: ielm,loop,ngauss,i
    character(18) :: inbuf18

    max_ngauss = 0

    if ( rank == 0 ) then

       unit_ps = 34
       do ielm=1,Nelement
          if ( ippform(ielm) == 2 ) open(unit_ps,file=file_ps(ielm),status="old")
          rewind unit_ps
          do loop=1,max_loop
             read(unit_ps,'(A)') inbuf18
             if ( inbuf18 == '### initial charge' ) then
                read(unit_ps,*) ngauss
                max_ngauss = max( ngauss, max_ngauss )
                exit
             end if
          end do ! loop
          if ( ippform(ielm) == 2 ) close(unit_ps)
          unit_ps = unit_ps + 1
       end do ! ielm

    end if

    call mpi_bcast(max_ngauss,1,mpi_integer,0,mpi_comm_world,ierr)

    if ( .not.allocated(cdd_coef) ) allocate( cdd_coef(3,max_ngauss,Nelement) )
    cdd_coef=0.0d0

    if ( rank == 0 ) then

       unit_ps = 34
       do ielm=1,Nelement
          if ( ippform(ielm) == 2 ) open(unit_ps,file=file_ps(ielm),status="old")
          rewind unit_ps
          do loop=1,max_loop
             read(unit_ps,'(A)') inbuf18
             if ( inbuf18 == '### initial charge' ) then
                read(unit_ps,*) ngauss
                write(*,'(1x,"ielm,ngauss=",2i5)') ielm,ngauss
                do i=1,ngauss
                   read(unit_ps,*) cdd_coef(1:3,i,ielm)
                   write(*,'(1x,i5,3g20.10)') i,cdd_coef(1:3,i,ielm)
                end do
                exit
             end if
          end do ! loop
          if ( ippform(ielm) == 2 ) close(unit_ps)
          unit_ps = unit_ps + 1
       end do! ielm

    end if

    i=size(cdd_coef)
    call mpi_bcast(cdd_coef,i,mpi_real8,0,mpi_comm_world,ierr)

  END SUBROUTINE read_coef_psv_initrho


END MODULE psv_initrho_module
