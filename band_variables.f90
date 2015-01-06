MODULE band_variables

  implicit none

  integer,parameter :: mnbk=100
  integer :: nbk,mb_band,mb2_band,maxiter_band
  integer :: nfki(mnbk),nskip_band
  real(8) :: ak(3,mnbk+1)
  real(8) :: esp_conv_tol=-1.d0

  integer,parameter :: unit_band_eigv=100
  integer,parameter :: unit_band_ovlp=110
  integer,parameter :: unit_band_dedk=120

CONTAINS

  SUBROUTINE read_band( rank, unit )
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i,ierr
    integer,parameter :: max_read=10000
    character(8) :: cbuf,ckey
    include 'mpif.h'
    esp_conv_tol = 1.d-5
    maxiter_band = 50
    nskip_band=0
    if ( rank == 0 ) then
       rewind unit
       do i=1,max_read
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:4) == 'BAND' ) then
             backspace(unit)
             read(unit,*) cbuf,nbk,mb_band,mb2_band,esp_conv_tol,maxiter_band
             read(unit,*) ak(1,1:nbk+1)
             read(unit,*) ak(2,1:nbk+1)
             read(unit,*) ak(3,1:nbk+1)
             read(unit,*) nfki(1:nbk)
          else if ( ckey(1:8) == 'SKIPBAND' ) then
             backspace(unit)
             read(unit,*) cbuf,nskip_band
          end if
       end do
999    continue
       write(*,*) "nbk     =",nbk
       write(*,*) "mb_band =",mb_band
       write(*,*) "mb2_band=",mb2_band
       if ( esp_conv_tol < 0.d0 ) esp_conv_tol=1.d-5
       write(*,*) "esp_von_tol=",esp_conv_tol
       write(*,*) "maxiter_band=",maxiter_band
       write(*,*) "nskip_band=",nskip_band
    end if
    call mpi_bcast(nbk,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(ak,3*mnbk,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(nfki,mnbk,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(mb_band,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(mb2_band,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(esp_conv_tol,1,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(maxiter_band,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(nskip_band,1,mpi_integer,0,mpi_comm_world,ierr)
  END SUBROUTINE read_band

END MODULE band_variables
