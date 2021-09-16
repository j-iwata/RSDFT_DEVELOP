module psv_initrho_module

  use pseudopot_module, only: cdd_coef, ippform
  use simc_module, only: fit_initrho_simc

  implicit none

  private
  public :: read_coef_psv_initrho

contains

  subroutine read_coef_psv_initrho
    implicit none
    integer :: unit_ps,ierr,Nelement,myrank
    integer :: ielm,ngauss,max_ngauss,i,n
    character(18) :: inbuf18
    logical :: flag
    real(8),allocatable :: r(:), d(:)
    real(8) :: pi4
    include 'mpif.h'

    call write_border( 0, " read_coef_psv_initrho(start)" )

    if ( all(ippform==2) ) then
      call write_border( 0, " read_coef_psv_initrho(return)" )
      return
    end if

    call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr )

    pi4 = 4.0d0*acos(-1.0d0)

    Nelement = size( ippform )

    max_ngauss = 0

    if ( myrank == 0 ) then

      unit_ps = 33

      do ielm=1,Nelement

        unit_ps = unit_ps + 1
        rewind unit_ps

        flag=.true.
        do
          read(unit_ps,'(A)',END=99) inbuf18
          if ( inbuf18 == '### initial charge' ) then
            read(unit_ps,*) ngauss
            max_ngauss = max( ngauss, max_ngauss )
            flag = .false.
            exit
          end if
        end do
        99 continue

        if ( flag ) max_ngauss = max( 4, max_ngauss )

      end do ! ielm

    end if

    call MPI_Bcast(max_ngauss,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    if ( .not.allocated(cdd_coef) ) allocate( cdd_coef(3,max_ngauss,Nelement) )
    cdd_coef=0.0d0

    if ( myrank == 0 ) then

      unit_ps = 33

      do ielm=1,Nelement

        unit_ps = unit_ps + 1
        rewind unit_ps

        flag = .true.
        do
          read(unit_ps,'(A)',END=999) inbuf18
          if ( inbuf18 == '### initial charge' ) then
            read(unit_ps,*) ngauss
            write(*,'(1x,"ielm,ngauss=",2i5)') ielm,ngauss
            do i=1,ngauss
              read(unit_ps,*) cdd_coef(1:3,i,ielm)
              write(*,'(1x,i5,3g20.10)') i,cdd_coef(1:3,i,ielm)
            end do
            flag = .false.
            exit
          end if
        end do
        999 continue

        if ( flag ) then
          rewind unit_ps
          read(unit_ps,*) n
          allocate( r(n) ); r=0.0d0
          allocate( d(n) ); d=0.0d0
          do i = 1, n
            read(unit_ps,*) r(i), d(i)
          end do
          d = r**2 * d * pi4
          call fit_initrho_simc( r, d, cdd_coef(:,:,ielm) )
          deallocate( d )
          deallocate( r )
        end if

      end do! ielm

    end if

    call MPI_Bcast( cdd_coef, size(cdd_coef), MPI_REAL8, 0, MPI_COMM_WORLD, ierr )

    call write_border( 0, " read_coef_psv_initrho(end)" )

  end subroutine read_coef_psv_initrho

  function psv_initrho( r, c )
    implicit none
    real(8) :: psv_initrho
    real(8),intent(in) :: r,c(:,:)
    real(8) :: pi4,r2
    integer :: j
    pi4=4.0d0*acos(-1.0d0)
    psv_initrho=0.0d0
    r2=r*r
    do j = 1, size(c,2)
      psv_initrho = psv_initrho + r2*( c(1,j) + c(2,j)*r2 )*exp( -c(3,j)*r2 )
    end do
    psv_initrho = pi4*psv_initrho
  end function psv_initrho

end module psv_initrho_module
