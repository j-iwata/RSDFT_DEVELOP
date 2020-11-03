module atomopt_io_module

  implicit none

  private
  public :: flag_continue_atomopt
  public :: read_atomopt_io
  public :: write_atomopt_io
  public :: read_diis_atomopt_io
  public :: write_diis_atomopt_io

  character(5),parameter :: version='ver01'
  character(8) :: file_name='wopt.dat'
  integer,parameter :: u=1

contains


  logical function flag_continue_atomopt( flag_in )
    implicit none
    logical,optional,intent(in) :: flag_in
    logical,save :: flag=.false.
    if ( present(flag_in) ) flag=flag_in
    flag_continue_atomopt=flag
  end function flag_continue_atomopt


  subroutine read_atomopt_io &
       ( loop_start, history, x, x0, g, g0, Hessian )
    implicit none
    integer,intent(inout) :: loop_start
    real(8),intent(inout) :: history(:,0:)
    real(8),intent(inout) :: x(:), x0(:), g(:), g0(:)
    real(8),intent(inout) :: Hessian(:,:)
    character(5) :: ver
    open(u,file=file_name,status='old',form='unformatted')
    read(u) ver
    if ( ver /= version ) call stop_program('wopt.dat: Invalid format !')
    read(u) loop_start
    read(u) history(:,0:loop_start-1)
    read(u) x
    read(u) x0
    read(u) g
    read(u) g0
    read(u) Hessian
    close(u)
  end subroutine read_atomopt_io


  subroutine write_atomopt_io &
       ( loop, history, x, x0, g, g0, Hessian )
    implicit none
    integer,intent(in) :: loop
    real(8),intent(in) :: history(:,0:)
    real(8),intent(in) :: x(:), x0(:), g(:), g0(:)
    real(8),intent(in) :: Hessian(:,:)
    open(u,file=file_name,form='unformatted')
    write(u) version
    write(u) loop
    write(u) history(:,:)
    write(u) x
    write(u) x0
    write(u) g
    write(u) g0
    write(u) Hessian
    close(u)
  end subroutine write_atomopt_io


  subroutine read_diis_atomopt_io &
       ( loop_start, history, ip, x, g, rho )
    implicit none
    integer,intent(inout) :: loop_start
    real(8),intent(inout) :: history(:,0:)
    integer,intent(inout) :: ip
    real(8),intent(inout) :: x(:,:,0:), g(:,:,0:)
    real(8),optional,intent(inout) :: rho(:,:,0:)
    character(5) :: ver
    integer :: jp,kp
    open(u,file=file_name,status='old',form='unformatted')
    read(u) ver
    if ( ver /= version ) call stop_program('wopt.dat: Invalid format !')
    read(u) loop_start
    read(u) history(:,0:loop_start-1)
    read(u) ip
    jp = ubound(x,3)
    if ( jp /= ip ) kp=min(ip,jp)
    read(u) x(:,:,0:kp)
    read(u) g(:,:,0:kp)
    if ( present(rho) ) read(u) rho(:,:,0:kp)
    close(u)
  end subroutine read_diis_atomopt_io


  subroutine write_diis_atomopt_io &
       ( loop, history, ip, x, g, rho )
    implicit none
    integer,intent(in) :: loop
    real(8),intent(in) :: history(:,0:)
    integer,intent(in) :: ip
    real(8),intent(in) :: x(:,:,0:), g(:,:,0:)
    real(8),optional,intent(in) :: rho(:,:,0:)
    open(u,file=file_name,form='unformatted')
    write(u) version
    write(u) loop
    write(u) history(:,:)
    write(u) ip
    write(u) x
    write(u) g
    if ( present(rho) ) write(u) rho
    close(u)
  end subroutine write_diis_atomopt_io


end module atomopt_io_module
