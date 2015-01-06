MODULE aa_module

  implicit none

  PRIVATE
  PUBLIC :: ax,aa,Va,read_aa,construct_aa,read_oldformat_aa &
           ,get_org_aa,set_org_aa

  real(8) :: ax,Va
  real(8) :: aa(3,3)

  real(8) :: ax_org
  real(8) :: aa_org(3,3)

CONTAINS

  SUBROUTINE read_aa(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i
    character(2) :: cbuf,ckey
    ax=0.d0
    aa(:,:)=0.d0
    aa(1,1)=1.d0
    aa(2,2)=1.d0
    aa(3,3)=1.d0
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey == "AX" ) then
             backspace(unit)
             read(unit,*) cbuf,ax
          else if ( ckey == "A1" ) then
             backspace(unit)
             read(unit,*) cbuf,aa(1:3,1)
          else if ( ckey == "A2" ) then
             backspace(unit)
             read(unit,*) cbuf,aa(1:3,2)
          else if ( ckey == "A3" ) then
             backspace(unit)
             read(unit,*) cbuf,aa(1:3,3)
          end if
       end do
999    continue
       write(*,*) "ax=",ax
       write(*,'(1x,"aa(1:3,1)=",3F20.15)') aa(:,1)
       write(*,'(1x,"aa(1:3,2)=",3F20.15)') aa(:,2)
       write(*,'(1x,"aa(1:3,3)=",3F20.15)') aa(:,3)
    end if
    call send_aa(0)
    call set_org_aa( ax, aa )
  END SUBROUTINE read_aa


  SUBROUTINE read_oldformat_aa(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    if ( rank == 0 ) then
       read(unit,*) ax
       read(unit,*) aa(1:3,1)
       read(unit,*) aa(1:3,2)
       read(unit,*) aa(1:3,3)
       write(*,*) "ax=",ax
       write(*,'(1x,"aa(1:3,1)=",3F20.15)') aa(:,1)
       write(*,'(1x,"aa(1:3,2)=",3F20.15)') aa(:,2)
       write(*,'(1x,"aa(1:3,3)=",3F20.15)') aa(:,3)
    end if
    call send_aa(0)
    ax_org=ax
    aa_org=aa(:,:)
  END SUBROUTINE read_oldformat_aa


  SUBROUTINE send_aa(rank)
    implicit none
    integer,intent(IN) :: rank
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(ax,1,MPI_REAL8,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(aa,9,MPI_REAL8,rank,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_aa


  SUBROUTINE construct_aa
    implicit none
    aa(:,:)=ax*aa(:,:)
    Va = aa(1,1)*aa(2,2)*aa(3,3)+aa(1,2)*aa(2,3)*aa(3,1) &
        +aa(1,3)*aa(2,1)*aa(3,2)-aa(1,3)*aa(2,2)*aa(3,1) &
        -aa(1,2)*aa(2,1)*aa(3,3)-aa(1,1)*aa(2,3)*aa(3,2)
  END SUBROUTINE construct_aa


  SUBROUTINE set_org_aa(ax_in,aa_in)
    implicit none
    real(8),intent(IN) :: ax_in,aa_in(3,3)
    ax_org=ax_in
    aa_org(:,:)=aa_in(:,:)
  END SUBROUTINE set_org_aa


  SUBROUTINE get_org_aa(ax_out,aa_out)
    implicit none
    real(8),intent(OUT) :: ax_out,aa_out(3,3)
    ax_out=ax_org
    aa_out(:,:)=aa_org(:,:)
  END SUBROUTINE get_org_aa


END MODULE aa_module
