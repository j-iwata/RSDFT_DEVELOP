module atom_module

  use lattice_module, only: lattice, read_lattice, get_inverse_lattice &
                           ,get_aa_lattice
  use cif_format_module
  use vasp_format_module

  implicit none

  private
  public :: atom
  public :: construct_atom
  public :: write_info_atom
  public :: read_atom
  public :: convert_to_aa_coordinates_atom
  public :: convert_to_xyz_coordinates_atom
  public :: write_coordinates_atom
  public :: shift_aa_coordinates_atom
  public :: write_xyz_atom
  public :: set_constraint_atom

  integer,parameter :: DP=kind(0.0d0)

  integer,PUBLIC :: Natom
  integer,PUBLIC :: Nelement
  integer,allocatable,PUBLIC :: ki_atom(:)
  integer,allocatable,PUBLIC :: zn_atom(:)
  integer,allocatable,PUBLIC :: md_atom(:)
  real(DP),allocatable,PUBLIC :: aa_atom(:,:)

  type atom
    integer :: natom, nelement
    integer,allocatable :: k(:)
    integer,allocatable :: z(:)
    real(DP),allocatable :: aaa(:,:)
    real(DP),allocatable :: xyz(:,:)
    real(DP),allocatable :: force(:,:)
  end type atom

  integer :: iformat=0
  integer :: iformat_org=0
  real(DP),parameter :: bohr=0.529177d0

  logical,allocatable :: md_flags(:,:)

contains


  subroutine read_atom( rank, unit, aa_obj )

    implicit none
    integer,intent(in) :: rank, unit
    type(lattice),intent(out) :: aa_obj
    character(8) :: cbuf, ckey, cdummy(2)
    character(80) :: line
    integer :: idummy(2),ierr
    include 'mpif.h'

    call write_border( 0, " read_atom(start)" )

    call read_lattice( aa_obj, unit )

    if ( rank == 0 ) then

      iformat=0

      rewind unit
      do
        read(unit,*,END=10) cbuf
        call convert_capital(cbuf,ckey)
        if ( cbuf == "AA" ) then
          iformat=1
          write(*,*) "AA format",iformat
          exit
        else if ( cbuf == "XYZ" ) then
          iformat=2
          write(*,*) "rsdft-XYZ format",iformat
          exit
        end if
      end do
      10 continue

      if ( iformat == 0 ) then
        call check_cif_format( unit, ierr )
        if ( ierr == 0 ) iformat=4
      end if

      if ( iformat == 0 ) then
        call check_vasp_format( unit, ierr )
        if ( ierr == 0 ) iformat=5
      end if

      if ( iformat == 0 ) then
        rewind unit
        idummy(:)=0
        read(unit,'(a)') line            ! Couting the number of
        read(line,*,END=11) idummy(1:2)  ! integer data in the first line.
        11 continue
        backspace( unit )
        if ( all(idummy/=0) ) then
          iformat=1
          write(*,*) "AA format is assumed",iformat
        else if ( idummy(1) /= 0 ) then
          iformat=3
          write(*,*) "XYZ format is assumed",iformat
        else
          write(*,*) "unknown format"
          call stop_program( "stop@check_format_atom" )
        end if
      end if

    end if ! rank == 0

    call MPI_Bcast(iformat,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    select case( iformat )
    case( 1,2 )
      call read_atom_rsdft( rank, unit )
    case( 3 )
      call read_atom_xyz( rank, unit )
    case( 4 )
      call read_atom_cif( rank, unit, aa_obj, aa_atom, ki_atom, md_atom, zn_atom )
      iformat=1
      Natom=size(ki_atom)
      Nelement=maxval(ki_atom)
    case( 5 )
      call read_atom_vasp( rank,unit,aa_obj,aa_atom,ki_atom,md_atom,zn_atom,md_flags )
      iformat=1
      Natom=size(ki_atom)
      Nelement=maxval(ki_atom)
    end select

    iformat_org = iformat

    call write_border( 0, " read_atom(end)" )

  end subroutine read_atom


  SUBROUTINE read_atom_rsdft( rank, unit )
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: idummy(99),i
    include 'mpif.h'

    if ( rank == 0 ) then
       idummy(:)=0
       read(unit,*) Nelement, Natom, idummy(1:Nelement)
       allocate( zn_atom(Nelement) ) ; zn_atom=0
       zn_atom(1:Nelement) = idummy(1:Nelement)
       write(*,*) "Nelment,Natom=",Nelement,Natom
       write(*,*) "zn_atom=",zn_atom(:)
    end if

    call mpi_bcast(Natom,1,MPI_INTEGER,0,MPI_COMM_WORLD,i)
    call mpi_bcast(Nelement,1,MPI_INTEGER,0,MPI_COMM_WORLD,i)

    allocate( aa_atom(3,Natom) ) ; aa_atom=0.d0
    allocate( ki_atom(Natom)   ) ; ki_atom=0
    allocate( md_atom(Natom)   ) ; md_atom=0
    if ( .not.allocated(zn_atom) ) then
       allocate( zn_atom(Nelement) ) ; zn_atom=0
    end if

    if ( rank == 0 ) then
       do i=1,Natom
          read(unit,*) ki_atom(i),aa_atom(1:3,i),md_atom(i)
       end do
       call write_read_atom_data( ki_atom, aa_atom, md_atom )
    end if

    call send_atom_2

  END SUBROUTINE read_atom_rsdft


  SUBROUTINE send_atom_2
    implicit none
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(ki_atom,Natom,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(aa_atom,3*Natom,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(zn_atom,Nelement,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(md_atom,Natom,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_atom_2


  SUBROUTINE construct_atom( x )
    implicit none
    type(atom) :: x
    x%natom = Natom
    x%nelement = Nelement
    allocate( x%k(x%natom)       ) ; x%k(:) = ki_atom(:)
    allocate( x%z(x%nelement)    ) ; x%z(:) = zn_atom(:)
    allocate( x%aaa(3,x%natom)   ) ; x%aaa(:,:) = aa_atom(:,:)
    allocate( x%xyz(3,x%natom)   ) ; x%xyz(:,:) = 0.0d0
    allocate( x%force(3,x%natom) ) ; x%force(:,:) = 0.0d0
  END SUBROUTINE construct_atom


  SUBROUTINE read_atom_xyz( rank, unit )
    implicit none
    integer,intent(IN) :: rank, unit
    integer :: i, j, idummy(99)
    character(2) :: cbuf
    include 'mpif.h'

    if ( rank == 0 ) then

       read(unit,*) Natom
       read(unit,*)

       allocate( aa_atom(3,Natom) ) ; aa_atom=0.0d0
       allocate( ki_atom(Natom)   ) ; ki_atom=0
       allocate( md_atom(Natom)   ) ; md_atom=0

       idummy(:)=0
       Nelement=0
       do i=1,Natom
          read(unit,*) cbuf, aa_atom(1:3,i)
          aa_atom(1:3,i)=aa_atom(1:3,i)/bohr
          call get_atomic_number( cbuf, ki_atom(i) )
          if ( Nelement==0 .or. all(idummy(1:Nelement)/=ki_atom(i)) ) then
             Nelement=Nelement+1
             idummy(Nelement)=ki_atom(i)
          end if
       end do

       allocate( zn_atom(Nelement) ) ; zn_atom=0
       zn_atom(1:Nelement) = idummy(1:Nelement)

       do i=1,Natom
          do j=1,Nelement
             if ( ki_atom(i) == zn_atom(j) ) then
                ki_atom(i)=j
                exit
             end if
          end do
       end do

       write(*,*) "Nelment,Natom=",Nelement,Natom
       write(*,*) "zn_atom=",zn_atom(:)
       call write_read_atom_data( ki_atom, aa_atom, md_atom )

    end if

    call mpi_bcast(Natom,1,MPI_INTEGER,0,MPI_COMM_WORLD,i)
    call mpi_bcast(Nelement,1,MPI_INTEGER,0,MPI_COMM_WORLD,i)

    if ( rank /= 0 ) then
       allocate( aa_atom(3,Natom)  ) ; aa_atom=0.d0
       allocate( ki_atom(Natom)    ) ; ki_atom=0
       allocate( md_atom(Natom)    ) ; md_atom=0
       allocate( zn_atom(Nelement) ) ; zn_atom=0
    end if

    call send_atom_2

  END SUBROUTINE read_atom_xyz


  SUBROUTINE write_read_atom_data( ki_atom, aa_atom, md_atom )
    implicit none
    integer,intent(IN) :: ki_atom(:), md_atom(:)
    real(8),intent(IN) :: aa_atom(:,:)
    integer :: Natom,i
    Natom = size( ki_atom )
    write(*,'(8x,a7,3a18,2x,a7)') &
         "ki_atom","aa_atom1","aa_atom2","aa_atom3","md_atom"
    if ( Natom <= 11 ) then
       do i=1,Natom
          write(*,'(1x,i5,2x,i7,3f18.12,4x,i5)') &
               i,ki_atom(i),aa_atom(:,i),md_atom(i)
       end do
    else
       do i=1,min(5,Natom)
          write(*,'(1x,i5,2x,i7,3f18.12,4x,i5)') &
               i,ki_atom(i),aa_atom(:,i),md_atom(i)
       end do
       write(*,'(1x,10x,".")')
       write(*,'(1x,10x,".")')
       write(*,'(1x,10x,".")')
       do i=Natom-5,Natom
          write(*,'(1x,i5,2x,i7,3f18.12,4x,i5)') &
               i,ki_atom(i),aa_atom(:,i),md_atom(i)
       end do
    end if
  END SUBROUTINE write_read_atom_data


  SUBROUTINE write_info_atom( zps, FilePS )
    implicit none
    real(8),intent(IN) :: zps(:)
    character(*),intent(IN) :: FilePS(:)
    integer :: i
    integer,allocatable :: num(:)
    logical :: disp_sw
    call write_border( 80, " write_info_atom(start)" )
    call check_disp_switch( disp_sw, 0 )
    if ( disp_sw ) then
       allocate( num(Nelement) ) ; num=0
       write(*,'(8x,3a8,4x,a)') "Zatom", "Zion", "Num", "FilePS"
       do i=1,Nelement
          num(i) = count( ki_atom == i )
          write(*,'(2i8,f8.2,i8,4x,a)') i,zn_atom(i),zps(i),num(i),FilePS(i)
       end do
       write(*,'(3x,"total",i8,f8.2,i8)') sum(zn_atom*num),sum(zps*num),Natom
       deallocate( num )
    end if
    call write_border( 80, " write_info_atom(end)" )
  END SUBROUTINE write_info_atom


  SUBROUTINE convert_to_aa_coordinates_atom( aa, aa_atom )
    implicit none
    type(lattice),intent(IN) :: aa
    real(8),intent(INOUT) :: aa_atom(:,:)
    real(8) :: aa_inverse(3,3), xyz_atom(3)
    integer :: i
    if ( iformat == 1 ) return
    iformat = 1
    call get_inverse_lattice( aa%LatticeVector, aa_inverse )
    do i=1,size( aa_atom, 2 )
       xyz_atom(1:3) = aa_atom(1:3,i)
       aa_atom(1:3,i) = matmul( aa_inverse(1:3,1:3), xyz_atom(1:3) )
    end do
  END SUBROUTINE convert_to_aa_coordinates_atom


  SUBROUTINE convert_to_xyz_coordinates_atom( aa, aa_atom )
    implicit none
    type(lattice),intent(IN) :: aa
    real(8),intent(INOUT) :: aa_atom(:,:)
    real(8) :: aa_inverse(3,3), xyz_atom(3)
    integer :: i
    if ( iformat == 2 .or. iformat == 3 ) return
    iformat = 2
    do i=1,size( aa_atom, 2 )
       xyz_atom(1:3) = matmul( aa%LatticeVector(1:3,1:3), aa_atom(1:3,i) )
       aa_atom(1:3,i) = xyz_atom(1:3)
    end do
  END SUBROUTINE convert_to_xyz_coordinates_atom


  SUBROUTINE get_dist( i, j, aa, aa_atom, d )
    implicit none
    integer,intent(IN)  :: i,j
    real(8),intent(IN)  :: aa(3,3),aa_atom(:,:)
    real(8),intent(OUT) :: d
    real(8) :: ri(3),rj(3)
    ri=matmul(aa,aa_atom(:,i))
    rj=matmul(aa,aa_atom(:,j))
    d=sqrt( sum((ri-rj)**2) )
  END SUBROUTINE get_dist


  SUBROUTINE shift_aa_coordinates_atom( aa_atom )
    implicit none
    real(8),intent(INOUT) :: aa_atom(:,:)
    integer :: i,j
!$OMP parallel do private( i,j )
    do i=1,size(aa_atom,2)
       do j=1,3
10        if ( aa_atom(j,i) < 0.0d0 ) then
             aa_atom(j,i)=aa_atom(j,i)+1.0d0
             goto 10
          else if ( aa_atom(j,i) >= 1.0d0 ) then
             aa_atom(j,i)=aa_atom(j,i)-1.0d0
             goto 10
          else
             cycle
          end if
       end do
    end do
!$OMP end parallel do
  END SUBROUTINE shift_aa_coordinates_atom


  SUBROUTINE write_coordinates_atom( unit, write_ctrl, file_name )
    implicit none
    integer,intent(IN) :: unit, write_ctrl
    character(*),optional,intent(IN) :: file_name
    integer :: a
    type(lattice) :: aa

    if ( present(file_name) ) open(unit,file=file_name)

    if ( write_ctrl == 1 .or. write_ctrl == 3 ) then
       call get_aa_lattice( aa )
       rewind unit
       write(unit,'("AX",1f20.15)') aa%LatticeConstant
       write(unit,'("A1",3f20.15)') aa%LatticeVector(1:3,1)/aa%LatticeConstant
       write(unit,'("A2",3f20.15)') aa%LatticeVector(1:3,2)/aa%LatticeConstant
       write(unit,'("A3",3f20.15)') aa%LatticeVector(1:3,3)/aa%LatticeConstant
       if ( iformat == 1 ) then
          write(unit,'("AA")')
       else if ( iformat == 2 .or. iformat == 3 ) then
          write(unit,'("XYZ")')
       end if
       write(unit,*) Nelement,Natom,zn_atom(1:Nelement), " /"
       if ( write_ctrl == 1 ) return
    end if

    if ( write_ctrl == 2 .or. write_ctrl == 3 ) then
       do a=1,Natom
          write(unit,'(1x,i5,3f20.12,i4)') &
               ki_atom(a),aa_atom(:,a),md_atom(a)
       end do
    end if

    if ( present(file_name) ) close(unit)

  END SUBROUTINE write_coordinates_atom


  subroutine write_xyz_atom( unit, xyz, comment )
    implicit none
    integer,intent(in) :: unit
    real(8),intent(in) :: xyz(:,:)
    character(*),intent(in) :: comment
    integer :: i
    character(2),external :: get_element_name
    write(unit,*) size(xyz,2)
    write(unit,*) trim(comment)
    do i=1,size(xyz,2)
      write(unit,'(a2,2x,3g24.15)') get_element_name(zn_atom(ki_atom(i))), xyz(:,i)
    end do
  end subroutine write_xyz_atom


  subroutine set_constraint_atom( tim, ierr )
    implicit none
    real(8),intent(inout) :: tim(:,:,:)
    integer,intent(out) :: ierr
    integer :: i_md,i
    call write_border( 0, ' set_constraint_atom(start)' )
    ierr = 0
    if ( .not.allocated(md_flags) ) then
      ierr = -1
      call write_border( 0, ' set_constraint_atom(return)' )
      return
    end if
    if ( size(tim,3) /= size(md_flags,2) ) call stop_program('@set_constraint_atom')
    do i_md = 1, size(tim,3)
      tim(:,:,i_md)=0.0d0
      do i = 1, 3
        if ( md_flags(i,i_md) ) tim(i,i,i_md)=1.0d0
      end do
    end do
    call write_border( 0, ' set_constraint_atom(end)' )
  end subroutine set_constraint_atom

end module atom_module
