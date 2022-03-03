module vasp_format_module

  implicit none

  private
  public :: check_vasp_format
  public :: read_atom_vasp

contains

  subroutine check_vasp_format( unit, ierr )
    implicit none
    integer,intent(in)  :: unit
    integer,intent(out) :: ierr
    character(100) :: cbuf
    character(1) :: c1
    integer :: i,j
    logical :: flag1,flag2(3),flag3,flag4,flag5,flag
    
    call write_border( 0, " check_vasp_format(start)" )

    ierr = -1

    rewind unit
    read(unit,*,END=99)
    read(unit,'(a)',END=99) cbuf
    !
    i = data_count( cbuf )
    if ( i == 1 ) then
      flag1 = check_numbers( cbuf, i )
    else
      flag1 = .false.
    end if
    !
    do j = 1, 3
      read(unit,'(a)',END=99) cbuf
      i = data_count( cbuf )
      if ( i == 3 ) then
        flag2(j) = check_numbers( cbuf, i )
      else
        flag2(j) = .false.
      end if
    end do
    !
    read(unit,'(a)',END=99) cbuf
    flag3 = check_capital( cbuf )
    !
    read(unit,'(a)',END=99) cbuf
    flag4 = check_numbers( cbuf, 1 )
    !
    read(unit,'(a)',END=99) cbuf
    c1 = adjustl(cbuf)
    flag5 = any( c1 == (/'s','S','d','D','c','C'/) )
    if ( .not.flag5 ) then
      read(unit,'(a)',END=99) cbuf
      c1 = adjustl(cbuf)
      flag5 = any( c1 == (/'d','D','c','C'/) )
    end if

    flag = flag1*all(flag2)*flag3*flag4*flag5
    if ( flag ) ierr=0 !means this is VASP format
write(*,*) ierr,flag1,flag2,flag3,flag4,flag5
    99 continue

    call write_border( 0, " check_vasp_format(end)" )

  end subroutine check_vasp_format

  integer function data_count( c )
    implicit none
    character(*),intent(in) :: c
    character(10) :: ctmp(10)
    ctmp=''
    read(c,*,END=9) ctmp
    9 continue
    data_count = count( ctmp /= '' )
  end function data_count

  logical function check_numbers( c, n )
    implicit none
    character(*),intent(in) :: c
    integer,intent(in) :: n
    character(10),allocatable :: ctmp(:)
    character(1) :: c1
    character(2) :: c2
    integer :: i,i1,i2
    logical :: flag
    allocate( ctmp(n) ); ctmp=''
    read(c,*) ctmp(1:n)
    flag = .true.
    do i = 1, n
      c1 = adjustl( ctmp(i) )
      i1 = iachar( c1 )
      if ( 48 <= i1 .and. i1 <= 57 ) then
        cycle
      else if ( i1 == 45 ) then
        if ( len_trim(adjustl(ctmp(i))) >= 2 ) then
          c2 = adjustl( ctmp(i) )
          i2 = iachar( c2(2:2) )
          if ( 48 <= i2 .and. i2 <= 57 ) cycle
        end if
      end if
      flag = .false.
    end do
    check_numbers = flag
  end function check_numbers

  logical function check_capital( c )
    implicit none
    character(*),intent(in) :: c
    character(1) :: c1
    integer :: i1
    logical :: flag
    c1 = adjustl( c )
    i1 = iachar( c1 )
    flag=.false.
    if ( 65 <= i1 .and. i1 <= 90 ) flag=.true.
    check_capital = flag
  end function check_capital

  subroutine read_atom_vasp &
       ( rank, unit, aa_obj, aa_atom, ki_atom, md_atom, zn_atom, md_flags )
    use rsdft_bcast_module, only: i_rsdft_bcast, d_rsdft_bcast, d_rsdft_bcast_tmp, &
                                  l_rsdft_bcast_tmp
    use lattice_module, only: lattice, get_inverse_lattice
    implicit none
    integer,intent(in) :: rank, unit
    type(lattice),intent(inout) :: aa_obj
    real(8),intent(inout),allocatable :: aa_atom(:,:)
    integer,intent(inout),allocatable :: ki_atom(:),md_atom(:),zn_atom(:)
    logical,allocatable,intent(inout) :: md_flags(:,:)
    logical :: flag_selective_dynamics
    integer :: i,j,n,natm,nelm,n_md,i_md
    character(40) :: cbuf
    character(1) :: flags(3), md_cflags(3,8)
    character(2) :: element_name(10)
    integer :: num_element(10)
    real(8),parameter :: bohr=0.529177d0
    real(8) :: ax,aa(3,3),aa_inv(3,3)
    real(8),allocatable :: rr_atom(:,:)

    call write_border( 0, " read_atom_vasp(start)" )

    if ( rank == 0 ) then

      rewind unit
      read(unit,'(a)') cbuf
      write(*,*) "-------------- Comment line"
      write(*,*) cbuf
      write(*,*) "---------- End Comment line"

      read(unit,*) ax
      aa_obj%LatticeConstant = ax / bohr

      write(*,'(1x,"LatticeConstant: ",f15.8,"(angs)",2x,f15.8,"(bohr)")') ax,ax/bohr

      read(unit,*) aa_obj%LatticeVector(:,1)
      read(unit,*) aa_obj%LatticeVector(:,2)
      read(unit,*) aa_obj%LatticeVector(:,3)

      where( abs(aa_obj%LatticeVector) < 1.d-5 )
        aa_obj%LatticeVector=0.0d0
      end where

      write(*,*) "LatticeVectors"
      write(*,'(1x,3f20.15)') aa_obj%LatticeVector(:,1)
      write(*,'(1x,3f20.15)') aa_obj%LatticeVector(:,2)
      write(*,'(1x,3f20.15)') aa_obj%LatticeVector(:,3)

      read(unit,'(a)') cbuf
      nelm = data_count( cbuf )

      backspace(unit)

      element_name = ''
      read(unit,*) (element_name(i),i=1,nelm)

      num_element = 0
      read(unit,*) (num_element(i),i=1,nelm)

      allocate( zn_atom(nelm) ); zn_atom=0
      do i = 1, nelm
        call get_atomic_number( element_name(i), zn_atom(i) )
      end do

      write(*,'(1x,"Element",2x,"AtomicNumber",2x,"Quantity")')
      do i = 1, nelm
        write(*,'(1x,"(",i1,")",a4,4x,i6,6x,i6)') i,element_name(i), zn_atom(i), num_element(i)
      end do

      natm = sum( num_element )
      write(*,*) "# of atoms =",natm

      read(unit,'(a)') cbuf

      flag_selective_dynamics = .false.
      if ( cbuf(1:1) == 's' .or. cbuf(1:1) == 'S' ) then
        flag_selective_dynamics = .true.
        write(*,*) trim(cbuf)
      else
        backspace(unit)
      end if

      allocate( aa_atom(3,natm) ); aa_atom=0.0d0
      allocate( ki_atom(natm)   ); ki_atom=0
      allocate( md_atom(natm)   ); md_atom=1

      read(unit,*) cbuf
      if ( .not.any( cbuf(1:1) == (/'c','C','d','D'/) ) ) read(unit,*) cbuf
      write(*,*) "Coordinates format: ",trim(cbuf)

      n_md=0
      if ( flag_selective_dynamics ) then
        md_cflags=''
        n = 0
        do j = 1, nelm
          do i = 1, num_element(j)
            n = n + 1
            ki_atom(n) = j
            read(unit,*) aa_atom(:,n), flags(1:3)
            call convert_to_capital_array( 3, flags )
            if ( n_md == 0 ) then
              n_md = n_md + 1
              md_cflags(:,n_md) = flags
              md_atom(n) = n_md
            else
              do i_md = 1, n_md
                if ( all(flags(:)==md_cflags(:,i_md)) ) then
                  md_atom(n) = i_md
                  exit
                else if ( i_md == n_md ) then
                  n_md = n_md + 1
                  md_cflags(:,n_md) = flags
                  md_atom(n) = n_md
                end if
              end do !i_md
            end if
          end do !i
        end do !j
        if ( n_md > 0 ) then
          allocate( md_flags(3,n_md) ); md_flags=.true.
          do i_md = 1, n_md
            do i = 1, 3
              md_flags(i,i_md) = ( md_cflags(i,i_md) == 'T' )
            end do
          end do
        end if
      else
        n = 0
        do j = 1, nelm
          do i = 1, num_element(j)
            n = n + 1
            read(unit,*) aa_atom(:,n)
            ki_atom(n) = j
          end do
        end do
      end if

      if ( cbuf(1:1) == 'd' .or. cbuf(1:1) == 'D' ) then
        continue
      else if ( cbuf(1:1) == 'c' .or. cbuf(1:1) == 'C' ) then
        allocate( rr_atom(3,natm) ); rr_atom=0.0d0
        rr_atom = aa_obj%LatticeConstant * aa_atom
        aa = aa_obj%LatticeConstant * aa_obj%LatticeVector
        call get_inverse_lattice( aa, aa_inv )
        aa_atom = matmul( aa_inv, rr_atom )
        deallocate( rr_atom )
      else
        call stop_program('Invalid keyword ? (stop@read_vasp)')
      end if

    end if ! rank == 0

    call i_rsdft_bcast( natm )
    call i_rsdft_bcast( nelm )
    if ( .not.allocated(aa_atom) ) then
      allocate( aa_atom(3,natm) ); aa_atom=0.0d0
      allocate( ki_atom(natm)   ); ki_atom=0
      allocate( md_atom(natm)   ); md_atom=0
      allocate( zn_atom(nelm)   ); zn_atom=0
    end if
    call d_rsdft_bcast_tmp( aa_atom, 3*natm )
    call i_rsdft_bcast( ki_atom )
    call i_rsdft_bcast( md_atom )
    call i_rsdft_bcast( zn_atom )

    call i_rsdft_bcast( n_md )
    if ( n_md > 0 ) then
      if ( .not.allocated(md_flags) ) then
        allocate( md_flags(3,n_md) ); md_flags=.true.
      end if
      call l_rsdft_bcast_tmp( md_flags, 3*n_md )
    end if

    call write_border( 0, " read_atom_vasp(end)" )

  end subroutine read_atom_vasp

end module vasp_format_module
