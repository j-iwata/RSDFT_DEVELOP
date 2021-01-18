module memory_module

  use rsdft_mpi_module, only: rsdft_allreduce

  implicit none

  private
  public :: check_memory

  interface check_memory
    module procedure check_memory_a, check_memory_b
  end interface check_memory

  real(8),parameter :: B2MB = 1.0d0/1024.0d0**2
  character(2),parameter :: indx_unit='MB'

contains

  subroutine check_memory_a( indx, n1, n2, n3, n4 ) 
    implicit none
    character(*),intent(in) :: indx
    integer,intent(in) :: n1
    integer,optional,intent(in) :: n2
    integer,optional,intent(in) :: n3
    integer,optional,intent(in) :: n4
    real(8) :: array_size, mem_unit_byte, mem_unit, mem
    logical :: disp_sw
    integer :: m1,m2,m3,m4
    call write_border(0,' check_memory(start)')
    array_size = n1
    m1=n1; call rsdft_allreduce( m1, op_in='max' )
    if ( present(n2) ) then
      array_size = array_size*n2
      m2=n2; call rsdft_allreduce( m2, op_in='max' )
    end if
    if ( present(n3) ) then
      array_size = array_size*n3
      m3=n3; call rsdft_allreduce( m3, op_in='max' )
    end if
    if ( present(n4) ) then
      array_size = array_size*n4
      m4=n4; call rsdft_allreduce( m4, op_in='max' )
    end if
    if ( indx == 'wf' .or. indx == 'dz' ) then
#ifdef _DRSDFT_
      mem_unit_byte =  8.0d0
#else
      mem_unit_byte = 16.0d0
#endif
    end if
    call check_disp_switch( disp_sw, 0 )
    if ( disp_sw ) then
      write(*,*) "indx, mem_unit_byte: ",indx, mem_unit_byte
      write(*,*) "array_size1(max)=",m1
      if ( present(n2) ) write(*,*) "array_size2(max)=",m2
      if ( present(n3) ) write(*,*) "array_size3(max)=",m3
      if ( present(n4) ) write(*,*) "array_size4(max)=",m4
    end if
    mem_unit = mem_unit_byte * B2MB
    mem = mem_unit * array_size
    call show_mem( mem )
    call write_border(0,' check_memory(end)')
  end subroutine check_memory_a


  subroutine check_memory_b( mem_unit_byte, n1, n2, n3, n4 ) 
    implicit none
    real(8),intent(in) :: mem_unit_byte
    integer,intent(in) :: n1
    integer,optional,intent(in) :: n2
    integer,optional,intent(in) :: n3
    integer,optional,intent(in) :: n4
    real(8) :: array_size, mem_unit, mem
    logical :: disp_sw
    integer :: m1,m2,m3,m4
    call write_border(0,' check_memory(start)')
    array_size = n1
    m1=n1; call rsdft_allreduce( m1, op_in='max' )
    if ( present(n2) ) then
      array_size = array_size*n2
      m2=n2; call rsdft_allreduce( m2, op_in='max' )
    end if
    if ( present(n3) ) then
      array_size = array_size*n3
      m3=n3; call rsdft_allreduce( m3, op_in='max' )
    end if
    if ( present(n4) ) then
      array_size = array_size*n4
      m4=n4; call rsdft_allreduce( m4, op_in='max' )
    end if
    call check_disp_switch( disp_sw, 0 )
    if ( disp_sw ) then
      write(*,*) "indx, mem_unit_byte: ","other", mem_unit_byte
      write(*,*) "array_size1(max)=",m1
      if ( present(n2) ) write(*,*) "array_size2(max)=",m2
      if ( present(n3) ) write(*,*) "array_size3(max)=",m3
      if ( present(n4) ) write(*,*) "array_size4(max)=",m4
    end if
    mem_unit = mem_unit_byte * B2MB
    mem = mem_unit * array_size
    call show_mem( mem )
    call write_border(0,' check_memory(end)')
  end subroutine check_memory_b

  subroutine show_mem( mem )
    implicit none
    real(8),intent(in) :: mem
    real(8) :: mem_tot, mem_min, mem_max
    logical :: disp_sw
    mem_tot=mem; call rsdft_allreduce( mem_tot )
    mem_max=mem; call rsdft_allreduce( mem_max, op_in='max' )
    mem_min=mem; call rsdft_allreduce( mem_min, op_in='min' )
    call check_disp_switch( disp_sw, 0 )
    if ( disp_sw ) then
      write(*,'( "total:", f12.3, 1x,"(",a2,")")') mem_tot,indx_unit
      write(*,'( "max  :", f12.3, 1x,"(",a2,")")') mem_max,indx_unit
      write(*,'( "min  :", f12.3, 1x,"(",a2,")")') mem_min,indx_unit
    end if 
  end subroutine show_mem

end module memory_module
