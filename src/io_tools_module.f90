module io_tools_module

  implicit none

  private
  public :: init_io_tools
  public :: IOTools_readStringKeyword
  public :: IOTools_readIntegerKeyword
  public :: IOTools_readReal8Keyword
  public :: IOTools_readIntegerString
  public :: IOTools_findKeyword
  public :: IOTools_parseInteger

  integer,parameter :: max_trial_read = 10000
  integer :: myrank = 0
  integer :: unit_default = 1
  integer :: unit_input_result = 110
  logical :: flag_init = .false.

  INTERFACE IOTools_readIntegerKeyword
     MODULE PROCEDURE IOTools_readIntegerKeyword_sca &
                     ,IOTools_readIntegerKeyword_vec
  END INTERFACE

  INTERFACE IOTools_readStringKeyword
     MODULE PROCEDURE IOTools_readStringKeyword_sca &
                     ,IOTools_readStringKeyword_vec
  END INTERFACE

  INTERFACE IOTools_readReal8Keyword
     MODULE PROCEDURE IOTools_readReal8Keyword_sca &
                     ,IOTools_readReal8Keyword_vec
  END INTERFACE

contains


  SUBROUTINE init_io_tools( myrank_in, unit_in )
    implicit none
    integer,intent(IN) :: myrank_in, unit_in
    myrank = myrank_in
    unit_default = unit_in
!    if ( myrank == 0 ) open(unit_input_result,file="input_result")
    flag_init = .true.
  END SUBROUTINE init_io_tools

  SUBROUTINE check_init
    implicit none
    integer :: ierr
#ifdef _NOMPI_
    flag_init = .true.
    return
#else
    include 'mpif.h'
    if ( .not.flag_init ) then
       call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )
    end if
    flag_init = .true.
#endif
  END SUBROUTINE check_init


  SUBROUTINE IOTools_readStringKeyword_sca( keyword, variable, unit_in )
    implicit none
    character(*),intent(IN) :: keyword
    character(*),intent(INOUT) :: variable
    integer,optional,intent(IN) :: unit_in
    character(20) :: cbuf,ckey
    integer :: i,unit
#ifndef _NOMPI_
    include 'mpif.h'
#endif
    call check_init
    unit=unit_default ; if ( present(unit_in) ) unit=unit_in
    if ( myrank == 0 ) then
       rewind unit
       do i=1,max_trial_read
          read(unit,*,END=999) cbuf
          call convertToCapital(cbuf,ckey)
          if ( ckey == keyword ) then
             backspace(unit)
             read(unit,*) cbuf, variable
             write(*,'(1x,A10," : ",A10)') keyword, variable
             exit
          end if
       end do ! i
999    continue
    end if
#ifndef _NOMPI_
    call MPI_BCAST( variable,len(variable),MPI_CHARACTER,0,MPI_COMM_WORLD,i )
#endif
  END SUBROUTINE IOTools_readStringKeyword_sca


  SUBROUTINE IOTools_readStringKeyword_vec( keyword, variables, unit_in )
    implicit none
    character(*),intent(IN) :: keyword
    character(*),intent(INOUT) :: variables(:)
    integer,optional,intent(IN) :: unit_in
    character(20) :: cbuf,ckey
    integer :: i,j,unit
#ifndef _NOMPI_
    include 'mpif.h'
#endif
    call check_init
    unit=unit_default ; if ( present(unit_in) ) unit=unit_in
    if ( myrank == 0 ) then
       rewind unit
       do i=1,max_trial_read
          read(unit,*,END=999) cbuf
          call convertToCapital(cbuf,ckey)
          if ( ckey == keyword ) then
             backspace(unit)
             read(unit,*) cbuf, variables(:)
             write(*,'(1x,A10," : ",A10)') keyword, variables(1)
             do j=2,size(variables)
                write(*,'(1x,10x,"   ",A10)') variables(j)
             end do
             exit
          end if
       end do ! i
999    continue
    end if
#ifndef _NOMPI_
    call MPI_BCAST( variables, len(variables)*size(variables), &
                    MPI_CHARACTER, 0, MPI_COMM_WORLD, i )
#endif
  END SUBROUTINE IOTools_readStringKeyword_vec


  SUBROUTINE IOTools_readIntegerKeyword_sca( keyword, variable, unit_in )
    implicit none
    character(*),intent(IN) :: keyword
    integer,intent(INOUT) :: variable
    integer,optional,intent(IN) :: unit_in
    character(20) :: cbuf,ckey
    integer :: i,unit
#ifndef _NOMPI_
    include 'mpif.h'
#endif
    call check_init
    unit=unit_default ; if ( present(unit_in) ) unit=unit_in
    if ( myrank == 0 ) then
       rewind unit
       do i=1,max_trial_read
          read(unit,*,END=999) cbuf
          call convertToCapital(cbuf,ckey)
          if ( ckey == keyword ) then
             backspace(unit)
             read(unit,*) cbuf,variable
             write(*,'(1x,A10," : ",3I10)') keyword,variable
             exit
          end if
       end do ! i
999    continue
       !write(unit_input_result,'(a10,3x,i10)') keyword,variable
    end if
#ifndef _NOMPI_
    call MPI_BCAST(variable,1,MPI_INTEGER,0,MPI_COMM_WORLD,i)
#endif
  END SUBROUTINE IOTools_readIntegerKeyword_sca


  subroutine IOTools_readIntegerKeyword_vec( keyword, variables, unit_in )
    implicit none
    character(*),intent(in) :: keyword
    integer,intent(inout) :: variables(:)
    integer,optional,intent(in) :: unit_in
    character(20) :: cbuf
    integer :: unit,i
#ifndef _NOMPI_
    include 'mpif.h'
#endif
    call write_border_f( 0, ' IOTools_readIntegerKeyword_vec(start)' )
    call check_init
    unit=unit_default; if ( present(unit_in) ) unit=unit_in
    if ( myrank == 0 ) then
      rewind unit
      if ( findKey(keyword,unit) ) then
        backspace(unit)
        read(unit,*) cbuf, variables(:)
        write(*,'(1x,A10," : ",7I10)') keyword,variables(:)
      else
        write(*,'(1x,A10," : ","not found")') keyword
      end if
    end if
#ifndef _NOMPI_
    call MPI_Bcast(variables,size(variables),MPI_INTEGER,0,MPI_COMM_WORLD,i)
#endif
    call write_border_f( 0, ' IOTools_readIntegerKeyword_vec(end)' )
  end subroutine IOTools_readIntegerKeyword_vec


  SUBROUTINE IOTools_readReal8Keyword_sca( keyword, variable, unit_in )
    implicit none
    character(*),intent(IN) :: keyword
    real(8),intent(INOUT) :: variable
    integer,optional,intent(IN) :: unit_in
    character(20) :: cbuf,ckey
    integer :: i,unit
#ifndef _NOMPI_
    include 'mpif.h'
#endif
    call check_init
    unit=unit_default ; if ( present(unit_in) ) unit=unit_in
    if ( myrank == 0 ) then
       rewind unit
       do i=1,max_trial_read
          read(unit,*,END=999) cbuf
          call convertToCapital(cbuf,ckey)
          if ( ckey == keyword ) then
             backspace(unit)
             read(unit,*) cbuf,variable
             if ( abs(variable) < 1.d-2 .or. abs(variable) > 1.d3 ) then
                write(*,'(1x,A10," : ",ES15.7)') keyword,variable
             else
                write(*,'(1x,A10," : ",F15.10)') keyword,variable
             end if
             exit
          end if
       end do ! i
999    continue
    end if
#ifndef _NOMPI_
    call MPI_BCAST(variable,1,MPI_REAL8,0,MPI_COMM_WORLD,i)
#endif
  END SUBROUTINE IOTools_readReal8Keyword_sca


  SUBROUTINE IOTools_readReal8Keyword_vec( keyword, variables, unit_in )
    implicit none
    character(*),intent(IN) :: keyword
    real(8),intent(INOUT) :: variables(:)
    integer,optional,intent(IN) :: unit_in
    character(20) :: cbuf,ckey
    integer :: i,unit
#ifndef _NOMPI_
    include 'mpif.h'
#endif
    call check_init
    unit=unit_default ; if ( present(unit_in) ) unit=unit_in
    if ( myrank == 0 ) then
       rewind unit
       do i=1,max_trial_read
          read(unit,*,END=999) cbuf
          call convertToCapital(cbuf,ckey)
          if ( ckey == keyword ) then
             backspace(unit)
             read(unit,*) cbuf,variables(:)
             if ( any(abs(variables) < 1.d-2 .or. abs(variables) > 1.d3) ) then
                write(*,'(1x,A10," : ",10(3ES15.7,/))') keyword,variables(:)
             else
                write(*,'(1x,A10," : ",10(3F15.10,/))') keyword,variables(:)
             end if
             exit
          end if
       end do ! i
999    continue
    end if
#ifndef _NOMPI_
    call MPI_BCAST(variables,size(variables),MPI_REAL8,0,MPI_COMM_WORLD,i)
#endif
  END SUBROUTINE IOTools_readReal8Keyword_vec


  SUBROUTINE IOTools_readIntegerString( keyword, variable1, variable2, u, norewind )
    implicit none
    character(*),intent(IN) :: keyword
    integer,intent(INOUT) :: variable1
    character(*),intent(INOUT) :: variable2
    integer,optional,intent(IN) :: u
    logical,optional,intent(IN) :: norewind
    character(20) :: cbuf,ckey
    integer :: i,unit
#ifndef _NOMPI_
    include 'mpif.h'
#endif
    call check_init
    unit=unit_default ; if ( present(u) ) unit=u
    if ( myrank == 0 ) then
       if ( .not.present(norewind) ) rewind unit
       do i=1,max_trial_read
          read(unit,*,END=999) cbuf
          call convertToCapital(cbuf,ckey)
          if ( ckey == keyword ) then
             backspace(unit)
             read(unit,*) cbuf,variable1, variable2
             write(*,'(1x,A10," : ",i4,2x,a20)') keyword,variable1,variable2
             exit
          end if
       end do ! i
999    continue
    end if
#ifndef _NOMPI_
    call MPI_BCAST(variable1,1,MPI_INTEGER,0,MPI_COMM_WORLD,i)
    call MPI_BCAST(variable2,len(variable2),MPI_CHARACTER,0,MPI_COMM_WORLD,i)
#endif
  END SUBROUTINE IOTools_readIntegerString


  SUBROUTINE IOTools_findKeyword( keyword, hasKeyword, unit_out, flag_bcast )
    implicit none
    character(*),intent(IN) :: keyword
    logical,intent(OUT) :: hasKeyword
    integer,optional,intent(OUT) :: unit_out
    logical,optional,intent(IN) :: flag_bcast
    integer :: i,unit
    character(20) :: cbuf,ckey
#ifndef _NOMPI_
    include 'mpif.h'
#endif
    call check_init
    unit=unit_default
    if ( present(unit_out) ) unit_out=unit
    hasKeyword = .false.
    if ( myrank == 0 ) then
       rewind unit
       do i=1,max_trial_read
          read(unit,*,END=999) cbuf
          call convertToCapital(cbuf,ckey)
          if ( ckey == keyword ) then
             hasKeyword = .true.
             exit
          end if
       end do ! i
       if ( hasKeyword ) write(*,'(1x,A10)') keyword
999    continue
    end if
    if ( present(flag_bcast) ) then
#ifndef _NOMPI_
       call MPI_BCAST( hasKeyword, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, i )
#endif
    end if
  END SUBROUTINE IOTools_findKeyword


  subroutine IOTools_parseInteger( keyword, values, unit_in )
    use rsdft_bcast_module, only: l_rsdft_bcast, i_rsdft_bcast
    implicit none
    character(*),intent(in) :: keyword
    integer,allocatable,intent(inout) :: values(:)
    integer,optional,intent(in) :: unit_in
    integer :: unit,m1,m2,m3,m4,n1,n2,n,i,j,itmp(10)
    character(127) :: str,ctmp
    character(7) :: cint
    logical :: flag

    call write_border_f( 0, ' IOTools_parseInteger(start)' )

    call check_init
    unit=unit_default; if ( present(unit_in) ) unit=unit_in

    if ( allocated(values) ) deallocate(values)

    flag = .false.
    if ( myrank == 0 ) then

      rewind unit
      ! do
      !   read(unit,'(a)',END=999) str
      !   call convertToCapital( str )
      !   if ( index(str,keyword) > 0 ) exit
      ! end do

      flag = findKey( keyword, unit )
      if ( .not.flag ) then
        write(*,'(1x,a)') keyword//" : not found"
        goto 999
      end if

      backspace(unit)
      read(unit,'(a)') str

      m1 = len_trim( keyword )
      m2 = index( str, keyword(1:m1) )
      n1 = m1 + m2
      n2 = len_trim(str)
      ! write(*,*) str(1:n2)

      itmp=0
      n=0
      do
        ctmp = str(n1:n2)
        cint = ''
        read(ctmp,*) cint
        if ( not_integer(cint) ) exit
        n = n + 1
        read(ctmp,*) itmp(n)
        m1 = len_trim( cint )
        m2 = index( ctmp, cint(1:m1) )
        n1 = n1 + m1 + m2 - 1
        if ( n1 > n2 ) exit
      end do

      write(*,'(1x,a,10i4)') keyword//" : ", itmp(1:n)
      write(*,*) "(# of Parameters = ",n,")"
      allocate( values(n) ); values=0
      values = itmp(1:n)

      999 continue
    end if

    call l_rsdft_bcast( flag, 1, 0 )
    if ( flag ) then
      call i_rsdft_bcast( n, 1, 0 )
      if ( .not.allocated(values) ) then
        allocate( values(n) ); values=0
      end if
      call i_rsdft_bcast( values, n, 0 )
    end if

    call write_border_f( 0, ' IOTools_parseInteger(end)' )

  end subroutine IOTools_parseInteger


  logical function findKey( keyword, unit )
    implicit none
    character(*),intent(in) :: keyword
    integer,intent(in) :: unit
    character(20) :: ckey
    findKey = .false.
    if ( myrank == 0 ) then
      rewind unit
      do
        read(unit,*,END=999) ckey
        call convertToCapital( ckey )
        if ( ckey == keyword ) then
          findKey = .true.
          return
        end if
      end do
      999 continue
    end if
  end function findKey


  subroutine convertToCapital( cbuf, CKEY )
    implicit none
    character(*),intent(inout) :: cbuf
    character(*),optional,intent(out) :: CKEY
    integer :: j,k,n
    n=len_trim(cbuf)
    if ( present(CKEY) ) CKEY=cbuf(1:n)
    do j=1,n
      k=iachar( cbuf(j:j) )
      if ( 97 <= k .and. k <= 122 ) k=k-32
      if ( present(CKEY) ) then
        CKEY(j:j) = achar(k)
      else
        cbuf(j:j) = achar(k)
      end if
    end do
  end subroutine convertToCapital


  logical function not_integer( c )
    implicit none
    character(*),intent(in) :: c
    integer :: n,i,j
    n = len_trim(c)
    if ( n == 0 ) then
      not_integer = .true.
      return
    end if
    do i = 1, n
      j = iachar( c(i:i) )
      if ( .not.(48 <= j .and. j <= 57) ) then
        not_integer = .true.
        return
      end if
    end do
    not_integer = .false.
  end function not_integer

end module io_tools_module
