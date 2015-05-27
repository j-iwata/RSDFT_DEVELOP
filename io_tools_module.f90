MODULE IOTools

  implicit none

  PRIVATE
!  PUBLIC :: IOTools_findKeyword
  PUBLIC :: IOTools_readStringKeyword
!  PUBLIC :: IOTools_readRealKeyword
!  PUBLIC :: IOTools_readRealVectorKeyword
!  PUBLIC :: IOTools_readIntegerKeyword
!  PUBLIC :: IOTools_readLogicalKeyword
!  PUBLIC :: IOTools_readIntegerVectorKeyword

CONTAINS
#ifdef TEST
  SUBROUTINE IOTools_findKeyword( keyword, unit_number, hasKeyword )
    implicit none
    character(*),intent(IN) :: keyword
    integer,intent(IN) :: unit_number
    logical,intent(OUT) :: hasKeyword
    integer :: i
    character(10) :: cbuf,ckey
    integer :: keyword_length
    keyword_length=len_trim(keyword)
    hasKeyword=.false.
    rewind unit_number
    do i=1,10000
       read(unit_number,*,END=999) cbuf
       call convertToCapital(cbuf,ckey)
       if ( ckey==keyword ) then
          hasKeyword=.true.
          exit
       endif
    enddo
999 continue
    if ( hasKeyword ) write(*,'(1x,A10)') keyword
  END SUBROUTINE IOTools_findKeyword
#endif
  SUBROUTINE IOTools_readStringKeyword( keyword, unit, variable )
    implicit none
    character(*),intent(IN) :: keyword
    integer,intent(IN) :: unit
    character(*),intent(OUT) :: variable
    logical :: hasKeyword=.false.
    integer :: i
    character(10) :: cbuf,ckey
    integer :: keyword_length
    keyword_length=len_trim(keyword)
    rewind unit
    do i=1,10000
      read(unit,*,END=999) cbuf
      call convertToCapital(cbuf,ckey)
      if ( ckey == keyword ) then
        backspace(unit)
        read(unit,*) cbuf, variable
        hasKeyword=.true.
        exit
      end if
    end do
999 continue
    if ( hasKeyWord ) write(*,'(1x,A10," : ",A10)') keyword,variable
  END SUBROUTINE IOTools_readStringKeyword
#ifdef TEST
  SUBROUTINE IOTools_readRealKeyword(keyword,unit_number,keyword_variable)
    implicit none
    character(*),intent(IN) :: keyword
    integer,intent(IN) :: unit_number
    real(8),intent(OUT) :: keyword_variable
    logical :: hasKeyword=.false.
    integer :: i
    character(10) :: cbuf,ckey
    integer :: keyword_length
    keyword_length=len_trim(keyword)
    rewind unit_number
    do i=1,10000
      read(unit_number,*,END=999) cbuf
      call convertToCapital(cbuf,ckey)
      if ( ckey==keyword ) then
        backspace(unit_number)
        read(unit_number,*) cbuf,keyword_variable
        hasKeyword=.true.
        exit
      endif
    enddo
999 continue
    if ( hasKeyword ) write(*,'(1x,A10," : ",F20.12)') keyword,keyword_variable
  END SUBROUTINE IOTools_readRealKeyword

  SUBROUTINE IOTools_readRealVectorKeyword(keyword,unit_number,keyword_variable)
    implicit none
    character(*),intent(IN) :: keyword
    integer,intent(IN) :: unit_number
    real(8),intent(OUT) :: keyword_variable(1:3)
    logical :: hasKeyword=.false.
    integer :: i
    character(10) :: cbuf,ckey
    integer :: keyword_length
    keyword_length=len_trim(keyword)
    rewind unit_number
    do i=1,10000
      read(unit_number,*,END=999) cbuf
      call convertToCapital(cbuf,ckey)
      if ( ckey==keyword ) then
        backspace(unit_number)
        read(unit_number,*) cbuf,keyword_variable
        hasKeyword=.true.
        exit
      endif
    enddo
999 continue
    if ( hasKeyword ) write(*,'(1x,A10," : ",3F20.12)') keyword,keyword_variable
  END SUBROUTINE IOTools_readRealVectorKeyword

  SUBROUTINE IOTools_readIntegerKeyword(keyword,unit_number,keyword_variable)
    implicit none
    character(*),intent(IN) :: keyword
    integer,intent(IN) :: unit_number
    integer,intent(OUT) :: keyword_variable
    logical :: hasKeyword=.false.
    integer :: i
    character(10) :: cbuf,ckey
    integer :: keyword_length
    keyword_length=len_trim(keyword)
    rewind unit_number
    do i=1,10000
      read(unit_number,*,END=999) cbuf
      call convertToCapital(cbuf,ckey)
      if ( ckey==keyword ) then
        backspace(unit_number)
        read(unit_number,*) cbuf,keyword_variable
        hasKeyword=.true.
        exit
      endif
    enddo
999 continue
    if ( hasKeyword ) write(*,'(1x,A10," : ",I10)') keyword,keyword_variable
  END SUBROUTINE IOTools_readIntegerKeyword

  SUBROUTINE IOTools_readIntegerVectorKeyword(keyword,unit_number,keyword_variable)
    implicit none
    character(*),intent(IN) :: keyword
    integer,intent(IN) :: unit_number
    integer,intent(OUT) :: keyword_variable(1:3)
    logical :: hasKeyword=.false.
    integer :: i
    character(10) :: cbuf,ckey
    integer :: keyword_length
    keyword_length=len_trim(keyword)
    rewind unit_number
    do i=1,10000
      read(unit_number,*,END=999) cbuf
      call convertToCapital(cbuf,ckey)
      if ( ckey==keyword ) then
        backspace(unit_number)
        read(unit_number,*) cbuf,keyword_variable
        hasKeyword=.true.
        exit
      endif
    enddo
999 continue
    if ( hasKeyword ) write(*,'(1x,A10," : ",3I10)') keyword,keyword_variable
  END SUBROUTINE IOTools_readIntegerVectorKeyword

  SUBROUTINE IOTools_readLogicalKeyword(keyword,unit_number,keyword_variable)
    implicit none
    character(*),intent(IN) :: keyword
    integer,intent(IN) :: unit_number
    logical,intent(OUT) :: keyword_variable
    logical :: hasKeyword=.false.
    integer :: i
    character(10) :: cbuf,ckey
    integer :: keyword_length
    keyword_length=len_trim(keyword)
    rewind unit_number
    do i=1,10000
      read(unit_number,*,END=999) cbuf
      call convertToCapital(cbuf,ckey)
      if ( ckey==keyword ) then
        backspace(unit_number)
        read(unit_number,*) cbuf,keyword_variable
        hasKeyword=.true.
        exit
      endif
    enddo
999 continue
    if ( hasKeyword ) write(*,'(1x,A10," : ",L10)') keyword,keyword_variable
  END SUBROUTINE IOTools_readLogicalKeyword
#endif
  SUBROUTINE convertToCapital(cbuf,CKEY)
    implicit none
    character(*),intent(IN)  :: cbuf
    character(*),intent(OUT) :: CKEY
    integer :: j,k,n
    n=len_trim(cbuf)
    CKEY=""
    do j=1,n
      k=iachar( cbuf(j:j) )
      if ( k >= 97 ) k=k-32
      CKEY(j:j) = achar(k)
    end do
  END SUBROUTINE convertToCapital

ENDMODULE IOTools
