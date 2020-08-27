SUBROUTINE convert_capital(cbuf,CKEY)

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

END SUBROUTINE convert_capital


SUBROUTINE convert_to_capital( CKEY )

  implicit none

  character(*),intent(INOUT) :: CKEY
  integer :: j,k,n

  n=len_trim(CKEY)
  do j=1,n
     k=iachar( CKEY(j:j) )
     if ( k >= 97 ) k=k-32
     CKEY(j:j) = achar(k)
  end do

END SUBROUTINE convert_to_capital
