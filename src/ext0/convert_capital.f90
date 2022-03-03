SUBROUTINE convert_capital(cbuf,CKEY)

  implicit none

  character(*),intent(IN)  :: cbuf
  character(*),intent(OUT) :: CKEY
  integer :: j,k,n

  n=len_trim(cbuf)
  CKEY=""
  do j=1,n
    k=iachar( cbuf(j:j) )
    if ( 97 <= k .and. k <= 122 ) k=k-32
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


subroutine convert_to_capital_array( ndata, CKEY )

  implicit none
  integer,intent(in) :: ndata
  character(*),intent(inout) :: CKEY(ndata)
  integer :: j,k,n,idata

  do idata = 1, ndata

    n = len_trim( CKEY(idata) )
    do j = 1, n
      k = iachar( CKEY(idata)(j:j) )
      if ( k >= 97 ) k=k-32
      CKEY(idata)(j:j) = achar(k)
    end do

  end do

end subroutine convert_to_capital_array
