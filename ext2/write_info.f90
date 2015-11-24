SUBROUTINE write_border( n, indx )

  implicit none
  integer,intent(IN) :: n
  character(*),intent(IN) :: indx
  character(80) :: axx
  logical :: disp
  integer :: m=80
  integer,save :: u0=6, u1=60

  write(axx,'(i2)') m-len(indx)
  axx=adjustl(axx)
  axx="(a"//axx(1:len_trim(axx))//",a)"

  call check_disp_switch( disp, 0 )
  if ( disp ) then
     if ( n == 0 ) then
        write(u0,axx) repeat("-",m),indx
     else if ( n == 1 ) then
!        open(u1,file="RSDFT_LOG",position="append")
!        write(u1,axx) repeat("-",m),indx
!        close(u1)
     else
        write(u0,axx) repeat("-",m),indx
     end if
     if ( u0 /= u1 ) then
        open(u1,file="RSDFT_LOG",position="append")
        write(u1,axx) repeat("-",m),indx
        close(u1)
     end if
  end if

END SUBROUTINE write_border


SUBROUTINE check_disp_switch( disp_switch, i )
  implicit none
  logical,intent(INOUT) :: disp_switch
  integer,intent(IN)    :: i
  logical,save :: disp_switch_save=.false.
  if ( i == 0 ) then
     disp_switch = disp_switch_save
  else if ( i == 1 ) then
     disp_switch_save = disp_switch
  end if
END SUBROUTINE check_disp_switch
