SUBROUTINE write_border(n,indx)

  implicit none
  integer,intent(IN) :: n
  character(*),intent(IN) :: indx
  character(80) :: axx

  write(axx,'(i2)') n
  axx=adjustl(axx)
  axx="(a"//axx(1:len_trim(axx))//",a)"

  write(*,axx) repeat("-",n),indx

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
