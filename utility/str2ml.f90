PROGRAM str2ml

  implicit none

  integer,parameter :: u1=11
  integer :: i,MI,MKI,MIP
  real(8) :: aa(3,3),ax
  character(4) :: cbuf
  character(20) :: file_name
  real(8) :: qf
  real(8), parameter :: PI=atan(1.d0)*4.d0
  real(8) :: c1,c2,c3
  real(8) :: d1,d2,d3
  real(8) :: dd1,dd2,dd3
  real(8) :: e1,e2,e3
  real(8) :: f1,f2,f3
  real(8) :: g1,g2,g3

! --- input

  write(*,'(A)') "-------------------------------------------------------------"
  write(*,*) "file_name(970format)=?"
  read(*,*) file_name
  write(*,*) 'Ecut[Ry]=qf=?'
  read(*,*) qf

  write(*,'(A)') "-------------------------------------------------------------"
  write(*,'(2A)')       "file_name= ",file_name
  write(*,'(A,F25.15)') "Ecut[Ry]= ",qf
  write(*,*) 
 
! --- read data ---

  open(u1,file=file_name)
  read(u1,*) cbuf,ax
  do i=1,3
     read(u1,*) cbuf,aa(1:3,i)
  end do
  close(u1)
  aa(:,:)=aa(:,:)*ax

! --- estimate

  c1=sqrt( sum(aa(:,1)**2) )
  c2=sqrt( sum(aa(:,2)**2) )
  c3=sqrt( sum(aa(:,3)**2) )
      
  d1=int(c1/PI*sqrt(qf))+1
  d2=int(c2/PI*sqrt(qf))+1
  d3=int(c3/PI*sqrt(qf))+1

  Write(*,*) ' '
  write(*,'(A,F25.15)') "-- to use Ecut[Ry] of",qf
  Write(*,*) "ML1=",c1/PI*sqrt(qf)
  Write(*,*) "ML2=",c2/PI*sqrt(qf)
  Write(*,*) "ML3=",c3/PI*sqrt(qf)
  write(*,'(A)') "However, these are not integers."
  write(*,*) 

  E1=(PI*d1/c1)**2;  F1=(PI*(d1-1)/c1)**2
  E2=(PI*d2/c2)**2;  F2=(PI*(d2-1)/c2)**2
  E3=(PI*d3/c3)**2;  F3=(PI*(d3-1)/c3)**2

  write(*,'(A)') "-------------------------------------------------------------"
  Write(*,'(A)') "       ML1           Ecut[Ry]                   diff[Ry]" 
  write(*,'(A)') "-------------------------------------------------------------"
  Write(*,'(I10,F25.10,F25.10)') int(d1),    E1, E1-qf
  Write(*,'(I10,F25.10,F25.10)') int(d1-1),  F1, F1-qf
  Write(*,*) 

  write(*,'(A)') "-------------------------------------------------------------"
  Write(*,'(A)') "       ML2           Ecut[Ry]                   diff[Ry]" 
  write(*,'(A)') "-------------------------------------------------------------"
  Write(*,'(I10,F25.10,F25.10)') int(d2),    E2, E2-qf
  Write(*,'(I10,F25.10,F25.10)') int(d2-1),  F2, F2-qf
  Write(*,*) 

  write(*,'(A)') "-------------------------------------------------------------"
  Write(*,'(A)') "       ML3           Ecut[Ry]                   diff[Ry]" 
  write(*,'(A)') "-------------------------------------------------------------"
  Write(*,'(I10,F25.10,F25.10)') int(d3),    E3, E3-qf
  Write(*,'(I10,F25.10,F25.10)') int(d3-1),  F3, F3-qf
  Write(*,*) 

 
  IF ( abs(E1-qf) .le. abs(F1-qf) ) then
      dd1=d1;   G1=E1
  else
      dd1=d1-1; G1=F1
  end IF

  IF ( abs(E2-qf) .le. abs(F2-qf) ) then
      dd2=d2;   G2=E2
  else
      dd2=d2-1; G2=F2
  end IF

  IF ( abs(E3-qf) .le. abs(F3-qf) ) then
      dd3=d3;   G3=E3
  else
      dd3=d3-1; G3=F3
  end IF

  write(*,'(A)') "-------------------------------------------------------------"
  Write(*,'(A)') "     RECOMMEND_value         Ecut[Ry]"
  write(*,'(A)') "-------------------------------------------------------------"
  Write(*,'(A,I10,F28.15)') "ML1       ", int(dd1), G1
  Write(*,'(A,I10,F28.15)') "ML2       ", int(dd2), G2
  Write(*,'(A,I10,F28.15)') "ML3       ", int(dd3), G3

END PROGRAM  str2ml
