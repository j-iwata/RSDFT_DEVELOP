PROGRAM str2xyz

  implicit none

  integer,parameter :: u1=11,iu=10
  real(8),parameter :: ab=0.529177d0
  integer :: i,j,MI,MKI,MIP
  real(8) :: aa(3,3),ax
  integer :: zatom(300)
  character(100) :: cbuf
  integer,allocatable :: Kion(:)
  real(8),allocatable :: asi(:,:),rsi(:,:)
  real(8),allocatable :: asi2(:,:)
  character(10) :: dummy1,dummy2,dummy3,dummy4
  character(20) :: file_name
  integer :: icy, itlin
  character(2) :: ATOM_NAME(300)
  integer :: p1,   p2,   p3,   p
  real(8) :: p1_b, p2_b, p3_b
  real(8) :: p1_e, p2_e, p3_e
  integer :: q1_b, q2_b, q3_b
  integer :: q1_e, q2_e, q3_e
  integer :: i1,i2,i3
  integer :: flag
  integer :: MIP_count, MI1, MI2
  integer :: LC=0, il
  real(8) :: qq1(3),qq2(3)

! --- input

  write(*,*) "file_name(970format)=?"
  read(*,*) file_name
  write(*,*) 'P1_BEGIN, P1_END=? [e.g. -0.5,0.5]'
  read(*,*) p1_b, p1_e
  write(*,*) 'P2_BEGIN, P1_END=? [e.g. -0.5,0.5]'
  read(*,*) p2_b, p2_e
  write(*,*) 'P3_BEGIN, P3_END=? [e.g. -0.5,0.5]'
  read(*,*) p3_b, p3_e
  write(*,*) 

  write(*,'(A)') "--------------------------- FILE ----------------------------"
  write(*,'(2A)') "file_name= ",file_name
  write(*,*) 
  open(u1,file=file_name)

  write(*,'(A)') "--------------------------- RANGE ---------------------------"
  write(*,'(A,2(3X,F25.15))') 'P1_[ BEGIN, END ]=',p1_b,p1_e
  write(*,'(A,2(3X,F25.15))') 'P2_[ BEGIN, END ]=',p2_b,p2_e
  write(*,'(A,2(3X,F25.15))') 'P3_[ BEGIN, END ]=',p3_b,p3_e
  write(*,*) 

  q1_b=int(p1_b)-1; q1_e=int(p1_e)+1
  q2_b=int(p2_b)-1; q2_e=int(p2_e)+1
  q3_b=int(p3_b)-1; q3_e=int(p3_e)+1
!  write(*,'(A)') "-------------------------ENOUGH  BOX-------------------------"
!  write(*,'(A,2(3X,I5))') 'Q1_[ BEGIN, END ]=',q1_b,q1_e
!  write(*,'(A,2(3X,I5))') 'Q2_[ BEGIN, END ]=',q2_b,q2_e
!  write(*,'(A,2(3X,I5))') 'Q3_[ BEGIN, END ]=',q3_b,q3_e
!  write(*,*) 

  p1=q1_e-q1_b
  p2=q2_e-q2_b
  p3=q3_e-q3_b
  p=p1*p2*p3
!  write(*,'(A)') "------------------------- BOX  SIZE -------------------------"
!  write(*,'(A,4(3X,I5))') 'P1, P2, P3, P',p1,p2,p3,p
!  write(*,*)
!  write(*,*)


! --- length_log
  open(unit=100,file="length.dat")
  write(*,*) "MEASURE LENGTH? (0=NO, 1=YES)"
  read(*,*) il
  if (il==1) then
   write(*,*) "ATOM1="
   read(*,*) MI1
   write(*,*) "ATOM2="
   read(*,*) MI2
  else
   MI1=0
   MI2=0
  end if


! --- target file ---

  open(iu,file='str.xyz')

! --- periodic table --
  ATOM_NAME="XX"
  ATOM_NAME (1:2)  =(/ "H ",                                     "He" /)
  ATOM_NAME (3:10) =(/ "Li", "Be", "B ", "C ", "N ", "O ", "F ", "Ne" /)
  ATOM_NAME(11:18) =(/ "Na", "Mg", "Al", "Si", "P ", "S ", "Cl", "Ar" /)
  ATOM_NAME(19:20) =(/ "K ", "Ca" /)
  ATOM_NAME(21:30) =(/ "Sc", "Ti", "V ", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn" /)
  ATOM_NAME(31:36) =(/ "Ga", "Ge", "As", "Se", "Br", "Kr"/)
 
! --- read input data ---

  read(u1,*) cbuf,ax
  do i=1,3
     read(u1,*) cbuf,aa(1:3,i)
  end do
  aa=ax*aa*ab

  read(u1,*)
  read(u1,*) MKI,MI,zatom(1:MKI)

  allocate( asi(3,MI),asi2(3,MI),Kion(MI),rsi(3,MI) )

!---

  read(u1,*,END=999) cbuf
  flag=0; if (cbuf(1:9)=="#_STRLOG_") flag=1
  backspace(u1)

!---
777 continue  

  if (flag==1) read(u1,*,END=999) dummy1, dummy2, dummy3, dummy4, icy, itlin

  do i=1,MI
     read(u1,*,END=999) Kion(i),asi(:,i)
  end do

!--
  do j=1,MI
   do i=1,3

    do while ( asi(i,j) >= 1.0 ) 
     asi(i,j)=asi(i,j)-1.0
    end do
    do while ( asi(i,j) <  0.0 ) 
     asi(i,j)=asi(i,j)+1.0
    end do

   end do !i
  end do !j


! --- xyz file

!-- count

  MIP=0
  do i1=q1_b,q1_e
  do i2=q2_b,q2_e
  do i3=q3_b,q3_e

  asi2=asi
  asi2(1,:)=asi2(1,:)+dble(i1)
  asi2(2,:)=asi2(2,:)+dble(i2)
  asi2(3,:)=asi2(3,:)+dble(i3)
  do i=1,MI
   if ( (p1_b <= asi2(1,i)) .AND. (asi2(1,i) < p1_e) .AND. & 
        (p2_b <= asi2(2,i)) .AND. (asi2(2,i) < p2_e) .AND. & 
        (p3_b <= asi2(3,i)) .AND. (asi2(3,i) < p3_e) ) MIP=MIP+1 
  end do !i 
  end do !i3
  end do !i2
  end do !i1

  write(*,*) "atom/cell, atom/range, icy, itlin  =", MI, MIP, icy, itlin


!--write

  write(iu,'(I5)') MIP
  if (flag==1) then
    write(iu,'(A,2(2X,I6))')  "icy,itlin=",icy,itlin
  else
    write(iu,'(A)') "comment"     
  end if

  MIP_count=0

  do i1=q1_b,q1_e
  do i2=q2_b,q2_e
  do i3=q3_b,q3_e

  asi2=asi
  asi2(1,:)=asi2(1,:)+dble(i1)
  asi2(2,:)=asi2(2,:)+dble(i2)
  asi2(3,:)=asi2(3,:)+dble(i3)

  rsi=matmul( aa,asi2 )

  do i=1,MI   
   if ( (p1_b <= asi2(1,i)) .AND. ( asi2(1,i) < p1_e) .AND. & 
        (p2_b <= asi2(2,i)) .AND. ( asi2(2,i) < p2_e) .AND. & 
        (p3_b <= asi2(3,i)) .AND. ( asi2(3,i) < p3_e) ) then
      MIP_count=MIP_count+1
      write(iu,'(A2,3(2X,F25.15))') ATOM_NAME(zatom(Kion(i))),rsi(:,i)
      if (MIP_count==MI1) qq1(:)=rsi(:,i) 
      if (MIP_count==MI2) qq2(:)=rsi(:,i) 
   end if
  end do
 
  end do !i3
  end do !i2
  end do !i1
  
  LC=LC+1
  if (il==1) Write(100,'(I6,1X,F25.15,2(1X,I6))') LC, sqrt( sum(qq1(:)-qq2(:))**2 ), MI1,MI2

  if (flag==1) goto 777

!----
999 continue

  deallocate( rsi,Kion,asi )
  close(iu)
  close(u1)

END PROGRAM  str2xyz
