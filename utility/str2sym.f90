PROGRAM str2sym

  implicit none

  integer,parameter :: u1=11,iu=10
  real(8),parameter :: ab=0.529177d0
  integer :: i,j,k,MI,MKI,MIP,nsym
  real(8) :: aa(3,3),ax,RR0(3),RR1(3),RR2(3),RR3(3)
  integer :: zatom(300)
  character(100) :: cbuf
  integer,allocatable :: Kion(:)
  real(8),allocatable :: asi(:,:),XX(:,:)
  real(8),allocatable :: asi2(:,:)
  character(10) :: dummy1,dummy2,dummy3,dummy4
  character(20) :: file_name
  integer :: icy, itlin
  character(2) :: ATOM_NAME(300)
  integer :: i1,i2,i3
  integer :: gcounter, flag
  real(8) :: R1s,R2s,R3s, RRR
  real(8) :: R1p,R2p,R3p
  real(8) :: R1q,R2q,R3q
  Integer :: A11,A12,A13
  Integer :: A21,A22,A23
  Integer :: A31,A32,A33
  Integer :: AAA(9,999)

! --- input
  write(*,'(A)') "-------------------------------------------------------------"
  write(*,'(2A)') "file_name(970format)=?"
  Read(*,*) file_name

  write(*,'(A)') "--------------------------- FILE ----------------------------"
  write(*,'(2A)') "file_name= ",file_name
  write(*,*) 
  open(u1,file=file_name)

! --- target file ---

  open(iu,file='sym.dat')

! --- read input data ---

  read(u1,*) cbuf,ax
  do i=1,3
     read(u1,*) cbuf,aa(1:3,i)
  end do
  aa=ax*aa*ab

  read(u1,*)
  read(u1,*) MKI,MI,zatom(1:MKI)

  allocate( asi(3,MI),asi2(3,MI),Kion(MI) )

!---
666 continue

  read(u1,*) cbuf
  if (cbuf(1:1)=="#") goto 666
  backspace(u1)

!---
  do i=1,MI
     read(u1,*) Kion(i),asi(:,i)
  end do

  allocate( XX(3,MI) ) ; XX=0.0d0

  do i=1,MI
      XX(1:3,i)=asi(1:3,i)
      do j=1,3
        do  while ( XX(j,i) >  0.5d0 )
             XX(j,i)=XX(j,i)-1.0d0
        end do
        do  while ( XX(j,i) <= -0.5d0 )
            XX(j,i)=XX(j,i)+1.0d0
        end do
      end do
  end do ! i

  nsym=0 

          do A11=1,-1,-1
          do A12=1,-1,-1
          do A13=1,-1,-1
          do A21=1,-1,-1
          do A22=1,-1,-1
          do A23=1,-1,-1
          do A31=1,-1,-1
          do A32=1,-1,-1
          do A33=1,-1,-1

             gcounter=0

             do i=1,MI

                RR0(1:3) = XX(1:3,i)

                flag = 0

                do j=1,MI

                   if ( Kion(i) /= Kion(j) ) cycle

                   RR1(1:3) = XX(1:3,j)

                   RR2(1) = dble(A11)*RR1(1) + dble(A12)*RR1(2) + dble(A13)*RR1(3)
                   RR2(2) = dble(A21)*RR1(1) + dble(A22)*RR1(2) + dble(A23)*RR1(3)
                   RR2(3) = dble(A31)*RR1(1) + dble(A32)*RR1(2) + dble(A33)*RR1(3)

                   do k=1,3
                      do  while ( RR2(k) >  0.5d0 ) 
                         RR2(k)=RR2(k)-1.0d0
                      end do
                      do while ( RR2(k) <= -0.5d0 ) 
                         RR2(k)=RR2(k)+1.0d0
                      end do
                   end do ! k

                   RRR = sqrt( (RR2(1)-RR0(1))**2.0 &
                             + (RR2(2)-RR0(2))**2.0 + (RR2(3)-RR0(3))**2.0 )
                   if ( RRR < 1.d-6 ) then
                      flag=1
                   end if

                end do ! j

                gcounter = gcounter + flag

             end do ! i
           

             if ( gcounter == MI ) then

                R1q = sqrt( AA(1,1)**2.0 + AA(2,1)**2.0 + AA(3,1)**2.0 )
                R2q = sqrt( AA(1,2)**2.0 + AA(2,2)**2.0 + AA(3,2)**2.0 )
                R3q = sqrt( AA(1,3)**2.0 + AA(2,3)**2.0 + AA(3,3)**2.0 )

                R1s = dble(A11)*1.0 + dble(A12)*0.0 + dble(A13)*0.0
                R2s = dble(A21)*1.0 + dble(A22)*0.0 + dble(A23)*0.0
                R3s = dble(A31)*1.0 + dble(A32)*0.0 + dble(A33)*0.0
                R1p = sqrt( (AA(1,1)*R1s + AA(1,2)*R2s + AA(1,3)*R3s)**2.0 &
                          + (AA(2,1)*R1s + AA(2,2)*R2s + AA(2,3)*R3s)**2.0 &
                          + (AA(3,1)*R1s + AA(3,2)*R2s + AA(3,3)*R3s)**2.0   )

                R1s = dble(A11)*0.0 + dble(A12)*1.0 + dble(A13)*0.0
                R2s = dble(A21)*0.0 + dble(A22)*1.0 + dble(A23)*0.0
                R3s = dble(A31)*0.0 + dble(A32)*1.0 + dble(A33)*0.0
                R2p = sqrt( (AA(1,1)*R1s + AA(1,2)*R2s + AA(1,3)*R3s)**2.0 &
                          + (AA(2,1)*R1s + AA(2,2)*R2s + AA(2,3)*R3s)**2.0 &
                          + (AA(3,1)*R1s + AA(3,2)*R2s + AA(3,3)*R3s)**2.0  )

                R1s = dble(A11)*0.0 + dble(A12)*0.0 + dble(A13)*1.0
                R2s = dble(A21)*0.0 + dble(A22)*0.0 + dble(A23)*1.0
                R3s = dble(A31)*0.0 + dble(A32)*0.0 + dble(A33)*1.0
                R3p = sqrt( (AA(1,1)*R1s + AA(1,2)*R2s + AA(1,3)*R3s)**2.0 &
                          + (AA(2,1)*R1s + AA(2,2)*R2s + AA(2,3)*R3s)**2.0 &
                          + (AA(3,1)*R1s + AA(3,2)*R2s + AA(3,3)*R3s)**2.0   )

                if ( abs(R1q-R1p) <= 1.d-6 .and. &
                     abs(R2q-R2p) <= 1.d-6 .and. &
                     abs(R3q-R3p) <= 1.d-6        ) then
                   nsym=nsym+1
                   AAA(1,nsym)=A11
                   AAA(2,nsym)=A12
                   AAA(3,nsym)=A13
                   AAA(4,nsym)=A21
                   AAA(5,nsym)=A22
                   AAA(6,nsym)=A23
                   AAA(7,nsym)=A31
                   AAA(8,nsym)=A32
                   AAA(9,nsym)=A33
                end if

             end if

          end do
          end do
          end do
          end do
          end do
          end do
          end do
          end do
          end do

          deallocate( XX )

          write(*,*) "------------------------------------------------------"
          write(*,'(I0,42X,I0)')  nsym, 1
          do i=1,nsym
             write(*,'(3(3I3,3X),6X,3I3)')  AAA(1:9,i),0,0,0
          end do !i
          write(iu,'(I0,42X,I0)')  nsym, 1
          do i=1,nsym
             write(iu,'(3(3I3,3X),6X,3I3)')  AAA(1:9,i),0,0,0
          end do !i
          write(*,*) "------------------------------------------------------"
          write(*,*) "BE CAREFUL!! Bubun heishin & Rasen are NOT considered."
          write(*,*) "BE CAREFUL!! SPIN is not considered here."
          
          close(u1)
          close(iu)
  END Program str2sym
