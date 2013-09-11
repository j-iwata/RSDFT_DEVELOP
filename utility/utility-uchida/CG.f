       Program CG
      !-----------------------------------
       Implicit NONE

       Integer :: j1,m1,j2,m2,j3,m3
       Integer :: s,z,N
       Real(8) :: cg_coef
       Real(8) :: v
       Real(8) :: A1
       Real(8) :: A2
       Real(8) :: A3
       Real(8) :: A4
       Real(8) :: A5
       Real(8) :: A6
       Real(8) :: A7
       Real(8) :: A8
       Real(8) :: A9
       Real(8) :: A10
       Real(8) :: A11

       Real(8) :: F1
       Real(8) :: F2
       Real(8) :: F3

       Real(8) :: C1
       Real(8) :: C2
       Real(8) :: C3
       Real(8) :: C4
       Real(8) :: C5
       Real(8) :: C6

 
       Integer :: I1
       Integer :: I2
       Integer :: I3
       Integer :: I4
       Integer :: I5

       Integer :: counter

       counter=0

       Write(*,'(A)') 'Counter, j1,m1,j2,m2,j3,m3,  :  ,cg_coef'

       Do j1=0,2
        Do m1=-j1,j1
         Do j2=1,1
          Do m2=-j2,j2
           Do j3=abs(j1-j2),j1+j2
            Do m3=-j3,j3

     
!--
       cg_coef=0.d0

!--test
!       Do N=1,10
!         call factorial(v,N)
!         write(*,*) 'N, N!=', N, v
!       End Do !N

 
       if ( (m1+m2)==m3 ) cycle
!-------------------------------------------------- 
 
       s=j1+j2+j3 

       A1=2*j3+1
       call factorial(A2,s-2*j3)
       call factorial(A3,s-2*j2)
       call factorial(A4,s-2*j1)
       call factorial(A5,s+1)

       call factorial(A6,j1+m1)
       call factorial(A7,j1-m1)
       call factorial(A8,j2+m2)
       call factorial(A9,j2-m2)
       call factorial(A10,j3+m3)
       call factorial(A11,j3-m3)      
       
       F1=A1*sqrt(A1*A2*A3*A4/A5*A6*A7*A8*A9*A10*A11)    

!--
       F2=0.d0
       z=1
!--

 700   continue

       I1=j1+j2-j3-z     
       I2=j1-m1-z     
       I3=j2+m2-z   
       I4=j3-j2+m1+z     
       I5=j3-j1-m2+z
 
       if (        I1<0 
     &         .or.I2<0 
     &         .or.I3<0 
     &         .or.I4<0 
     &         .or.I5<0 ) goto 800

 
      call factorial(C1,I1)
      call factorial(C2,I2)
      call factorial(C3,I3)
      call factorial(C4,I4)
      call factorial(C5,I5)
      call factorial(C6,z)

      F3=((-1)**z)/(C6*C1*C2*C3*C4*C5) 

      F2=F2+F3
      z=z+1
      goto 700

!------------------------------------------------- 
 800   continue

       cg_coef=F1*F2
  
!-------------------------------------------------- 
 900   continue
 
        
       if (abs(cg_coef)>1.0D-10) then
           counter=counter+1
           Write(*,'(7I,A,F25.20)') 
     &     Counter, j1,m1,j2,m2,j3,m3,' :  ',cg_coef
       end if

       End Do
       End Do
       End Do
       End Do
       End Do
       End Do

       End Program CG
!--------------------------------
       Subroutine factorial(v,N)
      !------------------------- 
       Implicit NONE
     
       real(8) :: v   
       integer :: N,i

       v=1.d0
 
       if (N<0) then
         write(*,*) 'N<0',N
         stop
       end if

       Do i=1,N
        v=V*dble(i)
       End Do

       Return
       End Subroutine factorial
!--------------------------------
