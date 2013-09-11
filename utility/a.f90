
implicit none
integer :: i,j
real(8) :: pi,d,t,aL,z,a1,b1,c1,cs,a
real(8) :: c,b,a0,aa(3,3),v(3)

pi=4.d0*atan(1.d0)
d=57.233d0
z=0.46724d0
a0=8.4845d0

t=pi*d/180.d0
write(*,*) "t,cos(t)",t,cos(t)

a=5.996841785723262d0
b=0.2500348855797364d0
aa(1:3,1:3)=a
aa(1,1)=b ; aa(2,2)=b ; aa(3,3)=b
write(*,*) "a=",a
write(*,*) "b,b/a=",b,b/a
 
aL=sqrt( 3.d0*(2.d0*a+b)**2 )

v(1:3)=aa(1:3,1)+aa(1:3,2)+aa(1:3,3)
write(*,*) "|a1+a2+a3|",aL,sqrt(sum(v(1:3)**2))

c=(1.d0-z)*aL
write(*,*) c
write(*,*) c/aL,c/aL*0.5d0

do i=1,3
   do j=1,i
      write(*,*) i,j,sum( aa(1:3,i)*aa(1:3,j) )/a0**2
   end do
end do   


a1=9.d0
b1=-2.d0*a0**2*cos(t)-4.d0*a0**2
c1=a0**4*cos(t)**2

write(*,*) b1**2-4.d0*a1*c1
write(*,*) (-b1+sqrt(b1**2-4.d0*a1*c1) )/(2.d0*a1)
write(*,*) (-b1-sqrt(b1**2-4.d0*a1*c1) )/(2.d0*a1)

v(1)= (-b1+sqrt(b1**2-4.d0*a1*c1) )/(2.d0*a1)
v(2)= (-b1-sqrt(b1**2-4.d0*a1*c1) )/(2.d0*a1)
write(*,*) sqrt(v(1)),sqrt(v(2))
write(*,*) sqrt(a0**2-2.d0*v(1)),sqrt(a0**2-2.d0*v(2))

write(*,*) ( a0**2*cos(t)-v(1) )/(2.d0*sqrt(v(1))),( a0**2*cos(t)-v(2) )/(2.d0*sqrt(v(2)))

stop
end

