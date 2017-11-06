PROGRAM grid_estimate

  implicit none

  integer,parameter :: u0=970, u1=5
  integer,parameter :: nmax=1000
  integer :: i,n,ngg,i1,i2,i3,ngg0,m1,n1,m2,n2,m3,n3
  integer :: Ngrid(3),NGgrid(3)
  integer,allocatable :: lll(:,:,:)
  real(8),parameter :: pi=3.141592653589793d0
  real(8) :: ax,aa(3,3),al(3),bb(3,3),bl(3)
  real(8) :: ecut,ggmax,gg,gx,gy,gz,Hgrid(3)
  real(8) :: ecut_rsdft
  character(len=2) :: cbuf

  read(u0,*) cbuf, ax
  read(u0,*) cbuf, aa(:,1)
  read(u0,*) cbuf, aa(:,2)
  read(u0,*) cbuf, aa(:,3)

  aa=ax*aa

  call calc_length( aa, al )
!  do i=1,3
!     write(*,'(1x,"al(",i1,")",f20.15)') i,al(i)
!  end do

  call calc_bb( aa, bb )

  call calc_length( bb, bl )
!  do i=1,3
!     write(*,'(1x,"bl(",i1,")",f20.15)') i,bl(i)
!  end do

  write(*,*) "Input Ecut(Ry) = ?"
  read(u1,*) ecut



  ngg=1
  ggmax=0.0d0

  do n=1,nmax
     ngg0=ngg
     do i3=-n,n
     do i2=-n,n
     do i1=-n,n
        if ( abs(i1) < n .and. abs(i2) < n .and. abs(i3) < n ) cycle
        gx=i1*bb(1,1)+i2*bb(1,2)+i3*bb(1,3)
        gy=i1*bb(2,1)+i2*bb(2,2)+i3*bb(2,3)
        gz=i1*bb(3,1)+i2*bb(3,2)+i3*bb(3,3)
        gg=gx*gx+gy*gy+gz*gz
        if ( gg <= ecut ) then
           ngg=ngg+1
           ggmax=max(gg,ggmax)
        end if
     end do
     end do
     end do
     if ( ngg == ngg0 ) exit
  end do

  write(*,*) "# of G vectors=",ngg
  write(*,*) "GGmax =",ggmax

  allocate( lll(-n:n,-n:n,-n:n) ) ; lll=0

  do i3=-n,n
  do i2=-n,n
  do i1=-n,n
     gx=i1*bb(1,1)+i2*bb(1,2)+i3*bb(1,3)
     gy=i1*bb(2,1)+i2*bb(2,2)+i3*bb(2,3)
     gz=i1*bb(3,1)+i2*bb(3,2)+i3*bb(3,3)
     gg=gx*gx+gy*gy+gz*gz
     if ( gg <= ecut ) then
        lll(i1,i2,i3)=1
        lll(-i1,-i2,-i3)=1
     end if
  end do
  end do
  end do

!  write(*,*) count(lll>0)

  do i=-n,n
     if ( any(lll(i,:,:)/=0) ) then
        m1=i
        exit
     end if
  end do
  do i=n,-n,-1
     if ( any(lll(i,:,:)/=0) ) then
        n1=i
        exit
     end if
  end do
  do i=-n,n
     if ( any(lll(:,i,:)/=0) ) then
        m2=i
        exit
     end if
  end do
  do i=n,-n,-1
     if ( any(lll(:,i,:)/=0) ) then
        n2=i
        exit
     end if
  end do
  do i=-n,n
     if ( any(lll(:,:,i)/=0) ) then
        m3=i
        exit
     end if
  end do
  do i=n,-n,-1
     if ( any(lll(:,:,i)/=0) ) then
        n3=i
        exit
     end if
  end do

!  write(*,*) m1,n1
!  write(*,*) m2,n2
!  write(*,*) m3,n3

  Ngrid(1)=2*n1
  Ngrid(2)=2*n2
  Ngrid(3)=2*n3

  Hgrid(:)=al(:)/Ngrid(:)

  write(*,'(1x,"Ngrid",2(2x,a20))') "Hgrid(bohr)","Ecut(Ry)"
  do i=1,3
     write(*,'(1x,i5,2(2x,f20.15))') Ngrid(i),Hgrid(i),(pi/Hgrid(i))**2
  end do
  write(*,*) "Total # of grid points =",Ngrid(1)*Ngrid(2)*Ngrid(3)

  NGgrid=(Ngrid-1)/2
  call calc_rsdft_ecut( bb, NGgrid, ecut_rsdft ) 
  write(*,*) "Ecut_RSDFT(Ry)=",ecut_rsdft

  deallocate( lll )

END PROGRAM grid_estimate


SUBROUTINE calc_bb(aa,bb_out)
  implicit none
  real(8),intent(IN)  :: aa(3,3)
  real(8),intent(OUT) :: bb_out(3,3)
  real(8) :: Vaa,PI2
  PI2 = 2.d0*acos(-1.d0)
  Vaa = aa(1,1)*aa(2,2)*aa(3,3)+aa(1,2)*aa(2,3)*aa(3,1) &
       +aa(1,3)*aa(2,1)*aa(3,2)-aa(1,3)*aa(2,2)*aa(3,1) &
       -aa(1,2)*aa(2,1)*aa(3,3)-aa(1,1)*aa(2,3)*aa(3,2)
  bb_out(1,1) = aa(2,2)*aa(3,3) - aa(3,2)*aa(2,3)
  bb_out(2,1) = aa(3,2)*aa(1,3) - aa(1,2)*aa(3,3)
  bb_out(3,1) = aa(1,2)*aa(2,3) - aa(2,2)*aa(1,3)
  bb_out(1,2) = aa(2,3)*aa(3,1) - aa(3,3)*aa(2,1)
  bb_out(2,2) = aa(3,3)*aa(1,1) - aa(1,3)*aa(3,1)
  bb_out(3,2) = aa(1,3)*aa(2,1) - aa(2,3)*aa(1,1)
  bb_out(1,3) = aa(2,1)*aa(3,2) - aa(3,1)*aa(2,2)
  bb_out(2,3) = aa(3,1)*aa(1,2) - aa(1,1)*aa(3,2)
  bb_out(3,3) = aa(1,1)*aa(2,2) - aa(2,1)*aa(1,2)
  bb_out(:,:)=bb_out(:,:)*PI2/Vaa
END SUBROUTINE calc_bb


SUBROUTINE calc_length( v, vl )
  implicit none
  real(8),intent(IN) :: v(3,3)
  real(8),intent(OUT) :: vl(3)
  integer :: i
  do i=1,3
     vl(i) = sqrt( sum(v(:,i)**2) )
  end do
END SUBROUTINE calc_length


SUBROUTINE calc_rsdft_ecut( bb, NGgrid, Ecut )
  implicit none
  real(8),intent(IN) :: bb(3,3)
  integer,intent(IN) :: NGgrid(3)
  real(8),intent(OUT) :: Ecut
  real(8) :: b1,b2,b3,b12,c,s,r1,r2,r3,b0(3),Gcut
  b1=sqrt(sum(bb(1:3,1)**2))*(2*NGgrid(1)+1)
  b2=sqrt(sum(bb(1:3,2)**2))*(2*NGgrid(2)+1)
  b3=sqrt(sum(bb(1:3,3)**2))*(2*NGgrid(3)+1)
  b12=sum(bb(1:3,1)*bb(1:3,2))
  c=b12/(b1*b2)*(2*NGgrid(1)+1)*(2*NGgrid(2)+1)
  s=sqrt(1.d0-c**2)
  r1=s*b1
  r2=s*b2
  b0(1)=bb(2,1)*bb(3,2)-bb(3,1)*bb(2,2)
  b0(2)=bb(3,1)*bb(1,2)-bb(1,1)*bb(3,2)
  b0(3)=bb(1,1)*bb(2,2)-bb(2,1)*bb(1,2)
  r3=abs( sum(bb(1:3,3)*b0(1:3))*(2*NGgrid(3)+1)/sqrt(sum(b0**2)) )
  Gcut=min(r1,r2,r3)*0.5d0
  Ecut=Gcut**2
END SUBROUTINE calc_rsdft_ecut
