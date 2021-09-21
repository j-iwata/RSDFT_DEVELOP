program plot_vrho

  implicit none

  integer,parameter :: u0=80, u1=1, u2=970, iu=2
  integer :: ML,ML1,ML2,ML3,zatom(9)
  integer :: i1,i2,i3,i,s,MI,MKI
  integer :: m1,m2,m3,n1,n2,n3
  integer,allocatable :: LL(:,:),Kion(:)
  real(8),allocatable :: rho(:,:,:,:),vloc(:,:,:,:)
  real(8),allocatable :: vh(:,:,:),vxc(:,:,:),rtmp(:)
  real(8),allocatable :: asi(:,:),rsi(:,:)
  real(8) :: ax,aa(3,3),Va,dV,efermi
  character(64) :: cbuf
  logical :: flag_versioned=.false.

! ---

  open(u0,file='vrho.dat1',status='old',form='unformatted')

  read(u0) cbuf
  if ( cbuf(1:7) == "version" ) then
     write(*,*) "file format version: "//cbuf
     flag_versioned=.true.
  else
     write(*,*) "non-versioned (old) file format"
     flag_versioned=.false.
     rewind(u0)
 end if

  read(u0) ML,ML1,ML2,ML3
  write(*,*) "ML,ML1,ML2,ML3=",ML,ML1,ML2,ML3

  allocate( LL(3,ML) ) ; LL=0.0d0

  read(u0) LL

  if ( flag_versioned ) then

     read(u0) cbuf ; write(*,*) cbuf
     read(u0) aa
     read(u0) cbuf ; write(*,*) cbuf
     read(u0) MKI,MI
     allocate( asi(3,MI) ) ; asi=0.0d0
     allocate( rsi(3,MI) ) ; rsi=0.0d0
     allocate( Kion(MI)  ) ; Kion=0
     read(u0) asi
     read(u0) Kion
     read(u0) zatom(1:MKI)
     read(u0) cbuf ; write(*,*) cbuf
     read(u0) efermi
     write(*,*) "efrmi",efermi

  else

     read(u2,*) cbuf,ax
     read(u2,*) cbuf,aa(1:3,1)
     read(u2,*) cbuf,aa(1:3,2)
     read(u2,*) cbuf,aa(1:3,3)
     read(u2,*)
     read(u2,*) MKI,MI,zatom(1:MKI)

     allocate( asi(3,MI) ) ; asi=0.0d0
     allocate( rsi(3,MI) ) ; asi=0.0d0
     allocate( Kion(MI)  ) ; Kion=0

     do i=1,MI
        read(u2,*) Kion(i),asi(1:3,i)
     end do

     aa(:,:)=ax*aa(:,:)

  end if

  Va = aa(1,1)*aa(2,2)*aa(3,3)+aa(1,2)*aa(2,3)*aa(3,1) &
      +aa(1,3)*aa(2,1)*aa(3,2)-aa(1,3)*aa(2,2)*aa(3,1) &
      -aa(1,2)*aa(2,1)*aa(3,3)-aa(1,1)*aa(2,3)*aa(3,2)

! ---

!  do i=1,MI
!     asi(3,i) = asi(3,i) + 0.5d0
!     if ( asi(3,i) > 1.0d0 ) asi(3,i) = asi(3,i) - 1.0d0
!  end do

! ---

  rsi(:,:)=matmul( aa, asi )

! ---

  m1 = minval( LL(1,:) )
  m2 = minval( LL(2,:) )
  m3 = minval( LL(3,:) )
  n1 = maxval( LL(1,:) )
  n2 = maxval( LL(2,:) )
  n3 = maxval( LL(3,:) )

  allocate( rtmp(ML) ) ; rtmp=0.0d0
  allocate( rho(m1:n1,m2:n2,m3:n3,2)  ) ; rho=0.0d0
  allocate( vloc(m1:n1,m2:n2,m3:n3,2) ) ; vloc=0.0d0
  allocate( vh(m1:n1,m2:n2,m3:n3)     ) ; vh=0.0d0
  allocate( vxc(m1:n1,m2:n2,m3:n3)    ) ; vxc=0.0d0

  dV=abs(Va)/dble(ML)

  write(*,*) "ML=",ML

  do s=1,2

     read(u0,END=10) rtmp(:) ; write(*,*) "rho",sum(rtmp**2)
     do i=1,ML
        rho(LL(1,i),LL(2,i),LL(3,i),s)=rtmp(i)
     end do
     write(*,*) "for SPIN = ", s, " : sum rho  = ", sum( rho(:,:,:,s) )*dV

!     call gen_cube( s, rho(:,:,:,s), "rho" )

     read(u0) rtmp(:) ; write(*,*) "vloc",sum(rtmp**2)
     do i=1,ML
        vloc(LL(1,i),LL(2,i),LL(3,i),s)=rtmp(i)
     end do

!     call gen_cube( s, vloc(:,:,:,s), "vloc" )
     call planer_averaged_plot( s, vloc(:,:,:,s), "vloc" )

     if ( s == 2 ) then
        rho(:,:,:,1) = rho(:,:,:,1) - rho(:,:,:,2)
!        call gen_cube( 1, rho(:,:,:,1), "dspin" )
        write(*,*) "dspin(max,min)=",maxval(rho(:,:,:,1)),minval(rho(:,:,:,1))
        write(*,*) "sum  dspin  = ", sum( rho(:,:,:,1) )*dV
        write(*,*) "sum |dspin| = ", sum( abs(rho(:,:,:,1)) )*dV
     end if

     if ( s == 2 ) exit
write(*,*) "h"
     read(u0) rtmp(:) ; write(*,*) "vh",sum(rtmp**2)
     do i=1,ML
        vh(LL(1,i),LL(2,i),LL(3,i))=rtmp(i)
     end do
write(*,*) "g"
     call planer_averaged_plot( s, vh(:,:,:), "vh" )
write(*,*) "hh"
     read(u0,END=10) rtmp
     do i=1,ML
        vxc(LL(1,i),LL(2,i),LL(3,i))=rtmp(i)
     end do
     call planer_averaged_plot( s, vxc(:,:,:), "vxc" )
write(*,*) "hhh"
     vh = vloc(:,:,:,s) - vxc
     call planer_averaged_plot( s, vh(:,:,:), "ves" ) ! es: electro static

  end do ! s
10 continue

  close(u0)

contains


  subroutine gen_cube( s, f, fn0 )
    implicit none
    integer,intent(IN) :: s
    real(8),intent(IN) :: f(0:,0:,0:)
    character(*),intent(IN) :: fn0
    integer :: i,m3,i1,i2,i3,j3
    character(30) :: name, file_name
    real(8) :: r0(3), aa_del(3,3)
    character(1) :: cs
    real(8),allocatable :: fsft(:,:,:)

    aa_del(:,1)=aa(:,1)/ML1
    aa_del(:,2)=aa(:,2)/ML2
    aa_del(:,3)=aa(:,3)/ML3

    name  = 'SYS1'
    r0(:) = 0.0d0

    write(cs,'(i1.1)') s

    file_name = fn0//"_"//cs//".cub"

!    allocate( fsft(0:ML1-1,0:ML2-1,0:ML3-1) ) ; fsft=0.0d0
!    m3=ML3/2
!    do i3=0,ML3-1
!       j3=mod(i3+m3,ML3)
!       do i2=0,ML2-1
!       do i1=0,ML1-1
!          fsft(i1,i2,j3) = f(i1,i2,i3)
!       end do
!       end do
!    end do

    open(iu,file=file_name)
    write(iu,100) name
    write(iu,100) name
    write(iu,110) MI,r0
    write(iu,110) ML1,aa_del(:,1)
    write(iu,110) ML2,aa_del(:,2)
    write(iu,110) ML3,aa_del(:,3)
    do i=1,MI
       write(iu,110) zatom(Kion(i)),real(zatom(Kion(i)),8),rsi(:,i)
    end do
    do i1=0,ML1-1
       do i2=0,ML2-1
!          write(iu,*) fsft(i1,i2,:)
          write(iu,*) f(i1,i2,:)
       end do
    end do
    close(iu)

100 format(a4)
110 format(i5,4f12.6)

!    deallocate( fsft )

  end subroutine gen_cube


  subroutine planer_averaged_plot( s, f, fn0 )
    implicit none
    integer,intent(IN) :: s
    real(8),intent(IN) :: f(0:,0:,0:)
    character(*),intent(IN) :: fn0
    integer :: i,m
    character(30) :: file_name
    real(8) :: Haa(3),c
    character(1) :: cs
    real(8),allocatable :: fav(:,:)

    Haa(1)=sqrt( sum(aa(:,1)**2) )/ML1
    Haa(2)=sqrt( sum(aa(:,2)**2) )/ML2
    Haa(3)=sqrt( sum(aa(:,3)**2) )/ML3

    write(cs,'(i1.1)') s

    file_name = fn0//"_"//cs//".gp"

    m = max( ML1, ML2, ML3 )
    allocate( fav(0:m-1,3) ); fav=0.0d0

    c=1.0d0/(ML2*ML3)
    do i=0,ML1-1
       fav(i,1) = sum( f(i,:,:) )*c
    end do
    c=1.0d0/(ML3*ML1)
    do i=0,ML2-1
       fav(i,2) = sum( f(:,i,:) )*c
    end do
    c=1.0d0/(ML1*ML2)
    do i=0,ML3-1
       fav(i,3) = sum( f(:,:,i) )*c
    end do

    open(iu,file=file_name//"x")
    do i=0,ML1-1
       write(iu,'(1x,4f20.10)') Haa(1)*i, fav(i,1)
    end do
    close(iu)
    open(iu,file=file_name//"y")
    do i=0,ML2-1
       write(iu,'(1x,4f20.10)') Haa(2)*i, fav(i,2)
    end do
    close(iu)
    open(iu,file=file_name//"z")
    do i=0,ML3-1
       write(iu,'(1x,4f20.10)') Haa(3)*i, fav(i,3)
    end do
    close(iu)

    deallocate( fav )

  end subroutine planer_averaged_plot


end program plot_vrho
