PROGRAM cubegen_simple

  implicit none

  integer,parameter :: u1=5, u2=970, iu=10
  real(8),parameter :: ab=0.529177d0
  integer :: idum,MB,MBZ,i,j,k,n,ML,ML1,ML2,ML3,MI,MKI
  integer :: NBLK1,NBLK2,NBLK,np1,np2,np3,i1,i2,i3,n1,n2,n3
  integer :: MBwr1,MBwr2,MBZ0,MB0
  integer :: zatom(5)
  real(8) :: aa(3,3),ax,r0(3),aL(3),aa_del(3,3),va,fac
  real(8) :: c,d,a,b
  integer,allocatable :: Kion(:),mkd(:)
  real(8),allocatable :: asi(:,:),rsi(:,:),rho(:,:,:),occ(:,:)
  real(8),allocatable :: rtmp(:)
  complex(8) ,parameter :: zi=(0.d0,1.d0)
  complex(8),allocatable :: unk(:,:,:,:,:),utmp(:)
  character(4) :: name,cbuf
  integer,allocatable :: ip1(:),ip2(:),ip3(:)
  integer,allocatable :: nip1(:),nip2(:),nip3(:)
  integer,allocatable :: LL2(:,:),LLL2(:,:,:)

! --- read input data ---

  read(u1,*)
  read(u1,*)
  read(u1,*) cbuf,ax
  do i=1,3
     read(u1,*) cbuf,aa(1:3,i)
  end do
  read(u2,*) MKI,MI,zatom(1:MKI)

  if ( MKI > 5 ) then
     write(*,*) "MKI=",MKI,size(zatom)
     stop
  end if

  allocate( asi(3,MI),Kion(MI),rsi(3,MI) )

  do i=1,MI
     read(u2,*) Kion(i),asi(:,i)
  end do

  aa=ax*aa
  do i=1,3
     aL(i)=sqrt( sum(aa(:,i)**2) )
  end do

  rsi=matmul( aa,asi )

  ML1=1.d0
  ML2=1.d0
  ML3=1.d0

  do i=1,3
     aa_del(1,i)=aa(1,i)/ML1
     aa_del(2,i)=aa(2,i)/ML2
     aa_del(3,i)=aa(3,i)/ML3
  end do

! --- cube file

  name='SYS1'
  r0(:)=0.d0

  open(iu,file='tmp.cube')
  write(iu,100) name
  write(iu,100) name
  write(iu,110) MI,r0
  write(iu,110) ML1,aa_del(:,1)
  write(iu,110) ML2,aa_del(:,2)
  write(iu,110) ML3,aa_del(:,3)
  do i=1,MI
     write(iu,110) zatom(Kion(i)),real(zatom(Kion(i)),8),rsi(:,i)
  end do
  close(iu)

100 format(a4)
110 format(i5,4f12.6)

  deallocate( rsi,Kion,asi )

END PROGRAM cubegen_simple
