PROGRAM cubegen_vrho
  implicit none

  integer,parameter :: u1=1,u2=970,iu=2
  integer :: ML,ML1,ML2,ML3,zatom(9)
  integer :: i1,i2,i3,i,s,MI,MKI
  integer,allocatable :: LL(:,:),Kion(:)
  real(8),allocatable :: rho(:,:,:,:),vloc(:,:,:,:)
  real(8),allocatable :: vh(:,:,:),vxc(:,:,:),rtmp(:)
  real(8),allocatable :: asi(:,:),rsi(:,:)
  real(8) :: ax,aa(3,3)
  character(6) :: cbuf

! ---

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

  rsi(:,:)=matmul( aa,asi )

! ---

  open(u1,file='vrho.dat1',status='old',form='unformatted')

  read(u1) ML,ML1,ML2,ML3

  allocate( LL(3,ML) ) ; LL=0.0d0
  allocate( rtmp(ML) ) ; rtmp=0.0d0
  allocate( rho(0:ML1-1,0:ML2-1,0:ML3-1,2)  ) ; rho=0.0d0
  allocate( vloc(0:ML1-1,0:ML2-1,0:ML3-1,2) ) ; vloc=0.0d0
  allocate( vh(0:ML1-1,0:ML2-1,0:ML3-1)     ) ; vh=0.0d0
  allocate( vxc(0:ML1-1,0:ML2-1,0:ML3-1)    ) ; vxc=0.0d0

  read(u1) LL

  do s=1,2

     read(u1,END=10) rtmp
     do i=1,ML
        rho(LL(1,i),LL(2,i),LL(3,i),s)=rtmp(i)
     end do

     call gen_cube( s, rho(:,:,:,s), "rho" )

     read(u1,END=10) rtmp
     do i=1,ML
        vloc(LL(1,i),LL(2,i),LL(3,i),s)=rtmp(i)
     end do

     call gen_cube( s, vloc(:,:,:,s), "vloc" )

     if ( s == 2 ) exit

     read(u1,END=10) rtmp
     do i=1,ML
        vh(LL(1,i),LL(2,i),LL(3,i))=rtmp(i)
     end do

     read(u1,END=10) rtmp
     do i=1,ML
        vxc(LL(1,i),LL(2,i),LL(3,i))=rtmp(i)
     end do

  end do ! s
10 continue

  close(u1)

CONTAINS


  SUBROUTINE gen_cube( s, f, fn0 )
    implicit none
    integer,intent(IN) :: s
    real(8),intent(IN) :: f(0:,0:,0:)
    character(*),intent(IN) :: fn0
    integer :: i
    character(30) :: name, file_name
    real(8) :: r0(3), aa_del(3,3)
    character(1) :: cs

    aa_del(:,1)=aa(:,1)/ML1
    aa_del(:,2)=aa(:,2)/ML2
    aa_del(:,3)=aa(:,3)/ML3

    name  = 'SYS1'
    r0(:) = 0.0d0

    write(cs,'(i1.1)') s

    file_name = fn0//"_"//cs//".cub"

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
          write(iu,*) f(i1,i2,:)
       end do
    end do
    close(iu)

100 format(a4)
110 format(i5,4f12.6)

  END SUBROUTINE gen_cube


END PROGRAM cubegen_vrho
