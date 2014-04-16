PROGRAM band_plot

  implicit none

  integer,parameter :: u1=100,max_loop=100000
  integer,allocatable :: mbi(:)
  integer :: idummy,idummy2,mb_min,mb_max,i,j,loop
  integer :: nbk,mb,mbv,mbc,mb1,mb2
  real(8),parameter :: HT=27.2116d0
  real(8),allocatable :: eval(:,:),dk_bz(:,:)
  real(8) :: evb,ecb,esft,dummy,rdummy(3),bb(3,3)
  real(8) :: kx,ky,kz,d1,d2,d3

  open(u1,file="band_eigv",status='old')
  read(u1,*) bb(1:3,1)
  read(u1,*) bb(1:3,2)
  read(u1,*) bb(1:3,3)
  read(u1,*) idummy,mb,rdummy(1:3),mbv
  write(*,*) "mb,mbv  =",mb,mbv
  do loop=1,max_loop
     do i=1,mb
        read(u1,*)
     end do
     read(u1,*,END=10)
  end do
  stop "stop!!! too much data"
10 continue
  nbk=loop

  mb_max=mb
  mb_min=mb

  write(*,*) "mb_max =",mb_max
  write(*,*) "nbk    =",nbk

  allocate( mbi(nbk)         ) ; mbi=0
  allocate( dk_bz(4,nbk)     ) ; dk_bz=0.d0
  allocate( eval(mb_max,nbk) ) ; eval=0.d0

  mbi(:)=mb

  rewind u1
  read(u1,*)
  read(u1,*)
  read(u1,*)
  do i=1,nbk
     read(u1,*) idummy,idummy2,dk_bz(1:3,i),dummy
     do j=1,mbi(i)
        read(u1,*) idummy,eval(j,i),dummy
     end do
     if ( i == 1 ) cycle
     d1=dk_bz(1,i)-dk_bz(1,i-1)
     d2=dk_bz(2,i)-dk_bz(2,i-1)
     d3=dk_bz(3,i)-dk_bz(3,i-1)
     kx=bb(1,1)*d1+bb(1,2)*d2+bb(1,3)*d3
     ky=bb(2,1)*d1+bb(2,2)*d2+bb(2,3)*d3
     kz=bb(3,1)*d1+bb(3,2)*d2+bb(3,3)*d3
     dk_bz(4,i)=dk_bz(4,i-1)+sqrt(kx*kx+ky*ky+kz*kz)
  end do

  close(u1)

!------------------------------

  mbc = mbv + 1
  mb1 = max(1,mbv-200)
  mb2 = min(mbv+1000,mb)
  evb = maxval( eval(1:mbv,1:nbk) )
  ecb = minval( eval(mbc:mb,1:nbk) )

  write(*,*) "eg(HT,eV)=",ecb-evb,(ecb-evb)*HT

!------------------------------

  write(*,*) "evb=",evb
  write(*,*) "ecb=",ecb
  write(*,*) "Which energy is used as origin ? [ 0:0.0, 1:evb, 2:ecb ]"
  read(*,*) i
  select case(i)
  case default
  case(1)
     esft=evb
  case(2)
     esft=ecb
  end select

  call plot_eval(mb1,mb2,esft,HT)

!------------------------------

CONTAINS

  SUBROUTINE plot_eval(mb1,mb2,e0,ht)

    implicit none
    integer,intent(IN) :: mb1,mb2
    real(8),intent(IN) :: e0,ht
    integer :: i,i0,i1,k,unit0,unit

    unit0=2000
    unit=unit0
    do i0=mb1,mb2,50
       i1=min(i0+50-1,mb2)
       unit=unit+1
       rewind unit
       do k=1,nbk
          write(unit,'(1x,f10.5,1x,50f15.8)') dk_bz(4,k),( (eval(i,k)-e0)*ht,i=i0,i1 )
       end do
    end do

    open(unit0,file='plot_band',status='new')
    unit=unit0
    write(unit0,'(1x,"set nokey")')
    do i0=mb1,mb2,50
       i1=min(i0+50-1,mb2)
       unit=unit+1
       do i=i0,i1
          if ( i == mb1 ) then
             write(unit0,'(1x,a11,i4,a6,i3,a7)') 'plot "fort.',unit,'" u 1:',i-i0+2,' w lp \'
          else
             write(unit0,'(1x,4x,a7,i4,a6,i3,a7)') ',"fort.',unit,'" u 1:',i-i0+2,' w lp \'
          end if
       end do
    end do
    close(unit0)

  END SUBROUTINE plot_eval

END PROGRAM band_plot
