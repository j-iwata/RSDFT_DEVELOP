MODULE esm_rgrid_module

  use parallel_module

  implicit none

  PRIVATE
  PUBLIC :: Rsize1,Rsize2,LL,KK &
       ,ML_ESM,MK_ESM,ML0_ESM,ML1_ESM,MK0_ESM,MK1_ESM,LL_ESM &
       ,Nshift_ESM &
       ,Read_RgridESM, Init_RgridESM, ParallelInit_RgridESM

  real(8) :: Rsize1,Rsize2
  integer,allocatable :: LL(:,:),KK(:,:),LL_ESM(:,:)
  integer :: ML_ESM,MK_ESM,ML0_ESM,ML1_ESM,MK0_ESM,MK1_ESM
  integer :: Nshift_ESM(3)

  real(8) :: aa(3,3)
  integer :: Ngrid(0:3)
  integer :: Md
  integer :: Igrid(2,0:3)

  real(8),parameter :: ep=1.d-10
  integer,allocatable :: ichk(:,:)

CONTAINS


  SUBROUTINE Read_RgridESM(rank,unit)
    implicit none
    integer,intent(IN)  :: rank,unit
    integer :: i
    character(4) :: cbuf,ckey
    Rsize1=0.0d0
    Rsize2=0.0d0
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:4) == "RESM" ) then
             backspace(unit)
             read(unit,*) cbuf,Rsize1,Rsize2
          end if
       end do
999    continue
       write(*,*) "Rsize1=",Rsize1
       write(*,*) "Rsize2=",Rsize2
    end if
    call Send_RgridESM(rank)
  END SUBROUTINE Read_RgridESM


  SUBROUTINE Send_RgridESM(rank)
    implicit none
    integer,intent(IN)  :: rank
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(Rsize1,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Rsize2,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  END SUBROUTINE Send_RgridESM


  SUBROUTINE Init_RgridESM( aa_in, Ngrid_inout, Md_in )
    implicit none
    real(8),intent(IN)    :: aa_in(3,3)
    integer,intent(INOUT) :: Ngrid_inout(0:3)
    integer,intent(IN)    :: Md_in
    integer :: i1,i2,i3,j1,j2
    real(8) :: c1,c2,c3,x,y,rrxy

    aa(:,:)  = aa_in(:,:)
    Ngrid(:) = Ngrid_inout(:)
    Md       = Md_in

    Nshift_ESM(1) = Ngrid(1)/2
    Nshift_ESM(2) = Ngrid(2)/2
    Nshift_ESM(3) = Ngrid(3)/2

    c1=1.d0/Ngrid(1)
    c2=1.d0/Ngrid(2)
    c3=1.d0/Ngrid(3)

    i1=Ngrid(1)/2
    i2=Ngrid(2)/2
    allocate( ichk(-i1-Md:i1-1+Md,-i2-Md:i2-1+Md) ) ; ichk=0

    do i2=-Ngrid(2)/2,Ngrid(2)/2-1
    do i1=-Ngrid(1)/2,Ngrid(1)/2-1

       x=aa(1,1)*c1*i1+aa(1,2)*c2*i2
       y=aa(2,1)*c1*i1+aa(2,2)*c2*i2
       rrxy=x*x+y*y

       if ( rrxy - ep < Rsize1**2 ) then

          ichk(i1,i2)=1
          do j1=i1-Md,i1+Md
             x=aa(1,1)*c1*j1+aa(1,2)*c2*i2
             y=aa(2,1)*c1*j1+aa(2,2)*c2*i2
             rrxy=x*x+y*y
             if ( Rsize1**2 <= rrxy - ep ) then
                ichk(j1,i2)=2
             end if
          end do
          do j2=i2-Md,i2+Md
             x=aa(1,1)*c1*i1+aa(1,2)*c2*i2
             y=aa(2,1)*c1*i1+aa(2,2)*c2*j2
             rrxy=x*x+y*y
             if ( Rsize1**2 <= rrxy - ep ) then
                ichk(i1,j2)=2
             end if
          end do

       end if

    end do ! i1
    end do ! i2

    ML_ESM=count(ichk==1)*Ngrid(3)
    MK_ESM=count(ichk==2)*Ngrid(3)

    if ( disp_switch_parallel ) then
       write(*,*) "ML_ESM ( < Rsize1 )",ML_ESM
       write(*,*) "MK_ESM             ",MK_ESM
    end if

    Ngrid(0) = ML_ESM
    Ngrid_inout(0) = ML_ESM

  END SUBROUTINE Init_RgridESM


  SUBROUTINE ParallelInit_RgridESM( np, pinfo_grid, Nshift )
    implicit none
    integer,intent(IN)    :: np
    integer,intent(INOUT) :: pinfo_grid(8,0:np-1) 
    integer,intent(OUT)   :: Nshift(3)
    integer :: i1,i2,i3,m,n,ierr
    integer,allocatable :: irc(:,:)

    Nshift(1:3) = Nshift_ESM(1:3)

    Igrid(1,1) = pinfo_grid(1,myrank_g) + 1
    Igrid(2,1) = pinfo_grid(1,myrank_g) + pinfo_grid(2,myrank_g)
    Igrid(1,2) = pinfo_grid(3,myrank_g) + 1
    Igrid(2,2) = pinfo_grid(3,myrank_g) + pinfo_grid(4,myrank_g)
    Igrid(1,3) = pinfo_grid(5,myrank_g) + 1
    Igrid(2,3) = pinfo_grid(5,myrank_g) + pinfo_grid(6,myrank_g)

    Igrid(1,1:3) = Igrid(1,1:3) - Nshift_ESM(1:3)
    Igrid(2,1:3) = Igrid(2,1:3) - Nshift_ESM(1:3)

    m=0
    do i3=Igrid(1,3),Igrid(2,3)
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
       if ( ichk(i1,i2) == 1 ) m=m+1
    end do
    end do
    end do
    n=0
    do i3=Igrid(1,3),Igrid(2,3)
    do i2=Igrid(1,2)-Md,Igrid(2,2)+Md
    do i1=Igrid(1,1)-Md,Igrid(2,1)+Md
       if ( ichk(i1,i2) == 2 ) n=n+1
    end do
    end do
    end do

!---

    allocate( irc(2,0:np_grid-1) ) ; irc=0

    irc(1,myrank_g) = m
    irc(2,myrank_g) = n
    call mpi_allgather(irc(1,myrank_g),2,mpi_integer,irc,2,mpi_integer,comm_grid,ierr)

    do n=0,np_grid-1
       pinfo_grid(7,n) = sum( irc(1,0:n) ) - irc(1,n)
       pinfo_grid(8,n) = irc(1,n)
    end do

    ML0_ESM = sum( irc(1,0:myrank_g) ) - irc(1,myrank_g) + 1
    ML1_ESM = sum( irc(1,0:myrank_g) )
    MK0_ESM = sum( irc(2,0:myrank_g) ) - irc(2,myrank_g) + 1
    MK1_ESM = sum( irc(2,0:myrank_g) )

    deallocate( irc )

!---

    allocate( LL_ESM(3,ML0_ESM:ML1_ESM) ) ; LL_ESM=0
    allocate( LL(3,ML0_ESM:ML1_ESM)     ) ; LL=0
    allocate( KK(3,MK0_ESM:MK1_ESM)     ) ; KK=0

    m = ML1_ESM - ML0_ESM + 1
    n = MK1_ESM - MK0_ESM + 1
    call construct_esm_rgrid(m,LL,n,KK)

    LL_ESM(:,:) = LL(:,:)

    deallocate( ichk )

  END SUBROUTINE ParallelInit_RgridESM

  SUBROUTINE construct_esm_rgrid(ml,map1,mk,map2)
    implicit none
    integer,intent(IN)  :: ml
    integer,intent(OUT) :: map1(3,ml)
    integer,intent(IN)  :: mk
    integer,intent(OUT) :: map2(3,mk)
    integer :: i1,i2,i3,m,n
    map1(:,:)=0
    map2(:,:)=0
    m=0
    do i3=Igrid(1,3),Igrid(2,3)
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
       if ( ichk(i1,i2) == 1 ) then
          m=m+1
          map1(1,m)=i1
          map1(2,m)=i2
          map1(3,m)=i3
       end if
    end do
    end do
    end do
    n=0
    do i3=Igrid(1,3),Igrid(2,3)
    do i2=Igrid(1,2)-Md,Igrid(2,2)+Md
    do i1=Igrid(1,1)-Md,Igrid(2,1)+Md
       if ( ichk(i1,i2) == 2 ) then
          n=n+1
          map2(1,n)=i1
          map2(2,n)=i2
          map2(3,n)=i3
       end if
    end do
    end do
    end do
    if ( m /= ml ) then
       write(*,*) "ml/=m! : ml,m=",ml,m
       stop "stop@construct_esm_rgrid"
    end if
    if ( n /= mk ) then
       write(*,*) "mk/=n! : mk,n=",mk,n
       stop "stop@construct_esm_rgrid"
    end if
  END SUBROUTINE construct_esm_rgrid

END MODULE esm_rgrid_module
