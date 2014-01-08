MODULE esm_rgrid_module

  use rgrid_module
  use aa_module
  use parallel_module

  implicit none

  PRIVATE
  PUBLIC :: Rsize1,Rsize2,LL,KK,read_esm_rgrid &
       ,prep_esm_rgrid,construct0_esm_rgrid,construct_esm_rgrid &
       ,ML_ESM,MK_ESM,ML0_ESM,ML1_ESM,MK0_ESM,MK1_ESM,LL_ESM &
       ,Nshift_ESM

  real(8) :: Rsize1,Rsize2
  integer,allocatable :: LL(:,:),KK(:,:),LL_ESM(:,:)
  integer :: ML_ESM,MK_ESM,ML0_ESM,ML1_ESM,MK0_ESM,MK1_ESM
  integer :: Nshift_ESM(3)
  integer :: Md

  real(8),parameter :: ep=1.d-10
  integer,allocatable :: ichk(:,:)

CONTAINS


  SUBROUTINE read_esm_rgrid(rank,unit)
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
    call send_esm_rgrid(rank)
  END SUBROUTINE read_esm_rgrid


  SUBROUTINE send_esm_rgrid(rank)
    implicit none
    integer,intent(IN)  :: rank
    integer :: ierr
    call mpi_bcast(Rsize1,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Rsize2,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_esm_rgrid


  SUBROUTINE prep_esm_rgrid(Md_in)
    implicit none
    integer,intent(IN) :: Md_in
    Md=Md_in
    Nshift_ESM(1) = Ngrid(1)/2
    Nshift_ESM(2) = Ngrid(2)/2
    Nshift_ESM(3) = Ngrid(3)/2
    call flush(6)
    write(*,'(1x,"before",2i12,6i6)') Igrid
    Igrid(1,1:3) = Igrid(1,1:3) - Nshift_ESM(1:3)
    Igrid(2,1:3) = Igrid(2,1:3) - Nshift_ESM(1:3)
    write(*,'(1x,"after ",2i12,6i6)') Igrid
    call flush(6)
  END SUBROUTINE prep_esm_rgrid


  SUBROUTINE construct0_esm_rgrid
    implicit none
    integer :: i1,i2,i3,m,n,is(2),ir(2),ierr,j1,j2
    real(8) :: x,y,z,rrxy,c1,c2,c3
    integer,allocatable :: ichk1(:,:,:),ichk2(:,:,:),irc(:,:)

    c1=1.d0/Ngrid(1)
    c2=1.d0/Ngrid(2)
    c3=1.d0/Ngrid(3)

    i1=Ngrid(1)/2 ; i2=Ngrid(2)/2
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
    end do
    end do
!    do i1=-Ngrid(1)/2-Md,Ngrid(1)/2-1+Md
!       write(*,'(1x,60i1)') (ichk(i1,i2),i2=-Ngrid(1)/2-Md,Ngrid(1)/2-1+Md)
!    end do
!    stop

    ML_ESM=count(ichk==1)*Ngrid(3)
    MK_ESM=count(ichk==2)*Ngrid(3)

    call flush(6)
    write(*,*) "ML_ESM ( < Rsize1 )",ML_ESM
    write(*,*) "MK_ESM             ",MK_ESM
    call flush(6)

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

    allocate( irc(2,0:np_grid-1) ) ; irc=0
    irc(1,myrank_g) = m
    irc(2,myrank_g) = n
    call mpi_allgather(irc(1,myrank_g),2,mpi_integer,irc,2,mpi_integer,comm_grid,ierr)
    ML0_ESM = sum( irc(1,0:myrank_g) ) - irc(1,myrank_g) + 1
    ML1_ESM = sum( irc(1,0:myrank_g) )
    MK0_ESM = sum( irc(2,0:myrank_g) ) - irc(2,myrank_g) + 1
    MK1_ESM = sum( irc(2,0:myrank_g) )
    deallocate( irc )

    allocate( LL_ESM(3,ML0_ESM:ML1_ESM) ) ; LL_ESM=0
    allocate( LL(3,ML0_ESM:ML1_ESM) ) ; LL=0
    allocate( KK(3,MK0_ESM:MK1_ESM) ) ; KK=0
    m=ML1_ESM-ML0_ESM+1
    n=MK1_ESM-MK0_ESM+1
    call construct_esm_rgrid(m,LL,n,KK)
    LL_ESM(:,:)=LL(:,:)

  END SUBROUTINE construct0_esm_rgrid

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
