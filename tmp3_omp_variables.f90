MODULE omp_variables

!$ use omp_lib

  implicit none

  PRIVATE
  PUBLIC :: init_omp, Igrid_omp,Ngrid_omp

  integer,allocatable :: Igrid_omp(:,:,:)
  integer,allocatable :: Ngrid_omp(:,:)
  integer :: nthreads

CONTAINS


  SUBROUTINE init_omp(a1b,b1b,a2b,b2b,a3b,b3b,n1,n2,disp_switch)

    implicit none

    integer,intent(IN) :: a1b,b1b,a2b,b2b,a3b,b3b,n1,n2
    logical,intent(IN) :: disp_switch
    integer :: i1,i2,i3,m,n,k,j,i,ab1,ab2,ab3,nn(3),np(3)
    integer,allocatable :: ntmp(:,:)

    if ( disp_switch ) write(*,'(a60," init_omp")') repeat("-",60)

    ab1 = b1b-a1b+1
    ab2 = b2b-a2b+1
    ab3 = b3b-a3b+1

    m=1
    n=1

!$OMP parallel private(m,n,k,j,i)

!$  n=omp_get_num_threads()

!$OMP single

    nthreads = n
    allocate( Ngrid_omp(0:3,0:n-1)   ) ; Ngrid_omp=0
    allocate( Igrid_omp(2,0:3,0:n-1) ) ; Igrid_omp=0

    m = n

    k     = gcd(ab3,m)
    m     = m/k
    nn(3) = ab3/k

    if ( m == 1 ) then

       nn(2) = ab2
       nn(1) = ab1

    else

       k     = gcd(ab2,m)
       m     = m/k
       nn(2) = ab2/k

       if ( m == 1 ) then

          nn(1) = ab1

       else

          k     = gcd(ab1,m)
          m     = m/k
          nn(1) = ab1/k

          if ( m == 1 ) then
          else
             nn(1) = ab1
          end if

       end if

    end if

    if ( disp_switch ) then
       write(*,*) "nn(1:3)=",nn(1:3)
       write(*,*) "m      =",m
    end if

    if ( m == 1 ) then

       np(1)=ab1/nn(1)
       np(2)=ab2/nn(2)
       np(3)=ab3/nn(3)

       if ( disp_switch ) write(*,*) "np=",np

       n=maxval(np)
       allocate( ntmp(3,n) ) ; ntmp=0

       do i=1,ab1
          n=mod(i-1,np(1))+1
          ntmp(1,n)=ntmp(1,n)+1
       end do
       do i=1,ab2
          n=mod(i-1,np(2))+1
          ntmp(2,n)=ntmp(2,n)+1
       end do
       do i=1,ab3
          n=mod(i-1,np(3))+1
          ntmp(3,n)=ntmp(3,n)+1
       end do

       n=-1
       do k=1,np(3)
       do j=1,np(2)
       do i=1,np(1)
          n=n+1
          Ngrid_omp(1,n)=ntmp(1,i)
          Ngrid_omp(2,n)=ntmp(2,j)
          Ngrid_omp(3,n)=ntmp(3,k)
          if ( disp_switch ) write(*,'(1x,"Ngrid_omp",3i5)') Ngrid_omp(1:3,n)
       end do
       end do
       end do

    else

       np(1)=m
       np(2)=ab2/nn(2)
       np(3)=ab3/nn(3)

       if ( disp_switch ) write(*,*) "np_=",np

       n=maxval(np)
       allocate( ntmp(3,n) ) ; ntmp=0

       do i=1,ab1
          n=mod(i-1,np(1))+1
          ntmp(1,n)=ntmp(1,n)+1
       end do
       do i=1,ab2
          n=mod(i-1,np(2))+1
          ntmp(2,n)=ntmp(2,n)+1
       end do
       do i=1,ab3
          n=mod(i-1,np(3))+1
          ntmp(3,n)=ntmp(3,n)+1
       end do

       n=-1
       do k=1,np(3)
       do j=1,np(2)
       do i=1,np(1)
          n=n+1
          Ngrid_omp(1,n)=ntmp(1,i)
          Ngrid_omp(2,n)=ntmp(2,j)
          Ngrid_omp(3,n)=ntmp(3,k)
          if ( disp_switch ) write(*,'(1x,"Ngrid_omp_",3i5)') Ngrid_omp(1:3,n)
       end do
       end do
       end do

    end if

    do i=0,nthreads-1
       Ngrid_omp(0,i) = Ngrid_omp(1,i)*Ngrid_omp(2,i)*Ngrid_omp(3,i)
    end do

    i=-1
    do i3=1,np(3)
    do i2=1,np(2)
    do i1=1,np(1)
       i=i+1
       Igrid_omp(1,0,i) = sum( Ngrid_omp(0,0:i) ) - Ngrid_omp(0,i) + n1
       Igrid_omp(2,0,i) = Igrid_omp(1,0,i) + Ngrid_omp(0,i) - 1
       Igrid_omp(1,1,i) = sum(ntmp(1,1:i1))-ntmp(1,i1) + a1b
       Igrid_omp(2,1,i) = Igrid_omp(1,1,i) + ntmp(1,i1) - 1
       Igrid_omp(1,2,i) = sum(ntmp(2,1:i2))-ntmp(2,i2) + a2b
       Igrid_omp(2,2,i) = Igrid_omp(1,2,i) + ntmp(2,i2) - 1
       Igrid_omp(1,3,i) = sum(ntmp(3,1:i3))-ntmp(3,i3) + a3b
       Igrid_omp(2,3,i) = Igrid_omp(1,3,i) + ntmp(3,i3) - 1
       if ( disp_switch ) then
          write(*,'(1x,"IGIG",1x,9i8)') Igrid_omp(:,:,i),Ngrid_omp(0,i)
       end if
    end do
    end do
    end do

    if ( disp_switch ) write(*,*) "omp_get_num_threads=",nthreads

    deallocate( ntmp )

!$OMP end single

!$OMP end parallel

  END SUBROUTINE init_omp


  FUNCTION gcd(m0,n0)
    implicit none
    integer :: gcd,m0,n0
    integer :: m,n,mtmp,loop

    if ( m0 >= n0 ) then
       m=m0
       n=n0
    else
       m=n0
       n=m0
    end if

    do loop=1,10000
       if ( n == 0 ) exit
       mtmp = n
       n = mod(m,n)
       m = mtmp
    end do

    gcd = m

  END FUNCTION gcd


END MODULE omp_variables
