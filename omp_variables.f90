MODULE omp_variables

!$ use omp_lib

  implicit none

  PRIVATE
  PUBLIC :: a3b_omp,b3b_omp,n1_omp,n2_omp,init_omp

  integer,allocatable :: a3b_omp(:),b3b_omp(:),n1_omp(:),n2_omp(:)

CONTAINS


  SUBROUTINE init_omp(a1b,b1b,a2b,b2b,a3b,b3b,n1,n2,disp_switch)

    implicit none

    integer,intent(IN) :: a1b,b1b,a2b,b2b,a3b,b3b,n1,n2
    logical,intent(IN) :: disp_switch
    integer :: m,n,i
    integer,allocatable :: ic(:)

    m=0
    n=1

!$OMP parallel private(m)

!$  n=omp_get_num_threads()

!$OMP single
    allocate( a3b_omp(0:n-1) )
    allocate( b3b_omp(0:n-1) )
    allocate( n1_omp(0:n-1)  )
    allocate( n2_omp(0:n-1)  )
    allocate( ic(0:n-1) )
    ic(:)=(b3b-a3b+1)/n
    m=(b3b-a3b+1)-sum(ic)
    do i=0,m-1
       ic(i)=ic(i)+1
    end do
!$OMP end single

!$  m=omp_get_thread_num()

    a3b_omp(m) = a3b+sum(ic(0:m))-ic(m)
    b3b_omp(m) = a3b_omp(m)+ic(m)-1
    n1_omp(m)  = n1+(a3b_omp(m)-a3b)*(b2b-a2b+1)*(b1b-a1b+1)
    n2_omp(m)  = n1_omp(m)+(b3b_omp(m)-a3b_omp(m)+1)*(b2b-a2b+1)*(b1b-a1b+1)-1

    if ( disp_switch ) then
       write(*,'(1x,"omp_threads_num, omp_num_trhreads=",2i4,2i6,2i8)') &
            ,m,n,a3b_omp(m),b3b_omp(m),n1_omp(m),n2_omp(m)
    end if

!$OMP barrier

!$OMP single
    deallocate( ic )
!$OMP end single

!$OMP end parallel

  END SUBROUTINE init_omp

END MODULE omp_variables
