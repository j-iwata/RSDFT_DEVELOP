MODULE kinetic_sol_simple_module

!$  use omp_lib
  use rgrid_module
  use omp_variables
  use bc_module, only: www, bcset_1
  use kinetic_variables, only: coef_lap0, coef_lap, zcoef_kin, coef_nab &
       ,flag_nab, flag_n12, flag_n23, flag_n31, const_k2, ggg, Md
  use watch_module, only: watchb_omp, time_kine, time_bcfd

  implicit none

  PRIVATE
  PUBLIC :: op_kinetic_sol_simple

CONTAINS


  SUBROUTINE op_kinetic_sol_simple(k,tpsi,htpsi,n1,n2,ib1,ib2)

    implicit none
    integer,intent(IN) :: k,n1,n2,ib1,ib2
#ifdef _DRSDFT_
    real(8),intent(IN)    ::  tpsi(n1:n2,ib1:ib2)
    real(8),intent(INOUT) :: htpsi(n1:n2,ib1:ib2)
    real(8),allocatable :: wk(:,:,:,:)
    real(8),parameter :: zero=0.d0
#else
    complex(8),intent(IN)    ::  tpsi(n1:n2,ib1:ib2)
    complex(8),intent(INOUT) :: htpsi(n1:n2,ib1:ib2)
    complex(8),allocatable :: wk(:,:,:,:)
    complex(8),parameter :: zero=(0.d0,0.d0)
#endif
    integer :: i,ib,i1,i2,i3,nb,m,n,j
    integer :: a1,a2,a3,b1,b2,b3,p,mm,nn
    integer :: a1b,b1b,a2b,b2b,a3b,b3b,ab1,ab12
    real(8) :: c,d,et0,et1, ttmp(2)
    integer :: a1b_omp,b1b_omp,a2b_omp,b2b_omp,a3b_omp,b3b_omp
    integer :: n1_omp,n2_omp,ib1_omp,ib2_omp,nb_omp

!    call watchb_omp( ttmp )
!$OMP master
    time_bcfd(:,:)=0.0d0
!$OMP end master

    a1b = Igrid(1,1)
    b1b = Igrid(2,1)
    a2b = Igrid(1,2)
    b2b = Igrid(2,2)
    a3b = Igrid(1,3)
    b3b = Igrid(2,3)
    ab1 = (b1b-a1b+1)
    ab12= (b1b-a1b+1)*(b2b-a2b+1)

    nb = ib2-ib1+1

    mm=0
!$  mm=omp_get_thread_num()

    n1_omp = Igrid_omp(1,0,mm)
    n2_omp = Igrid_omp(2,0,mm)

    a1b_omp = Igrid_omp(1,1,mm)
    b1b_omp = Igrid_omp(2,1,mm)
    a2b_omp = Igrid_omp(1,2,mm)
    b2b_omp = Igrid_omp(2,2,mm)
    a3b_omp = Igrid_omp(1,3,mm)
    b3b_omp = Igrid_omp(2,3,mm)

    c=coef_lap0+const_k2(k)
    do ib=ib1,ib2
       do i=n1_omp,n2_omp
          htpsi(i,ib) = c*tpsi(i,ib)
       end do
    end do

!    call watchb_omp( ttmp, time_kine(1,1) )

    do ib=ib1,ib2
       do i3=a3b_omp,b3b_omp
       do i2=a2b_omp,b2b_omp
       do i1=a1b_omp,b1b_omp
          j=n1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
          www(i1,i2,i3,ib-ib1+1) = tpsi(j,ib)
       end do
       end do
       end do
    end do

!$OMP barrier

!    call watchb_omp( ttmp, time_kine(1,2) )

    call bcset_1(1,nb,Md,0)

!$OMP barrier

!    call watchb_omp( ttmp, time_kine(1,3) )

    do ib=ib1,ib2

       p=ib-ib1+1

       do i3=a3b_omp,b3b_omp
       do i2=a2b_omp,b2b_omp
       do i1=a1b_omp,b1b_omp
          j=n1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
          htpsi(j,ib)=htpsi(j,ib) &
               +coef_lap(1,6)*www(i1-6,i2,i3,p) &
               +coef_lap(1,5)*www(i1-5,i2,i3,p) &
               +coef_lap(1,4)*www(i1-4,i2,i3,p) &
               +coef_lap(1,3)*www(i1-3,i2,i3,p) &
               +coef_lap(1,2)*www(i1-2,i2,i3,p) &
               +coef_lap(1,1)*www(i1-1,i2,i3,p) &
               +coef_lap(1,1)*www(i1+1,i2,i3,p) &
               +coef_lap(1,2)*www(i1+2,i2,i3,p) &
               +coef_lap(1,3)*www(i1+3,i2,i3,p) &
               +coef_lap(1,4)*www(i1+4,i2,i3,p) &
               +coef_lap(1,5)*www(i1+5,i2,i3,p) &
               +coef_lap(1,6)*www(i1+6,i2,i3,p) &
               +coef_lap(2,6)*www(i1,i2-6,i3,p) &
               +coef_lap(2,5)*www(i1,i2-5,i3,p) &
               +coef_lap(2,4)*www(i1,i2-4,i3,p) &
               +coef_lap(2,3)*www(i1,i2-3,i3,p) &
               +coef_lap(2,2)*www(i1,i2-2,i3,p) &
               +coef_lap(2,1)*www(i1,i2-1,i3,p) &
               +coef_lap(2,1)*www(i1,i2+1,i3,p) &
               +coef_lap(2,2)*www(i1,i2+2,i3,p) &
               +coef_lap(2,3)*www(i1,i2+3,i3,p) &
               +coef_lap(2,4)*www(i1,i2+4,i3,p) &
               +coef_lap(2,5)*www(i1,i2+5,i3,p) &
               +coef_lap(2,6)*www(i1,i2+6,i3,p) &
               +coef_lap(3,6)*www(i1,i2,i3-6,p) &
               +coef_lap(3,5)*www(i1,i2,i3-5,p) &
               +coef_lap(3,4)*www(i1,i2,i3-4,p) &
               +coef_lap(3,3)*www(i1,i2,i3-3,p) &
               +coef_lap(3,2)*www(i1,i2,i3-2,p) &
               +coef_lap(3,1)*www(i1,i2,i3-1,p) &
               +coef_lap(3,1)*www(i1,i2,i3+1,p) &
               +coef_lap(3,2)*www(i1,i2,i3+2,p) &
               +coef_lap(3,3)*www(i1,i2,i3+3,p) &
               +coef_lap(3,4)*www(i1,i2,i3+4,p) &
               +coef_lap(3,5)*www(i1,i2,i3+5,p) &
               +coef_lap(3,6)*www(i1,i2,i3+6,p)
       end do
       end do
       end do

    end do ! ib

!$OMP barrier

!    call watchb_omp( ttmp, time_kine(1,4) )
 
  END SUBROUTINE op_kinetic_sol_simple


END MODULE kinetic_sol_simple_module
