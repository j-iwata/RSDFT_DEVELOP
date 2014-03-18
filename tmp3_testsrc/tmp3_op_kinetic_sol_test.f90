MODULE op_kinetic_sol_test_module

  use rgrid_module
  use omp_variables
  use bc_variables
  use bc_test_module
!  use bc_module

  implicit none

  PRIVATE
  PUBLIC :: op_kinetic_sol_test1, init_op_kinetic_sol_test &
       ,flag_first_time

!  integer :: Md
  real(8) :: coef_lap0,ggg(6)
  real(8),allocatable :: coef_lap(:,:),coef_nab(:,:)
  real(8),allocatable :: coef_nabk(:,:,:),const_k2(:)
  complex(8),allocatable :: zcoef_kin(:,:,:)
  logical :: flag_nab,flag_n12,flag_n23,flag_n31
  logical :: flag_first_time=.true.

CONTAINS


  SUBROUTINE init_op_kinetic_sol_test(md_in,mk_in,coef_lap_in,coef_nab_in &
       ,coef_lap0_in,coef_nabk_in,const_k2_in,zcoef_kin_in,ggg_in &
       ,flag_nab_in,flag_n12_in,flag_n23_in,flag_n31_in)
    implicit none
    integer,intent(IN) :: md_in,mk_in
    real(8),intent(IN) :: coef_lap_in(3,md_in),coef_nab_in(3,md_in)
    real(8),intent(IN) :: coef_lap0_in,coef_nabk_in(3,md_in,mk_in)
    real(8),intent(IN) :: const_k2_in(0:mk_in),ggg_in(6)
    complex(8),intent(IN) :: zcoef_kin_in(3,-md_in:md_in,mk_in)
    logical,intent(IN) :: flag_nab_in,flag_n12_in,flag_n23_in,flag_n31_in

    Md = md_in

    ggg(:) = ggg_in(:)

    coef_lap0 = coef_lap0_in

    allocate( coef_lap(3,Md) ) ; coef_lap=0.0d0
    allocate( coef_nab(3,Md) ) ; coef_nab=0.0d0
    coef_lap(:,:)=coef_lap_in(:,:)
    coef_nab(:,:)=coef_nab_in(:,:)

    allocate( coef_nabk(3,Md,mk_in)     ) ; coef_nabk=0.0d0
    allocate( const_k2(0:mk_in)         ) ; const_k2 =0.0d0
    allocate( zcoef_kin(3,-Md:Md,mk_in) ) ; zcoef_kin=(0.0d0,0.0d0)
    coef_nabk(:,:,:)=coef_nabk_in(:,:,:)
    const_k2(:)=const_k2_in(:)
    zcoef_kin(:,:,:)=zcoef_kin_in(:,:,:)

    flag_nab = flag_nab_in
    flag_n12 = flag_n12_in
    flag_n23 = flag_n23_in
    flag_n31 = flag_n31_in

    flag_first_time = .false.

  END SUBROUTINE init_op_kinetic_sol_test


  SUBROUTINE op_kinetic_sol_test1(k,tpsi,htpsi,n1,n2,ib1,ib2)
!$  use omp_lib
    implicit none
    integer,intent(IN) :: k,n1,n2,ib1,ib2
#ifdef _DRSDFT_
    real(8),intent(IN)  :: tpsi(n1:n2,ib1:ib2)
    real(8),intent(INOUT) :: htpsi(n1:n2,ib1:ib2)
    real(8),allocatable :: wk(:,:,:,:)
    real(8),parameter :: zero=0.d0
#else
    complex(8),intent(IN)  :: tpsi(n1:n2,ib1:ib2)
    complex(8),intent(INOUT) :: htpsi(n1:n2,ib1:ib2)
    complex(8),allocatable :: wk(:,:,:,:)
    complex(8),parameter :: zero=(0.d0,0.d0)
#endif
    integer :: i,ib,i1,i2,i3,nb,m,n,j
    integer :: a1,a2,a3,b1,b2,b3,p,mm,nn
    integer :: a1b,b1b,a2b,b2b,a3b,b3b,ab1,ab12
    real(8) :: c,d
    integer,allocatable :: ic(:)
    integer :: a1b_omp,b1b_omp,a2b_omp,b2b_omp,a3b_omp,b3b_omp,n1_omp,n2_omp
    integer :: ib1_omp,ib2_omp,nb_omp

    a1b = Igrid(1,1)
    b1b = Igrid(2,1)
    a2b = Igrid(1,2)
    b2b = Igrid(2,2)
    a3b = Igrid(1,3)
    b3b = Igrid(2,3)
    ab1 = (b1b-a1b+1)
    ab12= (b1b-a1b+1)*(b2b-a2b+1)

    nb = ib2-ib1+1

!$OMP parallel private(a3b_omp,b3b_omp,a2b_omp,b2b_omp,a1b_omp,b1b_omp &
!$OMP                 ,n1_omp,n2_omp,j,p,d,mm,c)

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

    do ib=ib1,ib2
!       j=n1_omp-1
       do i3=a3b_omp,b3b_omp
       do i2=a2b_omp,b2b_omp
       do i1=a1b_omp,b1b_omp
!          j=j+1
          j=n1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
          www(i1,i2,i3,ib-ib1+1) = tpsi(j,ib)
       end do
       end do
       end do
    end do

!$OMP barrier
    call bcset_test1(1,nb,Md,0)
!    call bcset(1,nb,Md,0)
!$OMP barrier

    if ( flag_nab ) then

       do ib=ib1,ib2

          p = ib-ib1+1

          do m=1,Md

             !j=n1_omp-1
             do i3=a3b_omp,b3b_omp
             do i2=a2b_omp,b2b_omp
             do i1=a1b_omp,b1b_omp
                !j=j+1
                j=n1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
                htpsi(j,ib)=htpsi(j,ib) &
                     +zcoef_kin(1,m,k) *www(i1+m,i2,i3,p) &
               +conjg(zcoef_kin(1,m,k))*www(i1-m,i2,i3,p) &
                     +zcoef_kin(2,m,k) *www(i1,i2+m,i3,p) &
               +conjg(zcoef_kin(2,m,k))*www(i1,i2-m,i3,p) &
                     +zcoef_kin(3,m,k) *www(i1,i2,i3+m,p) &
               +conjg(zcoef_kin(3,m,k))*www(i1,i2,i3-m,p)   
             end do
             end do
             end do

          end do ! m

       end do ! ib

    else

       do ib=ib1,ib2

          p=ib-ib1+1

          do m=1,Md

             !j=n1_omp-1
             do i3=a3b_omp,b3b_omp
             do i2=a2b_omp,b2b_omp
             do i1=a1b_omp,b1b_omp
                !j=j+1
                j=n1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
                htpsi(j,ib)=htpsi(j,ib) &
                  +coef_lap(1,m)*( www(i1+m,i2,i3,p)+www(i1-m,i2,i3,p) ) &
                  +coef_lap(2,m)*( www(i1,i2+m,i3,p)+www(i1,i2-m,i3,p) ) &
                  +coef_lap(3,m)*( www(i1,i2,i3+m,p)+www(i1,i2,i3-m,p) )
             end do
             end do
             end do

          end do ! m

       end do ! ib

    end if

!$OMP barrier

    if ( flag_n12 .or. flag_n23 .or. flag_n31 ) then

!$OMP single
       a1=a1b-Md ; b1=b1b+Md
       a2=a2b-Md ; b2=b2b+Md
       a3=a3b-Md ; b3=b3b+Md
       allocate( wk(a1:b1,a2:b2,a3:b3,nb) )
       wk=www
!$OMP end single

!       do n=1,nb
!          do i3=a3b_omp,b3b_omp
!          do i2=a2b_omp,b2b_omp
!          do i1=a1b_omp,b1b_omp
!             wk(i1,i2,i3,n)=www(i1,i2,i3,n)
!          end do
!          end do
!          end do
!       end do

       if ( flag_n12 ) then

          do n=1,nb
             !j=n1_omp-1
             do i3=a3b_omp,b3b_omp
             do i2=a2b_omp,b2b_omp
             do i1=a1b_omp,b1b_omp
                www(i1,i2,i3,n)=zero
             end do
             end do
             end do
          end do
          do n=1,nb
             do m=1,Md
                d=coef_nab(1,m)
                do i3=a3b_omp,b3b_omp
                do i2=a2b_omp,b2b_omp
                do i1=a1b_omp,b1b_omp
                   www(i1,i2,i3,n)=www(i1,i2,i3,n) &
                        -d*( wk(i1-m,i2,i3,n)-wk(i1+m,i2,i3,n) )
                end do
                end do
                end do
             end do
          end do

!$OMP barrier
          call bcset_test1(1,nb,Md,3)
!          call bcset(1,nb,Md,3)
!$OMP barrier

          do ib=ib1,ib2
             p=ib-ib1+1
             do m=1,Md
                d=-ggg(4)*coef_nab(2,m)
                !j=n1_omp-1
                do i3=a3b_omp,b3b_omp
                do i2=a2b_omp,b2b_omp
                do i1=a1b_omp,b1b_omp
                   !j=j+1
                   j=n1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
                   htpsi(j,ib)=htpsi(j,ib) &
                        -d*(www(i1,i2-m,i3,p)-www(i1,i2+m,i3,p))
                end do
                end do
                end do
             end do
          end do

       end if

       if ( flag_n23 ) then

          do n=1,nb
             do i3=a3b_omp,b3b_omp
             do i2=a2b_omp,b2b_omp
             do i1=a1b_omp,b1b_omp
                www(i1,i2,i3,n)=zero
             end do
             end do
             end do
          end do

          do n=1,nb
             do m=1,Md
                d=coef_nab(2,m)
                do i3=a3b_omp,b3b_omp
                do i2=a2b_omp,b2b_omp
                do i1=a1b_omp,b1b_omp
                   www(i1,i2,i3,n)=www(i1,i2,i3,n) &
                        -d*( wk(i1,i2-m,i3,n)-wk(i1,i2+m,i3,n) )
                end do
                end do
                end do
             end do
          end do

!$OMP barrier
          call bcset_test1(1,nb,Md,5)
!          call bcset(1,nb,Md,5)
!$OMP barrier

          do ib=ib1,ib2
             p=ib-ib1+1
             do m=1,Md
                d=-ggg(5)*coef_nab(3,m)
                !j=n1_omp-1
                do i3=a3b_omp,b3b_omp
                do i2=a2b_omp,b2b_omp
                do i1=a1b_omp,b1b_omp
                   !j=j+1
                   j=n1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
                   htpsi(j,ib)=htpsi(j,ib) &
                        -d*(www(i1,i2,i3-m,p)-www(i1,i2,i3+m,p))
                end do
                end do
                end do
             end do
          end do

       end if

       if ( flag_n31 ) then

          do n=1,nb
             do i3=a3b_omp,b3b_omp
             do i2=a2b_omp,b2b_omp
             do i1=a1b_omp,b1b_omp
                www(i1,i2,i3,n)=zero
             end do
             end do
             end do
          end do

          do n=1,nb
             do m=1,Md
                d=coef_nab(3,m)
                do i3=a3b_omp,b3b_omp
                do i2=a2b_omp,b2b_omp
                do i1=a1b_omp,b1b_omp
                   www(i1,i2,i3,n)=www(i1,i2,i3,n) &
                        -d*( wk(i1,i2,i3-m,n)-wk(i1,i2,i3+m,n) )
                end do
                end do
                end do
             end do
          end do

!$OMP barrier
          call bcset_test1(1,nb,Md,1)
!          call bcset(1,nb,Md,1)
!$OMP barrier

          do ib=ib1,ib2
             p=ib-ib1+1
             do m=1,Md
                d=-ggg(6)*coef_nab(1,m)
                !j=n1_omp-1
                do i3=a3b_omp,b3b_omp
                do i2=a2b_omp,b2b_omp
                do i1=a1b_omp,b1b_omp
                   !j=j+1
                   j=n1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
                   htpsi(j,ib)=htpsi(j,ib) &
                        -d*(www(i1-m,i2,i3,p)-www(i1+m,i2,i3,p))
                end do
                end do
                end do
             end do
          end do

       end if

!$OMP single
       deallocate( wk )
!$OMP end single

    end if
!$OMP end parallel
 
  END SUBROUTINE op_kinetic_sol_test1

END MODULE op_kinetic_sol_test_module
