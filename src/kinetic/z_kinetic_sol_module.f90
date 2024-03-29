module z_kinetic_sol_module

!$  use omp_lib
  use rgrid_module
  use omp_variables
  use bc_module, only: www, bcset_1, bcset_3
  use kinetic_variables, only: coef_lap0, coef_lap, zcoef_kin, coef_nab &
       ,flag_nab, flag_n12, flag_n23, flag_n31, const_k2, ggg, Md, wk
  use watch_module, only: watchb_omp, time_kine, time_bcfd

  implicit none

  private
  public :: z_op_kinetic_sol
  public :: z_construct_matrix_kinetic_sol

  integer,allocatable :: ijk(:,:)
  complex(8),allocatable :: zc1a(:),zc2a(:),zc3a(:)

  logical :: flag_clap=.true.
  logical :: flag_clap_setval=.true.
  real(8) :: clap11,clap12,clap13,clap14,clap15,clap16
  real(8) :: clap21,clap22,clap23,clap24,clap25,clap26
  real(8) :: clap31,clap32,clap33,clap34,clap35,clap36

contains


  subroutine z_op_kinetic_sol( tpsi, htpsi, k_in, vloc )
    implicit none
    complex(8),intent(in)    ::  tpsi(:,:)
    complex(8),intent(inout) :: htpsi(:,:)
    real(8),parameter :: zero=0.0d0
    integer,optional,intent(in) :: k_in
    real(8),optional,intent(in) :: vloc(:)

    integer :: i,ib,i1,i2,i3,nb,m,n,j,k,i0
    integer :: a1,a2,a3,b1,b2,b3,p,mm,nn
    integer :: a1b,b1b,a2b,b2b,a3b,b3b,ab1,ab12
    real(8) :: c,d,et0,et1, ttmp(2)
    integer,allocatable :: ic(:)
    integer :: a1b_omp,b1b_omp,a2b_omp,b2b_omp,a3b_omp,b3b_omp,n1_omp,n2_omp
    integer :: ib1_omp,ib2_omp,nb_omp
    logical :: disp

    !call watchb_omp( ttmp )
!!$OMP master
!    time_bcfd(:,:)=0.0d0
!!$OMP end master

    k=1 ; if ( present(k_in) ) k=k_in

    a1b = Igrid(1,1)
    b1b = Igrid(2,1)
    a2b = Igrid(1,2)
    b2b = Igrid(2,2)
    a3b = Igrid(1,3)
    b3b = Igrid(2,3)
    ab1 = (b1b-a1b+1)
    ab12= (b1b-a1b+1)*(b2b-a2b+1)

    nb = size( tpsi, 2 )
    i0 = Igrid(1,0)

!$omp single
    if ( flag_clap_setval ) then
      call check_disp_switch( disp, 0 )
      ! if ( disp ) write(*,*) "flag_clap=",flag_clap
      if ( Md >= 4 ) then
        clap11=coef_lap(1,1)
        clap21=coef_lap(2,1)
        clap31=coef_lap(3,1)
        clap12=coef_lap(1,2)
        clap22=coef_lap(2,2)
        clap32=coef_lap(3,2)
        clap13=coef_lap(1,3)
        clap23=coef_lap(2,3)
        clap33=coef_lap(3,3)
        clap14=coef_lap(1,4)
        clap24=coef_lap(2,4)
        clap34=coef_lap(3,4)
      end if
      if ( Md == 6 ) then
        clap15=coef_lap(1,5)
        clap25=coef_lap(2,5)
        clap35=coef_lap(3,5)
        clap16=coef_lap(1,6)
        clap26=coef_lap(2,6)
        clap36=coef_lap(3,6)
      end if
      flag_clap_setval=.false.
    end if
!$omp end single

! ---

    mm=0
!$  mm=omp_get_thread_num()

    n1_omp = Igrid_omp(1,0,mm) - i0 + 1
    n2_omp = Igrid_omp(2,0,mm) - i0 + 1

    a1b_omp = Igrid_omp(1,1,mm)
    b1b_omp = Igrid_omp(2,1,mm)
    a2b_omp = Igrid_omp(1,2,mm)
    b2b_omp = Igrid_omp(2,2,mm)
    a3b_omp = Igrid_omp(1,3,mm)
    b3b_omp = Igrid_omp(2,3,mm)

! ---

    c=coef_lap0+const_k2(k)

    if( .not.present(vloc) )then
      do ib=1,nb
      do i=n1_omp,n2_omp
        htpsi(i,ib) = c*tpsi(i,ib)
      end do
      end do
!$omp barrier
    end if

    !call watchb_omp( ttmp, time_kine(1,1) )

    do ib=1,nb
      do i3=a3b_omp,b3b_omp
      do i2=a2b_omp,b2b_omp
      do i1=a1b_omp,b1b_omp
        j=1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
        www(i1,i2,i3,ib) = tpsi(j,ib)
      end do
      end do
      end do
    end do

!$omp barrier
    !call watchb_omp( ttmp, time_kine(1,2) )

    call bcset_3(1,nb,Md,0)

!$omp barrier
    !call watchb_omp( ttmp, time_kine(1,3) )

    if ( flag_nab ) then

      if( present(vloc) )then

        do ib=1,nb
          !do j=n1_omp,n2_omp
          !  htpsi(j,ib) = (c+vloc(j))*tpsi(j,ib)
          !end do
          !!$omp barrier
          do i3=a3b_omp,b3b_omp
          do i2=a2b_omp,b2b_omp
          do i1=a1b_omp,b1b_omp
            j=1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
            htpsi(j,ib) = (c+vloc(j))*tpsi(j,ib)
          end do
          end do
          end do
          do m = 1, Md
            do i3=a3b_omp,b3b_omp
            do i2=a2b_omp,b2b_omp
            do i1=a1b_omp,b1b_omp
              j=1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
              htpsi(j,ib)=htpsi(j,ib) &
                           +zcoef_kin(1,m,k) *www(i1+m,i2,i3,ib) &
                     +conjg(zcoef_kin(1,m,k))*www(i1-m,i2,i3,ib) &
                           +zcoef_kin(2,m,k) *www(i1,i2+m,i3,ib) &
                     +conjg(zcoef_kin(2,m,k))*www(i1,i2-m,i3,ib) &
                           +zcoef_kin(3,m,k) *www(i1,i2,i3+m,ib) &
                     +conjg(zcoef_kin(3,m,k))*www(i1,i2,i3-m,ib)
            end do
            end do
            end do
          end do !m
        end do !ib
       
      else !present(vloc)

        do ib=1,nb
          do m = 1, Md
            do i3=a3b_omp,b3b_omp
            do i2=a2b_omp,b2b_omp
            do i1=a1b_omp,b1b_omp
              j=1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
              htpsi(j,ib)=htpsi(j,ib) &
                           +zcoef_kin(1,m,k) *www(i1+m,i2,i3,ib) &
                     +conjg(zcoef_kin(1,m,k))*www(i1-m,i2,i3,ib) &
                           +zcoef_kin(2,m,k) *www(i1,i2+m,i3,ib) &
                     +conjg(zcoef_kin(2,m,k))*www(i1,i2-m,i3,ib) &
                           +zcoef_kin(3,m,k) *www(i1,i2,i3+m,ib) &
                     +conjg(zcoef_kin(3,m,k))*www(i1,i2,i3-m,ib)
            end do
            end do
            end do
          end do !m
        end do !ib

      end if !present(vloc)

    else !flag_nab

      if ( present(vloc) ) then

        if ( flag_clap .and. Md == 4 ) then

          do ib=1,nb
            do i3=a3b_omp,b3b_omp
            do i2=a2b_omp,b2b_omp
            do i1=a1b_omp,b1b_omp
              j=1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
              htpsi(j,ib)= &
                   +clap14*www(i1-4,i2,i3,ib) &
                   +clap13*www(i1-3,i2,i3,ib) &
                   +clap12*www(i1-2,i2,i3,ib) &
                   +clap11*www(i1-1,i2,i3,ib) &
                   +(c+vloc(j))*www(i1,i2,i3,ib) &
                   +clap11*www(i1+1,i2,i3,ib) &
                   +clap12*www(i1+2,i2,i3,ib) &
                   +clap13*www(i1+3,i2,i3,ib) &
                   +clap14*www(i1+4,i2,i3,ib) &
                   +clap24*www(i1,i2-4,i3,ib) &
                   +clap23*www(i1,i2-3,i3,ib) &
                   +clap22*www(i1,i2-2,i3,ib) &
                   +clap21*www(i1,i2-1,i3,ib) &
                   +clap21*www(i1,i2+1,i3,ib) &
                   +clap22*www(i1,i2+2,i3,ib) &
                   +clap23*www(i1,i2+3,i3,ib) &
                   +clap24*www(i1,i2+4,i3,ib) &
                   +clap34*www(i1,i2,i3-4,ib) &
                   +clap33*www(i1,i2,i3-3,ib) &
                   +clap32*www(i1,i2,i3-2,ib) &
                   +clap31*www(i1,i2,i3-1,ib) &
                   +clap31*www(i1,i2,i3+1,ib) &
                   +clap32*www(i1,i2,i3+2,ib) &
                   +clap33*www(i1,i2,i3+3,ib) &
                   +clap34*www(i1,i2,i3+4,ib)
            end do
            end do
            end do
          end do !ib
         
        else if ( flag_clap .and. Md == 6 ) then

          do ib=1,nb
            do i3=a3b_omp,b3b_omp
            do i2=a2b_omp,b2b_omp
            do i1=a1b_omp,b1b_omp
              j=1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
              htpsi(j,ib)= &
                   +clap16*www(i1-6,i2,i3,ib) &
                   +clap15*www(i1-5,i2,i3,ib) &
                   +clap14*www(i1-4,i2,i3,ib) &
                   +clap13*www(i1-3,i2,i3,ib) &
                   +clap12*www(i1-2,i2,i3,ib) &
                   +clap11*www(i1-1,i2,i3,ib) &
                   +(c+vloc(j))*www(i1,i2,i3,ib) &
                   +clap11*www(i1+1,i2,i3,ib) &
                   +clap12*www(i1+2,i2,i3,ib) &
                   +clap13*www(i1+3,i2,i3,ib) &
                   +clap14*www(i1+4,i2,i3,ib) &
                   +clap15*www(i1+5,i2,i3,ib) &
                   +clap16*www(i1+6,i2,i3,ib) &
                   +clap26*www(i1,i2-6,i3,ib) &
                   +clap25*www(i1,i2-5,i3,ib) &
                   +clap24*www(i1,i2-4,i3,ib) &
                   +clap23*www(i1,i2-3,i3,ib) &
                   +clap22*www(i1,i2-2,i3,ib) &
                   +clap21*www(i1,i2-1,i3,ib) &
                   +clap21*www(i1,i2+1,i3,ib) &
                   +clap22*www(i1,i2+2,i3,ib) &
                   +clap23*www(i1,i2+3,i3,ib) &
                   +clap24*www(i1,i2+4,i3,ib) &
                   +clap25*www(i1,i2+5,i3,ib) &
                   +clap26*www(i1,i2+6,i3,ib) &
                   +clap36*www(i1,i2,i3-6,ib) &
                   +clap35*www(i1,i2,i3-5,ib) &
                   +clap34*www(i1,i2,i3-4,ib) &
                   +clap33*www(i1,i2,i3-3,ib) &
                   +clap32*www(i1,i2,i3-2,ib) &
                   +clap31*www(i1,i2,i3-1,ib) &
                   +clap31*www(i1,i2,i3+1,ib) &
                   +clap32*www(i1,i2,i3+2,ib) &
                   +clap33*www(i1,i2,i3+3,ib) &
                   +clap34*www(i1,i2,i3+4,ib) &
                   +clap35*www(i1,i2,i3+5,ib) &
                   +clap36*www(i1,i2,i3+6,ib)
            end do
            end do
            end do
          end do !ib

        else ! flag_clap .and. (Md==4.or.Md==6)

          do ib=1,nb
            !do j=n1_omp,n2_omp
            !  htpsi(j,ib) = (c+vloc(j))*tpsi(j,ib)
            !end do
            !!$omp barrier
            do i3=a3b_omp,b3b_omp
            do i2=a2b_omp,b2b_omp
            do i1=a1b_omp,b1b_omp
              j=1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
              htpsi(j,ib) = (c+vloc(j))*tpsi(j,ib)
            end do
            end do
            end do
            do m=1,Md
              do i3=a3b_omp,b3b_omp
              do i2=a2b_omp,b2b_omp
              do i1=a1b_omp,b1b_omp
                j=1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
                htpsi(j,ib)=htpsi(j,ib) &
                  +coef_lap(1,m)*( www(i1-m,i2,i3,ib)+www(i1+m,i2,i3,ib) ) &
                  +coef_lap(2,m)*( www(i1,i2-m,i3,ib)+www(i1,i2+m,i3,ib) ) &
                  +coef_lap(3,m)*( www(i1,i2,i3-m,ib)+www(i1,i2,i3+m,ib) )
              end do
              end do
              end do
            end do ! m
          end do !ib
         
        end if ! flag_clap .and. Md==4

      else !present(vloc)

        do ib=1,nb
          do m=1,Md
            do i3=a3b_omp,b3b_omp
            do i2=a2b_omp,b2b_omp
            do i1=a1b_omp,b1b_omp
              j=1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
              htpsi(j,ib)=htpsi(j,ib) &
                   +coef_lap(1,m)*( www(i1-m,i2,i3,ib)+www(i1+m,i2,i3,ib) ) &
                   +coef_lap(2,m)*( www(i1,i2-m,i3,ib)+www(i1,i2+m,i3,ib) ) &
                   +coef_lap(3,m)*( www(i1,i2,i3-m,ib)+www(i1,i2,i3+m,ib) )
            end do
            end do
            end do
          end do ! m
        end do ! ib

      end if !present(vloc)

    end if !flag_nabla

!$omp barrier
    !call watchb_omp( ttmp, time_kine(1,4) )

    if ( flag_n12 .or. flag_n23 .or. flag_n31 ) then

       !call watchb_omp( ttmp )

!$OMP workshare
       wk=www
!$OMP end workshare

       !call watchb_omp( ttmp, time_kine(1,2) )

       if ( flag_n12 ) then

          !call watchb_omp( ttmp )

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

          !call watchb_omp( ttmp, time_kine(1,4) )

!$OMP barrier
          call bcset_1(1,nb,Md,3)
!$OMP barrier

          !call watchb_omp( ttmp, time_kine(1,3) )

          do ib=1,nb
             do m=1,Md
                d=-ggg(4)*coef_nab(2,m)
                do i3=a3b_omp,b3b_omp
                do i2=a2b_omp,b2b_omp
                do i1=a1b_omp,b1b_omp
                   j=1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
                   htpsi(j,ib)=htpsi(j,ib) &
                        -d*(www(i1,i2-m,i3,ib)-www(i1,i2+m,i3,ib))
                end do
                end do
                end do
             end do
          end do

          !call watchb_omp( ttmp, time_kine(1,4) )

!$OMP barrier

       end if

       if ( flag_n23 ) then

          !call watchb_omp( ttmp )

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

          !call watchb_omp( ttmp, time_kine(1,4) )

!$OMP barrier
          call bcset_1(1,nb,Md,5)
!$OMP barrier

          !call watchb_omp( ttmp, time_kine(1,3) )

          do ib=1,nb
             do m=1,Md
                d=-ggg(5)*coef_nab(3,m)
                do i3=a3b_omp,b3b_omp
                do i2=a2b_omp,b2b_omp
                do i1=a1b_omp,b1b_omp
                   j=1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
                   htpsi(j,ib)=htpsi(j,ib) &
                        -d*(www(i1,i2,i3-m,ib)-www(i1,i2,i3+m,ib))
                end do
                end do
                end do
             end do
          end do

          !call watchb_omp( ttmp, time_kine(1,4) )

!$OMP barrier

       end if

       if ( flag_n31 ) then

          !call watchb_omp( ttmp )

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

          !call watchb_omp( ttmp, time_kine(1,4) )

!$OMP barrier
          call bcset_1(1,nb,Md,1)
!$OMP barrier

          !call watchb_omp( ttmp, time_kine(1,3) )

          do ib=1,nb
             do m=1,Md
                d=-ggg(6)*coef_nab(1,m)
                do i3=a3b_omp,b3b_omp
                do i2=a2b_omp,b2b_omp
                do i1=a1b_omp,b1b_omp
                   j=1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
                   htpsi(j,ib)=htpsi(j,ib) &
                        -d*(www(i1-m,i2,i3,ib)-www(i1+m,i2,i3,ib))
                end do
                end do
                end do
             end do
          end do

          !call watchb_omp( ttmp, time_kine(1,4) )

       end if

    end if

    !call watchb_omp( ttmp, time_kine(1,5) )

!$OMP master
    time_kine(1:2,6:11) = time_kine(1:2,6:11) + time_bcfd(1:2,1:6)
!$OMP end master

  end subroutine z_op_kinetic_sol


  subroutine z_construct_matrix_kinetic_sol( k, ML, Hmat )
    implicit none
    integer,intent(in) :: k,ML
    complex(8),intent(inout) :: Hmat(ML,ML)
    integer :: i,i1,i2,i3,m,n,j,j1,j2,j3
    integer :: ML1,ML2,ML3
    real(8) :: d

    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

    d = coef_lap0 + const_k2(k)

    do i=1,ML
       Hmat(i,i) = Hmat(i,i) + d
    end do

    if ( flag_nab ) then

       do i3=0,ML3-1
       do i2=0,ML2-1
       do i1=0,ML1-1

          i = 1 + i1 + i2*ML1 + i3*ML1*ML2

          do m=1,Md

             j1 = mod(i1+m+ML1,ML1)
             j  = 1 + j1 + i2*ML1 + i3*ML1*ML2
             Hmat(i,j) = Hmat(i,j) + zcoef_kin(1,m,k)

             j1 = mod(i1-m+ML1,ML1)
             j  = 1 + j1 + i2*ML1 + i3*ML1*ML2
             Hmat(i,j) = Hmat(i,j) + conjg( zcoef_kin(1,m,k) )

             j2 = mod(i2+m+ML2,ML2)
             j  = 1 + i1 + j2*ML1 + i3*ML1*ML2
             Hmat(i,j) = Hmat(i,j) + zcoef_kin(2,m,k)

             j2 = mod(i2-m+ML2,ML2)
             j  = 1 + i1 + j2*ML1 + i3*ML1*ML2
             Hmat(i,j) = Hmat(i,j) + conjg( zcoef_kin(2,m,k) )

             j3 = mod(i3+m+ML3,ML3)
             j  = 1 + i1 + i2*ML1 + j3*ML1*ML2
             Hmat(i,j) = Hmat(i,j) + zcoef_kin(3,m,k)

             j3 = mod(i3-m+ML3,ML3)
             j  = 1 + i1 + i2*ML1 + j3*ML1*ML2
             Hmat(i,j) = Hmat(i,j) + conjg( zcoef_kin(3,m,k) )

          end do ! m

       end do ! i1
       end do ! i2
       end do ! i3

    else

       do i3=0,ML3-1
       do i2=0,ML2-1
       do i1=0,ML1-1

          i = 1 + i1 + i2*ML1 + i3*ML1*ML2

          do m=1,Md

             j1 = mod(i1+m+ML1,ML1)
             j  = 1 + j1 + i2*ML1 + i3*ML1*ML2
             Hmat(i,j) = Hmat(i,j) + coef_lap(1,m)

             j1 = mod(i1-m+ML1,ML1)
             j  = 1 + j1 + i2*ML1 + i3*ML1*ML2
             Hmat(i,j) = Hmat(i,j) + coef_lap(1,m)

             j2 = mod(i2+m+ML2,ML2)
             j  = 1 + i1 + j2*ML1 + i3*ML1*ML2
             Hmat(i,j) = Hmat(i,j) + coef_lap(2,m)

             j2 = mod(i2-m+ML2,ML2)
             j  = 1 + i1 + j2*ML1 + i3*ML1*ML2
             Hmat(i,j) = Hmat(i,j) + coef_lap(2,m)

             j3 = mod(i3+m+ML3,ML3)
             j  = 1 + i1 + i2*ML1 + j3*ML1*ML2
             Hmat(i,j) = Hmat(i,j) + coef_lap(3,m)

             j3 = mod(i3-m+ML3,ML3)
             j  = 1 + i1 + i2*ML1 + j3*ML1*ML2
             Hmat(i,j) = Hmat(i,j) + coef_lap(3,m)

          end do ! m

       end do ! i1
       end do ! i2
       end do ! i3

    end if

    if ( flag_n12 .or. flag_n23 .or. flag_n31 ) then

       if ( flag_n12 ) then

          do i3=0,ML3-1
          do i2=0,ML2-1
          do i1=0,ML1-1

             i = 1 + i1 + i2*ML1 + i3*ML1*ML2

             do n=1,Md
             do m=1,Md

                d = -ggg(4)*coef_nab(1,m)*coef_nab(2,n)

                j1 = mod(i1+m+ML1,ML1)
                j2 = mod(i2+n+ML2,ML2)
                j  = 1 + j1 + j2*ML1 + i3*ML1*ML2
                Hmat(i,j) = Hmat(i,j) + d

                j1 = mod(i1-m+ML1,ML1)
                j2 = mod(i2+n+ML2,ML2)
                j  = 1 + j1 + j2*ML1 + i3*ML1*ML2
                Hmat(i,j) = Hmat(i,j) - d

                j1 = mod(i1+m+ML1,ML1)
                j2 = mod(i2-n+ML2,ML2)
                j  = 1 + j1 + j2*ML1 + i3*ML1*ML2
                Hmat(i,j) = Hmat(i,j) - d

                j1 = mod(i1-m+ML1,ML1)
                j2 = mod(i2-n+ML2,ML2)
                j  = 1 + j1 + j2*ML1 + i3*ML1*ML2
                Hmat(i,j) = Hmat(i,j) + d

             end do ! m
             end do ! n

          end do ! i1
          end do ! i2
          end do ! i3

       end if

       if ( flag_n23 ) then

          do i3=0,ML3-1
          do i2=0,ML2-1
          do i1=0,ML1-1

             i = 1 + i1 + i2*ML1 + i3*ML1*ML2

             do n=1,Md
             do m=1,Md

                d = -ggg(5)*coef_nab(2,m)*coef_nab(3,n)

                j2 = mod(i2+m+ML2,ML2)
                j3 = mod(i3+n+ML3,ML3)
                j  = 1 + i1 + j2*ML1 + j3*ML1*ML2
                Hmat(i,j) = Hmat(i,j) + d

                j2 = mod(i2-m+ML2,ML2)
                j3 = mod(i3+n+ML3,ML3)
                j  = 1 + i1 + j2*ML1 + j3*ML1*ML2
                Hmat(i,j) = Hmat(i,j) - d

                j2 = mod(i2+m+ML2,ML2)
                j3 = mod(i3-n+ML3,ML3)
                j  = 1 + i1 + j2*ML1 + j3*ML1*ML2
                Hmat(i,j) = Hmat(i,j) - d

                j2 = mod(i2-m+ML2,ML2)
                j3 = mod(i3-n+ML3,ML3)
                j  = 1 + i1 + j2*ML1 + j3*ML1*ML2
                Hmat(i,j) = Hmat(i,j) + d

             end do ! m
             end do ! n

          end do ! i1
          end do ! i2
          end do ! i3

       end if

       if ( flag_n31 ) then

          do i3=0,ML3-1
          do i2=0,ML2-1
          do i1=0,ML1-1

             i = 1 + i1 + i2*ML1 + i3*ML1*ML2

             do n=1,Md
             do m=1,Md

                d = -ggg(6)*coef_nab(1,m)*coef_nab(3,n)

                j1 = mod(i1+m+ML1,ML1)
                j3 = mod(i3+n+ML3,ML3)
                j  = 1 + j1 + i2*ML1 + j3*ML1*ML2
                Hmat(i,j) = Hmat(i,j) + d

                j1 = mod(i1-m+ML1,ML1)
                j3 = mod(i3+n+ML3,ML3)
                j  = 1 + j1 + i2*ML1 + j3*ML1*ML2
                Hmat(i,j) = Hmat(i,j) - d

                j1 = mod(i1+m+ML1,ML1)
                j3 = mod(i3-n+ML3,ML3)
                j  = 1 + j1 + i2*ML1 + j3*ML1*ML2
                Hmat(i,j) = Hmat(i,j) - d

                j1 = mod(i1-m+ML1,ML1)
                j3 = mod(i3-n+ML3,ML3)
                j  = 1 + j1 + i2*ML1 + j3*ML1*ML2
                Hmat(i,j) = Hmat(i,j) + d

             end do ! m
             end do ! n

          end do ! i1
          end do ! i2
          end do ! i3

       end if

    end if
 
  end subroutine z_construct_matrix_kinetic_sol


end module z_kinetic_sol_module
