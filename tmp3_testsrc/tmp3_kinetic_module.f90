MODULE kinetic_module

  use rgrid_module, only: Hgrid,Igrid
  use kinetic_mol_module
  use esm_kinetic_module
  use op_kinetic_sol_test_module

  implicit none

  PRIVATE
  PUBLIC :: Md, ggg, get_ggg_kinetic, get_coef_kinetic &
           ,read_kinetic, op_kinetic &
           ,coef_lap0,coef_lap,SYStype,coef_nab,a2x_nab,a2y_nab,a2z_nab &
           ,read_oldformat_kinetic

  integer :: Md
  real(8) :: ggg(6)
  real(8) :: coef_lap0
  real(8),allocatable :: coef_lap(:,:),coef_nab(:,:)
  real(8),allocatable :: coef_nabk(:,:,:),const_k2(:)
  complex(8),allocatable :: zcoef_kin(:,:,:)
  logical :: flag_nab,flag_n12,flag_n23,flag_n31
  integer :: SYStype=0
  real(8) :: a2x_nab(3),a2y_nab(3),a2z_nab(3)

CONTAINS

  SUBROUTINE read_kinetic(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i
    character(7) :: cbuf,ckey
    Md = 6
    SYStype = 0
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:2) == "MD" ) then
             backspace(unit)
             read(unit,*) cbuf,Md
          else if ( ckey(1:7) == "SYSTYPE" ) then
             backspace(unit)
             read(unit,*) cbuf,SYStype
          end if
       end do
999    continue
       write(*,*) "Md =",Md
       write(*,*) "SYStype =",SYStype
    end if
    call send_kinetic(0)
  END SUBROUTINE read_kinetic


  SUBROUTINE read_oldformat_kinetic(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    if ( rank == 0 ) then
       read(unit,*) Md, SYStype
       write(*,*) "Md =",Md
       write(*,*) "SYStype =",SYStype
    end if
    call send_kinetic(0)
  END SUBROUTINE read_oldformat_kinetic


  SUBROUTINE send_kinetic(rank)
    implicit none
    integer,intent(IN) :: rank
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(Md,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(SYStype,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_kinetic


  SUBROUTINE get_ggg_kinetic(aa,bb)
    implicit none
    real(8),intent(IN) :: aa(3,3),bb(3,3)
    real(8) :: const,a1,a2,a3
    const=1.d0/(4.d0*acos(-1.d0)**2)
    a1 = sqrt( sum( aa(:,1)**2 ) )
    a2 = sqrt( sum( aa(:,2)**2 ) )
    a3 = sqrt( sum( aa(:,3)**2 ) )
    ggg(1) = a1*a1*sum(bb(:,1)*bb(:,1))*const
    ggg(2) = a2*a2*sum(bb(:,2)*bb(:,2))*const
    ggg(3) = a3*a3*sum(bb(:,3)*bb(:,3))*const
    ggg(4) = a1*a2*sum(bb(:,1)*bb(:,2))*const
    ggg(5) = a2*a3*sum(bb(:,2)*bb(:,3))*const
    ggg(6) = a3*a1*sum(bb(:,3)*bb(:,1))*const
  END SUBROUTINE get_ggg_kinetic


  SUBROUTINE get_coef_kinetic(aa,bb,MBZ,kbb,disp_switch,SYStype_in)
    use fd_module
    implicit none
    real(8),intent(IN) :: aa(3,3),bb(3,3)
    integer,intent(IN) :: MBZ
    real(8),intent(IN) :: kbb(3,MBZ)
    logical,intent(IN) :: disp_switch
    integer,optional,intent(IN) :: SYStype_in
    integer :: m,n,k,is,i
    real(8) :: c1,c2,c3,kx,ky,kz,pi2
    real(8) :: a1,a2,a3,H1,H2,H3
    complex(8),parameter :: zi=(0.d0,1.d0)
    real(8),allocatable :: nab(:),lap(:)
    logical :: first_time = .true.

    if ( present(SYStype_in) ) then
       SYStype=SYStype_in
    end if

    pi2 = 2.d0*acos(-1.d0)
    a1  = sqrt(sum(aa(1:3,1)**2))/pi2
    a2  = sqrt(sum(aa(1:3,2)**2))/pi2
    a3  = sqrt(sum(aa(1:3,3)**2))/pi2

    if ( first_time ) then

       first_time = .false.

       allocate( coef_lap(3,Md) )
       allocate( coef_nab(3,Md) )

       allocate( lap(-Md:Md) ) ; lap=0.d0
       allocate( nab(-Md:Md) ) ; nab=0.d0

       call get_coef_laplacian_fd(Md,lap)
       call get_coef_nabla_fd(Md,nab)

       if ( disp_switch ) then
          do i=0,Md
             write(*,'(1x,2f15.10,2x,2f15.10)') lap(i),lap(-i),nab(i),nab(-i)
          end do
       end if

       call get_ggg_kinetic(aa,bb)

       if ( disp_switch ) write(*,'(1x,"ggg=",6f10.5)') ggg

       flag_n12 = .false.
       flag_n23 = .false.
       flag_n31 = .false.
       if ( ggg(4) /= 0.d0 ) flag_n12 = .true.
       if ( ggg(5) /= 0.d0 ) flag_n23 = .true.
       if ( ggg(6) /= 0.d0 ) flag_n31 = .true.

       H1  = Hgrid(1)
       H2  = Hgrid(2)
       H3  = Hgrid(3)

       c1 = -0.5d0*ggg(1)/H1**2
       c2 = -0.5d0*ggg(2)/H2**2
       c3 = -0.5d0*ggg(3)/H3**2

       coef_lap0 = lap(0)*(c1+c2+c3)
       do n=1,Md
          coef_lap(1,n)=lap(n)*c1
          coef_lap(2,n)=lap(n)*c2
          coef_lap(3,n)=lap(n)*c3
       end do

       do n=1,Md
          coef_nab(1,n)=nab(n)/H1
          coef_nab(2,n)=nab(n)/H2
          coef_nab(3,n)=nab(n)/H3
       end do

       if ( disp_switch ) then
          write(*,'(1x,3x,3x,a20,3x,a20)') "lap","coef_lap"
          write(*,'(1x,i3,3x,f20.15,3x,f20.15)') 0,lap(0),coef_lap0
          do n=1,Md
             write(*,'(1x,i3,3x,f20.15,3x,3f20.15)') n,lap(n),coef_lap(1:3,n)
          end do
          write(*,'(1x,3x,3x,a20,3x,a20)') "nab","coef_nab"
          do n=1,Md
             write(*,'(1x,i3,3x,f20.15,3x,3f20.15)') n,nab(n),coef_nab(1:3,n)
          end do
       end if

       deallocate( nab,lap )

       a2x_nab(1) = bb(1,1)*a1
       a2x_nab(2) = bb(1,2)*a2
       a2x_nab(3) = bb(1,3)*a3
       a2y_nab(1) = bb(2,1)*a1
       a2y_nab(2) = bb(2,2)*a2
       a2y_nab(3) = bb(2,3)*a3
       a2z_nab(1) = bb(3,1)*a1
       a2z_nab(2) = bb(3,2)*a2
       a2z_nab(3) = bb(3,3)*a3

    end if

    if ( allocated(const_k2)   ) deallocate( const_k2 )
    if ( allocated(zcoef_kin) ) deallocate( zcoef_kin )
    if ( allocated(coef_nabk) ) deallocate( coef_nabk )
    allocate( coef_nabk(3,Md,MBZ) )
    allocate( zcoef_kin(3,-Md:Md,MBZ) )
    allocate( const_k2(0:MBZ) ) ; const_k2=0.d0

    flag_nab = .false.
    do k=1,MBZ
       kx=bb(1,1)*kbb(1,k)+bb(1,2)*kbb(2,k)+bb(1,3)*kbb(3,k)
       ky=bb(2,1)*kbb(1,k)+bb(2,2)*kbb(2,k)+bb(2,3)*kbb(3,k)
       kz=bb(3,1)*kbb(1,k)+bb(3,2)*kbb(2,k)+bb(3,3)*kbb(3,k)
       c1=a1*( bb(1,1)*kx+bb(2,1)*ky+bb(3,1)*kz )
       c2=a2*( bb(1,2)*kx+bb(2,2)*ky+bb(3,2)*kz )
       c3=a3*( bb(1,3)*kx+bb(2,3)*ky+bb(3,3)*kz )
       if ( c1/=0.d0 .or. c2/=0.d0 .or. c3/=0.d0 ) flag_nab=.true.
       do n=1,Md
          coef_nabk(1,n,k)=coef_nab(1,n)*c1
          coef_nabk(2,n,k)=coef_nab(2,n)*c2
          coef_nabk(3,n,k)=coef_nab(3,n)*c3
       end do
       const_k2(k) = 0.5d0*( kx*kx + ky*ky + kz*kz )
    end do

    do k=1,MBZ
       do n=1,Md
          zcoef_kin(1:3,-n,k)=coef_lap(1:3,n)+zi*coef_nabk(1:3,n,k)
          zcoef_kin(1:3, n,k)=coef_lap(1:3,n)-zi*coef_nabk(1:3,n,k)
       end do
    end do

  END SUBROUTINE get_coef_kinetic


  SUBROUTINE op_kinetic(k,tpsi,htpsi,n1,n2,ib1,ib2)
    use bc_module
    use bc_variables, only: init_bc_test,allocate_bc_test
    use parallel_module, only: MB_d
    implicit none
    integer,intent(IN) :: k,n1,n2,ib1,ib2
#ifdef _DRSDFT_
    real(8),intent(IN)  :: tpsi(n1:n2,ib1:ib2)
    real(8),intent(INOUT) :: htpsi(n1:n2,ib1:ib2)
#else
    complex(8),intent(IN)  :: tpsi(n1:n2,ib1:ib2)
    complex(8),intent(INOUT) :: htpsi(n1:n2,ib1:ib2)
#endif
    select case(SYStype)
    case default
!       call op_kinetic_sol(k,tpsi,htpsi,n1,n2,ib1,ib2)
       if ( flag_first_time ) then
          call init_op_kinetic_sol_test(Md,size(coef_nabk,3) &
               ,coef_lap,coef_nab &
               ,coef_lap0,coef_nabk,const_k2,zcoef_kin,ggg &
               ,flag_nab,flag_n12,flag_n23,flag_n31)
          call init_bc_test(maxval(n_neighbor),fdinfo_send,fdinfo_recv &
               ,n_neighbor)
          call allocate_bc_test(MB_d,Igrid)
       end if
       call op_kinetic_sol_test1(k,tpsi,htpsi,n1,n2,ib1,ib2)
    case(1)
       call op_kinetic_mol(n1,n2,ib1,ib2,tpsi,htpsi)
    case(3)
       call op_esm_kinetic(k,n1,n2,ib1,ib2,tpsi,htpsi)
    end select
  END SUBROUTINE op_kinetic


  SUBROUTINE op_kinetic_sol(k,tpsi,htpsi,n1,n2,ib1,ib2)
    use bc_module, only: www,bcset
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
    integer :: a1b,b1b,a2b,b2b,a3b,b3b
    real(8) :: c,d
    integer,allocatable :: ic(:)
    integer :: a3b_omp,b3b_omp,n1_omp,n2_omp

    a1b = Igrid(1,1)
    b1b = Igrid(2,1)
    a2b = Igrid(1,2)
    b2b = Igrid(2,2)
    a3b = Igrid(1,3)
    b3b = Igrid(2,3)

    nb = ib2-ib1+1

!$OMP parallel private(a3b_omp,b3b_omp,n1_omp,n2_omp,j,p,d,mm,nn)

    nn=1
!$  nn=omp_get_num_threads()
!$OMP single
    allocate( ic(0:nn-1) )
    ic(:)=(b3b-a3b+1)/nn
    mm=(b3b-a3b+1)-sum(ic)
    do i=0,mm-1
       ic(i)=ic(i)+1
    end do
!$OMP end single
    mm=0
!$  mm=omp_get_thread_num()
    a3b_omp=a3b+sum(ic(0:mm))-ic(mm)
    b3b_omp=a3b_omp+ic(mm)-1
    n1_omp=n1+(a3b_omp-a3b)*(b2b-a2b+1)*(b1b-a1b+1)
    n2_omp=n1_omp+(b3b_omp-a3b_omp+1)*(b2b-a2b+1)*(b1b-a1b+1)-1
!$OMP barrier
!$OMP single
    deallocate( ic )
!$OMP end single

    do ib=ib1,ib2
       j=n1_omp-1
       do i3=a3b_omp,b3b_omp
       do i2=a2b,b2b
       do i1=a1b,b1b
          j=j+1
          www(i1,i2,i3,ib-ib1+1)=tpsi(j,ib)
       end do
       end do
       end do
    end do
!$OMP barrier

!$OMP single
    call bcset(1,nb,Md,0)
    c=coef_lap0+const_k2(k)
!$OMP end single

    do ib=ib1,ib2
       do i=n1_omp,n2_omp
          htpsi(i,ib)=htpsi(i,ib)+c*tpsi(i,ib)
       end do
    end do

    if ( flag_nab ) then
       do ib=ib1,ib2
          do m=1,Md
             j=n1_omp-1
             do i3=a3b_omp,b3b_omp
             do i2=a2b,b2b
             do i1=a1b,b1b
                j=j+1
                htpsi(j,ib)=htpsi(j,ib) &
                     +zcoef_kin(1,m,k) *www(i1+m,i2,i3,ib-ib1+1) &
               +conjg(zcoef_kin(1,m,k))*www(i1-m,i2,i3,ib-ib1+1) &
                     +zcoef_kin(2,m,k) *www(i1,i2+m,i3,ib-ib1+1) &
               +conjg(zcoef_kin(2,m,k))*www(i1,i2-m,i3,ib-ib1+1) &
                     +zcoef_kin(3,m,k) *www(i1,i2,i3+m,ib-ib1+1) &
               +conjg(zcoef_kin(3,m,k))*www(i1,i2,i3-m,ib-ib1+1)   
             end do
             end do
             end do
          end do
       end do
    else
       do ib=ib1,ib2
          p=ib-ib1+1
          do m=1,Md
             j=n1_omp-1
             do i3=a3b_omp,b3b_omp
             do i2=a2b,b2b
             do i1=a1b,b1b
                j=j+1
                htpsi(j,ib)=htpsi(j,ib) &
                  +coef_lap(1,m)*( www(i1+m,i2,i3,p)+www(i1-m,i2,i3,p) ) &
                  +coef_lap(2,m)*( www(i1,i2+m,i3,p)+www(i1,i2-m,i3,p) ) &
                  +coef_lap(3,m)*( www(i1,i2,i3+m,p)+www(i1,i2,i3-m,p) )
             end do
             end do
             end do
          end do
       end do
    end if
!$OMP barrier

    if ( flag_n12 .or. flag_n23 .or. flag_n31 ) then

!$OMP single
       a1=a1b-Md ; b1=b1b+Md
       a2=a2b-Md ; b2=b2b+Md
       a3=a3b-Md ; b3=b3b+Md
       allocate( wk(a1:b1,a2:b2,a3:b3,nb) )
!$OMP end single

!$OMP workshare
       wk(:,:,:,:)=www(:,:,:,1:nb)
!$OMP end workshare

       if ( flag_n12 ) then
!$OMP workshare
          www(:,:,:,:)=zero
!$OMP end workshare
          do n=1,nb
             do m=1,Md
                d=coef_nab(1,m)
                do i3=a3b_omp,b3b_omp
                do i2=a2b,b2b
                do i1=a1b,b1b
                   www(i1,i2,i3,n)=www(i1,i2,i3,n) &
                        -d*( wk(i1-m,i2,i3,n)-wk(i1+m,i2,i3,n) )
                end do
                end do
                end do
             end do
          end do
!$OMP barrier
!$OMP single
          call bcset(1,nb,Md,3)
!$OMP end single
          do ib=ib1,ib2
             p=ib-ib1+1
             do m=1,Md
                d=-ggg(4)*coef_nab(2,m)
                j=n1_omp-1
                do i3=a3b_omp,b3b_omp
                do i2=a2b,b2b
                do i1=a1b,b1b
                   j=j+1
                   htpsi(j,ib)=htpsi(j,ib) &
                        -d*(www(i1,i2-m,i3,p)-www(i1,i2+m,i3,p))
                end do
                end do
                end do
             end do
          end do
!$OMP barrier
       end if

       if ( flag_n23 ) then
!$OMP workshare
          www(:,:,:,:)=zero
!$OMP end workshare
          do n=1,nb
             do m=1,Md
                d=coef_nab(2,m)
                do i3=a3b_omp,b3b_omp
                do i2=a2b,b2b
                do i1=a1b,b1b
                   www(i1,i2,i3,n)=www(i1,i2,i3,n) &
                        -d*( wk(i1,i2-m,i3,n)-wk(i1,i2+m,i3,n) )
                end do
                end do
                end do
             end do
          end do
!$OMP barrier
!$OMP single
          call bcset(1,nb,Md,5)
!$OMP end single
          do ib=ib1,ib2
             p=ib-ib1+1
             do m=1,Md
                d=-ggg(5)*coef_nab(3,m)
                j=n1_omp-1
                do i3=a3b_omp,b3b_omp
                do i2=a2b,b2b
                do i1=a1b,b1b
                   j=j+1
                   htpsi(j,ib)=htpsi(j,ib) &
                        -d*(www(i1,i2,i3-m,p)-www(i1,i2,i3+m,p))
                end do
                end do
                end do
             end do
          end do
!$OMP barrier
       end if

       if ( flag_n31 ) then
!$OMP workshare
          www(:,:,:,:)=zero
!$OMP end workshare
          do n=1,nb
             do m=1,Md
                d=coef_nab(3,m)
                do i3=a3b_omp,b3b_omp
                do i2=a2b,b2b
                do i1=a1b,b1b
                   www(i1,i2,i3,n)=www(i1,i2,i3,n) &
                        -d*( wk(i1,i2,i3-m,n)-wk(i1,i2,i3+m,n) )
                end do
                end do
                end do
             end do
          end do
!$OMP barrier
!$OMP single
          call bcset(1,nb,Md,1)
!$OMP end single
          do ib=ib1,ib2
             p=ib-ib1+1
             do m=1,Md
                d=-ggg(6)*coef_nab(1,m)
                j=n1_omp-1
                do i3=a3b_omp,b3b_omp
                do i2=a2b,b2b
                do i1=a1b,b1b
                   j=j+1
                   htpsi(j,ib)=htpsi(j,ib) &
                        -d*(www(i1-m,i2,i3,p)-www(i1+m,i2,i3,p))
                end do
                end do
                end do
             end do
          end do
!$OMP barrier
       end if

!$OMP single
       deallocate( wk )
!$OMP end single

    end if

!$OMP end parallel
 
  END SUBROUTINE op_kinetic_sol

END MODULE kinetic_module
