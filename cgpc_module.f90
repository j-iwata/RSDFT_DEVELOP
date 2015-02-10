MODULE cgpc_module

!$  use omp_lib
  use rgrid_module
  use parallel_module
  use kinetic_module, only: SYStype
  use bc_module
  use kinetic_variables, only: Md, ggg
  use array_bound_module, only: ML_0,ML_1

  use cgpc_diag_module
  use cgpc_seitsonen_module
  use cgpc_2_module
  use cgpc_gausseidel_module

  implicit none

  PRIVATE
  PUBLIC :: preconditioning,read_cgpc,read_oldformat_cgpc,init_cgpc

  integer :: mloop
  integer :: iswitch_cgpc

#ifdef _DRSDFT_
  real(8),allocatable :: ftmp2(:,:),gtmp2(:,:)
  real(8),parameter :: zero=0.d0
#else
  complex(8),allocatable :: ftmp2(:,:),gtmp2(:,:)
  complex(8),parameter :: zero=(0.d0,0.d0)
#endif


CONTAINS


  SUBROUTINE read_cgpc(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i
    character(5) :: cbuf,ckey
    mloop=3
    iswitch_cgpc=1
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:5) == "MLOOP" ) then
             backspace(unit)
             read(unit,*) cbuf,mloop
          else if ( ckey(1:4) == "IPC" ) then
             backspace(unit)
             read(unit,*) cbuf,iswitch_cgpc
          end if
       end do
999    continue
       write(*,*) "mloop=",mloop
       write(*,*) "iswitch_cgpc=",iswitch_cgpc
    end if
    call send_cgpc(0)
  END SUBROUTINE read_cgpc


  SUBROUTINE read_oldformat_cgpc(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    if ( rank == 0 ) then
       read(unit,*) mloop
       write(*,*) "mloop=",mloop
    end if
    call send_cgpc(0)
  END SUBROUTINE read_oldformat_cgpc


  SUBROUTINE send_cgpc(rank)
    implicit none
    integer,intent(IN) :: rank
    integer :: ierr
    call mpi_bcast(mloop,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(iswitch_cgpc,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_cgpc


  SUBROUTINE init_cgpc(n1,n2,k,s,dV_in)
    implicit none
    integer,intent(IN) :: n1,n2,k,s
    real(8),intent(IN) :: dV_in
    select case( iswitch_cgpc )
    case( 0 )
    case( 1 )
!       if ( disp_switch_parallel ) write(*,*) "--- cgpc_1 ---"
       if ( mloop == 0 ) mloop=3
    case( 2 )
!       if ( disp_switch_parallel ) write(*,*) "--- cgpc_2 ---"
       call init_cgpc_2( mloop, SYStype )
    case( 3 )
       if ( disp_switch_parallel ) write(*,*) "--- cgpc_gausseidel ---"
       call init_cgpc_gausseidel(Igrid,Ngrid,Hgrid,ggg,mloop)
    case( 10 )
       if ( disp_switch_parallel ) write(*,*) "--- cgpc_diag ---"
       call init_cgpc_diag(n1,n2,k,s,dV_in)
    case( 20 )
       if ( disp_switch_parallel ) write(*,*) "--- cgpc_seitsonen ---"
       call init_cgpc_seitsonen(n1,n2,comm_grid,dV_in)
    end select
  END SUBROUTINE init_cgpc


  SUBROUTINE preconditioning(E,k,s,nn,ML0,xk,gk,Pgk)
    implicit none
    integer,intent(IN) :: k,s,nn,ML0
    real(8),intent(IN) :: E(nn)
#ifdef _DRSDFT_
    real(8),intent(INOUT) :: gk(ML0,nn),Pgk(ML0,nn)
    real(8),intent(IN) :: xk(ML0,nn)
#else
    complex(8),intent(INOUT) :: gk(ML0,nn),Pgk(ML0,nn)
    complex(8),intent(IN) :: xk(ML0,nn)
#endif

    select case( iswitch_cgpc )
    case( 0 )
       return
    case( 1 )
       call preconditioning_1(E,k,s,nn,ML0,gk,Pgk)
    case( 2 )
       call preconditioning_2(E,k,s,nn,ML0,gk,Pgk)
    case( 3 )
       call cgpc_gausseidel(ML0,nn,Pgk)
    case( 10 )
       call cgpc_diag(ML0,nn,Pgk)
    case( 20 )
       call cgpc_seitsonen(k,s,ML0,nn,E,xk,Pgk)
    end select

  END SUBROUTINE preconditioning


  SUBROUTINE preconditioning_1(E,k,s,nn,ML0,gk,Pgk)
    implicit none
    integer,intent(IN)    :: k,s,nn,ML0
    real(8),intent(IN)    :: E(nn)
#ifdef _DRSDFT_
    real(8),intent(INOUT) :: gk(ML0,nn),Pgk(ML0,nn)
    real(8),allocatable :: rk_pc(:,:),pk_pc(:,:)
#else
    complex(8),intent(INOUT) :: gk(ML0,nn),Pgk(ML0,nn)
    complex(8),allocatable :: rk_pc(:,:),pk_pc(:,:)
#endif
    real(8),parameter :: ep=1.d-24
    real(8) :: rr0(nn),rr1(nn),sb(nn),pAp(nn),a(nn),E0(nn),b,c
    integer :: i,n,ierr,iloop !,n1,n2

    E0(:) = 0.d0

    if ( mloop<=0 ) return

    allocate( rk_pc(ML0,nn),pk_pc(ML0,nn) )
    allocate( gtmp2(ML0,nn),ftmp2(ML0,nn) )

!$OMP parallel workshare
    gtmp2(1:ML0,1:nn)=gk(1:ML0,1:nn)
!$OMP end parallel workshare

    call precond_cg_mat(E0,k,s,ML0,nn)

    do n=1,nn
!$OMP parallel do
       do i=1,ML0
          rk_pc(i,n)=gk(i,n)-ftmp2(i,n)
          pk_pc(i,n)=rk_pc(i,n)
       end do
!$OMP end parallel do
    end do

    do n=1,nn
       c=0.d0
!$OMP parallel do reduction(+:c)
       do i=1,ML0
          c=c+abs(rk_pc(i,n))**2
       end do
!$OMP end parallel do
       sb(n)=c
    end do

    call mpi_allreduce(sb,rr0,nn,mpi_real8,mpi_sum,comm_grid,ierr)
    if ( all(rr0(1:nn) < ep) ) then
       deallocate( ftmp2,gtmp2,pk_pc,rk_pc )
       return
    end if

    do iloop=1,mloop+1

!$OMP parallel workshare
       gtmp2(1:ML0,1:nn)=pk_pc(1:ML0,1:nn)
!$OMP end parallel workshare

       call precond_cg_mat(E0,k,s,ML0,nn)

!$OMP parallel

       do n=1,nn
!$OMP single
          c=0.d0
!$OMP end single
!$OMP do reduction(+:c)
          do i=1,ML0
#ifdef _DRSDFT_
             c=c+pk_pc(i,n)*ftmp2(i,n)
#else
             c=c+conjg(pk_pc(i,n))*ftmp2(i,n)
#endif
          end do
!$OMP end do
!$OMP single
          sb(n)=c
!$OMP end single
       end do

!$OMP single
       call mpi_allreduce(sb,pAp,nn,mpi_real8,mpi_sum,comm_grid,ierr)
!$OMP end single

       do n=1,nn
!$OMP single
          a(n)=rr0(n)/pAp(n)
!$OMP end single
!$OMP do
          do i=1,ML0
             rk_pc(i,n)=rk_pc(i,n)-a(n)*ftmp2(i,n)
          end do
!$OMP end do
       end do

       do n=1,nn
!$OMP single
          c=0.d0
!$OMP end single
!$OMP do reduction(+:c)
          do i=1,ML0
             c=c+abs(rk_pc(i,n))**2
          end do
!$OMP end do
!$OMP single
          sb(n)=c
!$OMP end single
       end do

!$OMP single
       call mpi_allreduce(sb,rr1,nn,mpi_real8,mpi_sum,comm_grid,ierr)
!$OMP end single

!$OMP end parallel

       if ( iloop==mloop+1 ) then
          exit
       end if

       do n=1,nn
          b=rr1(n)/rr0(n)
          rr0(n)=rr1(n)
!$OMP parallel do
          do i=1,ML0
             Pgk(i,n)=Pgk(i,n)+a(n)*pk_pc(i,n)
             pk_pc(i,n)=rk_pc(i,n)+b*pk_pc(i,n)
          end do
!$OMP end parallel do
       end do

    end do ! iloop

    deallocate( ftmp2,gtmp2,pk_pc,rk_pc )

    return

  END SUBROUTINE preconditioning_1


  SUBROUTINE precond_cg_mat(E,k,s,mm,nn)
    implicit none
    integer,intent(IN) :: k,s,mm,nn
    real(8),intent(INOUT) :: E(nn)
    select case(SYStype)
    case(0)
       call precond_cg_mat_sol(E,k,s,mm,nn)
    case(1)
       call precond_cg_mat_mol(E,k,s,mm,nn)
    case(3)
       call precond_cg_mat_esm(E,k,s,mm,nn)
    end select
  END SUBROUTINE precond_cg_mat


  SUBROUTINE precond_cg_mat_sol(E,k,s,mm,nn)
    implicit none
    integer,intent(IN) :: k,s,mm,nn
    real(8),intent(INOUT) :: E(nn)
    real(8) :: c,c1,c2,c3,d
    integer :: m,n,i,i1,i2,i3,j !,n1,n2
    integer :: a1b,b1b,a2b,b2b,a3b,b3b
    integer,allocatable :: ic(:)
    integer :: a3b_omp,b3b_omp,n1_omp,mt,nt

!$OMP parallel private(a3b_omp,b3b_omp,n1_omp,d,i,mt,nt)

    nt = 1
!$  nt = omp_get_num_threads()

!$OMP single

    a1b = Igrid(1,1)
    b1b = Igrid(2,1)
    a2b = Igrid(1,2)
    b2b = Igrid(2,2)
    a3b = Igrid(1,3)
    b3b = Igrid(2,3)

    c  =  ggg(1)/Hgrid(1)**2+ggg(2)/Hgrid(2)**2+ggg(3)/Hgrid(3)**2
    c1 = -0.5d0/Hgrid(1)**2*ggg(1)
    c2 = -0.5d0/Hgrid(2)**2*ggg(2)
    c3 = -0.5d0/Hgrid(3)**2*ggg(3)

    allocate( ic(0:nt-1) )
    ic(:)=(b3b-a3b+1)/nt
    mt=(b3b-a3b+1)-sum(ic)
    do i=0,mt-1
       ic(i)=ic(i)+1
    end do
!$OMP end single

    mt=0
!$  mt=omp_get_thread_num()
    a3b_omp=a3b+sum(ic(0:mt))-ic(mt)
    b3b_omp=a3b_omp+ic(mt)-1
    n1_omp=1+(a3b_omp-a3b)*(b2b-a2b+1)*(b1b-a1b+1)
!$OMP barrier

!$OMP single
    deallocate( ic )
!$OMP end single

    do n=1,nn
       d=c-E(n)
!$OMP do
       do i=1,mm
          ftmp2(i,n)=d*gtmp2(i,n)
       end do
!$OMP end do
    end do

!$OMP workshare
    www(:,:,:,:)=zero
!$OMP end workshare
    do n=1,nn
       i=n1_omp-1
       do i3=a3b_omp,b3b_omp
       do i2=a2b,b2b
       do i1=a1b,b1b
          i=i+1 ; www(i1,i2,i3,n)=gtmp2(i,n)
       end do
       end do
       end do
    end do
!$OMP barrier

!$OMP single
    call bcset(1,nn,1,0)
!$OMP end single

    do n=1,nn
       i=n1_omp-1
       do i3=a3b_omp,b3b_omp
       do i2=a2b,b2b
       do i1=a1b,b1b
          i=i+1
          ftmp2(i,n)=ftmp2(i,n)+c1*( www(i1-1,i2,i3,n)+www(i1+1,i2,i3,n) ) &
                               +c2*( www(i1,i2-1,i3,n)+www(i1,i2+1,i3,n) ) &
                               +c3*( www(i1,i2,i3-1,n)+www(i1,i2,i3+1,n) )
       end do
       end do
       end do
    end do
!$OMP barrier

!$OMP end parallel

    return
  END SUBROUTINE precond_cg_mat_sol


  SUBROUTINE precond_cg_mat_mol(E,k,s,mm,nn)
    use rgrid_mol_module
    implicit none
    integer,intent(IN) :: k,s,mm,nn
    real(8),intent(INOUT) :: E(nn)
    real(8) :: c,c1,c2,c3,d
    integer :: m,n,i,i1,i2,i3,j

    c  =  3.d0/Hsize**2
    c1 = -0.5d0/Hsize**2
    c2 = -0.5d0/Hsize**2
    c3 = -0.5d0/Hsize**2

!$OMP parallel

    do n=1,nn
       d=c-E(n)
!$OMP do
       do i=1,mm
          ftmp2(i,n)=d*gtmp2(i,n)
       end do
!$OMP end do
    end do

!$OMP workshare
    www(:,:,:,:)=zero
!$OMP end workshare
    do n=1,nn
       do i=ML_0,ML_1
          www( LL(1,i),LL(2,i),LL(3,i),n ) = gtmp2(i-ML_0+1,n)
       end do
    end do
!$OMP barrier

!$OMP single
    call bcset(1,nn,1,0)
!$OMP end single

    do n=1,nn
!$OMP do
       do i=ML_0,ML_1
          i1=LL(1,i)
          i2=LL(2,i)
          i3=LL(3,i)
          j=i-ML_0+1
          ftmp2(j,n)=ftmp2(j,n)+c1*( www(i1-1,i2,i3,n)+www(i1+1,i2,i3,n) ) &
                               +c2*( www(i1,i2-1,i3,n)+www(i1,i2+1,i3,n) ) &
                               +c3*( www(i1,i2,i3-1,n)+www(i1,i2,i3+1,n) )
       end do
!$OMP end do
    end do
!$OMP barrier

!$OMP end parallel

    return
  END SUBROUTINE precond_cg_mat_mol


  SUBROUTINE precond_cg_mat_esm(E,k,s,mm,nn)
    use esm_rgrid_module
    implicit none
    integer,intent(IN) :: k,s,mm,nn
    real(8),intent(INOUT) :: E(nn)
    real(8) :: c,c1,c2,c3,d
    integer :: m,n,i,i1,i2,i3,j,m1,m2,m3

    c  =  1.d0/Hgrid(1)**2 + 1.d0/Hgrid(2)**2 + 1.d0/Hgrid(3)**2
    c1 = -0.5d0/Hgrid(1)**2
    c2 = -0.5d0/Hgrid(2)**2
    c3 = -0.5d0/Hgrid(3)**2
    m1 = Nshift_ESM(1)
    m2 = Nshift_ESM(2)
    m3 = Nshift_ESM(3)

!$OMP parallel

    do n=1,nn
       d=c-E(n)
!$OMP do
       do i=1,mm
          ftmp2(i,n)=d*gtmp2(i,n)
       end do
!$OMP end do
    end do
!$OMP workshare
    www(:,:,:,:)=zero
!$OMP end workshare
    do n=1,nn
       do i=ML0_ESM,ML1_ESM
          i1=LL_ESM(1,i)+m1
          i2=LL_ESM(2,i)+m2
          i3=LL_ESM(3,i)+m3
          www(i1,i2,i3,n) = gtmp2(i-ML0_ESM+1,n)
       end do
    end do
!$OMP barrier

!$OMP single
    call bcset(1,nn,Md,0)
!$OMP end single
    do n=1,nn
       do i=MK0_ESM,MK1_ESM
          i1=KK(1,i)+m1
          i2=KK(2,i)+m2
          i3=KK(3,i)+m3
          www(i1,i2,i3,n) = zero
       end do
    end do
!$OMP barrier

    do n=1,nn
!$OMP do
       do i=ML0_ESM,ML1_ESM
          i1=LL_ESM(1,i)+m1
          i2=LL_ESM(2,i)+m2
          i3=LL_ESM(3,i)+m3
          j=i-ML0_ESM+1
          ftmp2(j,n)=ftmp2(j,n)+c1*( www(i1-1,i2,i3,n)+www(i1+1,i2,i3,n) ) &
                               +c2*( www(i1,i2-1,i3,n)+www(i1,i2+1,i3,n) ) &
                               +c3*( www(i1,i2,i3-1,n)+www(i1,i2,i3+1,n) )
       end do
!$OMP end do
    end do
!$OMP barrier

!$OMP end parallel

    return
  END SUBROUTINE precond_cg_mat_esm


END MODULE cgpc_module
