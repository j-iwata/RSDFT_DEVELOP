MODULE nabla_sym_module

  use grid_module
  use rgrid_module
  use fd_module
  use parallel_module, only: comm_grid, ir_grid, id_grid
  use symmetry_module, only: get_mat_symmetry, isymmetry, &
                             construct_matrix_symmetry
  use kinetic_sym_ini_module, only: get_mat_kinetic_sym_ini
  use enlarge_array_module
  use rsdft_mpi_module

  implicit none

  PRIVATE
  PUBLIC :: init_nabla_sym
  PUBLIC :: op_nabla_sym

  integer,PUBLIC :: isymmetry_nabla=0

  integer,allocatable :: rga(:,:,:)
  real(8),allocatable :: rgb(:,:,:)

  integer :: n_fd_points
  integer,allocatable :: n_fd_pt(:)
  integer,allocatable :: fd_neighbor_pt(:,:,:)
  real(8),allocatable :: fd_coef(:,:)
  integer,allocatable :: LLL(:,:,:)

  real(8),allocatable :: work(:)

  logical :: init_flag=.false.
  real(8),allocatable :: coef_nab(:,:)

CONTAINS


  SUBROUTINE init_nabla_sym

    implicit none
    type( fd ) :: nabla
    integer :: Md, nsym
    integer :: isym,m,j(3),Rj(3),i1,i2,i3,k
    real(8),allocatable :: tmp(:,:,:,:), tmp0(:,:,:,:)
    logical :: disp

    if ( init_flag ) return

    call write_border( 0, " init_nabla_sym(start)" )
    call check_disp_switch( disp, 0 )

    nsym=0
    if ( isymmetry /= 0 ) then
       call get_mat_symmetry( rga, rgb )
    else
       call get_mat_kinetic_sym_ini( rga, rgb )
    end if
    if ( allocated(rga) ) nsym=size(rga,3)
    if ( disp ) write(*,*) "nsym=",nsym
    if ( nsym == 0 ) goto 900

    isymmetry_nabla = 1

    call construct_nabla_fd( nabla )
    Md = nabla%md

    allocate( coef_nab(3,Md) ) ; coef_nab=0.0d0
    coef_nab(1,:)=nabla%coef(1:Md)
    coef_nab(2,:)=nabla%coef(1:Md)
    coef_nab(3,:)=nabla%coef(1:Md)

    call construct_map_3d_to_1d_grid( Ngrid, Igrid, comm_grid, LLL )

    allocate( tmp(-2*Md:2*Md,-2*Md:2*Md,-2*Md:2*Md,3)  ) ; tmp=0.0d0
    allocate( tmp0(-2*Md:2*Md,-2*Md:2*Md,-2*Md:2*Md,3) ) ; tmp0=0.0d0

    i1=0
    i2=0
    i3=0

    do m=1,Md

       do isym=1,nsym

          tmp=0.0d0

          j(1) = i1 + m
          j(2) = i2
          j(3) = i3
          Rj(:) = matmul( rga(:,:,isym), j(:) )
          tmp(Rj(1),Rj(2),Rj(3),1) = tmp(Rj(1),Rj(2),Rj(3),1) + coef_nab(1,m)
          j(1) = i1 - m
          j(2) = i2
          j(3) = i3
          Rj(:) = matmul( rga(:,:,isym), j(:) )
          tmp(Rj(1),Rj(2),Rj(3),1) = tmp(Rj(1),Rj(2),Rj(3),1) - coef_nab(1,m)

          j(1) = i1
          j(2) = i2 + m
          j(3) = i3
          Rj(:) = matmul( rga(:,:,isym), j(:) )
          tmp(Rj(1),Rj(2),Rj(3),2) = tmp(Rj(1),Rj(2),Rj(3),2) + coef_nab(2,m)
          j(1) = i1
          j(2) = i2 - m
          j(3) = i3
          Rj(:) = matmul( rga(:,:,isym), j(:) )
          tmp(Rj(1),Rj(2),Rj(3),2) = tmp(Rj(1),Rj(2),Rj(3),2) - coef_nab(2,m)

          j(1) = i1
          j(2) = i2
          j(3) = i3 + m
          Rj(:) = matmul( rga(:,:,isym), j(:) )
          tmp(Rj(1),Rj(2),Rj(3),3) = tmp(Rj(1),Rj(2),Rj(3),3) + coef_nab(3,m)
          j(1) = i1
          j(2) = i2
          j(3) = i3 - m
          Rj(:) = matmul( rga(:,:,isym), j(:) )
          tmp(Rj(1),Rj(2),Rj(3),3) = tmp(Rj(1),Rj(2),Rj(3),3) - coef_nab(3,m)

          do k=1,3
             tmp0(:,:,:,k) = tmp0(:,:,:,k) + tmp(:,:,:,1)*rgb(k,1,isym) &
                                           + tmp(:,:,:,2)*rgb(k,2,isym) &
                                           + tmp(:,:,:,3)*rgb(k,3,isym)
          end do

       end do ! isym

    end do ! Md

    tmp0=tmp0/dble(nsym)

    do m=1,3
       call init2( 2*Md, m, tmp0(:,:,:,m) )
    end do

    deallocate( tmp0 )
    deallocate( tmp  )

    allocate( work(Ngrid(0)) ) ; work=0.0d0

    init_flag = .true.

900 call write_border( 0, " init_nabla_sym(end)" )

  END SUBROUTINE init_nabla_sym


  SUBROUTINE init2( n, k, tmp )
    implicit none
    integer,intent(IN) :: n, k
    real(8),intent(IN) :: tmp(-n:n,-n:n,-n:n)
    integer :: i1,i2,i3,m
    call write_border( 0, " init2(in nabla_sym, start)" )
    m = count( abs(tmp) > 1.d-10 )
    if ( n_fd_points == 0 ) then
       n_fd_points=m
       allocate( n_fd_pt(3)                      ) ; n_fd_pt=0
       allocate( fd_neighbor_pt(3,n_fd_points,3) ) ; fd_neighbor_pt=0
       allocate( fd_coef(n_fd_points,3)          ) ; fd_coef=0.0d0
    else if ( n_fd_points < m ) then
       n_fd_points=m
       call enlarge_array( fd_neighbor_pt, 3, m, 3 )
       call enlarge_array( fd_coef, m, 3 )
    end if
    m=0
    do i3=-n,n
    do i2=-n,n
    do i1=-n,n
       if ( abs(tmp(i1,i2,i3)) > 1.d-10 ) then
          m=m+1
          fd_neighbor_pt(1,m,k) = i1
          fd_neighbor_pt(2,m,k) = i2
          fd_neighbor_pt(3,m,k) = i3
          fd_coef(m,k) = tmp(i1,i2,i3)
       end if
    end do
    end do
    end do
    n_fd_pt(k)=m
    call write_border( 0, " init2(in nabla_sym, end)" )
  END SUBROUTINE init2


  SUBROUTINE construct_matrix_nabla( k, f, Df )

    use bc_module
    use lattice_module
    !use symmetry_2_module

    implicit none
    integer,intent(IN)  :: k
    real(8),intent(IN)  :: f(:)
    real(8),intent(OUT) :: Df(:)
    integer :: isym,i,i1,i2,i3,m,j1,j2,j3,j,n1,n2,Md
    integer :: mmin,mmax,msum,nsym,loop
    real(8) :: diag_min,diag_max,d,diag_sum
    real(8),allocatable :: Gmat(:,:),Gmat2(:,:),Gmat3(:,:)
    type( grid ) :: rgrid
    type( fd   ) :: nabla
    integer,allocatable :: LLL(:,:,:)
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    complex(8),allocatable :: Gf(:),Gf2(:),Gf3(:),Gmatf(:)
    real(8),allocatable :: SymMat(:,:), w1(:,:), w2(:,:), w3(:,:)
    type(lattice) :: aa,bb
    real(8) :: pi2,b(3,3)

    call write_border( 0, "construct_matrix_nabla(start)" )

    nsym=0
    if ( isymmetry /= 0 ) then
       call get_mat_symmetry( rga, rgb )
    else
       call get_mat_kinetic_sym_ini( rga, rgb )
    end if
    if ( allocated(rga) ) nsym=size(rga,3)

    if ( nsym == 0 ) goto 900

    call construct_nabla_fd( nabla )
    Md = nabla%md

    call get_range_rgrid( rgrid )

    call construct_map_3d_to_1d_grid( Ngrid, Igrid, comm_grid, LLL )

! ---

    call get_aa_lattice( aa )
    call get_reciprocal_lattice( aa, bb )
    pi2      = 2.0d0*acos(-1.0d0)
    b(1:3,1) = aa%Length(1)*bb%LatticeVector(1:3,1)/( pi2*rgrid%spacing(1) )
    b(1:3,2) = aa%Length(2)*bb%LatticeVector(1:3,2)/( pi2*rgrid%spacing(2) )
    b(1:3,3) = aa%Length(3)*bb%LatticeVector(1:3,3)/( pi2*rgrid%spacing(3) )

! ---

    n1=Igrid(1,0)
    n2=Igrid(2,0)

    allocate( Gf(n1:n2)  ) ; Gf =(0.0d0,0.0d0)
    allocate( Gf2(n1:n2) ) ; Gf2=(0.0d0,0.0d0)
    allocate( Gf3(n1:n2) ) ; Gf3=(0.0d0,0.0d0)
    allocate( Gmatf(n1:n2) ) ; Gmatf=(0.0d0,0.0d0)
    allocate( Gmat(rgrid%g1%size_global,rgrid%g1%size_global)  ) ; Gmat=0.0d0
    allocate( Gmat2(rgrid%g1%size_global,rgrid%g1%size_global) ) ; Gmat2=0.0d0
    allocate( Gmat3(rgrid%g1%size_global,rgrid%g1%size_global) ) ; Gmat3=0.0d0
    allocate( SymMat(rgrid%g1%size_global,rgrid%g1%size_global) ) ; SymMat=0.0d0
    allocate( w1(rgrid%g1%size_global,rgrid%g1%size_global) ) ; w1=0.0d0
    allocate( w2(rgrid%g1%size_global,rgrid%g1%size_global) ) ; w2=0.0d0
    allocate( w3(rgrid%g1%size_global,rgrid%g1%size_global) ) ; w3=0.0d0

    www=(0.0d0,0.0d0)
    i=n1-1
    do i3=Igrid(1,3),Igrid(2,3)
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
       i=i+1
       www(i1,i2,i3,1) = f(i-n1+1)
    end do
    end do
    end do
    call bcset_3(1,1,Md,0)

    do m=1,Md

       i=rgrid%g1%head-1
       do i3=Igrid(1,3),Igrid(2,3)
       do i2=Igrid(1,2),Igrid(2,2)
       do i1=Igrid(1,1),Igrid(2,1)
          i=i+1

          Gf(i) =Gf(i) +coef_nab(1,m)*( www(i1+m,i2,i3,1)-www(i1-m,i2,i3,1) )
          Gf2(i)=Gf2(i)+coef_nab(2,m)*( www(i1,i2+m,i3,1)-www(i1,i2-m,i3,1) )
          Gf3(i)=Gf3(i)+coef_nab(3,m)*( www(i1,i2,i3+m,1)-www(i1,i2,i3-m,1) )

          j1 = mod( i1+m+Ngrid(1), Ngrid(1) )
          j2 = i2
          j3 = i3
          j  = LLL(j1,j2,j3)
          Gmat(i,j) = Gmat(i,j) + coef_nab(1,m)
          j1 = mod( i1-m+Ngrid(1), Ngrid(1) )
          j2 = i2
          j3 = i3
          j  = LLL(j1,j2,j3)
          Gmat(i,j) = Gmat(i,j) - coef_nab(1,m)

          j1 = i1
          j2 = mod( i2+m+Ngrid(2), Ngrid(2) )
          j3 = i3
          j  = LLL(j1,j2,j3)
          Gmat2(i,j) = Gmat2(i,j) + coef_nab(2,m)
          j1 = i1
          j2 = mod( i2-m+Ngrid(2), Ngrid(2) )
          j3 = i3
          j  = LLL(j1,j2,j3)
          Gmat2(i,j) = Gmat2(i,j) - coef_nab(2,m)

          j1 = i1
          j2 = i2
          j3 = mod( i3+m+Ngrid(3), Ngrid(3) )
          j  = LLL(j1,j2,j3)
          Gmat3(i,j) = Gmat3(i,j) + coef_nab(3,m)
          j1 = i1
          j2 = i2
          j3 = mod( i3-m+Ngrid(3), Ngrid(3) )
          j  = LLL(j1,j2,j3)
          Gmat3(i,j) = Gmat3(i,j) - coef_nab(3,m)

       end do
       end do
       end do
    end do ! m

    Gmatf=zero
    do i=Igrid(1,0),Igrid(2,0)
       do j=1,Ngrid(0)
          Gmatf(i) = Gmatf(i) + Gmat(i,j)*f(j)
       end do
    end do
    write(*,*) count(abs(Gmat)>1.d-8)
    write(*,*) sum(abs(Gmatf)),sum(abs(Gf)),sum(abs(f))
    write(*,*) sum(abs(Gmatf-Gf))

    Gmatf=zero
    do i=Igrid(1,0),Igrid(2,0)
       do j=1,Ngrid(0)
          Gmatf(i) = Gmatf(i) + Gmat2(i,j)*f(j)
       end do
    end do
    write(*,*) count(abs(Gmat2)>1.d-8)
    write(*,*) sum(abs(Gmatf)),sum(abs(Gf2)),sum(abs(f))
    write(*,*) sum(abs(Gmatf-Gf2))

    Gmatf=zero
    do i=Igrid(1,0),Igrid(2,0)
       do j=1,Ngrid(0)
          Gmatf(i) = Gmatf(i) + Gmat3(i,j)*f(j)
       end do
    end do
    write(*,*) count(abs(Gmat3)>1.d-8)
    write(*,*) sum(abs(Gmatf)),sum(abs(Gf3)),sum(abs(f))
    write(*,*) sum(abs(Gmatf-Gf3))

    write(*,*) "k=",k,sum(abs(rga)),sum(abs(rgb))

    if ( nsym /= 0 ) then

       do loop=1,1

!       w3(:,:)=Gmat(:,:) ! backup
          w3=0.0d0
       do isym=1,nsym
          write(*,'(1x,"isym/nsym=",i3," /",i3)') isym,nsym
          call construct_matrix_symmetry( isym, Ngrid(0), SymMat )
!          call construct_matrix_symmetry_2 &
!               ( rga(:,:,isym), (/0.0d0,0.0d0,0.0d0/), SymMat )

!          w2(:,:)=b(k,1)*Gmat(:,:)+b(k,2)*Gmat2(:,:)+b(k,3)*Gmat3(:,:)

          w2(:,:) = 0.0d0
          w1(:,:) = matmul( Gmat, transpose(SymMat) )
          w2(:,:) = w2(:,:) + matmul( SymMat, w1 )*rgb(k,1,isym)
          w1(:,:) = matmul( Gmat2, transpose(SymMat) )
          w2(:,:) = w2(:,:) + matmul( SymMat, w1 )*rgb(k,2,isym)
          w1(:,:) = matmul( Gmat3, transpose(SymMat) )
          w2(:,:) = w2(:,:) + matmul( SymMat, w1 )*rgb(k,3,isym)

          w3(:,:)=w3(:,:)+w2(:,:)

!          Gmat(:,:) = w2(:,:)
!          select case(k)
!          case(2)
!             Gmat(:,:)=Gmat2(:,:)
!          case(3)
!             Gmat(:,:)=Gmat3(:,:)
!          end select

          Gmatf=0.0d0
          do i=Igrid(1,0),Igrid(2,0)
          do j=1,Ngrid(0)
             Gmatf(i) = Gmatf(i) + w2(i,j)*f(j)
          end do
          end do

          if ( k==1 ) write(*,*) sum(abs(Gmatf)),sum(abs(Gf)),sum(abs(Gmatf-Gf))
          if ( k==2 ) write(*,*) sum(abs(Gmatf)),sum(abs(Gf2)),sum(abs(Gmatf-Gf2))
          if ( k==3 ) write(*,*) sum(abs(Gmatf)),sum(abs(Gf3)),sum(abs(Gmatf-Gf3))
       end do ! isym
       Gmat=w3/nsym

       end do ! loop

    end if

    mmin=Ngrid(0) ; mmax=0 ; msum=0
    diag_min=1.d100 ; diag_max=-1.d100 ; diag_sum=0.0d0
    Gmatf=0.0d0
    do i=Igrid(1,0),Igrid(2,0)
       do j=1,Ngrid(0)
          Gmatf(i) = Gmatf(i) + Gmat(i,j)*f(j)
       end do
       m=count(abs(Gmat(i,:))>1.d-8)
       mmin=min(m,mmin)
       mmax=max(m,mmax)
       msum=msum+m
       d=Gmat(i,i)
       diag_min=min(d,diag_min)
       diag_max=min(d,diag_max)
       diag_sum=diag_sum+d
    end do
    write(*,*) count(abs(Gmat)>1.d-8),mmin,mmax,msum/Ngrid(0)
    write(*,*) diag_min,diag_max,diag_sum/Ngrid(0)
    if ( k==1 ) write(*,*) sum(abs(Gmatf)),sum(abs(Gf)),sum(abs(Gmatf-Gf))
    if ( k==2 ) write(*,*) sum(abs(Gmatf)),sum(abs(Gf2)),sum(abs(Gmatf-Gf2))
    if ( k==3 ) write(*,*) sum(abs(Gmatf)),sum(abs(Gf3)),sum(abs(Gmatf-Gf3))

    Df(:) = Gmatf(:)
!    select case(k)
!    case(1)
!       Df=Gf
!    case(2)
!       Df=Gf2
!    case(3)
!       Df=Gf3
!    end select

! ---

    deallocate( SymMat, w1, w2, w3 )
    deallocate( Gmat, Gmat2, Gmat3 )
    deallocate( Gf, Gf2, Gf3 )
    deallocate( Gmatf )

900 call write_border( 0, "construct_matrix_nabla(end)" )

  END SUBROUTINE construct_matrix_nabla


  SUBROUTINE op_nabla_sym( k, f, Df )

    implicit none
    integer,intent(IN)  :: k
    real(8),intent(IN)  :: f(:)
    real(8),intent(OUT) :: Df(:)
    integer :: n1,n2,i1,i2,i3,i,j1,j2,j3,j,k1,k2,k3

    n1=Igrid(1,0)
    n2=Igrid(2,0)

    work(n1:n2)=f(:)
    call rsdft_allgatherv( work(n1:n2), work, ir_grid, id_grid, comm_grid )

    Df(:)=0.0d0

    i=0
    do i3=Igrid(1,3),Igrid(2,3)
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
       i=i+1

       do j=1,n_fd_pt(k)

          j1 = i1 + fd_neighbor_pt(1,j,k)
          j2 = i2 + fd_neighbor_pt(2,j,k)
          j3 = i3 + fd_neighbor_pt(3,j,k)

          k1 = mod( j1+Ngrid(1),Ngrid(1) )
          k2 = mod( j2+Ngrid(2),Ngrid(2) )
          k3 = mod( j3+Ngrid(3),Ngrid(3) )

          Df(i) = Df(i) + fd_coef(j,k)*work( LLL(k1,k2,k3) )

       end do ! j

    end do ! i1
    end do ! i2
    end do ! i3

  END SUBROUTINE op_nabla_sym


END MODULE nabla_sym_module
