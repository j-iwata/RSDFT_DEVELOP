MODULE esm_cylindrical_test

  use aa_module
  use rgrid_module, only: Ngrid,Hgrid,Igrid,dV
  use modified_bessel_module
  use density_module
  use array_bound_module
  use parallel_module
  use ps_local_module
  use ps_local_rs_module
  use hartree_module

  implicit none

  PRIVATE
  PUBLIC :: esm_test1,esm_test2,esm_test3

  real(8) :: const_rrigt
  real(8) :: radius_esm1, radius_esm2

  integer :: ngrd_radial
  integer,allocatable :: index_xy_to_radial(:,:)
  real(8),allocatable :: xy_grid(:,:,:) !,rho_ps(:)
  real(8),allocatable :: radial_grid(:)

  complex(8),allocatable :: v_outer(:,:,:)

CONTAINS

  SUBROUTINE esm_test1

    integer :: io(2,4)
    integer :: ngrd_xy
    integer,allocatable :: isort(:)
    real(8),allocatable :: rr_grid(:)
    integer :: i,i1,i2,j
    real(8) :: r,rr,rmin,x,y,xmin,ymin

    radius_esm1 =  5.0d0 ! density
    radius_esm2 = 10.0d0 ! ESM

    io(1:2,1) = (/ 0, 0 /)
    io(1:2,2) = (/ Ngrid(1), 0 /)
    io(1:2,3) = (/ 0, Ngrid(2) /)
    io(1:2,4) = (/ Ngrid(1), Ngrid(2) /)

    allocate( xy_grid(2,0:Ngrid(1)-1,0:Ngrid(2)-1) )
    allocate( index_xy_to_radial(0:Ngrid(1)-1,0:Ngrid(2)-1) )

    ngrd_xy=Ngrid(1)*Ngrid(2)
    allocate( rr_grid(ngrd_xy) ) ; rr_grid=0.d0
    allocate( isort(ngrd_xy) ) ; isort=0

    ngrd_xy=0
    do i2=0,Ngrid(2)-1
    do i1=0,Ngrid(1)-1
       ngrd_xy=ngrd_xy+1
       rmin=1.d10
       xmin=1.d10
       ymin=1.d10
       do j=1,4
          x=( i1 - io(1,j) )*Hgrid(1)
          y=( i2 - io(2,j) )*Hgrid(2)
          r=x*x+y*y
          if ( r < rmin ) then
             rmin=r
             xmin=x
             ymin=y
          end if
       end do
       rr_grid(ngrd_xy) = rmin
       xy_grid(1,i1,i2) = xmin
       xy_grid(2,i1,i2) = ymin
    end do
    end do

    call indexx(ngrd_xy,rr_grid,isort)

    rmin=-1.d10
    ngrd_radial=0
    do i=1,ngrd_xy
       r=rr_grid(isort(i))
       if ( r-rmin > 1.d-10 ) then
          ngrd_radial=ngrd_radial+1
          rmin=r
       end if
    end do

    allocate( radial_grid(ngrd_radial) )

    rmin=-1.d10
    ngrd_radial=0
    do i=1,ngrd_xy
       r=rr_grid(isort(i))
       if ( r-rmin > 1.d-10 ) then
          ngrd_radial=ngrd_radial+1
          rmin=r
          radial_grid(ngrd_radial)=sqrt(rmin)
       end if
    end do

    do i2=0,Ngrid(2)-1
    do i1=0,Ngrid(1)-1
       x=xy_grid(1,i1,i2)
       y=xy_grid(2,i1,i2)
       r=sqrt(x*x+y*y)
       do j=1,ngrd_radial
          rr=radial_grid(j)
          if ( abs(r-rr) < 1.d-10 ) then
             index_xy_to_radial(i1,i2) = j
             exit
          end if
       end do
    end do
    end do

    deallocate( isort )
    deallocate( rr_grid )

  END SUBROUTINE esm_test1


  SUBROUTINE esm_test2
    use kinetic_module
    use pseudopot_module
    integer :: ikz,i,i1,i2,i3,j,j1,j2,ma,ima,ierr,n,s
    integer :: a1,a2,a3,b1,b2,b3
    real(8) :: pi4,pi2,r,r1,r2,kz,phi,x,y,z,c,pi
    real(8) :: p1,p2,p3,p4
    real(8),allocatable :: rhot(:,:,:),rwork(:)
    real(8),allocatable :: work(:,:,:),vtmp(:)
    complex(8),parameter :: zi=(0.0d0,1.0d0),zero=(0.d0,0.d0)
    complex(8) :: phase_k,phase_m
    complex(8),allocatable :: fmk_inner(:,:)
    integer :: ikz_max,ima_max
    real(8),allocatable :: Rmk_inner(:,:,:),Rmk_outer(:,:,:)

    pi  = acos(-1.d0)
    pi2 = 2.d0*acos(-1.d0)
    pi4 = 4.d0*acos(-1.d0)

    a1 = Igrid(1,1)
    b1 = Igrid(2,1)
    a2 = Igrid(1,2)
    b2 = Igrid(2,2)
    a3 = Igrid(1,3)
    b3 = Igrid(2,3)

    write(*,*) "sum(rho_ps)=",sum(rho_ps)*dV
    write(*,*) "minval(rho_ps)=",minval(rho_ps)
    write(*,*) "maxval(rho_ps)=",maxval(rho_ps)

    allocate( v_outer(a1-Md:b1+Md,a2-Md:b2+Md,a3-Md:b3+Md) ) ; v_outer=zero
    allocate( rhot(a1:b1,a2:b2,a3:b3) ) ; rhot=0.d0

    c = sum(rho_ps)/Ngrid(0)

    i=ML_0-1
    do i3=Igrid(1,3),Igrid(2,3)
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
       i=i+1
       rhot(i1,i2,i3)=rho_ps(i)-c
    end do
    end do
    end do

    write(*,*) "sum(rhot)=",sum(rhot)*dV
    write(*,*) "minval(rhot)=",minval(rhot)
    write(*,*) "maxval(rhot)=",maxval(rhot)

    c = -2.d0*dV/aa(3,3)
    ikz_max = 8
    ima_max = 4

    allocate( Rmk_inner(ngrd_radial,-ima_max:ima_max,-ikz_max:ikz_max) )
    allocate( Rmk_outer(ngrd_radial,-ima_max:ima_max,-ikz_max:ikz_max) )
    allocate( fmk_inner(-ima_max:ima_max,-ikz_max:ikz_max) )
    Rmk_inner = 0.d0
    Rmk_outer = 0.d0
    fmk_inner = zero

    do ikz=-ikz_max,ikz_max
       if ( ikz == 0 ) cycle
       kz=pi2/aa(3,3)*ikz
    do ima=-ima_max,ima_max
       ma=abs(ima)

       x=abs(kz)*radius_esm2
       const_rrigt = bessk(ma,x)/bessi(ma,x)

       do j=1,ngrd_radial
          r=radial_grid(j)
          x=abs(kz)*r
          if ( r <= radius_esm1 ) then
             Rmk_inner(j,ima,ikz) = bessi(ma,x)
          else
!             Rmk_outer(j,ima,ikz) = -bessk(ma,x)
             Rmk_outer(j,ima,ikz) = Rrigt(ma,x)
          end if
       end do

       do i3=Igrid(1,3),Igrid(2,3)
          z=Hgrid(3)*i3 !-0.5d0*aa(3,3)
          phase_k = dcmplx( cos(kz*z),-sin(kz*z) )
          do i2=Igrid(1,2),Igrid(2,2)
          do i1=Igrid(1,1),Igrid(2,1)
             j1=index_xy_to_radial(i1,i2)
             r1=radial_grid(j1)
             if ( radius_esm1 < r1 ) cycle
             x = xy_grid(1,i1,i2)
             y = xy_grid(2,i1,i2)
             if ( r1 == 0.d0 ) then
                phi = 0.d0
             else
                phi = sign( acos(x/r1),y )
             end if
             phase_m = dcmplx( cos(ima*phi),-sin(ima*phi) )
             fmk_inner(ima,ikz) = fmk_inner(ima,ikz) &
                  + phase_k*phase_m*Rmk_inner(j1,ima,ikz)*rhot(i1,i2,i3)
          end do ! i1
          end do ! i2
       end do ! i3

       do i3=Igrid(1,3),Igrid(2,3)
          z=Hgrid(3)*i3 !-0.5d0*aa(3,3)
          phase_k = dcmplx( cos(kz*z),sin(kz*z) )
          do i2=Igrid(1,2),Igrid(2,2)
          do i1=Igrid(1,1),Igrid(2,1)
             j1=index_xy_to_radial(i1,i2)
             r1=radial_grid(j1)
             if ( r1 <= radius_esm1 .or. radius_esm2 < r1 ) cycle
             x = xy_grid(1,i1,i2)
             y = xy_grid(2,i1,i2)
             if ( r1 == 0.d0 ) then
                phi = 0.d0
             else
                phi = sign( acos(x/r1),y )
             end if
             phase_m = dcmplx( cos(ima*phi),sin(ima*phi) )
             v_outer(i1,i2,i3) = v_outer(i1,i2,i3) &
                  +phase_k*phase_m*Rmk_outer(j1,ima,ikz) &
                  *fmk_inner(ima,ikz)*c
          end do ! i1
          end do ! i2
       end do ! i3

    end do ! ima
    end do ! ikz

    v_outer = v_outer + 0.5d0

!    p1=parloc(1,1) ; p2=sqrt(parloc(2,1))
!    p3=parloc(3,1) ; p4=sqrt(parloc(4,1))
!    v_outer(:,:,:)=0.d0
!    do i3=Igrid(1,3)-Md,Igrid(2,3)+Md
!       z=i3*Hgrid(3)-aa(3,3)*0.5d0
!    do i2=Igrid(1,2)-Md,Igrid(2,2)+Md
!    do i1=Igrid(1,1)-Md,Igrid(2,1)+Md
!       x=i1*Hgrid(1)-aa(1,1)*0.5d0
!       y=i2*Hgrid(2)-aa(2,2)*0.5d0!       r=sqrt( x**2 + y**2 + z**2 )
!       if ( r < 1.d-9 ) then
!          v_outer(i1,i2,i3)=-2.d0*Zps(1)/sqrt(pi)*(p1*p2+p3*p4)
!       else
!          v_outer(i1,i2,i3)=-Zps(1)/r*( p1*erf(p2*r)+p3*erf(p4*r) )
!       end if
!    end do
!    end do
!    end do
!    i=0
!    do i3=Igrid(1,3)-Md,Igrid(2,3)+Md
!       z=i3*Hgrid(3)-aa(3,3)*0.5d0
!       do i2=Igrid(1,2),Igrid(2,2)
!       do i1=Igrid(1,1),Igrid(2,1)
!          j1=index_xy_to_radial(i1,i2)
!          r1=radial_grid(j1)
!          if ( radius_esm2 < r1 ) cycle
!          i=i+1
!          r=sqrt( r1**2 + z**2 )
!          if ( r < 1.d-9 ) then
!             v_outer(i1,i2,i3)=-2.d0*Zps(1)/sqrt(pi)*(p1*p2+p3*p4)
!          else
!             v_outer(i1,i2,i3)=-Zps(1)/r*( p1*erf(p2*r)+p3*erf(p4*r) )
!          end if
!       end do
!       end do
!    end do

    write(*,*) i,count( v_outer(a1:b1,a2:b2,a3:b3)/=0.d0 )
    write(*,*) i,count( v_outer(a1:b1,a2:b2,a3-Md:b3+Md)/=0.d0 )

    write(*,*) "sum(v_outer)=",sum(v_outer)

    rewind 10
    do i3=a3-Md,b3+Md
    do i2=a2-Md,b2+Md
    do i1=a1-Md,b1+Md
       if ( i2==0 .and. i3==0 ) then
          write(10,'(1x,i6,2f20.10,2g20.10)') &
               i1,i1*Hgrid(1),v_outer(i1,i2,i3)
       end if
    end do
    end do
    end do
!    stop
    if ( myrank == 0 ) write(*,*) "------- esm_test3"
    allocate( vtmp(ML_0:ML_1) ) ; vtmp=0.d0
    call random_number(vtmp) ; vtmp=vtmp-0.5d0
    write(*,*) sum(rhot)*dV,sum(rho_ps)*dV
    call esm_test3(ML_0,ML_1,rhot,vtmp)

    deallocate( fmk_inner )
    deallocate( Rmk_inner )
    deallocate( rhot )
    deallocate( v_outer )

  END SUBROUTINE esm_test2

  FUNCTION Rleft(n,x)
    integer :: n
    real(8) :: x,Rleft
    Rleft=bessi(n,x)
  END FUNCTION Rleft

  FUNCTION Rrigt(n,x)
    integer :: n
    real(8) :: x,Rrigt
    Rrigt = const_rrigt*bessi(n,x) - bessk(n,x)
  END FUNCTION Rrigt

  SUBROUTINE esm_test3(n1,n2,rho_in,v_inout)
    use kinetic_module
    use bc_module
    use fd_module
    implicit none
    integer,intent(IN)  :: n1,n2
    real(8),intent(IN)  :: rho_in(n1:n2)
    real(8),intent(INOUT) :: v_inout(n1:n2)
    integer,parameter :: max_iter=2000
    integer :: a1,a2,a3,b1,b2,b3,i,i1,i2,i3,j1,j
    integer :: m,mm,nn,nnn,ierr,iter,mp
    integer,allocatable :: grdmap(:,:)
    real(8),parameter :: tol=1.d-24
    real(8) :: pi4,r1,c1,c2,c3,sum0,sum1,ee,e0,sb(3),rb(3)
    real(8) :: ak,ck,sum2,gg_0,gg_1
    real(8),allocatable :: rk(:),qk(:),v0(:),xk(:)
    real(8),allocatable :: pk(:),b(:),lap(:)
    complex(8),allocatable :: wk(:,:,:)
    complex(8),parameter :: z0=(0.d0,0.d0)
    logical,allocatable :: mat(:,:)

    a1=Igrid(1,1)
    b1=Igrid(2,1)
    a2=Igrid(1,2)
    b2=Igrid(2,2)
    a3=Igrid(1,3)
    b3=Igrid(2,3)

    pi4 = 4.d0*acos(-1.d0)
    c1  = 1.d0/Hgrid(1)**2
    c2  = 1.d0/Hgrid(2)**2
    c3  = 1.d0/Hgrid(3)**2

    allocate( lap(-Md:Md) ) ; lap=0.d0
    call get_coef_laplacian_fd(Md,lap)
    lap(0)=0.5d0*lap(0)

    mp = 0 !Ngrid(3)/2

    nn =0
    mm =0
    nnn=0
    do i3=a3,b3
       do i2=a2,b2
       do i1=a1,b1
          j1=index_xy_to_radial(i1,i2)
          r1=radial_grid(j1)
          if ( r1 <= radius_esm1 ) then
             mm=mm+1
          else if ( r1 <= radius_esm2 ) then
             nn=nn+1
          else
             nnn=nnn+1
          end if
       end do
       end do
    end do
    write(*,*) "# of grid points (<= rad1)         ",mm
    write(*,*) "# of grid points (rad1 < r <= rad2)",nn
    write(*,*) "# of grid points (rad2+rad1)       ",nn+mm
    write(*,*) "# of grid points (other)",nnn

    allocate( grdmap(0:3,mm) ) ; grdmap=0
    i=n1-1
    j=0
    do i3=a3,b3
       do i2=a2,b2
       do i1=a1,b1
          i=i+1
          j1=index_xy_to_radial(i1,i2)
          r1=radial_grid(j1)
          if ( r1 <= radius_esm1 ) then
             j=j+1
             grdmap(0,j)=i
             grdmap(1,j)=i1
             grdmap(2,j)=i2
             grdmap(3,j)=i3
          end if
       end do
       end do
    end do
    write(*,*) "check: j,mm=",j,mm

    allocate( wk(a1-Md:b1+Md,a2-Md:b2+Md,a3-Md:b3+Md) ) ; wk=0.d0
    allocate(  b(mm) )
    allocate( xk(mm) )

    do j=1,mm
       b(j) = -pi4*rho_in( grdmap(0,j) )
    end do

    www=z0
    do j=1,mm
       i1=grdmap(1,j)
       i2=grdmap(2,j)
       i3=grdmap(3,j)
       www(i1,i2,i3,1)=b(j)
    end do

    do j=1,mm
       xk(j) = v_inout( grdmap(0,j) )
!       i1=grdmap(1,j)
!       i2=grdmap(2,j)
!       i3=grdmap(3,j)
!       xk(j)=v_outer(i1,i2,i3)
    end do

    allocate( rk(mm) )
    allocate( pk(mm) )
    allocate( qk(mm) )
    allocate( v0(mm) )

    www(:,:,:,:)=z0
    do j=1,mm
       i1=grdmap(1,j)
       i2=grdmap(2,j)
       i3=grdmap(3,j)
       www(i1,i2,i3,1) = xk(j)
    end do
    do i3=a3,b3
       do i2=a2,b2
       do i1=a1,b1
          j1=index_xy_to_radial(i1,i2)
          r1=radial_grid(j1)
          if ( r1 <= radius_esm1 .or. radius_esm2 < r1 ) cycle
          www(i1,i2,i3,1)=v_outer(i1,i2,i3)
       end do
       end do
    end do
    call bcset(1,1,Md,0)

    rewind 11
    rewind 12
    rewind 13
    do i3=a3-Md,b3+Md
    do i2=a2-Md,b2+Md
    do i1=a1-Md,b1+Md
       if ( i2==0 .and. i3==mp ) then
          write(11,'(1x,i6,f20.10,2g20.10)') i1,i1*Hgrid(1),www(i1,i2,i3,1)
       end if
       if ( i1==0 .and. i3==mp ) then
          write(12,'(1x,i6,f20.10,2g20.10)') i2,i2*Hgrid(2),www(i1,i2,i3,1)
       end if
       if ( i2==0 .and. i1==0 ) then
          write(13,'(1x,i6,f20.10,2g20.10)') i3,i3*Hgrid(3),www(i1,i2,i3,1)
       end if
    end do
    end do
    end do

    qk(:)=0.d0
    do j=1,mm
       i1=grdmap(1,j)
       i2=grdmap(2,j)
       i3=grdmap(3,j)
       do m=0,Md
          qk(j) = qk(j) &
               + c1*lap(m)*( www(i1+m,i2,i3,1) + www(i1-m,i2,i3,1) ) &
               + c2*lap(m)*( www(i1,i2+m,i3,1) + www(i1,i2-m,i3,1) ) &
               + c3*lap(m)*( www(i1,i2,i3+m,1) + www(i1,i2,i3-m,1) )
       end do
    end do

    rk(:) = b(:) - qk(:)
    pk(:) = rk(:)

    ee = 0.d0
    gg_1 = 1.d10
    sum1 = 1.d10
    sum2 = 1.d10

    do iter=1,max_iter

       gg_0 = gg_1

       sum0 = sum( rk(:)*rk(:) )*dV
       call mpi_allreduce(sum0,gg_1,1,mpi_real8,mpi_sum,comm_grid,ierr)

       if ( sum1 < tol ) exit

       if ( iter > 1 ) then
          ck = gg_1/gg_0
          pk(:) = rk(:) + ck*pk(:)
       end if

       www(:,:,:,:)=z0
       do j=1,mm
          i1=grdmap(1,j)
          i2=grdmap(2,j)
          i3=grdmap(3,j)
          www(i1,i2,i3,1)=pk(j)
       end do
       do i3=a3,b3
          do i2=a2,b2
          do i1=a1,b1
             j1=index_xy_to_radial(i1,i2)
             r1=radial_grid(j1)
             if ( r1 <= radius_esm1 .or. radius_esm2 < r1 ) cycle
             www(i1,i2,i3,1)=z0
          end do
          end do
       end do
       call bcset(1,1,Md,0)

       qk(:)=0.d0
       do j=1,mm
          i1=grdmap(1,j)
          i2=grdmap(2,j)
          i3=grdmap(3,j)
          do m=0,Md
             qk(j) = qk(j) &
                  + c1*lap(m)*( www(i1+m,i2,i3,1) + www(i1-m,i2,i3,1) ) &
                  + c2*lap(m)*( www(i1,i2+m,i3,1) + www(i1,i2-m,i3,1) ) &
                  + c3*lap(m)*( www(i1,i2,i3+m,1) + www(i1,i2,i3-m,1) )
          end do
       end do

       sum0 = sum( pk(:)*qk(:) )*dV
       call mpi_allreduce(sum0,ak,1,mpi_real8,mpi_sum,comm_grid,ierr)
       ak=gg_1/ak

       v0(:) = xk(:)

       xk(:) = xk(:) + ak*pk(:)
       rk(:) = rk(:) - ak*qk(:)

       e0 = ee

       sb(1) = sum( rk(:)*rk(:) )*dV
       sb(2) = sum( xk(:)*b(:) )*dV *0.5d0/(-pi4)
       sb(3) = sum( (v0(:)-xk(:))**2 )/Ngrid(0)
       call mpi_allreduce(sb,rb,3,mpi_real8,mpi_sum,comm_grid,ierr) 
       sum1=sb(3)
       write(*,'(1x,i6,3g20.10)') iter,rb(1:3)

    end do ! iter

    www=z0
    do j=1,mm
       i1=grdmap(1,j)
       i2=grdmap(2,j)
       i3=grdmap(3,j)
       www(i1,i2,i3,1) = xk(j)
    end do
    do i3=a3,b3
    do i2=a2,b2
    do i1=a1,b1
       j1=index_xy_to_radial(i1,i2)
       r1=radial_grid(j1)
       if ( r1 <= radius_esm1 .or. radius_esm2 < r1 ) cycle
       www(i1,i2,i3,1) = v_outer(i1,i2,i3)
    end do
    end do
    end do
    call bcset(1,1,Md,0)

    rewind 21
    rewind 22
    rewind 23
    do i3=a3-Md,b3+Md
    do i2=a2-Md,b2+Md
    do i1=a1-Md,b1+Md
       if ( i3 == mp .and. i2 == 0 ) then
          write(21,'(1x,i6,f20.10,2g20.10)') i1,i1*Hgrid(1),www(i1,i2,i3,1)
       end if
       if ( i1 == 0 .and. i3 == mp ) then
          write(22,'(1x,i6,f20.10,2g20.10)') i2,i2*Hgrid(2),www(i1,i2,i3,1)
       end if
       if ( i1 == 0 .and. i2 == 0 ) then
          write(23,'(1x,i6,f20.10,2g20.10)') i3,i3*Hgrid(3),www(i1,i2,i3,1)
       end if
    end do
    end do
    end do

    www=z0
    do j=1,mm
       i1=grdmap(1,j)
       i2=grdmap(2,j)
       i3=grdmap(3,j)
       www(i1,i2,i3,1) = xk(j)
    end do
    do i3=a3,b3
       do i2=a2,b2
       do i1=a1,b1
          j1=index_xy_to_radial(i1,i2)
          r1=radial_grid(j1)
          if ( r1 <= radius_esm1 .or. radius_esm2 < r1 ) cycle
          www(i1,i2,i3,1)=v_outer(i1,i2,i3)
       end do
       end do
    end do
    call bcset(1,1,Md,0)

    wk=z0
    do j=1,mm
       i1=grdmap(1,j)
       i2=grdmap(2,j)
       i3=grdmap(3,j)
       do m=0,Md
          wk(i1,i2,i3)=wk(i1,i2,i3) &
               + c1*lap(m)*( www(i1+m,i2,i3,1) + www(i1-m,i2,i3,1) ) &
               + c2*lap(m)*( www(i1,i2+m,i3,1) + www(i1,i2-m,i3,1) ) &
               + c3*lap(m)*( www(i1,i2,i3+m,1) + www(i1,i2,i3-m,1) )
       end do
    end do

    rewind 14
    www=0.d0
    do j=1,mm
       i1=grdmap(1,j)
       i2=grdmap(2,j)
       i3=grdmap(3,j)
       www(i1,i2,i3,1)=b(j)
    end do
    do i1=a1,b1
       write(14,'(1x,i6,f20.10,6g20.10)') i1,i1*Hgrid(1) &
      ,real( wk(i1,0,mp) ),real( wk(0,i1,mp) ),real( wk(0,0,i1) ) &
      ,real( www(i1,0,mp,1) ),real( www(0,i1,mp,1) ),real( www(0,0,i1,1) )
    end do

    deallocate( v0 )
    deallocate( qk ) 
    deallocate( pk ) 
    deallocate( rk )
    deallocate( xk )
    deallocate( b )

    deallocate( grdmap )
    deallocate( wk )
    deallocate( lap )

  END SUBROUTINE esm_test3

  SUBROUTINE esm_test3_test(n1,n2,rho_in,v_inout)
    use kinetic_module
    use bc_module
    use fd_module
    implicit none
    integer,intent(IN)  :: n1,n2
    real(8),intent(IN)  :: rho_in(n1:n2)
    real(8),intent(INOUT) :: v_inout(n1:n2)
    integer,parameter :: max_iter=2000
    integer :: a1,a2,a3,b1,b2,b3,i,i1,i2,i3,j1,j
    integer :: m,mm,nn,nnn,ierr,iter
    integer,allocatable :: grdmap(:,:)
    real(8),parameter :: tol=1.d-24
    real(8) :: pi4,r1,c1,c2,c3,sum0,sum1,ee,e0,sb(3),rb(3)
    real(8) :: d1,d2,d3
    real(8) :: ak,ck,sum2,gg_0,gg_1
    real(8),allocatable :: rk(:),qk(:),v0(:),xk(:)
    real(8),allocatable :: pk(:),b(:),lap(:),nab(:)
    real(8),allocatable :: wr(:,:,:)
    complex(8),allocatable :: wk(:,:,:)
    complex(8),parameter :: z0=(0.d0,0.d0)
    logical,allocatable :: mat(:,:)

    a1=Igrid(1,1)
    b1=Igrid(2,1)
    a2=Igrid(1,2)
    b2=Igrid(2,2)
    a3=Igrid(1,3)
    b3=Igrid(2,3)

    allocate( lap(-Md:Md) ) ; lap=0.d0
    call get_coef_laplacian_fd(Md,lap)
    lap(0)=0.5d0*lap(0)

    c1=1.d0/Hgrid(1)**2
    c2=1.d0/Hgrid(2)**2
    c3=1.d0/Hgrid(3)**2

    pi4 = 4.d0*acos(-1.d0)

    allocate( wk(a1-Md:b1+Md,a2-Md:b2+Md,a3-Md:b3+Md) ) ; wk=z0
    allocate( wr(a1-Md:b1+Md,a2-Md:b2+Md,a3-Md:b3+Md) ) ; wr=0.d0

    allocate( b(n1:n2) )
    allocate( xk(n1:n2) )
    allocate( rk(n1:n2) )
    allocate( pk(n1:n2) )
    allocate( qk(n1:n2) )
    allocate( v0(n1:n2) )

    b(:)  = -pi4*rho_in(:)
    xk(:) = v_inout(:)

    wr(:,:,:)=0.d0
    i=n1-1
    do i3=a3,b3
    do i2=a2,b2
    do i1=a1,b1
       i=i+1
       wr(i1,i2,i3)=b(i)
    end do
    end do
    end do
    www(:,:,:,1) = v_outer(:,:,:)
    wk(:,:,:)=z0
    do i3=a3,b3
    do i2=a2,b2
    do i1=a1,b1
       do m=0,Md
          wk(i1,i2,i3)=wk(i1,i2,i3) &
               + c1*lap(m)*( www(i1+m,i2,i3,1) + www(i1-m,i2,i3,1) ) &
               + c2*lap(m)*( www(i1,i2+m,i3,1) + www(i1,i2-m,i3,1) ) &
               + c3*lap(m)*( www(i1,i2,i3+m,1) + www(i1,i2,i3-m,1) )
       end do
    end do
    end do
    end do
    rewind 21
    rewind 22
    rewind 23
    m=Ngrid(1)/2
    do i3=a3-Md,b3+Md
    do i2=a2-Md,b2+Md
    do i1=a1-Md,b1+Md
       if ( i3==m .and. i2==m ) then
          write(21,'(1x,i6,f20.10,3g20.10)') i1,i1*Hgrid(1) &
               ,real(www(i1,i2,i3,1)),real(wk(i1,i2,i3)),wr(i1,i2,i3)
       end if
       if ( i3==m .and. i1==m ) then
          write(22,'(1x,i6,f20.10,3g20.10)') i2,i2*Hgrid(2) &
               ,real(www(i1,i2,i3,1)),real(wk(i1,i2,i3)),wr(i1,i2,i3)
       end if
       if ( i1==m .and. i2==m ) then
          write(23,'(1x,i6,f20.10,3g20.10)') i3,i3*Hgrid(3) &
               ,real(www(i1,i2,i3,1)),real(wk(i1,i2,i3)),wr(i1,i2,i3)
       end if
    end do
    end do
    end do
    write(*,*) "sum(wr),sum(wk)",sum(wr)*dV/pi4,sum(wk)*dV/pi4


    www(:,:,:,1) = v_outer(:,:,:)
    i=n1-1
    do i3=a3,b3
    do i2=a2,b2
    do i1=a1,b1
       i=i+1
       www(i1,i2,i3,1) = xk(i)
    end do
    end do
    end do
    qk(:)=0.d0
    i=n1-1
    do i3=a3,b3
    do i2=a2,b2
    do i1=a1,b1
       i=i+1
       do m=0,Md
          qk(i)=qk(i) &
               + c1*lap(m)*( www(i1+m,i2,i3,1) + www(i1-m,i2,i3,1) ) &
               + c2*lap(m)*( www(i1,i2+m,i3,1) + www(i1,i2-m,i3,1) ) &
               + c3*lap(m)*( www(i1,i2,i3+m,1) + www(i1,i2,i3-m,1) )
       end do
    end do
    end do
    end do

    rk(:) = b(:) - qk(:)
    pk(:) = rk(:)

    ee   = 0.d0
    gg_1 = 1.d10
    sum1 = 1.d10
    sum2 = 1.d10

    do iter=1,max_iter

       gg_0 = gg_1

       sum0 = sum( rk(:)*rk(:) )*dV
       call mpi_allreduce(sum0,gg_1,1,mpi_real8,mpi_sum,comm_grid,ierr)

       if ( sum1 < tol ) exit

       if ( iter > 1 ) then
          ck = gg_1/gg_0
          pk(:) = rk(:) + ck*pk(:)
       end if

       www(:,:,:,:)=z0
       i=n1-1
       do i3=a3,b3
       do i2=a2,b2
       do i1=a1,b1
          i=i+1
          www(i1,i2,i3,1)=pk(i)
       end do
       end do
       end do
       qk(:)=0.d0
       i=n1-1
       do i3=a3,b3
       do i2=a2,b2
       do i1=a1,b1
          i=i+1
          do m=0,Md
             qk(i) = qk(i) &
                  + c1*lap(m)*( www(i1+m,i2,i3,1) + www(i1-m,i2,i3,1) ) &
                  + c2*lap(m)*( www(i1,i2+m,i3,1) + www(i1,i2-m,i3,1) ) &
                  + c3*lap(m)*( www(i1,i2,i3+m,1) + www(i1,i2,i3-m,1) )
          end do
       end do
       end do
       end do

       sum0 = sum( pk(:)*qk(:) )*dV
       call mpi_allreduce(sum0,ak,1,mpi_real8,mpi_sum,comm_grid,ierr)
       ak=gg_1/ak

       v0(:) = xk(:)

       xk(:) = xk(:) + ak*pk(:)
       rk(:) = rk(:) - ak*qk(:)

       e0 = ee

       sb(1) = sum( rk(:)*rk(:) )*dV
       sb(2) = sum( xk(:)*b(:) )*dV *0.5d0/(-pi4)
       sb(3) = sum( (v0(:)-xk(:))**2 )/Ngrid(0)
       call mpi_allreduce(sb,rb,3,mpi_real8,mpi_sum,comm_grid,ierr) 
       sum1=sb(3)
       write(*,'(1x,i6,3g20.10)') iter,rb(1:3)

    end do ! iter

    www(:,:,:,1)=v_outer(:,:,:)
    i=n1-1
    do i3=a3,b3
    do i2=a2,b2
    do i1=a1,b1
       i=i+1
       www(i1,i2,i3,1)=xk(i)
    end do
    end do
    end do

    rewind 30
    rewind 31
    rewind 32
    do i3=a3-Md,b3+Md
    do i2=a2-Md,b2+Md
    do i1=a1-Md,b1+Md
       if ( i3==20 .and. i2==20 ) then
          write(30,'(1x,i6,f20.10,2g20.10)') i1,i1*Hgrid(1),www(i1,i2,i3,1)
       end if
       if ( i3==20 .and. i1==20 ) then
          write(31,'(1x,i6,f20.10,2g20.10)') i2,i2*Hgrid(2),www(i1,i2,i3,1)
       end if
       if ( i1==20 .and. i2==20 ) then
          write(32,'(1x,i6,f20.10,2g20.10)') i3,i3*Hgrid(3),www(i1,i2,i3,1)
       end if
    end do
    end do
    end do

    wk=0.d0
    do i3=a3,b3
    do i2=a2,b2
    do i1=a1,b1
       do m=0,Md
          wk(i1,i2,i3)=wk(i1,i2,i3) &
                  + c1*lap(m)*( www(i1+m,i2,i3,1) + www(i1-m,i2,i3,1) ) &
                  + c2*lap(m)*( www(i1,i2+m,i3,1) + www(i1,i2-m,i3,1) ) &
                  + c3*lap(m)*( www(i1,i2,i3+m,1) + www(i1,i2,i3-m,1) )
       end do
    end do
    end do
    end do

    rewind 11
    rewind 12
    rewind 13
    do i3=a3-Md,b3+Md
    do i2=a2-Md,b2+Md
    do i1=a1-Md,b1+Md
       if ( i3==20 .and. i2==20 ) then
          write(11,'(1x,i6,f20.10,2g20.10)') i1,i1*Hgrid(1),wk(i1,i2,i3)
       end if
       if ( i3==20 .and. i1==20 ) then
          write(12,'(1x,i6,f20.10,2g20.10)') i2,i2*Hgrid(2),wk(i1,i2,i3)
       end if
       if ( i1==20 .and. i2==20 ) then
          write(13,'(1x,i6,f20.10,2g20.10)') i3,i3*Hgrid(3),wk(i1,i2,i3)
       end if
    end do
    end do
    end do

    deallocate( v0 )
    deallocate( qk ) 
    deallocate( pk ) 
    deallocate( rk )
    deallocate( xk )
    deallocate( b )

    deallocate( wk )
    deallocate( lap )

  END SUBROUTINE esm_test3_test

END MODULE esm_cylindrical_test
