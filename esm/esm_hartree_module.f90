MODULE esm_hartree_module

  use aa_module
  use rgrid_module
  use esm_rgrid_module
  use ps_local_rs_module
  use modified_bessel_module
  use parallel_module
  use esm_genpot_module
  use bc_module
  use fd_module

  implicit none

  PRIVATE
  PUBLIC :: calc_esm_hartree, init_esm_hartree

  complex(8),allocatable :: v_outer(:)
  integer :: Md
  logical :: flag_init=.true.

CONTAINS


  SUBROUTINE init_esm_hartree(Md_in)
    implicit none
    integer,intent(IN) :: Md_in
    Md = Md_in
    flag_init = .false.
  END SUBROUTINE init_esm_hartree


  SUBROUTINE calc_esm_hartree(n1,n2,nsp,rho_in,v_inout,e_out)
    implicit none
    integer,intent(IN)    :: n1,n2,nsp
    real(8),intent(IN)    :: rho_in(n1:n2,nsp)
    real(8),intent(INOUT) :: v_inout(n1:n2)
    real(8),intent(OUT)   :: e_out
    integer :: ikz,i,i1,i2,i3,j,j1,j2,ma,ima,ierr,n,s
    integer :: a1,a2,a3,b1,b2,b3
    real(8),parameter :: ep=1.d-10
    real(8) :: pi4,pi2,r,r1,r2,kz,phi,x,y,z,c,c0,pi,kr
    real(8) :: p1,p2,p3,p4,const,c1,c2,c3
    real(8),allocatable :: rhot(:),vtmp(:),vtmp0(:)
    complex(8),parameter :: zero=(0.d0,0.d0)
    complex(8) :: phase_k,phase_m
    complex(8) :: fmk_inner,ztmp

    if ( flag_init ) then
       write(*,*) "init_hartree should be called first"
       stop "stop@calc_esm_hartree"
    end if

    e_out = 0.d0

    pi  = acos(-1.d0)
    pi2 = 2.d0*acos(-1.d0)
    pi4 = 4.d0*acos(-1.d0)

    allocate( v_outer(MK0_ESM:MK1_ESM) ) ; v_outer=zero
    allocate( rhot(ML0_ESM:ML1_ESM)    ) ; rhot=0.d0

    do s=1,nsp
       rhot(:) = rhot(:) + rho_in(:,s)
    end do

!--- neutralize
!
    c0 = sum(rhot)*dV
    call mpi_allreduce(c0,c,1,mpi_real8,mpi_sum,comm_grid,ierr)
    do i=ML0_ESM,ML1_ESM
       if ( LL_ESM(1,i) == 0 .and. LL_ESM(2,i) == 0 ) then
          rhot(i) = rhot(i) - c/(dV*Ngrid(3))
       end if
    end do

    const = -2.d0*dV/aa(3,3)

    c1 = 1.d0/Ngrid(1)
    c2 = 1.d0/Ngrid(2)
    c3 = 1.d0/Ngrid(3)

    fmk_inner = zero

    do ikz=-ikz_max,ikz_max
       if ( ikz == 0 ) cycle
       kz=pi2/aa(3,3)*ikz
    do ima=-ima_max,ima_max
       ma=abs(ima)

       x=abs(kz)*Rsize2
       const_rrigt = bessk(ma,x)/bessi(ma,x)

       fmk_inner=zero
       do i=ML0_ESM,ML1_ESM
          i1=LL(1,i)
          i2=LL(2,i)
          i3=LL(3,i)
          x = aa(1,1)*i1*c1 + aa(1,2)*i2*c2 + aa(1,3)*i3*c3
          y = aa(2,1)*i1*c1 + aa(2,2)*i2*c2 + aa(2,3)*i3*c3
          z = aa(3,1)*i1*c1 + aa(3,2)*i2*c2 + aa(3,3)*i3*c3
          r = sqrt(x*x+y*y)
          phase_k = dcmplx( cos(kz*z),-sin(kz*z) )
          if ( r == 0.d0 ) then
             phi = 0.d0
          else
             phi = sign( acos(x/r),y )
          end if
          phase_m = dcmplx( cos(ima*phi),-sin(ima*phi) )
          kr=abs(kz)*r
          fmk_inner = fmk_inner + phase_k*phase_m*bessi(ma,kr)*rhot(i)
       end do
       fmk_inner=fmk_inner*const

       ztmp=fmk_inner
       call mpi_allreduce(ztmp,fmk_inner,1,mpi_complex16,mpi_sum,comm_grid,ierr)
!       write(*,'(1x,i4,2i6,2g20.10)') myrank,ikz,ima,fmk_inner
!       call flush(6)

       do i=MK0_ESM,MK1_ESM
          i1=KK(1,i)
          i2=KK(2,i)
          i3=KK(3,i)
          x=aa(1,1)*i1*c1+aa(1,2)*i2*c2+aa(1,3)*i3*c3
          y=aa(2,1)*i1*c1+aa(2,2)*i2*c2+aa(2,3)*i3*c3
          z=aa(3,1)*i1*c1+aa(3,2)*i2*c2+aa(3,3)*i3*c3
          r=sqrt(x*x+y*y)
          if ( iswitch_bc == 1 .and. Rsize2 <= r + ep ) cycle
          phase_k = dcmplx( cos(kz*z),sin(kz*z) )
          if ( r == 0.d0 ) then
             phi = 0.d0
          else
             phi = sign( acos(x/r),y )
          end if
          phase_m = dcmplx( cos(ima*phi),sin(ima*phi) )
          kr=abs(kz)*r
          v_outer(i) = v_outer(i) &
               +phase_k*phase_m*Rrigt(ma,kr)*fmk_inner
       end do

    end do ! ima
    end do ! ikz

    v_outer(:) = v_outer(:) + v_const

    do i=MK0_ESM,MK1_ESM
       i1=KK(1,i)
       i2=KK(2,i)
       i3=KK(3,i)
       x=aa(1,1)*i1*c1+aa(1,2)*i2*c2+aa(1,3)*i3*c3
       y=aa(2,1)*i1*c1+aa(2,2)*i2*c2+aa(2,3)*i3*c3
       z=aa(3,1)*i1*c1+aa(3,2)*i2*c2+aa(3,3)*i3*c3
       r=sqrt(x*x+y*y)
       if ( iswitch_bc == 1 ) then
          if ( Rsize2 <= r + ep ) cycle
          v_outer(i) = v_outer(i) - 2.d0*c/aa(3,3)*log(r/Rsize2)
       else
          v_outer(i) = v_outer(i) - 2.d0*c/aa(3,3)*log(r)
       end if
    end do
    rhot(:)=0.d0
    do s=1,nsp
       rhot(:) = rhot(:) + rho_in(:,s)
    end do

    if ( myrank == 0 ) write(*,*) "------- esm_test3"

    call esm_test3(ML0_ESM,ML1_ESM,rhot,v_inout,e_out)

    deallocate( rhot )
    deallocate( v_outer )

  END SUBROUTINE calc_esm_hartree


  SUBROUTINE esm_test3(n1,n2,rho_in,v_inout,e_out)
    implicit none
    integer,intent(IN)  :: n1,n2
    real(8),intent(IN)  :: rho_in(n1:n2)
    real(8),intent(INOUT) :: v_inout(n1:n2)
    real(8),intent(OUT) :: e_out
    integer,parameter :: max_iter=2000
    integer :: a1,a2,a3,b1,b2,b3,i,i1,i2,i3,j1,j
    integer :: m,mm,nn,nnn,ierr,iter,mp,m1,m2,m3
    integer :: u1,u2,u3
    integer,allocatable :: grdmap(:,:)
    real(8),parameter :: tol=1.d-24
    real(8) :: pi4,r1,c1,c2,c3,sum0,sum1,ee,e0,sb(3),rb(3)
    real(8) :: ak,ck,sum2,gg_0,gg_1
    real(8),allocatable :: rk(:),qk(:),v0(:),xk(:)
    real(8),allocatable :: pk(:),b(:),lap(:)
    complex(8),allocatable :: wk(:,:,:)
    complex(8),parameter :: z0=(0.d0,0.d0)
    logical,allocatable :: mat(:,:)

    e_out = 0.d0

    m1=Nshift_ESM(1)
    m2=Nshift_ESM(2)
    m3=Nshift_ESM(3)
    a1=Igrid(1,1)+m1
    b1=Igrid(2,1)+m1
    a2=Igrid(1,2)+m2
    b2=Igrid(2,2)+m2
    a3=Igrid(1,3)+m3
    b3=Igrid(2,3)+m3

    allocate( wk(-m1-Md:m1-1+Md,-m2-Md:m2-1+Md,-m3-Md:m3-1+Md) )
    wk=(0.d0,0.d0)

    mp = Ngrid(1)/2

    pi4 = 4.d0*acos(-1.d0)
    c1  = 1.d0/Hgrid(1)**2
    c2  = 1.d0/Hgrid(2)**2
    c3  = 1.d0/Hgrid(3)**2

    allocate( lap(-Md:Md) ) ; lap=0.0d0
    call get_coef_lapla_fd(Md,lap)
    lap(0)=0.5d0*lap(0)

    allocate(  b(n1:n2) )
    allocate( xk(n1:n2) )
    allocate( rk(n1:n2) )
    allocate( pk(n1:n2) )
    allocate( qk(n1:n2) )
    allocate( v0(n1:n2) )

    do j=n1,n2
       b(j) = -pi4*rho_in(j)
    end do

    do j=n1,n2
       xk(j) = v_inout(j)
    end do

    www(:,:,:,:)=z0
    do j=n1,n2
       i1=LL(1,j)+m1
       i2=LL(2,j)+m2
       i3=LL(3,j)+m3
       www(i1,i2,i3,1) = xk(j)
    end do
    call bcset(1,1,Md,0)
    do j=MK0_ESM,MK1_ESM
       i1=KK(1,j)+m1
       i2=KK(2,j)+m2
       i3=KK(3,j)+m3
       www(i1,i2,i3,1)=v_outer(j)
    end do

!    allocate( mat(a1-Md:b1+Md,a2-Md:b2+Md) ) ; mat=.false.
!    i3=-1
!    do i2=a2-Md,b2+Md
!    do i1=a1-Md,b1+Md
!       if ( www(i1,i2,i3,1) /=z0 ) mat(i1,i2)=.true.
!    end do
!    end do
!    rewind 100+myrank
!    do i2=a2-Md,b2+Md
!       write(100+myrank,'(1x,52l1)') ( mat(i1,i2),i1=a1-Md,b1+Md )
!    end do
!    deallocate( mat )

    goto 20
    u1=10*myrank+11
    u2=10*myrank+12
    u3=10*myrank+13
    rewind u1
    rewind u2
    rewind u3
    do i3=a3-Md,b3+Md
    do i2=a2-Md,b2+Md
    do i1=a1-Md,b1+Md
       if ( i2==m2 .and. i3==m3 ) then
!          write(u1,'(1x,i6,f20.10,2g20.10)') i1,i1*Hgrid(1),www(i1,i2,i3,1)
          write(u1,*) i1,i1*Hgrid(1),real(www(i1,i2,i3,1))
       end if
       if ( i1==m1 .and. i3==m3 ) then
!          write(u2,'(1x,i6,f20.10,2g20.10)') i2,i2*Hgrid(2),www(i1,i2,i3,1)
          write(u2,*) i2,i2*Hgrid(2),real(www(i1,i2,i3,1))
       end if
       if ( i2==m2 .and. i1==m1 ) then
!          write(u3,'(1x,i6,f20.10,2g20.10)') i3,i3*Hgrid(3),www(i1,i2,i3,1)
          write(u3,*) i3,i3*Hgrid(3),real(www(i1,i2,i3,1))
       end if
    end do
    end do
    end do
20 continue

    qk(:)=0.d0
    do j=ML0_ESM,ML1_ESM
       i1=LL(1,j)+m1
       i2=LL(2,j)+m2
       i3=LL(3,j)+m3
       do m=0,Md
          qk(j) = qk(j) &
               + c1*lap(m)*( www(i1+m,i2,i3,1) + www(i1-m,i2,i3,1) ) &
               + c2*lap(m)*( www(i1,i2+m,i3,1) + www(i1,i2-m,i3,1) ) &
               + c3*lap(m)*( www(i1,i2,i3+m,1) + www(i1,i2,i3-m,1) )
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

       if ( sum1 < tol ) then
          exit
          if ( myrank==0 ) write(*,'(1x,i6,3g20.10,i4)') iter,rb(1:3),myrank
       end if

       if ( iter > 1 ) then
          ck = gg_1/gg_0
          pk(:) = rk(:) + ck*pk(:)
       end if

       www(:,:,:,:)=z0
       do j=ML0_ESM,ML1_ESM
          i1=LL(1,j)+m1
          i2=LL(2,j)+m2
          i3=LL(3,j)+m3
          www(i1,i2,i3,1)=pk(j)
       end do
       call bcset(1,1,Md,0)
       do j=MK0_ESM,MK1_ESM
          i1=KK(1,j)+m1
          i2=KK(2,j)+m2
          i3=KK(3,j)+m3
          www(i1,i2,i3,1)=z0
       end do

       qk(:)=0.d0
       do j=ML0_ESM,ML1_ESM
          i1=LL(1,j)+m1
          i2=LL(2,j)+m2
          i3=LL(3,j)+m3
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
       e0 = ee

       xk(:) = xk(:) + ak*pk(:)
       rk(:) = rk(:) - ak*qk(:)

       sb(1) = sum( rk(:)*rk(:) )*dV
       sb(2) = sum( xk(:)*b(:) )*dV *0.5d0/(-pi4)
       sb(3) = sum( (v0(:)-xk(:))**2 )/Ngrid(0)
       call mpi_allreduce(sb,rb,3,mpi_real8,mpi_sum,comm_grid,ierr) 
       sum1=rb(3)
!       call flush(6)
!       write(*,'(1x,i6,3g20.10,i4)') iter,rb(1:3),myrank
!       call flush(6)

       e_out = rb(2)
       v_inout(:) = xk(:)

    end do ! iter


    goto 1
    www=z0
    do j=ML0_ESM,ML1_ESM
       i1=LL(1,j)+m1
       i2=LL(2,j)+m2
       i3=LL(3,j)+m3
       www(i1,i2,i3,1) = xk(j)
    end do
    call bcset(1,1,Md,0)
    do j=MK0_ESM,MK1_ESM
       i1=KK(1,j)+m1
       i2=KK(2,j)+m2
       i3=KK(3,j)+m3
       www(i1,i2,i3,1) = v_outer(j)
    end do
    if ( myrank == 0 ) then
    u1=211
    u2=212
    u3=213
    rewind u1
    rewind u2
    rewind u3
    do i3=a3-Md,b3+Md
    do i2=a2-Md,b2+Md
    do i1=a1-Md,b1+Md
       if ( i3 == m3 .and. i2 == m2 ) then
          write(u1,'(1x,i6,f20.10,2g20.10)') &
               i1-m1,(i1-m1)*Hgrid(1),www(i1,i2,i3,1)
       end if
       if ( i1 == m1 .and. i3 == m3 ) then
          write(u2,'(1x,i6,f20.10,2g20.10)') &
               i2-m2,(i2-m2)*Hgrid(2),www(i1,i2,i3,1)
       end if
       if ( i1 == m1 .and. i2 == m2 ) then
          write(u3,'(1x,i6,f20.10,2g20.10)') &
               i3-m3,(i3-m3)*Hgrid(3),www(i1,i2,i3,1)
       end if
    end do
    end do
    end do
    end if
!    call end_mpi_parallel ; stop
1 continue
    deallocate( v0 )
    deallocate( qk ) 
    deallocate( pk ) 
    deallocate( rk )
    deallocate( xk )
    deallocate( b )

    deallocate( wk )
    deallocate( lap )

    www=z0

  END SUBROUTINE esm_test3

END MODULE esm_hartree_module
