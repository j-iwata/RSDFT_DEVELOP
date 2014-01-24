MODULE esm_genpot_module

  use aa_module
  use rgrid_module
  use esm_rgrid_module
  use ps_local_rs_module
  use modified_bessel_module
  use parallel_module
  use kinetic_module

  implicit none

  PRIVATE
  PUBLIC :: read_esm_genpot,esm_genpot,ikz_max,ima_max &
           ,Rleft,Rrigt,v_const,iswitch_bc,const_rrigt

  integer :: ikz_max,ima_max,iswitch_bc
  complex(8),allocatable :: v_outer(:) 
  real(8) :: const_rrigt, v_const
  integer,allocatable :: map1_grid(:,:),map2_grid(:,:)

CONTAINS


  SUBROUTINE read_esm_genpot(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: ierr,i
    character(6) :: cbuf,ckey
    ikz_max=0
    ima_max=0
    v_const=0.d0
    iswitch_bc=0
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:5) == "KMESM" ) then
             backspace(unit)
             read(unit,*) cbuf,ikz_max,ima_max
          else if ( ckey(1:6) == "VCONST" ) then
             backspace(unit)
             read(unit,*) cbuf,v_const
          else if ( ckey(1:4) == "SWBC"   ) then
             backspace(unit)
             read(unit,*) cbuf,iswitch_bc
          end if
       end do
999    continue
       write(*,*) "ikz_max,ima_max=",ikz_max,ima_max
       write(*,*) "iswitch_bc     =",iswitch_bc
       if ( iswitch_bc == 0 ) v_const=0.d0
       write(*,*) "v_const        =",v_const
    end if
    call mpi_bcast(ikz_max,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(ima_max,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(v_const,1,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(iswitch_bc,1,mpi_integer,0,mpi_comm_world,ierr)
  END SUBROUTINE read_esm_genpot


  SUBROUTINE esm_genpot(v_inout)
    implicit none
    real(8),intent(INOUT) :: v_inout(*)
    integer :: ikz,i,i1,i2,i3,j,j1,j2,ma,ima,ierr,n,s
    integer :: a1,a2,a3,b1,b2,b3
    real(8),parameter :: ep=1.d-10
    real(8) :: pi4,pi2,r,r1,r2,kz,phi,x,y,z,c,c0,pi,kr
    real(8) :: p1,p2,p3,p4,const,c1,c2,c3
    real(8),allocatable :: rhot(:),vtmp(:),vtmp0(:)
    complex(8),parameter :: zero=(0.d0,0.d0)
    complex(8) :: phase_k,phase_m
    complex(8) :: fmk_inner,ztmp

    pi  = acos(-1.d0)
    pi2 = 2.d0*acos(-1.d0)
    pi4 = 4.d0*acos(-1.d0)

    allocate( v_outer(MK0_ESM:MK1_ESM) ) ; v_outer=zero
    allocate( rhot(ML0_ESM:ML1_ESM)    ) ; rhot=0.d0

    rhot(:)=rho_ps(:)

!--- neutralize ---
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

!    write(*,*) "maxval(v_vouter)",maxval(abs(v_outer)**2)

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
    rhot(:) = rho_ps(:)

    if ( myrank == 0 ) write(*,*) "------- esm_test3"

    call esm_test3(ML0_ESM,ML1_ESM,rhot,v_inout)

    deallocate( rhot )
    deallocate( v_outer )

  END SUBROUTINE esm_genpot

  FUNCTION Rleft(n,x)
    integer :: n
    real(8) :: x,Rleft
    Rleft = bessi(n,x)
  END FUNCTION Rleft

  FUNCTION Rrigt(n,x)
    integer :: n
    real(8) :: x,Rrigt
    select case(iswitch_bc)
    case default
       Rrigt = - bessk(n,x)
    case(1)
       Rrigt = const_rrigt*bessi(n,x) - bessk(n,x)
    end select
  END FUNCTION Rrigt


  SUBROUTINE esm_test3(n1,n2,rho_in,v_inout)
    use kinetic_module
    use bc_module
    use fd_module
    use pseudopot_module
    implicit none
    integer,intent(IN)  :: n1,n2
    real(8),intent(IN)  :: rho_in(n1:n2)
    real(8),intent(INOUT) :: v_inout(n1:n2)
    integer,parameter :: max_iter=2000
    integer :: a1,a2,a3,b1,b2,b3,i,i1,i2,i3,j1,j2,j3,j
    integer :: m,mm,nn,nnn,ierr,iter,mp,m1,m2,m3
    integer :: u1,u2,u3,ik
    integer,allocatable :: grdmap(:,:)
    real(8),parameter :: tol=1.d-24
    real(8) :: pi4,r1,c1,c2,c3,sum0,sum1,ee,e0,sb(3),rb(3)
    real(8) :: ak,ck,sum2,gg_0,gg_1
    real(8) :: r,x,y,z,p1,p2,p3,p4
    real(8),allocatable :: rk(:),qk(:),v0(:),xk(:)
    real(8),allocatable :: pk(:),b(:),lap(:)
    complex(8),allocatable :: wk(:,:,:)
    complex(8),parameter :: z0=(0.d0,0.d0)
    logical,allocatable :: mat(:,:)

    m1=Nshift_ESM(1)
    m2=Nshift_ESM(2)
    m3=Nshift_ESM(3)
    a1=Igrid(1,1)+m1
    b1=Igrid(2,1)+m1
    a2=Igrid(1,2)+m2
    b2=Igrid(2,2)+m2
    a3=Igrid(1,3)+m3
    b3=Igrid(2,3)+m3

    allocate( wk(a1-Md:b1+Md,a2-Md:b2+Md,a3-Md:b3+Md) )
    wk=(0.d0,0.d0)

    mp = Nshift_ESM(1)

    pi4 = 4.d0*acos(-1.d0)
    c1  = 1.d0/Hgrid(1)**2
    c2  = 1.d0/Hgrid(2)**2
    c3  = 1.d0/Hgrid(3)**2

    allocate( lap(-Md:Md) ) ; lap=0.d0
    call get_coef_laplacian_fd(Md,lap)
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

    goto 33
    wk(:,:,:) = www(:,:,:,1)
    ik=1
    p1=parloc(1,ik) ; p2=sqrt(parloc(2,ik))
    p3=parloc(3,ik) ; p4=sqrt(parloc(4,ik))
    www=z0
    do j=MK0_ESM,MK1_ESM
       i1=KK(1,j) !+m1
       i2=KK(2,j) !+m2
       i3=KK(3,j) !+m3
       x=i1*Hgrid(1)
       y=i2*Hgrid(2)
       z=i3*Hgrid(3)
       r=sqrt(x*x+y*y+z*z)
       www(i1+m1,i2+m2,i3+m3,1)=-Zps(ik)/r*( p1*erf(p2*r)+p3*erf(p4*r) )
    end do
    do i3=Igrid(1,3)-Md,Igrid(1,3)-1
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
       x=i1*Hgrid(1)
       y=i2*Hgrid(2)
       z=i3*Hgrid(3)
       r=sqrt(x*x+y*y+z*z)
       www(i1+m1,i2+m2,i3+m3,1)=-Zps(ik)/r*( p1*erf(p2*r)+p3*erf(p4*r) )
    end do
    end do
    end do
    do i3=Igrid(2,3)+1,Igrid(2,3)+Md
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
       x=i1*Hgrid(1)
       y=i2*Hgrid(2)
       z=i3*Hgrid(3)
       r=sqrt(x*x+y*y+z*z)
       www(i1+m1,i2+m2,i3+m3,1)=-Zps(ik)/r*( p1*erf(p2*r)+p3*erf(p4*r) )
    end do
    end do
    end do
33  continue
    

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
    if ( myrank == 0 ) then
    u1=21
    u2=22
    u3=23
    rewind u1
    rewind u2
    rewind u3
    do i3=a3-Md,b3+Md
    do i2=a2-Md,b2+Md
    do i1=a1-Md,b1+Md
       j1=i1-m1
       j2=i2-m2
       j3=i3-m3
       if ( i2==m2 .and. i3==m3 ) then
          write(u1,'(1x,i6,f20.10,2g20.10)') &
               j1,j1*Hgrid(1),real(www(i1,i2,i3,1)),real(wk(i1,i2,i3))
!          write(u1,'(1x,i6,f20.10,2g20.10)') &
!          i1,i1*Hgrid(1),real(www(i1,i2,i3,1)),real(wk(i1,i2,i3))
!          write(u1,*) i1-m1,(i1-m1)*Hgrid(1),real(www(i1,i2,i3,1))
       end if
       if ( i1==m1 .and. i3==m3 ) then
          write(u2,'(1x,i6,f20.10,2g20.10)') &
               j2,j2*Hgrid(2),real(www(i1,i2,i3,1)),real(wk(i1,i2,i3))
!          write(u2,*) i2-m2,(i2-m2)*Hgrid(2) &
!          ,real(www(i1,i2,i3,1)),real(wk(i1,i2,i3))
       end if
       if ( i2==m2 .and. i1==m1 ) then
          write(u3,'(1x,i6,f20.10,2g20.10)') &
               j3,j3*Hgrid(3),real(www(i1,i2,i3,1)),real(wk(i1,i2,i3))
!          write(u3,*) i3-m3,(i3-m3)*Hgrid(3) &
!          ,real(www(i1,i2,i3,1)),real(wk(i1,i2,i3))
       end if
    end do
    end do
    end do
    end if
    call end_mpi_parallel ; stop
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

       if ( sum1 < tol ) exit

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
       call flush(6)
       write(*,'(1x,i6,3g20.10,i4)') iter,rb(1:3),myrank
       call flush(6)

       v_inout(:) = xk(:)

    end do ! iter

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

    goto 22
    ik=1
    p1=parloc(1,ik) ; p2=sqrt(parloc(2,ik))
    p3=parloc(3,ik) ; p4=sqrt(parloc(4,ik))
    www=z0
    do j=MK0_ESM,MK1_ESM
       i1=KK(1,j) !+m1
       i2=KK(2,j) !+m2
       i3=KK(3,j) !+m3
       x=i1*Hgrid(1)
       y=i2*Hgrid(2)
       z=i3*Hgrid(3)
       r=sqrt(x*x+y*y+z*z)
       www(i1+m1,i2+m2,i3+m3,1)=-Zps(ik)/r*( p1*erf(p2*r)+p3*erf(p4*r) )
    end do
    do i3=Igrid(1,3)-Md,Igrid(1,3)-1
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
       x=i1*Hgrid(1)
       y=i2*Hgrid(2)
       z=i3*Hgrid(3)
       r=sqrt(x*x+y*y+z*z)
       www(i1+m1,i2+m2,i3+m3,1)=-Zps(ik)/r*( p1*erf(p2*r)+p3*erf(p4*r) )
    end do
    end do
    end do
    do i3=Igrid(2,3)+1,Igrid(2,3)+Md
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
       x=i1*Hgrid(1)
       y=i2*Hgrid(2)
       z=i3*Hgrid(3)
       r=sqrt(x*x+y*y+z*z)
       www(i1+m1,i2+m2,i3+m3,1)=-Zps(ik)/r*( p1*erf(p2*r)+p3*erf(p4*r) )
    end do
    end do
    end do
    do j=ML0_ESM,ML1_ESM
       i1=LL(1,j)+m1
       i2=LL(2,j)+m2
       i3=LL(3,j)+m3
       www(i1,i2,i3,1) = xk(j)
    end do
22  continue

    goto 1
    if ( myrank == 0 ) then
    u1=111
    u2=112
    u3=113
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
    call end_mpi_parallel ; stop
1   continue

    goto 2
    wk=z0
    do i=ML0_ESM,ML1_ESM
       i1=LL(1,i)+m1
       i2=LL(2,i)+m2
       i3=LL(3,i)+m3
       do m=0,Md
          wk(i1,i2,i3) = wk(i1,i2,i3) &
               + c1*lap(m)*( www(i1+m,i2,i3,1) + www(i1-m,i2,i3,1) ) &
               + c2*lap(m)*( www(i1,i2+m,i3,1) + www(i1,i2-m,i3,1) ) &
               + c3*lap(m)*( www(i1,i2,i3+m,1) + www(i1,i2,i3-m,1) )
       end do
    end do
    www=z0
    do i=ML0_ESM,ML1_ESM
       i1=LL(1,i)+m1
       i2=LL(2,i)+m2
       i3=LL(3,i)+m3
       www(i1,i2,i3,1) = b(i)
    end do
    u1=10*myrank+14
    rewind u1
    do i3=a3-Md,b3+Md
    do i2=a2-Md,b2+Md
    do i1=a1-Md,b1+Md
       if ( i3 == mp .and. i2 == mp ) then
          write(u1,'(1x,i6,f20.10,2g20.10)') &
               i1-m1,(i1-m1)*Hgrid(1),real(wk(i1,i2,i3)),real(www(i1,i2,i3,1))
       end if
    end do
    end do
    end do
2   continue

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

END MODULE esm_genpot_module
