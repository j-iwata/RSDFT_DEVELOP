MODULE kinetic_esm_module

  use esm_rgrid_module
  use rgrid_module
  use bc_module

  implicit none

  PRIVATE
  PUBLIC :: get_coef_kinetic_esm,op_kinetic_esm

  real(8) :: coef_lap0
  real(8),allocatable :: coef_lap(:,:),coef_nab(:,:)
  integer :: Md
  real(8),allocatable :: const_k2(:),coef_nabk(:,:,:)
  complex(8),allocatable :: zcoef_kin(:,:,:)
  logical :: flag_nab

CONTAINS


  SUBROUTINE get_coef_kinetic_esm(aa,bb,MBZ,kbb,Md_in)
    use fd_module
    implicit none
    real(8),intent(IN) :: aa(3,3),bb(3,3)
    integer,intent(IN) :: MBZ,Md_in
    real(8),intent(IN) :: kbb(3,MBZ)
    real(8),allocatable :: lap(:),nab(:)
    integer :: m,k,n
    real(8) :: a1,a2,a3,c1,c2,c3,pi2,kx,ky,kz
    complex(8),parameter :: zi=(0.d0,1.d0)
    Md=Md_in
    allocate( coef_lap(3,Md) ) ; coef_lap=0.d0
    allocate( coef_nab(3,Md) ) ; coef_nab=0.d0
    allocate( lap(-Md:Md) ) ; lap=0.d0
    allocate( nab(-Md:Md) ) ; nab=0.d0
    call get_coef_laplacian_fd(Md,lap)
    call get_coef_laplacian_fd(Md,nab)
    coef_lap0 = -0.5d0*lap(0)*(1.d0/Hgrid(1)**2+1.d0/Hgrid(2)**2+1.d0/Hgrid(3)**2)
    do m=1,Md
       coef_lap(1,m) = -0.5d0*lap(m)/Hgrid(1)**2
       coef_lap(2,m) = -0.5d0*lap(m)/Hgrid(2)**2
       coef_lap(3,m) = -0.5d0*lap(m)/Hgrid(3)**2
    end do
    do m=1,Md
       coef_nab(1,m) = nab(m)/Hgrid(1)
       coef_nab(2,m) = nab(m)/Hgrid(2)
       coef_nab(3,m) = nab(m)/Hgrid(3)
    end do
    deallocate( nab )
    deallocate( lap )

    allocate( coef_nabk(3,Md,MBZ) ) ; coef_nabk=0.d0
    allocate( zcoef_kin(3,-Md:Md,MBZ) ) ; zcoef_kin=(0.d0,0.d0)
    allocate( const_k2(MBZ) ) ; const_k2=0.d0
    pi2=2.d0*acos(-1.d0)
    a1  = sqrt(sum(aa(1:3,1)**2))/pi2
    a2  = sqrt(sum(aa(1:3,2)**2))/pi2
    a3  = sqrt(sum(aa(1:3,3)**2))/pi2
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
          zcoef_kin(1:3,-n,k)= zi*coef_nabk(1:3,n,k)
          zcoef_kin(1:3, n,k)=-zi*coef_nabk(1:3,n,k)
       end do
    end do
  END SUBROUTINE get_coef_kinetic_esm


  SUBROUTINE op_kinetic_esm(k,n1,n2,ib1,ib2,tpsi,htpsi)
    integer,intent(IN) :: k,n1,n2,ib1,ib2
#ifdef _DRSDFT_
    real(8),intent(IN)  :: tpsi(n1:n2,ib1:ib2)
    real(8),intent(INOUT) :: htpsi(n1:n2,ib1:ib2)
    real(8),parameter :: zero=0.d0
#else
    complex(8),intent(IN)  :: tpsi(n1:n2,ib1:ib2)
    complex(8),intent(INOUT) :: htpsi(n1:n2,ib1:ib2)
    complex(8),parameter :: zero=(0.d0,0.d0)
#endif
    integer :: i,ib,i1,i2,i3,m,n
    real(8) :: c

    do ib=ib1,ib2
       do i=n1,n2
          i1=LL_ESM(1,i)+Nshift_ESM(1)
          i2=LL_ESM(2,i)+Nshift_ESM(2)
          i3=LL_ESM(3,i)+Nshift_ESM(3)
          www(i1,i2,i3,ib-ib1+1) = tpsi(i,ib)
       end do
    end do
    call bcset(1,ib2-ib1+1,Md,0)
    do ib=ib1,ib2
       do i=MK0_ESM,MK1_ESM
          i1=KK(1,i)+Nshift_ESM(1)
          i2=KK(2,i)+Nshift_ESM(2)
          i3=KK(3,i)+Nshift_ESM(3)
          www(i1,i2,i3,ib-ib1+1) = zero
       end do
    end do

    c = coef_lap0 + const_k2(k)

    do ib=ib1,ib2
       do i=n1,n2
          htpsi(i,ib) = htpsi(i,ib) + c*tpsi(i,ib)
       end do
    end do

    do ib=ib1,ib2
       n=ib-ib1+1
       do m=1,Md
          do i=n1,n2
             i1=LL_ESM(1,i)+Nshift_ESM(1)
             i2=LL_ESM(2,i)+Nshift_ESM(2)
             i3=LL_ESM(3,i)+Nshift_ESM(3)
             htpsi(i,ib)=htpsi(i,ib) &
                  +coef_lap(1,m)*( www(i1+m,i2,i3,n)+www(i1-m,i2,i3,n) ) &
                  +coef_lap(2,m)*( www(i1,i2+m,i3,n)+www(i1,i2-m,i3,n) ) &
                  +coef_lap(3,m)*( www(i1,i2,i3+m,n)+www(i1,i2,i3-m,n) )
          end do
       end do
    end do

    if ( flag_nab ) then

       do ib=ib1,ib2
          n=ib-ib1+1
          do m=1,Md
             do i=n1,n2
                i1=LL_ESM(1,i)+Nshift_ESM(1)
                i2=LL_ESM(2,i)+Nshift_ESM(2)
                i3=LL_ESM(3,i)+Nshift_ESM(3)
                htpsi(i,ib)=htpsi(i,ib) &
                     +zcoef_kin(1,m,k) *www(i1+m,i2,i3,ib-ib1+1) &
               +conjg(zcoef_kin(1,m,k))*www(i1-m,i2,i3,ib-ib1+1) &
                     +zcoef_kin(2,m,k) *www(i1,i2+m,i3,ib-ib1+1) &
               +conjg(zcoef_kin(2,m,k))*www(i1,i2-m,i3,ib-ib1+1) &
                     +zcoef_kin(3,m,k) *www(i1,i2,i3+m,ib-ib1+1) &
               +conjg(zcoef_kin(3,m,k))*www(i1,i2,i3-m,ib-ib1+1)   
             end do
          end do
       end do

    end if

  END SUBROUTINE op_kinetic_esm


END MODULE kinetic_esm_module
