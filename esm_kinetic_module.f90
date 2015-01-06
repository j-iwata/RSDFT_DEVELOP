MODULE esm_kinetic_module

  use esm_rgrid_module
  use rgrid_module
  use bc_module
  use kinetic_variables, only: Md, coef_lap0, coef_lap &
                              ,const_k2, zcoef_kin, flag_nab

  implicit none

  PRIVATE
  PUBLIC :: op_esm_kinetic

CONTAINS

  SUBROUTINE op_esm_kinetic(k,n1,n2,ib1,ib2,tpsi,htpsi)
    implicit none
    integer,intent(IN) :: k,n1,n2,ib1,ib2
#ifdef _DRSDFT_
    real(8),intent(IN)    ::  tpsi(n1:n2,ib1:ib2)
    real(8),intent(INOUT) :: htpsi(n1:n2,ib1:ib2)
    real(8),parameter :: zero=0.d0
#else
    complex(8),intent(IN)    ::  tpsi(n1:n2,ib1:ib2)
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

  END SUBROUTINE op_esm_kinetic

END MODULE esm_kinetic_module
