module fermi_dirac_module

  implicit none

  private
  public :: fermi_dirac
  public :: entropy_fd

contains

  subroutine fermi_dirac( occ,Efermi,esp,ekbt,ZNel,wbz,dspn )
    implicit none
    real(8),intent(inout) :: occ(:,:,:)
    real(8),intent(out) :: Efermi
    real(8),intent(in) :: esp(:,:,:)
    real(8),intent(in) :: ekbt, ZNel
    real(8),intent(in) :: wbz(:)
    real(8),optional,intent(in) :: dspn
    real(8) :: ef1,ef2,ef
    integer :: ib,id,n,k,s,efconv
    real(8) :: ctime0,ctime1,etime0,etime1
    real(8) :: zne,octmp,ff,xx,ff0
    real(8),parameter :: eps=1.0d-10, big=36.9d0
    integer,parameter :: mxcycl=100
    integer,parameter :: kinteg=-1
    integer :: MB, MBZ, nspin, mb1, mb2
    logical :: flag, disp_switch
    real(8) :: diff

    call write_border( 1, "fermi_dirac(start)" )
    call check_disp_switch( disp_switch, 0 )

    MB = size( occ, 1 )
    MBZ = size( occ, 2 )
    nspin = size( occ, 3 )
    mb1 = 1
    mb2 = MB

    flag=.false.
    if ( present(dspn) ) then
       if ( nspin == 2 ) flag=.true.
    end if

    diff=abs(2*MB-ZNel)
    if ( diff < eps ) then
       occ(1:MB,1:MBZ,1:nspin)=1.0d0
       ef=maxval(esp)
       goto 100
    end if

    loop_mb1: do mb1 = 1, MB

! Set upper & lower boundarires of Fermi energy

       ef1=minval( esp(mb1:mb2,1:MBZ,1:nspin) )
       ef2=maxval( esp(mb1:mb2,1:MBZ,1:nspin) )
       if ( ef1==ef2 ) then
          ef1=ef1-0.01d0
          ef2=ef2+0.01d0
       end if

! Safety margin for highly degenerate systems & artificial fault

       ef2 = ef2 + min( ekbt*1.0d2, 0.1d0 )

       if ( .not.flag ) then
          zne=ZNel-2.0d0*(mb1-1)
          do id=1,mxcycl
             ef=0.5d0*(ef1+ef2)
             octmp=0.0d0
             do ib=mb1,mb2
             do s=1,nspin
             do k=1,MBZ
                select case (kinteg)
                !case default
                !   xx=(esp(ib,k,s)-ef)/ekbt
                !   ff=ff0(kinteg,xx)
                case (-1)
                   if ( ekbt<eps ) then
                      if ( (esp(ib,k,s)-ef )>0.0d0 ) then
                         ff=0.0d0
                      else if( ( esp(ib,k,s)-ef)<0.0d0 ) then
                         ff=1.0d0
                      else
                         ff=0.5d0
                      end if
                   else
                      xx=( esp(ib,k,s)-ef )/ekbt
                      if ( xx>big ) then
                         ff=0.0d0
                      else if ( xx<-big ) then
                         ff=1.0d0
                      else
                         ff=1.0d0/(exp(xx)+1.0d0)
                      end if
                   end if
                end select
                octmp=octmp+ff*wbz(k)*2.0d0/dble(nspin)
                occ(ib,k,s)=ff
             end do !k
             end do !s
             end do !ib
             if ( octmp-zne>eps ) then
                ef2=ef
             else if ( octmp-zne<-eps ) then
                ef1=ef
             else
                goto 100
             end if
          end do !id

       else !flag

          efconv=0
          if (disp_switch) write(*,*) "total spin density is fixed!!"
          do s=1,nspin
             ef1=minval( esp(mb1:mb2,1:MBZ,s) )
             ef2=maxval( esp(mb1:mb2,1:MBZ,s) )
             if ( ef1==ef2 ) then
                ef1=ef1-0.01d0
                ef2=ef2+0.01d0
             end if
             ef2 = ef2 + min( ekbt*1.d2, 0.1d0 )
             zne=(ZNel+dble(3-2*s)*dspn)/2.0d0-dble(mb1-1)
             do id=1,mxcycl
                ef=0.5d0*(ef1+ef2)
                octmp=0.0d0
                do ib=mb1,mb2
                do k=1,MBZ
                   select case(kinteg)
                   !case default
                   !   xx=(esp(ib,k,s)-ef)/ekbt
                   !   ff=ff0(kinteg,xx)
                   case(-1)
                      if ( ekbt<eps ) then
                         if ( (esp(ib,k,s)-ef )>0.0d0 ) then
                            ff=0.0d0
                         else if( ( esp(ib,k,s)-ef)<0.0d0 ) then
                            ff=1.0d0
                         else
                            ff=0.5d0
                         end if
                      else
                         xx=( esp(ib,k,s)-ef )/ekbt
                         if ( xx>big ) then
                            ff=0.0d0
                         else if ( xx<-big ) then
                            ff=1.0d0
                         else
                            ff=1.0d0/(exp(xx)+1.0d0)
                         end if
                      end if
                   end select
                   octmp=octmp+ff*wbz(k)
                   occ(ib,k,s)=ff
                end do !k
                end do !ib
                if ( octmp-zne>eps ) then
                   ef2=ef
                else if ( octmp-zne<-eps ) then
                   ef1=ef
                else
                   efconv=efconv+1
                   exit
                end if
             end do !id
          end do !s
          if ( efconv==2 ) goto 100

       end if !flag

       !if (DISP_SWITCH) then
       !   write(6,*)' EF IS NOT CONVERGED'
       !   write(6,*)' Check the # of electron, mb1, and mb2'
       !   write(6,*)' EF1 & EF2=',ef1,ef2
       !   write(6,*)' octmp,zne=',octmp,zne
       !   do s=1,nspin
       !   do k=1,MBZ
       !      write(*,*) "s,k =",s,k
       !      do ib=mb1,mb2
       !         write(*,*) ib,occ(ib,k,s),esp(ib,k,s)
       !      end do
       !   end do
       !   end do
       !end if

    end do loop_mb1

    call stop_program('FERMI(1)')

100 continue

    do s=1,nspin
    do k=1,MBZ
       do ib=1,mb1-1
          occ(ib,k,s)=2.0d0*wbz(k)/dble(nspin)
       end do
       do ib=mb1,mb2
          occ(ib,k,s)=2.0d0*occ(ib,k,s)*wbz(k)/dble(nspin)
       end do
       do ib=mb2+1,MB
          occ(ib,k,s)=0.0d0
       end do
    end do
    end do

    octmp=sum(occ(1:MB,1:MBZ,1:nspin))
    octmp=ZNel-octmp
!      if ( abs(octmp)/=0.d0 ) then
!         loop1 : do ib=mb2,mb1,-1
!            do s=1,nspin
!            do k=1,MBZ
!               if ( occ(ib,k,s)>0.d0 ) then
!                  n=ib
!                  exit loop1
!               end if
!            end do
!            end do
!         end do loop1
!         octmp=octmp/real(MBZ*nspin,8)
!         do s=1,nspin
!         do k=1,MBZ
!            occ(n,k,s)=occ(n,k,s)+octmp/dble(MBZ*nspin)
!         end do
!         end do
!      end if
!!SF120502
!!      if ( abs(octmp)/=0.d0 ) then
!!         do s=1,nspin
!!            loop1 : do ib=mb2,mb1,-1
!!            do k=1,MBZ
!!               if ( occ(ib,k,s)>0.d0 ) then
!!                  n=ib
!!                  exit loop1
!!               end if
!!            end do
!!            end do loop1
!!            do k=1,MBZ
!!               occ(n,k,s)=occ(n,k,s)+octmp/dble(MBZ*nspin)
!!            end do
!!         end do
!!      end if
    if ( abs(octmp)>0.0d0 ) then
       do s=1,nspin
          do k=1,MBZ
             loop1 : do ib=mb2,mb1,-1
                if ( occ(ib,k,s)>0.0d0 ) then
                   occ(ib,k,s)=occ(ib,k,s)+octmp/dble(MBZ*nspin)
                   exit loop1
                end if
             end do loop1
          end do
       end do
    end if
!!/SF120502
    !if ( nspin == 2 ) then
    !   octmp = sum(occ(1:MB,1:MBZ,1:nspin))
    !   do s=1,nspin
    !      Ntot(s) = sum(occ(1:MB,1:MBZ,s))*ZNel/octmp &
    !           & +ZNext/2.d0
    !   end do
    !end if

    Efermi=ef

    call write_border( 1, "fermi_dirac(end)" )

    return
  end subroutine fermi_dirac

  !FUNCTION ff0(n,x)
  !  implicit none
  !  integer :: n,i
  !  real(8) :: x,ff,sqrtPi,ff0
  !  INTERFACE
  !     FUNCTION bberf(x)
  !       real(8) :: bberf,x
  !     END FUNCTION bberf
  !     FUNCTION factorial(n)
  !       integer :: n
  !       real(8) :: factorial
  !     END FUNCTION factorial
  !     FUNCTION hermitepoly(n,x)
  !       integer :: n
  !       real(8) :: hermitepoly,x
  !     END FUNCTION hermitepoly
  !  END INTERFACE
  !  sqrtPi=sqrt(acos(-1.d0))
  !  ff0 = 0.5d0*(1.d0-bberf(x))
  !  if ( n<=0 ) return
  !  do i=1,n
  !     ff0 = ff0 +(-1.d0)**i/(factorial(i)*4.d0**i*sqrtPi) &
  !          & *HermitePoly(2*i-1,x) *exp(-x**2)
  !  end do
  !  return
  !END FUNCTION ff0

  !FUNCTION factorial(n)
  !  implicit none
  !  integer :: n,i,m
  !  real(8) :: factorial
  !  m=1
  !  do i=2,n
  !     m=m*i
  !  end do
  !  factorial=dble(m)
  !  return
  !END FUNCTION factorial

  !RECURSIVE FUNCTION HermitePoly(n,x) RESULT(HerPol)
  !  implicit none
  !  integer :: n
  !  real(8) :: x,HerPol
  !  if ( n <= 0 ) then
  !     HerPol = 1.d0
  !  else if ( n == 1 ) then
  !     HerPol = 2.d0*x
  !  else
  !     HerPol = 2.d0*x*HermitePoly(n-1,x) &
  !          & -2.d0*dble(n-1)*HermitePoly(n-2,x)
  !  end if
  !  return
  !END FUNCTION HermitePoly

  subroutine entropy_fd( Eentropy, ekbt, occ )
    implicit none
    real(8), intent(out) :: Eentropy
    real(8), intent(in)  :: ekbt
    real(8), intent(in)  :: occ(:,:,:)
    real(8) :: factor1,factor2,Sf,f,g
    integer :: n,k,s,nb,nk,ns
    nb=size(occ,1)
    nk=size(occ,2)
    ns=size(occ,3)
    factor1=ns/2.0d0  ! ns=1: factor1=0.5,  ns=2: factor1=1.0
    factor2=2.0d0/ns  ! ns=1: factor2=2.0,  ns=2: factor2=1.0
    Eentropy=0.0d0
    do s=1,ns
    do k=1,nk
    do n=1,nb
       Sf=0.0d0
       f=occ(n,k,s)*factor1
       if ( f > 0.0d0 ) Sf=Sf-f*log(f)
       g=1.0d0-f
       if ( g > 0.0d0 ) Sf=Sf-g*log(g)
       Eentropy=Eentropy+factor2*ekbt*Sf
    end do
    end do
    end do
  end subroutine entropy_fd

end module fermi_dirac_module
