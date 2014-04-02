MODULE ps_read_YB_module

  implicit none

  PRIVATE
  PUBLIC :: ps_read_YB, ps_yb

  integer,parameter :: nrmax = 5000
  integer,parameter :: lmax  = 3

  TYPE yb
     integer :: nrr
     integer :: norb
     integer :: NRc(lmax+1)
     integer :: inorm(lmax+1)
     integer :: lo(lmax+1)
     real(8) :: znuc
     real(8) :: vps(nrmax,lmax+1)
     real(8) :: ups(nrmax,lmax+1)
     real(8) :: vql(nrmax)
     real(8) :: cdc(nrmax)
     real(8) :: cdd(nrmax)
     real(8) :: rr(nrmax)
     real(8) :: rx(nrmax)
     real(8) :: Rc(lmax+1)
     real(8) :: anorm(lmax+1)
  END TYPE yb

  type(yb) :: ps_yb

CONTAINS

  SUBROUTINE ps_read_YB(g)
    implicit none
    integer,intent(IN) :: g
    real(8),parameter :: E2=14.39965d0, H2M=3.80998d0
    real(8),parameter :: a_B=0.529177249d0
    real(8),parameter :: Ry=13.6056981d0, Pi=3.141592653589793d0
    integer :: i,j,L,Mlps0,Lref,nrr,norb
    real(8) :: rPC,dr,c,znuc_tmp
    real(8),allocatable :: tmp(:,:)

    ps_yb%nrr     =0
    ps_yb%norb    =0
    ps_yb%znuc    =0.d0
    ps_yb%vps(:,:)=0.d0
    ps_yb%ups(:,:)=0.d0
    ps_yb%vql(:)  =0.d0
    ps_yb%cdc(:)  =0.d0
    ps_yb%cdd(:)  =0.d0
    ps_yb%rr(:)   =0.d0
    ps_yb%rx(:)   =0.d0
    ps_yb%Rc(:)   =0.d0
    ps_yb%NRc(:)  =0
    ps_yb%anorm(:)=0.d0
    ps_yb%inorm(:)=0
    ps_yb%lo(:)   =0

! Read

     read(g,*) nrr, dr, Mlps0, ps_yb%znuc

     nrr = nrr + 1

     if ( nrr > nrmax .or. Mlps0 > lmax ) then
        write(*,*) "nrr0,nrmax,Mlps0,lmax",nrr,nrmax,Mlps0,lmax
        stop "Array size is small (stop at KY_format)"
     end if

     read(g,*) rPC, ( ps_yb%Rc(L), L=1,Mlps0+1 )

     do i=1,nrr
        read(g,*) ps_yb%rr(i),ps_yb%cdc(i),( ps_yb%vps(i,L), L=1,Mlps0+1 )
     end do
     do i=1,nrr
        read(g,*) ps_yb%rr(i),( ps_yb%ups(i,L), L=1,Mlps0+1 )
     end do
!
     i=0
     do L=0,Mlps0-1
        i=i+1
        ps_yb%lo(i)=L
     end do
     norb=i
     Lref=Mlps0

     do j=1,norb
        do i=1,nrr
           if ( ps_yb%rr(i) >= ps_yb%Rc(j) ) then
              ps_yb%NRc(j)=i
              exit
           end if
        end do
        ps_yb%Rc(j) = ps_yb%rr(ps_yb%NRc(j))
     end do

     ps_yb%rx(1:nrr) = dr

     ps_yb%vql(1:nrr) = ps_yb%vps(1:nrr,Mlps0+1)

     allocate( tmp(nrr,norb) ) ; tmp=0.0d0
     tmp(:,:)=ps_yb%vps(1:nrr,1:norb)
     do i=1,norb
        ps_yb%vps(1:nrr,i)=(tmp(1:nrr,i)-ps_yb%vql(1:nrr))*ps_yb%ups(1:nrr,i)
        ps_yb%anorm(i)=sum( ps_yb%vps(1:nrr,i)*ps_yb%ups(1:nrr,i) )*dr
        if ( ps_yb%anorm(i) < 0.0d0 ) then
           ps_yb%inorm(i)=-1
           ps_yb%anorm(i)=abs( ps_yb%anorm(i) )
        else
           ps_yb%inorm(i)=1
        end if
        ps_yb%vps(1:nrr,i)=ps_yb%vps(1:nrr,i)/sqrt( ps_yb%anorm(i) )
     end do
     deallocate( tmp )

     ps_yb%rr(:)    = ps_yb%rr(:)/a_B
     ps_yb%rx(:)    = ps_yb%rx(:)/a_B
     ps_yb%Rc(:)    = ps_yb%Rc(:)/a_B
     ps_yb%vql(:)   = ps_yb%vql(:)/(2.d0*Ry)
     ps_yb%vps(:,:) = ps_yb%vps(:,:)/sqrt(2.d0*Ry/a_B)
     ps_yb%anorm(:) = ps_yb%anorm(:)/(2.d0*Ry)

     ps_yb%cdd(:)=0.d0
     znuc_tmp=0.d0
     do j=1,norb
        c=2.d0*(2*ps_yb%lo(j)+1)
        znuc_tmp=znuc_tmp+c
        write(*,*) j,ps_yb%lo(j),znuc_tmp,c
        if ( znuc_tmp <= ps_yb%znuc ) then
           ps_yb%cdd(1:nrr)=ps_yb%cdd(1:nrr)+c*ps_yb%ups(1:nrr,j)**2
        end if
     end do

     ps_yb%nrr  = nrr
     ps_yb%norb = norb

     write(*,*) "*** KY format ***"
     write(*,*) "Znuc=",ps_yb%znuc
     write(*,*) "# of radial mesh points =",nrr
     write(*,*) "# of orbitals =",norb
     write(*,*) "angular momentum =",ps_yb%lo(1:norb)
     write(*,*) "reference L =",Lref
     write(*,*) "cut off radius =",ps_yb%Rc(1:norb)
     write(*,*) "uVu integral (anorm) ="
     write(*,*) "sum(rhov)=",sum(ps_yb%cdd*ps_yb%rx)
     write(*,'(1x,8f10.5)') ( ps_yb%inorm(i)*ps_yb%anorm(i),i=1,norb )

     return
   END SUBROUTINE ps_read_YB

 END MODULE ps_read_YB_module
