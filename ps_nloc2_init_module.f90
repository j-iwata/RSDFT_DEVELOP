MODULE ps_nloc2_init_module

!  use pseudopot_module
  use VarPSMember

  use Filtering, only: opFiltering

  implicit none

  PRIVATE
  PUBLIC :: rad1,dviod,read_ps_nloc2_init,read_oldformat_ps_nloc2_init &
           ,ps_nloc2_init,ps_nloc2_init_derivative &
           ,rcfac,qcfac,etafac

  real(8),allocatable :: dviod(:,:,:)
  real(8) :: rcfac,qcfac,etafac

CONTAINS


  SUBROUTINE read_ps_nloc2_init(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i
    character(3) :: cbuf,ckey
    rcfac =1.5d0
    qcfac =1.0d0
    etafac=8.0d0
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:3) == "NL2" ) then
             backspace(unit)
             read(unit,*) cbuf,rcfac,qcfac,etafac
          end if
       end do
999    continue
       write(*,*) "rcfac =",rcfac
       write(*,*) "qcfac =",qcfac
       write(*,*) "etafac=",etafac
    end if
    call send_ps_nloc2_init(0)
  END SUBROUTINE read_ps_nloc2_init


  SUBROUTINE read_oldformat_ps_nloc2_init(rank,unit)
    integer,intent(IN) :: rank,unit
    if ( rank == 0 ) then
       read(unit,*) rcfac,qcfac,etafac
       write(*,*) "rcfac, qcfac =",rcfac,qcfac
       write(*,*) "etafac       =",etafac
    end if
    call send_ps_nloc2_init(0)
  END SUBROUTINE read_oldformat_ps_nloc2_init


  SUBROUTINE send_ps_nloc2_init(rank)
    implicit none
    integer,intent(IN) :: rank
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(rcfac  ,1,mpi_real8,rank,mpi_comm_world,ierr)
    call mpi_bcast(qcfac  ,1,mpi_real8,rank,mpi_comm_world,ierr)
    call mpi_bcast(etafac ,1,mpi_real8,rank,mpi_comm_world,ierr)
  END SUBROUTINE send_ps_nloc2_init


  SUBROUTINE ps_nloc2_init(qcut)
    use atom_module, only: Natom,ki_atom
    use maskf_module
use parallel_module
    implicit none
    real(8),intent(IN) :: qcut
    integer :: i,j,ik,iorb,L,m,m0,m1,m2,MMr,NRc,iloc(1)
    real(8),parameter :: dr=2.d-3
    real(8) :: qc,Rc,sum0,const
    real(8) :: x,y,y0,dy,dy0,maxerr
    real(8) :: r,r1,sb0x,sb0y,sb1x,sb1y
    real(8),allocatable :: vrad(:),tmp(:),wm(:,:,:),vtmp(:,:,:)

if (myrank==0) write(400+myrank,*) ">>>>> inside ps_nloc2_init"
    qc = qcut*qcfac
    if ( qc<=0.d0 ) qc=qcut
if (myrank==0) write(400+myrank,*) "before allocateRps"
    call allocateRps
if (myrank==0) write(400+myrank,*) "after allocateRps"

    do ik=1,Nelement_local
       MMr=Mr(ik)
       do iorb=1,norb(ik)
          Rc=Rps(iorb,ik)*rcfac
          iloc=minloc( abs(rad(1:MMr,ik)-Rc) )
          NRc=iloc(1) ; if ( rad(NRc,ik)<Rc ) NRc=NRc+1
          if ( NRc>MMr ) then
             write(*,*) "NRc,MMr=",NRc,MMr
             stop "rcfac is too large."
          end if
          NRps(iorb,ik)=NRc
          Rps(iorb,ik)=rad(NRc,ik)
       end do
    end do

if (myrank==0) write(400+myrank,*) "ps_nloc2_init 1"
    NRc=max_psgrd
    m=max_psorb
    allocate( wm(NRc,m,Nelement_local) )
if (myrank==0) write(400+myrank,*) "ps_nloc2_init 2"

    do ik=1,Nelement_local
       do iorb=1,norb(ik)
          NRc=NRps(iorb,ik)
          Rc=Rps(iorb,ik)
          call makemaskf(etafac)
          maxerr=0.d0
          do i=1,NRc
             x=rad(i,ik)/Rc
             if ( x<=dxm ) then
                y0=1.d0 ; dy0=0.d0
             else
                m0=int(x/dxm)
                dy0=1.d10
                do m=1,20
                   m1=max(m0-m,1) ; m2=min(m0+m,nmsk)
                   call polint(xm(m1),maskr(m1),m2-m1+1,x,y,dy)
                   if ( abs(dy)<dy0 ) then
                      y0=y ; dy0=abs(dy)
                   end if
                end do
             end if
             wm(i,iorb,ik)=y0
             maxerr=max(maxerr,dy0)
          end do
       end do
    end do

if (myrank==0) write(400+myrank,*) "ps_nloc2_init 3"
    do ik=1,Nelement_local
       do iorb=1,norb(ik)
          NRps(iorb,ik)=Rps(iorb,ik)/dr+1
          if ( (NRps(iorb,ik)-1)*dr < Rps(iorb,ik) ) then
             NRps(iorb,ik)=NRps(iorb,ik)+1
          end if
       end do
    end do
    MMr=max( maxval(Mr),maxval(NRps) )

if (myrank==0) write(400+myrank,*) "ps_nloc2_init 4"
    if ( MMr>maxval(Mr) ) then
       m0=size(viod,1)
       m1=size(viod,2)
       m2=size(viod,3)
       allocate( vtmp(m0,m1,m2) )
       vtmp=viod
       deallocate( viod )
       allocate( viod(MMr,m1,m2) ) ; viod=0.d0
       viod(1:m0,1:m1,1:m2) = vtmp(1:m0,1:m1,1:m2)
       deallocate( vtmp )
    end if

!    allocate( rad1(MMr,Nelement_local) ) ; rad1=0.d0
    do ik=1,Nelement_local
       do i=1,MMr
          rad1(i,ik)=(i-1)*dr
       end do
    end do

    NRc=maxval(NRps0)
    allocate( vrad(NRc),tmp(NRc) )

if (myrank==0) write(400+myrank,*) "ps_nloc2_init 5"
    const=2.d0/acos(-1.d0)

    do ik=1,Nelement_local
       do iorb=1,norb(ik)
          L=lo(iorb,ik)
          NRc=NRps0(iorb,ik)
          vrad(1:NRc)=rad(1:NRc,ik)*viod(1:NRc,iorb,ik) &
                     *rab(1:NRc,ik)/wm(1:NRc,iorb,ik)
          
          call opFiltering( qc,L,NRc,NRps(iorb,ik),rad(1,ik),rad1(1,ik),vrad,viod(1,iorb,ik) )

       end do ! iorb
    end do ! ik
    deallocate( vrad,tmp )

if (myrank==0) write(400+myrank,*) "ps_nloc2_init 6"
    do ik=1,Nelement_local
       do iorb=1,norb(ik)
          L=lo(iorb,ik)
          NRc=NRps(iorb,ik)
          Rc=Rps(iorb,ik)
          call makemaskf(etafac)
          maxerr=0.d0
          do i=1,NRc
             x=(i-1)*dr/Rc
             if ( x<=dxm ) then
                y0=1.d0 ; dy0=0.d0
             else
                m0=int(x/dxm)
                dy0=1.d10
                do m=1,20
                   m1=max(m0-m,1) ; m2=min(m0+m,nmsk)
                   call polint(xm(m1),maskr(m1),m2-m1+1,x,y,dy)
                   if ( abs(dy)<dy0 ) then
                      y0=y ; dy0=abs(dy)
                   end if
                end do
             end if
             if ( maxerr<dy0 ) maxerr=dy0
             viod(i,iorb,ik)=y0*viod(i,iorb,ik)
          end do
       end do
    end do

    deallocate( wm )
if (myrank==0) write(400+myrank,*) "ps_nloc2_init 7"

    return
  END SUBROUTINE ps_nloc2_init


  SUBROUTINE ps_nloc2_init_derivative
    implicit none
    integer :: ik,L,NRc,J,iorb,i,m,m1,m2,lm
    real(8) :: maxerr,y,dy,y0,dy0
    real(8) :: pi4,const
    real(8),allocatable :: dvrad(:,:,:)

    pi4 = 4.d0*acos(-1.d0)

    lm=0
    do ik=1,Nelement_local
       m=0
       do iorb=1,norb(ik)
          if ( lo(iorb,ik)==0 ) then
             m=m+1
          else
             m=m+3
          end if
       end do
       lm=max(m,lm)
    end do
    NRc=maxval(NRps)

    allocate( dviod(NRc,lm,Nelement_local) )
    dviod=0.d0

    do ik=1,Nelement_local
       do iorb=1,norb(ik)
          L=lo(iorb,ik)
          NRc=NRps(iorb,ik)
          maxerr=0.d0
          do i=1,NRc
             dy0=1.d10
             do m=1,20
                m1=max(i-m,1) ; m2=min(i+m,NRc)
                call dpolint( rad1(m1,ik),viod(m1,iorb,ik),m2-m1+1 &
                     ,rad1(i,ik),y,dy )
                if ( abs(dy)<dy0 ) then
                   y0=y ; dy0=abs(dy)
                end if
             end do
             dviod(i,iorb,ik)=y0
             maxerr=max(maxerr,dy0)
          end do
       end do
    end do

    NRc=maxval(NRps)
    allocate( dvrad(NRc,lm,Nelement_local) ) ; dvrad=0.d0

    do ik=1,Nelement_local
       lm=0
       do iorb=1,norb(ik)
          L=lo(iorb,ik)
          NRc=NRps(iorb,ik)
          do J=abs(L-1),L+1
             lm=lm+1
             const=0.5d0*(2.d0+L*(L+1)-J*(J+1))
             do i=1,NRc
                dvrad(i,lm,ik)=rad1(i,ik)**2*dviod(i,iorb,ik) &
                     +const*rad1(i,ik)*viod(i,iorb,ik)
             end do
          end do
       end do
    end do
    const=sqrt(pi4/3.d0)
    do ik=1,Nelement_local
       lm=0
       do iorb=1,norb(ik)
          L=lo(iorb,ik)
          NRc=NRps(iorb,ik)
          do J=abs(L-1),L+1
             lm=lm+1
             do i=1,NRc
                dviod(i,lm,ik)=const*dvrad(i,lm,ik)
             end do
          end do
       end do
    end do

    deallocate( dvrad )

  END SUBROUTINE ps_nloc2_init_derivative

END MODULE ps_nloc2_init_module
