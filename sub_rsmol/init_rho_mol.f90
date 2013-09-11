!--------1---------2---------3---------4---------5---------6---------7--

      SUBROUTINE init_rho_mol
      use global_variables
      implicit none

      integer :: a,ik,MMr,i,ir,m,m1,m2,NR0(10),iloc(1)
      integer :: n1,n2,ierr,ispin,i1,i2,i3
      logical :: flag_alloc
      real(8),parameter :: ep=1.d-9,eps=1.d-10
      real(8) :: x,y,z,r,v0,v,err0,err,errmax,sum_rho(2,2),pi4,Rc2,fac,r2
      real(8) :: ctime0,etime0,ctime1,etime1,mem,memax
      real(8),allocatable :: rho_tmp(:)

      call watch(ctime0,etime0)

      if ( DISP_SWITCH) then
         write(*,'(a60," init_rho_mol")') repeat("-",60) 
      end if

      if ( MKI>10 ) stop "init_rho"

      n1  = idisp(myrank)+1
      n2  = idisp(myrank)+ircnt(myrank)
      pi4 = 16.d0*atan(1.d0)
      rho(:,:)=0.d0
      flag_alloc=.false.
      mem=0.d0
      memax=0.d0

      if ( .not.allocated(LL) ) then

         flag_alloc=.true.

!- allocate -----------------------------------------------------------------
         allocate( LL(3,n1:n2) ) ; LL=0 ; mem=bsintg*(n2-n1+1)*3 ; memax=mem
!----------------------------------------------------------------------------

         Rc2=Rsize**2
         i=n1-1
         do i3=a3b,b3b
         do i2=a2b,b2b
         do i1=a1b,b1b
            select case(SYStype)
            case(1)
               r2=H1*H1*(i1*i1+i2*i2+i3*i3)
               if ( r2<Rc2+eps ) then
                  i=i+1
                  LL(1,i)=i1
                  LL(2,i)=i2
                  LL(3,i)=i3
               end if
            case(2)
               r2=H1*H1*(i1*i1+i2*i2)
               z=abs(i3*H1)
               if ( r2<Rc2+eps .and. z<Zsize+eps ) then
                  i=i+1
                  LL(1,i)=i1
                  LL(2,i)=i2
                  LL(3,i)=i3
               end if
            end select
         end do
         end do
         end do

      end if

      do ik=1,MKI
         do ir=2,Mr(ik)
            cdd(ir,ik)=cdd(ir,ik)/(pi4*rad(ir,ik)**2)
         end do
         if ( rad(1,ik)==0.d0 ) then
            cdd(1,ik)=cdd(2,ik)
         else
            cdd(1,ik)=cdd(1,ik)/(pi4*rad(1,ik)**2)
         end if
         iloc=maxloc(cdd(:,ik))
         do ir=iloc(1),Mr(ik)
            if ( cdd(ir,ik)<ep ) then
               NR0(ik)=ir
               exit
            end if
         end do
      end do

!- allocate ------------------------------------------------------------------------
      allocate( rho_tmp(n1:n2) ) ; mem=mem+bdreal*(n2-n1+1) ; memax=max(mem,memax)
!-----------------------------------------------------------------------------------

      errmax=0.d0
      do a=1,MI
         ik=Kion(a)
         MMr=NR0(ik)
         rho_tmp(:)=0.d0
         do i=n1,n2
            x=LL(1,i)*H1-asi(1,a)
            y=LL(2,i)*H1-asi(2,a)
            z=LL(3,i)*H1-asi(3,a)
            r=sqrt(x*x+y*y+z*z)
            if ( r<=rad(MMr,ik) ) then
               do ir=1,MMr
                  if ( rad(ir,ik)>=r ) exit
               end do
               ir=min(ir,MMr)
               err0=1.d10
               do m=2,2
                  m1=max(1,ir-m)
                  m2=min(ir+m,MMr)
                  call polint(rad(m1,ik),cdd(m1,ik),m2-m1+1,r,v,err)
!                  if ( abs(err)<err0 ) then
                   v0=v ; err0=abs(err)
!                  if ( err0<ep ) exit
!                  end if
               end do
            else
               v0=0.d0
               err0=0.d0
            end if
            rho_tmp(i)=rho_tmp(i)+v0
            errmax=max(errmax,err0)
         end do
         do ispin=1,nspin
            fac=(Zps(ik)-(-1)**ispin*dspin(a))/Zps(ik)/dble(nspin)
            rho(n1:n2,ispin) = rho(n1:n2,ispin) + rho_tmp(n1:n2)*fac
         end do
      end do

!- deallocate ---------------------------------------------
      mem=mem-bdreal*size(rho_tmp) ; deallocate( rho_tmp )
!----------------------------------------------------------

!- deallocate ---------------------------------------------
      if ( flag_alloc ) then
         mem=mem-bsintg*size(LL)
         deallocate( LL )
      end if
!----------------------------------------------------------

      if (DISP_SWITCH) write(*,*) "errmax(polint)=",errmax

      if ( isymmetry==1 ) then
         do ispin=1,nspin
            sum_rho(ispin,1)=sum(weight(n1:n2)*rho(n1:n2,ispin))*dV
         end do
      else
         do ispin=1,nspin
            sum_rho(ispin,1)=sum(rho(n1:n2,ispin))*dV
         end do
      end if

      call mpi_allreduce(sum_rho(1,1),sum_rho(1,2),nspin,mpi_real8,mpi_sum,comm_grid,ierr)

      if (DISP_SWITCH) then
         do ispin=1,nspin
            write(*,*) "sum_rho, Ntot =",sum_rho(ispin,2),Ntot(ispin)
         end do
      end if

      do ispin=1,nspin
         fac=Ntot(ispin)/sum_rho(ispin,2)
         rho(n1:n2,ispin)=fac*rho(n1:n2,ispin)
      end do

      call watch(ctime1,etime1)
      if(DISP_SWITCH)then
         write(*,*) "TIME(INIT_RHO_MOL)=",ctime1-ctime0,etime1-etime0
         write(*,*) "MEM(MB)",memax*B2MB,mem
      end if

      return
      END SUBROUTINE init_rho_mol
