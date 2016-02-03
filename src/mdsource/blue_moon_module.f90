MODULE blue_moon_module

  use parallel_module
  use cpmd_variables, only: ia,ib,ic,cval,DT2BYM,DTB2MI &
                           ,xlagr,ylagr,Rion0,Rion,anorm,fc,fv,Velocity &
                           ,ipvt,nodim,mcnstr,ityp,dt,amu,pmass &
                           ,index,cnpar,itime
  use atom_module, only: Natom, ki_atom, zn_atom
  use aa_module, only: aa

  implicit none

  PRIVATE
  PUBLIC :: cpmdshake,rattle,write_blue_data,read_blue,gettau,solvs &
           ,dealloc_blue,diffd,difft,funcd,funct,fillc,cnstfc,solve &
           ,bcast_blue_data,coornum,pbcdis

CONTAINS

!==================================================================
  SUBROUTINE cpmdshake
!==--------------------------------------------------------------==
!==  GENERAL CONSTRAINTS ON POSITIONS FOR VELOCITY VERLET        ==
!==--------------------------------------------------------------==
    implicit none
    real(8) :: velp(3,Natom)
!   Variables
    integer,parameter :: maxrat=5000, mrdiis=5
    real(8),parameter :: tolf=1.d-7, tolx=1.d-8
    real(8) :: vtest
    real(8) :: tau0(3,Natom),taup(3,Natom)
    real(8) :: v(mrdiis+1),diism(mrdiis+1,mrdiis+1),cnmax,fact,errf,errx,xval
    real(8) :: dasum,ddot,pm
    real(8) :: tscr(3,Natom),dx(nodim),asl(mcnstr,mcnstr)
    real(8) :: xlo(mcnstr,mrdiis),err(mcnstr,mrdiis)
    real(8) :: an0(nodim,mcnstr)
    integer :: i,j,is,iter,k,info,idiis,m,i_dim
    integer,save :: istart1
    data      istart1 /0/
!   SAVE      ISTART1,IP_TSCR,IP_DX,IP_ASL,IP_XLO,IP_ERR,IP_IPVT,IP_AN0
!   LOGICAL   OLDSTATUS
!==--------------------------------------------------------------==
!   change constraint value according to growth rate
!    do i=1,mcnstr
!       cnsval(i)=cnsval(i)+grate(i)*dt
!       if(cnsval_dest(i).ne.-999.d0) then
!          if(grate(i).gt.0.d0.and.cnsval(i).gt.cnsval_dest(i)) then
!             cnsval(i)=cnsval_dest(i) ! increase
!          else if(grate(i).lt.0.d0.and.cnsval(i).lt.cnsval_dest(i)) then
!             cnsval(i)=cnsval_dest(i) ! decrease
!          endif
!       endif
!    enddo

!GLOBAL: ylagr(alloc),xlagr(alloc),mcnstr,ipvt
!GLOBAL: ia(alloc),ib(alloc),cval(alloc),istart

!----------------------------------------------allocate and set !parameters
    if ( istart1 == 0 ) then
       allocate( DT2BYM(nodim) )
       allocate( DTB2MI(nodim) )
!       allocate( pm_dim(nodim) )
       do i=1,Natom
          pm = pmass( zn_atom(ki_atom(i)) )*amu
!          pm_dim(3*i-2)=pmass(iatom(i))*amu
!          pm_dim(3*i-1)=pmass(iatom(i))*amu
!          pm_dim(3*i)  =pmass(iatom(i))*amu
!          DT2BYM(3*i-2)=(dt*dt)/pm_dim(3*i-2)
!          DT2BYM(3*i-1)=(dt*dt)/pm_dim(3*i-1)
!          DT2BYM(3*i)  =(dt*dt)/pm_dim(3*i)
!          DTB2MI(3*i-2)=dt/(2.0d0*pm_dim(3*i-2))
!          DTB2MI(3*i-1)=dt/(2.0d0*pm_dim(3*i-1))
!          DTB2MI(3*i)  =dt/(2.0d0*pm_dim(3*i))
          DT2BYM(3*i-2)=(dt*dt)/pm
          DT2BYM(3*i-1)=(dt*dt)/pm
          DT2BYM(3*i)  =(dt*dt)/pm
          DTB2MI(3*i-2)=dt/(2.0d0*pm)
          DTB2MI(3*i-1)=dt/(2.0d0*pm)
          DTB2MI(3*i)  =dt/(2.0d0*pm)
       end do
    end if
    call dcopy(mcnstr,ylagr(1),1,xlagr(1),1) !point: YLAGR comes from Rattle
    if ( istart1 == 0 ) then
       xlagr(:)=0.0d0
       istart1=1
    endif
!--------------------------------------------------------------------------

    tau0(:,:)=Rion0(:,:)
    taup(:,:)=Rion(:,:)

    dx(:)=0.0d0
    call cnstfc(tau0,dx) !point: TAU0 is position unrenewed

    call dcopy(mcnstr*nodim,anorm(1,1),1,an0(1,1),1)
    xval=mcnstr
! First guess for the force
    dx(:)=0.0d0
    do i=1,mcnstr
       do j=1,nodim
          dx(j)=dx(j)+xlagr(i)*an0(j,i)
       end do
    end do
! Update the positions
    tscr(:,:)=0.0d0
    call gettau(tscr,dx)
    do is=1,Natom
       fact=-dt*DTB2MI(IS)
       tau0(1,is)=taup(1,is)+fact*tscr(1,is)
       tau0(2,is)=taup(2,is)+fact*tscr(2,is)
       tau0(3,is)=taup(3,is)+fact*tscr(3,is)
    end do

! Iterativ calculation of lambda

    do iter=1,maxrat

!      Calculate constraint function and forces
!      for the current value of lambda

       dx(:)=0.0d0
       call cnstfc(tau0,dx)

       call dscal(mcnstr,-1.0d0,fc(1),1)
       errf=dasum(mcnstr,fc(1),1)

       if ( myrank == 0 ) write(*,*) "ERRF--->", errf, iter
!---------------------------------------------------------------------
       if ( errf < tolf ) then
          Rion(:,:)=tau0(:,:)
          goto 100
       end if
!---------------------------------------------------------------------
!      PERFORM DIIS CALCULATION

!      Derivatives of sigma wrt lambda
       do i=1,mcnstr
          do j=1,mcnstr
             asl(i,j)=0.d0
             do k=1,nodim
                asl(i,j)=asl(i,j)-dt*DT2BYM(k)*anorm(k,i)*an0(k,j)
             end do
          end do
       end do

!      Solve for asl*cg=fc
!      LAPACK matrix solver
       if ( nint(xval) == mcnstr ) then
          info=0
          call DGESV(mcnstr,1,asl,mcnstr,ipvt,fc,mcnstr,info)
          if ( info /= 0 ) stop " ipmdshake|','error in dgesv"
       else
          call solvs(mcnstr,asl,fc)
       end if
       errx=dasum(mcnstr,fc(1),1)
!      DIIS!
       idiis=mod(iter-1,mrdiis)+1
       call dcopy(mcnstr,xlagr(1),1,xlo(1,idiis),1)
       call dcopy(mcnstr,fc(1),1,err(1,idiis),1)
       if ( iter > mrdiis ) then
          m=mrdiis+1
          diism(:,:)=0.0d0
          v(:)=0.0d0
          do i=1,mrdiis
             do j=1,mrdiis
                diism(i,j)=ddot(mcnstr,err(1,i),1,err(1,j),1)
             end do
             diism(m,i)=1.d0
             diism(i,m)=1.d0
          end do
          v(m)=1.d0
          call solve(diism,m,m,v)
          fc(:)=0.0d0
          xlagr(:)=0.0d0
          do i=1,mrdiis
             do j=1,mcnstr
                fc(j)=fc(j)+v(i)*err(j,i)
                xlagr(j)=xlagr(j)+v(i)*xlo(j,i)
             end do
          end do
       end if
       call DAXPY(mcnstr,1.0d0,fc(1),1,xlagr(1),1)
       if ( errx < tolx ) goto 100
!      Update forces
       dx(:)=0.0d0
       do i=1,mcnstr
          do j=1,nodim
             dx(j)=dx(j)+xlagr(i)*an0(j,i)
          end do
       end do
!      Update the positions
       tscr(:,:)=0.0d0
       call gettau(tscr,dx)
       do is=1,Natom
          FACT=-dt*DTB2MI(IS)
          tau0(1,is)=taup(1,is)+fact*tscr(1,is)
          tau0(2,is)=taup(2,is)+fact*tscr(2,is)
          tau0(3,is)=taup(3,is)+fact*tscr(3,is)
       end do
       Rion(:,:)=tau0(:,:)
    end do
    if ( myrank == 0 ) write(*,'(30a)') 'cpmdshake| did not converge!!!'
    stop

!   MODIFY THE VELOCITIIES

100 continue

    call dcopy(3*Natom,tau0(1,1),1,taup(1,1),1)
!   Update velocities
    dx(:)=0.0d0
    do i=1,mcnstr
       do j=1,nodim
          dx(j)=dx(j)+xlagr(i)*an0(j,i)
       end do
    end do
    tscr(:,:)=0.0d0
    call gettau(tscr,dx)
    velp(:,:)=Velocity(:,:)
!   FIXME: add OpenMP?? AK
    do is=1,Natom
       fact=-DTB2MI(is)
       velp(1,is)=velp(1,is)+fact*tscr(1,is)
       velp(2,is)=velp(2,is)+fact*tscr(2,is)
       velp(3,is)=velp(3,is)+fact*tscr(3,is)
    end do

    Velocity(:,:)=velp(:,:)
!==--------------------------------------------------------------==
    return
  END SUBROUTINE cpmdshake
!

!==================================================================
  SUBROUTINE rattle
!==--------------------------------------------------------------==
!==  GENERAL CONSTRAINTS ON VELOCITIES FOR VELOCITY VERLET       ==
!==--------------------------------------------------------------==
    implicit none
    real(8) :: tau0(3,Natom),velp(3,Natom)
    real(8) :: tols,tolx
    parameter (tols=1.d-12,tolx=1.d-8)
    real(8),allocatable :: aux(:) 
    real(8) :: errx,cnmax,xval
    real(8) :: dasum
    real(8) :: dx(nodim),dvv(nodim),asl(mcnstr,mcnstr)
    external  dasum
    integer :: naux,i,j,k,info,nrank
    logical :: oldstatus
!==--------------------------------------------------------------==

    if ( mcnstr == 0 ) return

! moved out of the if to keep in initialized (Joost) ...

    naux=6*nodim+mcnstr
    allocate( aux(naux) )

    velp(:,:)=Velocity(:,:)
    tau0(:,:)=Rion(:,:) 
      
    dx(:)=0.0d0
    call cnstfc(tau0,dx)
    call puttau(velp,dx)
    asl(:,:)=0.0d0
    do i=1,mcnstr
       dvv(i)=0.d0
       do j=1,nodim
          dvv(i)=dvv(i)+anorm(j,i)*dx(j)
       end do
       do j=1,mcnstr
          do k=1,nodim
             asl(i,j)=asl(i,j)+anorm(k,i)*anorm(k,j)*DT2BYM(k)/dt/2.d0
          end do
       end do
    end do

    call DGELSS &
         (mcnstr,mcnstr,1,asl,mcnstr,dvv,mcnstr,fc,1.d-12,nrank,aux,naux,info)
    if ( info /= 0 ) stop "rattle','dgelss ended with info.ne.0"

    do i=1,mcnstr
       ylagr(i)=dvv(i)
    end do
    do j=1,nodim
       dvv(j)=0.d0
       do i=1,mcnstr
          dvv(j)=dvv(j)-ylagr(i)*anorm(j,i)*DT2BYM(j)/dt/2.0d0
       end do
    end do
    do j=1,nodim
       dx(j)=dx(j)+dvv(j)
    end do
    call gettau(velp,dx)

    Velocity(:,:)=velp(:,:)

!   Check accuracy
    do i=1,mcnstr
       fc(i)=0.d0
       do j=1,nodim
          fc(i)=fc(i)+anorm(j,i)*dx(j)
       end do
    end do
    errx=dasum(mcnstr,fc(1),1)

    if ( errx > tolx ) stop ' rattle| CONSTRAINT CAN NOT BE FULFILLED'

    deallocate( aux )
!==--------------------------------------------------------------==
    return
  END SUBROUTINE rattle
!==================================================================


!==================================================================
  SUBROUTINE gettau(tau1,xpar)
    implicit none
    real(8) ::  tau1(3,Natom),xpar(nodim)
    integer :: j,k,i
!==--------------------------------------------------------------==
    k=0
    do i=1,Natom
       do j=1,3
          k=k+1
          tau1(j,i) = xpar(k)
       end do
    end do
!==--------------------------------------------------------------==
    return
  END SUBROUTINE gettau

!==================================================================
  SUBROUTINE puttau(tau1,xpar)
    implicit none
    real(8) ::  tau1(3,Natom),xpar(nodim)
    integer :: j,k,i
!==--------------------------------------------------------------==
    k=0
    do i=1,Natom
       do j=1,3
          k=k+1
          xpar(k) = tau1(j,i)
       end do
    end do
!==--------------------------------------------------------------==
    return
  END SUBROUTINE puttau


!==================================================================
  SUBROUTINE solvs(m,asl,fc_dum)
    implicit none
    integer :: m
    real(8) :: asl(m,m),fc_dum(m)
    integer,parameter :: mx=100
    integer il(mx),i,n,info
!==--------------------------------------------------------------==
    if ( m > mx ) stop "solvs','mx'"
    n=0
    do i=1,m
       if ( abs(asl(i,i)) > 1.d-12 ) then
          n=n+1
          il(i)=n
          if ( i /= n ) then
             call dcopy(m,asl(1,i),1,asl(1,n),1)
             call dcopy(m,asl(i,1),m,asl(n,1),m)
             fc_dum(n)=fc_dum(i)
          end if
       else
          il(i)=0
       end if
    end do
    info=0
    call DGESV(n,1,asl,m,ipvt,fc_dum,m,info)
    if ( info /= 0 ) stop "' solvs|','error in dgesv'"
    call dcopy(n,fc_dum,1,asl,1)
    do i=1,m
       if ( il(i) /= 0 ) then
          fc_dum(i)=asl(il(i),1)
       else
          fc_dum(i)=0.d0
       end if
    end do
!==--------------------------------------------------------------==
    return
  END SUBROUTINE solvs
!==================================================================
!==================================================================
  SUBROUTINE solve(b,ldb,ndim,v)
!==--------------------------------------------------------------==
    implicit none
!   Arguments
    integer :: ldb,ndim, info
    real(8) :: b(:,:),v(:)
!   Variables
    integer :: nr,lw
    integer,parameter :: maxdis=20
    integer,parameter ::mrdiis=5
    real(8) :: diism(mrdiis+1,mrdiis+1)
    real(8) :: toleig
    real(8) :: scr1(maxdis+1,maxdis+1),scr2(maxdis+1,maxdis+1)
!==--------------------------------------------------------------==
    if ( ndim > maxdis+1 ) then
       write(*,*) 'solve! ndim=',ndim,' maxdis+1=',maxdis+1
       write(*,*) 'SOLVE','MDIM GREATER THAN MAXDIS+1'
       stop
    end if
    toleig=0.2d-15
    lw=(maxdis+1)*(maxdis+1)
    call dgelss(ndim,ndim,1,b,ldb,v,ldb,scr1,toleig,nr,scr2,lw,info)
    if ( info /= 0 ) then
       write(*,*) 'solve! info=',info
       write(*,*) 'SOLVE','COULD NOT SOLVE DIIS EQUATION'
       stop
    end if
    return
  END SUBROUTINE solve
!------------------------------------------------------------------------
!------------------------------------------------------------------------

  SUBROUTINE cnstfc(tau0,dxpar)
    implicit none
    real(8) :: tau0(3,Natom)
    real(8) :: dxpar(nodim)
    real(8) :: x1(3),x2(3),x3(3),dx(9)
    real(8) :: fcstr(nodim)
    real(8) :: c_kappa,c_rc
    real(8),allocatable :: askel(:,:,:)
    integer :: kmax,k,n,i,j,iaa,ibb,icc
    integer :: lskptr(3,Natom)

    if ( mcnstr == 0 ) return

    anorm(:,:)=0.0d0
    dx(:)=0.0d0

    if ( ityp /= 4 .and. ityp /= 2 .and. ityp /= 6 ) then

       stop "THIS TYPE OF CONSTRAINT HAS NOT BEEN IMPLIMENTED YET !!!"

    else if ( ityp == 4 .or. ityp == 2 .or. ityp == 6 ) then

       if ( ityp == 4 ) then
          allocate( askel(nodim,mcnstr,6) )     !ASKEL for DIST!!!!!
          do i=1,mcnstr
!            ..distance
             call fillc(ia(i),tau0,x1)             !a(i) is index for R1
             call fillc(ib(i),tau0,x2)             !b(i) is index for R2
             call pbc_correction( cval(i),x1,x2 )  !x1&x2 are modified by pbc
             call funcd(fc(i),fv(i),cval(i),x1,x2) !FC and FV must be global
             call diffd(dx,x1,x2)                  !derivative of sigma
             kmax=6
             do k=1,kmax
                do n=1,nodim
                   askel(:,:,:)=0.0d0
                   iaa=(ia(i)-1)*3
                   ibb=(ib(i)-1)*3
                   askel(iaa+1,i,1)=1.0d0        !ASKEL for DIST!!!!
                   askel(iaa+2,i,2)=1.0d0
                   askel(iaa+3,i,3)=1.0d0
                   askel(ibb+1,i,4)=1.0d0
                   askel(ibb+2,i,5)=1.0d0
                   askel(ibb+3,i,6)=1.0d0
                   anorm(n,i)=anorm(n,i)+dx(k)*askel(n,i,k)  
                end do
             end do
          end do
          deallocate( askel )
       end if

       if ( ityp == 2 ) then
          allocate( askel(nodim,mcnstr,9) )       !ASKEL for angle!!!!!
          do i=1,mcnstr
!            ..angle
             call fillc(ia(i),tau0,x1)
             call fillc(ib(i),tau0,x2)
             call fillc(ic(i),tau0,x3)
             call funct(fc(i),fv(i),cval(i),x1,x2,x3)
             call difft(dx,x1,x2,x3)
             kmax=9
             do k=1,kmax
                do n=1,nodim
                   askel(:,:,:)=0.0d0
                   iaa=(ia(i)-1)*3
                   ibb=(ib(i)-1)*3
                   icc=(ic(i)-1)*3
                   askel(iaa+1,i,1)=1.0d0        !ASKEL for DIST!!!!
                   askel(iaa+2,i,2)=1.0d0
                   askel(iaa+3,i,3)=1.0d0
                   askel(ibb+1,i,4)=1.0d0
                   askel(ibb+2,i,5)=1.0d0
                   askel(ibb+3,i,6)=1.0d0
                   askel(icc+1,i,7)=1.0d0
                   askel(icc+2,i,8)=1.0d0
                   askel(icc+3,i,9)=1.0d0
                   anorm(n,i)=anorm(n,i)+dx(k)*askel(n,i,k)
                end do
             end do
          end do
          deallocate( askel )
       end if

       if ( ityp == 6 ) then
          allocate( askel(nodim,mcnstr,6) )
          lskptr(:,:)=0
          do i=1,Natom
             do k=1,3
                lskptr(k,i)=i*3-(3-k)
             end do
          end do
          do i=1,mcnstr
             c_kappa=cnpar(1,i)
             c_rc   =cnpar(2,i)
             call coornum(ia(i),c_kappa,c_rc,tau0 &
                         ,anorm(:,i),lskptr,fc(i),fv(i),cval(i))
          end do
          deallocate( askel )
       end if

       do i=1,mcnstr
          do j=1,nodim
             fcstr(j)=fcstr(j)+xlagr(i)*anorm(j,i)
          end do
       end do
       call daxpy(nodim,1.0d0,fcstr(1),1,dxpar(1),1)

    end if

    return
  END SUBROUTINE cnstfc
!=======================================================================

!==================================================================
  SUBROUTINE fillc(iat,tau,x)
    implicit none
    integer :: iat
    real(8) ::  tau(3,Natom),x(3)
    integer :: id,naa,ityp
    x(1)=tau(1,iat)
    x(2)=tau(2,iat)
    x(3)=tau(3,iat)
    return
  END SUBROUTINE fillc
!------------------------------------------------------------------------
!------------------------------------------------------------------------
  SUBROUTINE funcd(fd,d,R0,x,y)
!   FUNCTION: R - R0  !Distance between x and y,fd becomes FC
    implicit none
    real(8) :: d, fd, r0, x(3), y(3)
    real(8) :: t3, t6, t9
    t3 = (x(1)-y(1))**2
    t6 = (x(2)-y(2))**2
    t9 = (x(3)-y(3))**2
    d  = sqrt(t3+t6+t9)
    fd = d-R0
    return
  END SUBROUTINE funcd
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
  SUBROUTINE funct(ft,t,T0,x,y,z)
!   FUNCTION: COS(theta) - COS(theta0)
    implicit none
!   Argument variables
    real(8) :: ft, t, t0, x(3), y(3), z(3)
!   Local variables
    real(8),parameter :: epsilon=1.d-12
    real(8) :: ac, caa, cab, cbb, t1, t2, t3, t4, t5, t6
    t1 = (x(1)-y(1))
    t2 = (x(2)-y(2))
    t3 = (x(3)-y(3))
    t4 = (y(1)-z(1))
    t5 = (y(2)-z(2))
    t6 = (y(3)-z(3))
    caa = t1*t1+t2*t2+t3*t3
    cab = t1*t4+t2*t5+t3*t6
    cbb = t4*t4+t5*t5+t6*t6
    if ( caa < epsilon .or. cbb < epsilon ) then
!      2 points at the same place (T.D. 16/12/1999).
!      We give the desired value (T0) for satisfied constraint.
!      ac=-cos(T0)
       ft=0.D0
       t = T0
    else
       ac = cab/sqrt(caa*cbb)
       if ( ac < -1.D0 ) ac=-1.D0
       if ( ac >  1.D0 ) ac= 1.D0
       ft = ac+cos(T0)
       t = acos(-ac)
    end if
    return
  END SUBROUTINE funct
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
  SUBROUTINE diffd(dR,x,y)
    implicit none
    real(8) :: dR(6),x(3),y(3)
    real(8),parameter :: epsilon=1.d-12
    real(8) :: r, t3, t6, t9
    t3 = (x(1)-y(1))
    t6 = (x(2)-y(2))
    t9 = (x(3)-y(3))
    r  = sqrt(t3*t3+t6*t6+t9*t9)
    if ( r < epsilon ) then
       dR(1) = 0.D0
       dR(2) = 0.D0
       dR(3) = 0.D0
       dR(4) = 0.D0
       dR(5) = 0.D0
       dR(6) = 0.D0
    else
       dR(1) =  t3/r
       dR(2) =  t6/r
       dR(3) =  t9/r
       dR(4) = -t3/r
       dR(5) = -t6/r
       dR(6) = -t9/r
    end if
    return
  END SUBROUTINE diffd
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
  SUBROUTINE difft(dT,x,y,z)
    implicit none
    real(8) :: dT(9), x(3), y(3), z(3)
    real(8),parameter :: epsilon=1.d-12
    real(8) :: caa, cab, cbb, ccc, t1, t2, t3, t4, t5, t6
    integer :: i
    t1 = (x(1)-y(1))
    t2 = (x(2)-y(2))
    t3 = (x(3)-y(3))
    t4 = (y(1)-z(1))
    t5 = (y(2)-z(2))
    t6 = (y(3)-z(3))
    caa = t1*t1+t2*t2+t3*t3
    cab = t1*t4+t2*t5+t3*t6
    cbb = t4*t4+t5*t5+t6*t6
    if ( caa < epsilon .or. cbb < epsilon ) then
       do i=1,9
          dT(i)=0.D0
       end do
    else
       ccc = -1.D0/sqrt(caa*cbb)
       dT(1) =  ccc*(cab/caa*t1-t4)
       dT(2) =  ccc*(cab/caa*t2-t5)
       dT(3) =  ccc*(cab/caa*t3-t6)
       dT(4) = -ccc*(cab/caa*t1-cab/cbb*t4+t1-t4)
       dT(5) = -ccc*(cab/caa*t2-cab/cbb*t5+t2-t5)
       dT(6) = -ccc*(cab/caa*t3-cab/cbb*t6+t3-t6)
       dT(7) = -ccc*(cab/cbb*t4-t1)
       dT(8) = -ccc*(cab/cbb*t5-t2)
       dT(9) = -ccc*(cab/cbb*t6-t3)
    end if
    return
  END SUBROUTINE difft
!---------------------------------------------------------------------------
!--------1---------2--------3---------4---------5---------6---------7--
  SUBROUTINE pbc_correction( cnstr_dist, x, y )
    implicit none
    real(8),intent(IN) :: cnstr_dist
    real(8),intent(INOUT) :: x(3), y(3)
    integer :: i1,i2,i3,j1,j2,j3
    real(8) :: ti(3),tj(3),xtmp(3),ytmp(3),cc,dd,minerr
    cc=cnstr_dist**2
    dd=sum( (x(:)-y(:))**2 )
    minerr=abs(cc-dd)
    xtmp(:)=x(:)
    ytmp(:)=y(:)
    do i3=-1,1
    do i2=-1,1
    do i1=-1,1
       ti(1)=aa(1,1)*i1+aa(1,2)*i2+aa(1,3)*i3 + x(1)
       ti(2)=aa(2,1)*i1+aa(2,2)*i2+aa(2,3)*i3 + x(2)
       ti(3)=aa(3,1)*i1+aa(3,2)*i2+aa(3,3)*i3 + x(3)
       do j3=-1,1
       do j2=-1,1
       do j1=-1,1
          tj(1)=aa(1,1)*j1+aa(1,2)*j2+aa(1,3)*j3 + y(1)
          tj(2)=aa(2,1)*j1+aa(2,2)*j2+aa(2,3)*j3 + y(2)
          tj(3)=aa(3,1)*j1+aa(3,2)*j2+aa(3,3)*j3 + y(3)
          dd=sum( (ti(:)-tj(:))**2 )
          if ( abs(cc-dd) < minerr ) then
             xtmp(:)=ti(:)
             ytmp(:)=tj(:)
             minerr=abs(cc-dd)
          end if
       end do ! j1
       end do ! j2
       end do ! j3
    end do ! i1
    end do ! i2
    end do ! i3
    x(:)=xtmp(:)
    y(:)=ytmp(:)
  END SUBROUTINE pbc_correction
!----------------------------------------------------------------------
!---------------------------------------------------------------------------
  SUBROUTINE coornum(ia_dum,c1_kappa,c1_rc,tscr,an,lsk,fci,fvi,cval_dum)
!==--------------------------------------------------------------==
    implicit none
    integer :: ia_dum,lsk(:,:)
    real(8) ::  tscr(:,:),an(:),fci,fvi,cval_dum,c1_kappa,c1_rc,xx(3)
    real(8) ::  x0,y0,z0,dx,dy,dz,dd,ff,df,fff
    integer iat,k,l1,l2,l3
    integer naa, jj,kx,jx,jy,jz,j,is,i
    real(8),save :: cval_save
!==--------------------------------------------------------------==
    if ( ia_dum <= Natom ) then
       l1 = lsk(1,ia_dum)
       l2 = lsk(2,ia_dum)
       l3 = lsk(3,ia_dum)
    endif
    call fillc(ia_dum,tscr,xx)
    x0 = xx(1)
    y0 = xx(2)
    z0 = xx(3)
    fvi = 0.d0
    dd  = 0.0d0
    do iat=1,Natom
!-------------------------------------------------------
       do i=1,10
          if ( ki_atom(iat) == index(i) ) goto 101
       enddo
!-------------------------------------------------------
       call pbcdis(tscr(:,iat),x0,y0,z0,dx,dy,dz,dd)
!       dx=tscr(1,iat)-x0
!       dy=tscr(2,iat)-Y0
!       dz=tscr(3,iat)-z0
!       dd=dsqrt(dx*dx+dy*dy+dz*dz)
       if ( dd > 1.d-2 ) then
          df=c1_kappa*(dd-c1_rc)
          fvi=fvi+1.d0/(dexp(df)+1.d0)
          ff=-0.5d0*c1_kappa/(cosh(df)+1.d0)/dd
          if ( lsk(1,iat) /= 0 ) then
             k=lsk(1,iat)
             an(k)=an(k)+ff*dx
          endif
          if ( lsk(2,iat) /= 0 ) then
             k=lsk(2,iat)
             an(k)=an(k)+ff*dy
          endif
          if ( lsk(3,iat) /= 0 ) then
             k=lsk(3,iat)
             an(k)=an(k)+ff*dz
          end if
          if ( l1 /= 0 ) an(l1)=an(l1)-ff*dx
          if ( l2 /= 0 ) an(l2)=an(l2)-ff*dy
          if ( l3 /= 0 ) an(l3)=an(l3)-ff*dz
       end if
!------------------------------------------------------
101    continue
!------------------------------------------------------
    end do
    if ( cval_dum == 999.9d0 .and. itime == 1 ) then
       cval_save=fvi
    end if
    if ( cval_dum == 999.9d0 ) then
       fci=fvi-cval_save
    else
       fci=fvi-cval_dum
    end if
    return
  END SUBROUTINE coornum
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  SUBROUTINE pbcdis(target,x0,y0,z0,dx,dy,dz,dd)
    implicit none
    integer :: fact1(27),fact2(27),fact3(27)
    integer :: i
    real(8) :: xx(27),yy(27),zz(27)
    real(8) :: target(3,1)
    real(8) :: d(27)
    real(8) :: x0,y0,z0,dd,dx,dy,dz
!           1 2 3 4 5 6 7 8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27      
    fact1=(/0,1,0,0,1,1,0,1,-1, 0, 0,-1, 1, 1,-1, 0, 0,-1,-1, 0, 1, 1,-1, 1,-1,-1,-1/)        
    fact2=(/0,0,1,0,1,0,1,1, 0,-1, 0, 1,-1, 0, 0,-1, 1,-1, 0,-1, 1,-1, 1,-1, 1,-1,-1/)        
    fact3=(/0,0,0,1,0,1,1,1, 0, 0,-1, 0, 0,-1, 1, 1,-1, 0,-1,-1,-1, 1, 1,-1,-1, 1,-1/)    
    dd=0.0d0
    dx=0.0d0
    dy=0.0d0
    dz=0.0d0
!    write(*,*) "aa(1,1)--->", aa(1,1) 
    do i=1,27             
       xx(i)=target(1,1)-x0+fact1(i)*aa(1,1)+fact2(i)*aa(1,2)+fact3(i)*aa(1,3)
       yy(i)=target(2,1)-y0+fact1(i)*aa(2,1)+fact2(i)*aa(2,2)+fact3(i)*aa(2,3)
       zz(i)=target(3,1)-z0+fact1(i)*aa(3,1)+fact2(i)*aa(3,2)+fact3(i)*aa(3,3)
    end do

    do i=1,27
       d(i)=DSQRT(xx(i)*xx(i)+yy(i)*yy(i)+zz(i)*zz(i))
    enddo

    dd=minval(d)

    do i=1,27
       if ( dd == d(i) ) then
          dx=xx(i)
          dy=yy(i)
          dz=zz(i)
       end if
    end do

    return
  END SUBROUTINE pbcdis
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
  SUBROUTINE read_blue
    implicit none
    integer :: i,ispecies

    ia(:)=0
    ib(:)=0
    ic(:)=0
    index(:)=0
    cval(:)=0.0d0
    cnpar(:,:)=0.0d0

    open(888,file='blue_control.dat')
    read(888,*) ityp
    read(888,*) mcnstr
    if ( ityp == 4 ) then
       read(888,*) (ia(i),i=1,mcnstr)
       read(888,*) (ib(i),i=1,mcnstr)
       read(888,*) (cval(i),i=1,mcnstr)
    end if
    if ( ityp == 2 ) then
       read(888,*) (ia(i),i=1,mcnstr)
       read(888,*) (ib(i),i=1,mcnstr)
       read(888,*) (ic(i),i=1,mcnstr)
       read(888,*) (cval(i),i=1,mcnstr)
    end if
    if ( ityp == 6 ) then
       read(888,*) ispecies
       read(888,*) (ia(i)     ,i=1,mcnstr)
       read(888,*) (cnpar(1,i),i=1,mcnstr)  !c_kappa
       read(888,*) (cnpar(2,i),i=1,mcnstr)  !c_rc
       read(888,*) (cval(i),i=1,mcnstr)
       read(888,*) (index(i),i=1,ispecies)
    end if
    close(888)
    return
  END SUBROUTINE read_blue
!---------------------------------------------------------------------------------------
  SUBROUTINE bcast_blue_data 
    implicit none
    integer :: ierr

    call mpi_bcast(ityp,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(mcnstr,1,mpi_integer,0,mpi_comm_world,ierr)

!     Bcast and allocation
    if ( ityp == 4 ) then
       write(*,*) "ityp--->",ityp,myrank
       write(*,*) "mcnstr->",mcnstr,myrank
       call mpi_bcast(ia(1),10,mpi_integer,0,mpi_comm_world,ierr)
       call mpi_bcast(ib(1),10,mpi_integer,0,mpi_comm_world,ierr)
       call mpi_bcast(cval(1),10,mpi_real8,0,mpi_comm_world,ierr)
    end if


    if ( ityp == 2 ) then
       write(*,*) "ityp--->",ityp,myrank
       write(*,*) "mcnstr->",mcnstr,myrank
       call mpi_bcast(ia(1),10,mpi_integer,0,mpi_comm_world,ierr)
       call mpi_bcast(ib(1),10,mpi_integer,0,mpi_comm_world,ierr)
       call mpi_bcast(ic(1),10,mpi_integer,0,mpi_comm_world,ierr)
       call mpi_bcast(cval(1),10,mpi_real8,0,mpi_comm_world,ierr)
    end if

    if ( ityp == 6 ) then
       call mpi_bcast(ia(1),10,mpi_integer,0,mpi_comm_world,ierr)
       call mpi_bcast(index(1),10,mpi_integer,0,mpi_comm_world,ierr)
       call mpi_bcast(cnpar(1,1),20,mpi_integer,0,mpi_comm_world,ierr)
       call mpi_bcast(cval(1),10,mpi_real8,0,mpi_comm_world,ierr)
    end if

    allocate( ylagr(mcnstr) )
    allocate( xlagr(mcnstr) )
    ylagr=0.0d0
    xlagr=0.0d0
    nodim=3*Natom
    allocate(fc(mcnstr))
    allocate(fv(mcnstr))
    allocate(anorm(nodim,mcnstr))
    allocate(ipvt(mcnstr))

    write(*,*) "PASS THE MYRANK-->", myrank
    return
  END SUBROUTINE bcast_blue_data
!---------------------------------------------------------------------------------------
  SUBROUTINE dealloc_blue
    implicit none
    deallocate(ylagr)
    deallocate(xlagr)
    deallocate(ipvt)
    deallocate(DT2BYM,DTB2MI)
!    deallocate(pm_dim)
    return
  END SUBROUTINE dealloc_blue
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
  SUBROUTINE write_blue_data(nfi)
    implicit none
    integer :: nfi,j
    real(8) :: fval
    open(889,file='CONSTRAINT')
    if ( mod(nfi,1) == 0 ) then
       do j=1,mcnstr
          fval=fv(j)
          write(889,'(i7,2x,i4,5x,2(1pe20.10))') nfi,j,xlagr(j),fval
       enddo
    end if
    return
  END SUBROUTINE write_blue_data
!----------------------------------------------------------------------------------------
END MODULE blue_moon_module



