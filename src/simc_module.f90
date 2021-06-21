MODULE simc_module

  use bberf_module

  implicit none

  PRIVATE
  PUBLIC :: simc
  PUBLIC :: fit_initrho_simc

  real(8),allocatable :: rads(:),vins(:),wgt(:),finp(:),rr(:)
  real(8) :: zvs

CONTAINS


  SUBROUTINE simc(rad,vin,rc,zv,parloc,mesh)
    implicit none
    integer,intent(IN) :: mesh
    real(8),intent(IN) :: rad(mesh),vin(mesh),zv,rc
    real(8),intent(OUT) :: parloc(4)
    integer :: ierr
    call simc_0(rad,vin,rc,zv,parloc,mesh,ierr)
    if ( ierr < 0 ) then
      call simc_1(rad,vin,rc,zv,parloc,mesh)
    end if
  END SUBROUTINE simc


  SUBROUTINE simc_1( rad, vin, rc, zv, parloc, mesh )

    implicit none
    integer,intent(IN) :: mesh
    real(8),intent(IN) :: rad(mesh),vin(mesh),zv,rc
    real(8),intent(OUT) :: parloc(4)
    integer,parameter :: lwa0=15, nummin=6
    integer :: info,k,num,ipvt(3),maxnum, lwa
    real(8) :: x(3), nxtsmp, pi
    real(8),allocatable :: wa(:),fvec(:),fjac(:,:)
    real(8),parameter :: lambda=3.5d0
    real(8),parameter :: x1ini=1.0d0, x2ini=0.4d0, x3ini=0.6d0
    real(8),parameter :: tol=1.0d-5, smpstp=0.2d0
    real(8),parameter :: rmax=10.0d0, vmax=100.0d0, vrzmin=3.0d-6

    call write_border( 1, 'simc_1(start)' )

    allocate( rads(mesh)   ) ; rads=0.0d0
    allocate( vins(mesh)   ) ; vins=0.0d0
    allocate( wgt(mesh)    ) ; wgt=0.0d0
    allocate( fvec(mesh)   ) ; fvec=0.0d0
    allocate( fjac(mesh,3) ) ; fjac=0.0d0

    pi = 4.0d0*atan(1.0d0)

    num=0
    nxtsmp=0.0d0

    do k=1,mesh

       if ( rad(k) > nxtsmp ) then

          nxtsmp = nxtsmp + smpstp

          if ( abs(vin(k)) <= vmax ) then
             num = num+1
             rads(num) = rad(k)
             vins(num) = vin(k)
             wgt(num)  = 1.0d0-exp(-(1.2d0*rad(k)/rc)**lambda)
          end if
          if ( (abs(vin(k)*rad(k)+zv) < vrzmin .or. &
                rad(k) > rmax ) .and. num > nummin ) exit

       end if

    end do ! k

    maxnum = num
    lwa = lwa0 + maxnum

    allocate( wa(lwa) ) ; wa=0.0d0

    zvs  = zv
    x(1) = x1ini
    x(2) = x2ini
    x(3) = x3ini

    call levenberg_marquardt(num,3,x,fvec,fjac)

    if ( x(2) < 0.0d0 .or. x(3) < 0.0d0 ) then
       stop "simc: illegally converged."
    end if

    parloc(1) = x(1)
    parloc(2) = x(2)
    parloc(3) = 1.0d0 - x(1)
    parloc(4) = x(3)

    deallocate( wa   )
    deallocate( fjac )
    deallocate( fvec )
    deallocate( wgt  )
    deallocate( vins )
    deallocate( rads )

    call write_border( 1, 'simc_1(end)' )

    return
  END SUBROUTINE simc_1


  SUBROUTINE uscfit_1( m, n, x, fvec, fjac, ldfjac, iflag )
    implicit none
    integer,intent(IN) :: m,n,ldfjac,iflag
    real(8),intent(OUT) :: fvec(m),fjac(ldfjac,n)
    real(8),intent(INOUT) :: x(n)
    real(8) :: pi,xxxx
    integer :: i

    pi = 4.0d0*atan(1.0d0)
    xxxx = 1.0d0 - x(1)

    if ( x(2) < 0.0d0 ) x(2)=0.0d0
    if ( x(3) < 0.0d0 ) x(3)=0.0d0

    if ( iflag == 1 ) then

       do i=1,m
          fvec(i) = (- zvs/rads(i)*( x(1)*bberf(sqrt(x(2))*rads(i)) &
                                   + xxxx*bberf(sqrt(x(3))*rads(i)) &
                                   ) - vins(i) )*sqrt(wgt(i))
       end do

    else if ( iflag == 2 ) then

       do i=1,m
          fjac(i,1) = -zvs/rads(i)*( bberf(sqrt(x(2)*rads(i))) &
                                   - bberf(sqrt(x(3)*rads(i))) &
                                   )*sqrt(wgt(i))
          fjac(i,2) = -zvs/sqrt(pi*x(2))*x(1) &
                      *exp(-x(2)*rads(i)**2)*sqrt(wgt(i))
          fjac(i,3) = -zvs/sqrt(pi*x(3))*(1.0d0-x(1)) &
                      *exp(-x(3)*rads(i)**2)*sqrt(wgt(i))
       end do

    else
       stop "Error in vlfit: iflag must be 1 or 2."
    end if

    return
  END SUBROUTINE uscfit_1


  SUBROUTINE levenberg_marquardt(ndim,nvec,x,fvec,fjac)

    implicit none
    integer,intent(IN) :: ndim,nvec
    real(8),intent(INOUT) :: x(nvec)
    real(8),intent(INOUT) :: fvec(ndim),fjac(ndim,nvec)
    integer,parameter :: max_loop=10000, max_loop0=20
    integer,allocatable :: ipiv(:)
    integer :: m,n,ierr,loop,loop0,num_conv
    real(8),parameter :: delta=1.d-8
    real(8),allocatable :: Hes(:,:),dJ(:),du(:),xtmp(:),ftmp(:),xmin(:)
    real(8) :: c,J,Jtmp,err,Jmin,errmin
    character(80) :: error_message

    allocate( Hes(nvec,nvec) ) ; Hes=0.0d0
    allocate( dJ(nvec)       ) ; dJ=0.0d0
    allocate( du(nvec)       ) ; du=0.0d0
    allocate( ipiv(nvec)     ) ; ipiv=0
    allocate( xtmp(nvec)     ) ; xtmp=0.0d0
    allocate( ftmp(ndim)     ) ; ftmp=0.0d0
    allocate( xmin(ndim)     ) ; xmin=0.0d0

    xmin = x
    Jmin = 1.d100
    num_conv = 0
    errmin = 1.d100

    do loop0=1,max_loop0

       c=0.0001d0

       call uscfit_1( ndim, nvec, x, fvec, fjac, ndim, 1 )
       J=sum( fvec(:)**2 )
       call uscfit_1( ndim, nvec, x, fvec, fjac, ndim, 2 )

       do loop=1,max_loop

          do n=1,nvec
             do m=1,nvec
                Hes(m,n) = sum( fjac(:,m)*fjac(:,n) )
             end do
             Hes(n,n) = Hes(n,n) + c*Hes(n,n)
          end do
          do n=1,nvec
             dJ(n) = sum( fvec(:)*fjac(:,n) )
          end do

          du(:) = -dJ(:)
          call dgesv(nvec,1,Hes,nvec,ipiv,du,nvec,ierr)

          xtmp(:)=x(:)+du(:)
          call uscfit_1( ndim, nvec, xtmp, ftmp, fjac, ndim, 1 )
          Jtmp=sum( ftmp(:)**2 )

          err = sum( du(:)**2 )
          errmin = min( err, errmin )
          if ( err < delta ) then
             num_conv = num_conv + 1
             exit
          end if

          if ( Jtmp > J ) then
             c=10.0d0*c
          else if ( Jtmp <= J ) then
             c=0.1d0*c
             J=Jtmp
             x(:)=xtmp(:)
             fvec(:)=ftmp(:)
             call uscfit_1( ndim, nvec, x, fvec, fjac, ndim, 2 )
          end if

       end do ! loop

!       if ( loop <= max_loop ) then
          if ( Jtmp < Jmin ) then
             Jmin = Jtmp
             xmin = xtmp
          end if
!       end if

       call random_number(x)

    end do ! loop0

    if ( num_conv == 0 ) then
       write(error_message,*) "Jmin,errmin=",Jmin,errmin
       call write_string( error_message )
!       stop "error@levenberg_marquardt"
    end if

    x = xmin

! ---
    deallocate( ftmp )
    deallocate( xtmp )
    deallocate( ipiv )
    deallocate( du   )
    deallocate( dJ   )
    deallocate( Hes  )
  END SUBROUTINE levenberg_marquardt

!#ifdef _TEST_
  SUBROUTINE simc_0(rad,vin,rc,zv,parloc,mesh,ierr)

    implicit none
    integer :: mesh
    real(8) :: rad(mesh),vin(mesh),parloc(4),zv,rc
    integer,intent(out) :: ierr
    integer,parameter :: maxnum = 2000
    integer,parameter :: lwa=5*3+maxnum
    integer :: k,num,info,ipvt(3)
    real(8) :: pi,nxtsmp,x(3)
    real(8),allocatable :: fvec(:),fjac(:,:),wa(:)
    real(8),parameter :: lambda = 3.5d0
    real(8),parameter :: x1ini=1.0d0, x2ini=0.4d0, x3ini=0.6d0
    real(8),parameter :: tol=1.0d-5
    real(8),parameter :: rmax=10.0d0, vmax=100.0d0, vrzmin=3.0d-6
    real(8),parameter :: smpstp=0.2d0, nummin=6

    call write_border( 1, 'simc_0(start)' )

    ierr=0

    pi = 4.0d0*atan(1.0d0)

    allocate( fvec(maxnum)   ) ; fvec=0.0d0
    allocate( fjac(maxnum,3) ) ; fjac=0.0d0
    allocate( wa(lwa)        ) ; wa=0.0d0

    allocate( rads(maxnum)   ) ; rads=0.0d0
    allocate( vins(maxnum)   ) ; vins=0.0d0
    allocate( wgt(maxnum)    ) ; wgt=0.0d0

    num=0
    nxtsmp=0.0d0
    do k=1,mesh
       if ( rad(k) > nxtsmp ) then
          nxtsmp = nxtsmp + smpstp
          if ( abs(vin(k)) <= vmax ) then ![ 1.0*rc -- 1.5*rc in original ]
             num=num+1
             if ( num > maxnum ) then
                call stop_program( "simc_0: Too many sample points." )
             end if
             rads(num) = rad(k)
             vins(num) = vin(k)
             wgt(num)  = 1.0d0-exp(-(1.2d0*rad(k)/rc)**lambda)
          end if
          if ( ( abs(vin(k)*rad(k)+zv) < vrzmin .or. rad(k) > rmax ) &
               .and. num > nummin ) exit
       end if
    end do
    zvs = zv

    x(1:3) = (/ x1ini, x2ini, x3ini /)

    call lmder1(uscfit,num,3,x,fvec,fjac,maxnum,tol,info,ipvt,wa,lwa)

    if ( info==0 .or. info==4 .or. info==5 .or. info==6 .or. info==7 ) then
       write(*,*) 'simc_0: Not converged (info=',info,')'
       write(*,*) 'x(1) = ',x(1)
       write(*,*) 'x(2) = ',x(2)
       write(*,*) 'x(3) = ',x(3)
       write(*,'(A5,3A20)') 'k','rads(k)','vins(k)','wgt(k)'
       do k=1,num
          write(*,'(I5,3g20.12)')k,rads(k),vins(k),wgt(k)
       end do
       !call stop_program( "simc_0(2)" )
       ierr=-1
       goto 999
    end if

    if ( x(2) < 0.0d0 .or. x(3) < 0.0d0 ) then
       !call stop_program( "simc_0: illegally converged." )
       ierr=-1
       goto 999
    end if

    parloc(1) = x(1)
    parloc(2) = x(2)
    parloc(3) = 1.0d0 - x(1)
    parloc(4) = x(3)

999 deallocate( wgt  )
    deallocate( vins )
    deallocate( rads )
    deallocate( wa   )
    deallocate( fjac )
    deallocate( fvec )

    call write_border( 1, 'simc_0(end)' )

  END SUBROUTINE simc_0
!
!     fitting function for simc
!
  SUBROUTINE uscfit(m,n,x,fvec,fjac,ldfjac,iflag)
    implicit none
    integer :: m,n,ldfjac,iflag
    real(8) :: x(n),fvec(m),fjac(ldfjac,n)
    real(8) :: pi
    integer :: i
    pi = 4.0d0*atan(1.0d0)
    if ( x(2) < 0.0d0 ) x(2)=0.0d0
    if ( x(3) < 0.0d0 ) x(3)=0.0d0
    if ( iflag == 1 ) then
       do i=1,m
          fvec(i) = (- zvs/rads(i)*(x(1)*bberf(sqrt(x(2))*rads(i)) &
               + (1.0d0-x(1))*bberf(sqrt(x(3))*rads(i)))-vins(i))*sqrt(wgt(i))
       end do
    else if ( iflag == 2 ) then
       do i=1,m
          fjac(i,1) = -zvs/rads(i) &
               *(bberf(sqrt(x(2)*rads(i)))-bberf(sqrt(x(3)*rads(i)))) &
               *sqrt(wgt(i))
          fjac(i,2) = -zvs/sqrt(pi*x(2))*x(1) &
               *exp(-x(2)*rads(i)**2)*sqrt(wgt(i))
          fjac(i,3) = -zvs/sqrt(pi*x(3))*(1.0d0-x(1)) &
               *exp(-x(3)*rads(i)**2)*sqrt(wgt(i))
       end do
    else
       call stop_program( "Error in vlfit: iflag must be 1 or 2." )
    end if
  END SUBROUTINE uscfit
!#endif


  SUBROUTINE fit_initrho_simc( rad_grid, rho, x )
    implicit none
    real(8),intent(IN)  :: rad_grid(:), rho(:)
    real(8),intent(OUT) :: x(:,:)
    real(8),allocatable :: fvec(:),fjac(:,:),wa(:),y(:)
    real(8),parameter :: tol=1.0d-5
    integer,allocatable :: ipvt(:)
    integer :: m,n,info,lwa
    m=size(rho)
    n=size(x)
    lwa=5*n+m
    allocate( finp(m)   ) ; finp=rho
    allocate( rr(m)     ) ; rr=rad_grid(:)
    allocate( fvec(m)   ) ; fvec=0.0d0
    allocate( fjac(m,n) ) ; fjac=0.0d0
    allocate( wa(lwa)   ) ; wa=0.0d0
    allocate( ipvt(n)   ) ; ipvt=0
! --- initial values of parameters (Si_psv.dat) ---
    x=0.0d0
    x(1:3,1) = (/ 0.877033839192d-2, 0.126999133549d-4, 0.124643066149d0 /)
    x(1:3,2) = (/ 0.851848549949d-1,-0.252885386521d-6, 0.297264415878d0 /)
    x(1:3,3) = (/-0.227044323933d+1, 0.458882466862d0 , 0.106830656608d+1/)
    x(1:3,4) = (/ 0.220373089866d+1, 0.504784143365d0 , 0.150801442780d+1/)
! ---
    call lmder1( fit_func_sub,m,n,x,fvec,fjac,m,tol,info,ipvt,wa,lwa )
    deallocate( ipvt )
    deallocate( wa   )
    deallocate( fjac )
    deallocate( fvec )
    deallocate( rr   )
    deallocate( finp )
  END SUBROUTINE fit_initrho_simc

  SUBROUTINE fit_func_sub( m,n,x,fvec,fjac,ldfjac,iflag )
    implicit none
    integer :: m,n,ldfjac,iflag
    real(8) :: x(n),fvec(m),fjac(ldfjac,n)
    real(8),parameter :: pi4=12.56637061435917d0
    real(8) :: r2
    integer :: i,j
    select case( iflag )
    case( 1 )
       fvec=0.0d0
       do j=1,12,3
       do i=1,m
          r2=rr(i)**2
          fvec(i) = fvec(i) + pi4*r2*(x(j)+x(j+1)*r2)*exp(-x(j+2)*r2)
       end do
       end do
       fvec(:) = fvec(:) - finp(:)
    case( 2 )
       do j=1,12,3
       do i=1,m
          r2=rr(i)**2
          fjac(i,j  ) = pi4*r2*exp(-x(j+2)*r2)
          fjac(i,j+1) = pi4*r2*r2*exp(-x(j+2)*r2)
          fjac(i,j+2) = pi4*r2*(x(j)+x(j+1)*r2)*(-r2)*exp(-x(j+2)*r2)
       end do
       end do
    case default
       call stop_program("error@ fit_func_sub in simc_module")
    end select
  END SUBROUTINE fit_func_sub

END MODULE simc_module
