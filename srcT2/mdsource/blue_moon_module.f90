MODULE blue_moon_module

  use parallel_module
  use cpmd_variables, only: amu, pmass, dt, itime
  use atom_module, only: Natom, ki_atom, zn_atom
  use aa_module, only: aa

  implicit none

  PRIVATE
  PUBLIC :: shake, rattle, write_blue_data, read_blue

  real(8),allocatable :: xlagr(:), ylagr(:)
  real(8),allocatable :: Rion0(:,:)
  real(8),allocatable :: dtm3(:), dtm1(:)
  integer,allocatable :: ia(:),ib(:),ic(:)
  real(8),allocatable :: cval(:), cnpar(:,:)
  real(8),allocatable :: anorm(:,:)
  real(8),allocatable :: fc(:), fv(:)
  integer,allocatable :: ipvt(:)
  integer :: ityp, mcnstr, nodim
  integer :: index(10)

CONTAINS


  SUBROUTINE shake( Rion, Velocity )

    implicit none
    real(8),intent(INOUT) :: Rion(:,:)
    real(8),intent(INOUT) :: Velocity(:,:)
    integer,parameter :: maxrat=5000, mrdiis=5
    integer :: i,j,iter,k,info,idiis,m,i_dim
    integer,save :: istart1=0
    real(8),parameter :: tolf=1.d-7, tolx=1.d-8
    real(8) :: tau0(3,Natom),taup(3,Natom)
    real(8) :: v(mrdiis+1),diism(mrdiis+1,mrdiis+1),cnmax,fact,errf,errx
    real(8) :: pm
    real(8) :: tscr(3,Natom),dx(nodim),asl(mcnstr,mcnstr)
    real(8) :: xlo(mcnstr,mrdiis),err(mcnstr,mrdiis)
    real(8) :: an0(nodim,mcnstr)

    if ( mcnstr == 0 ) return

    if ( istart1 == 0 ) then
       allocate( Rion0(3,Natom) ) ; Rion0=0.0d0
       allocate( dtm3(nodim)  ) ; dtm3=0.0d0
       allocate( dtm1(Natom)  ) ; dtm1=0.0d0
       do i=1,Natom
          pm = pmass( zn_atom(ki_atom(i)) )*amu
          dtm3(3*i-2)=dt/(2.0d0*pm)
          dtm3(3*i-1)=dt/(2.0d0*pm)
          dtm3(3*i)  =dt/(2.0d0*pm)
          dtm1(i)    =dt/(2.0d0*pm)
       end do
       xlagr(:)=0.0d0
       istart1 =1
    else
       xlagr(:)=ylagr(:)
    end if

! ---

    Rion0(:,:) = Rion(:,:)
    Rion(:,:)  = Rion(:,:) + Velocity(:,:)*dt

! ---

    tau0(:,:) = Rion0(:,:)
    taup(:,:) = Rion(:,:)

    call calc_constraint( tau0, anorm, fc )

    an0(:,:) = anorm(:,:)

    dx(:)=0.0d0
    do i=1,mcnstr
       do j=1,nodim
          dx(j)=dx(j)+xlagr(i)*an0(j,i)
       end do
    end do

    call put_1d_to_2d( tscr, dx )
    do i=1,Natom
       fact=-dt*dtm1(i)
       tau0(1:3,i)=taup(1:3,i)+fact*tscr(1:3,i)
    end do

! --- Iteration ---

    do iter=1,maxrat

       call calc_constraint( tau0, anorm, fc )

       fc(:)=-fc(:)
       errf = sum( abs(fc) )

!       if ( myrank == 0 ) write(*,*) "ERRF--->", errf, iter

       if ( errf < tolf ) then
          Rion(:,:)=tau0(:,:)
          goto 100
       end if

! --- DIIS CALCULATION ---

       do i=1,mcnstr
          do j=1,mcnstr
             asl(i,j)=0.0d0
             do k=1,nodim
                asl(i,j)=asl(i,j)-2.0d0*dt*dt*dtm3(k)*anorm(k,i)*an0(k,j)
             end do
          end do
       end do

       call DGESV(mcnstr,1,asl,mcnstr,ipvt,fc,mcnstr,info)
       if ( info /= 0 ) call stop_program( "stop@shake(1)" )

       errx = sum( abs(fc) )

! --- DIIS!

       idiis=mod(iter-1,mrdiis)+1
       xlo(:,idiis)=xlagr(:)
       err(:,idiis)=fc(:)

       if ( iter > mrdiis ) then

          m=mrdiis+1
          diism(:,:)=0.0d0
          v(:)=0.0d0
          do i=1,mrdiis
             do j=1,mrdiis
                diism(i,j)=sum( err(:,i)*err(:,j) )
             end do
             diism(m,i)=1.0d0
             diism(i,m)=1.0d0
          end do
          v(m)=1.0d0
          call least_square(diism,m,m,v)
          fc(:)=0.0d0
          xlagr(:)=0.0d0
          do i=1,mrdiis
             do j=1,mcnstr
                fc(j)=fc(j)+v(i)*err(j,i)
                xlagr(j)=xlagr(j)+v(i)*xlo(j,i)
             end do
          end do

       end if

       xlagr(:) = xlagr(:) + fc(:)

       if ( errx < tolx ) goto 100

!--- Update forces

       dx(:)=0.0d0
       do i=1,mcnstr
          do j=1,nodim
             dx(j)=dx(j)+xlagr(i)*an0(j,i)
          end do
       end do

!--- Update the positions

       call put_1d_to_2d( tscr, dx )
       do i=1,Natom
          fact=-dt*dtm1(i)
          tau0(1:3,i)=taup(1:3,i)+fact*tscr(1:3,i)
       end do
       Rion(:,:)=tau0(:,:)

    end do ! iter

    if ( myrank == 0 ) write(*,'(30a)') 'shake| did not converge!!!'

    stop

!   MODIFY THE VELOCITIIES

100 continue

    dx(:)=0.0d0
    do i=1,mcnstr
       do j=1,nodim
          dx(j)=dx(j)+xlagr(i)*an0(j,i)
       end do
    end do

    call put_1d_to_2d( tscr, dx )

    do i=1,Natom
       Velocity(1:3,i) = Velocity(1:3,i) - dtm1(i)*tscr(1:3,i)
    end do

  END SUBROUTINE shake


  SUBROUTINE rattle( Rion, Velocity )

    implicit none
    real(8),intent(IN) :: Rion(:,:)
    real(8),intent(INOUT) :: Velocity(:,:)
    real(8) :: tau0(3,Natom)
    real(8),parameter :: tols=1.d-12, tolx=1.d-8
    real(8),allocatable :: aux(:) 
    real(8) :: errx,cnmax
    real(8) :: dx(nodim),dvv(nodim),asl(mcnstr,mcnstr)
    integer :: naux,i,j,k,info,nrank

    if ( mcnstr == 0 ) return

    tau0(:,:)=Rion(:,:) 
      
    call calc_constraint( tau0, anorm, fc )
    call put_2d_to_1d( Velocity, dx )
    asl(:,:)=0.0d0
    do i=1,mcnstr
       dvv(i)=0.0d0
       do j=1,nodim
          dvv(i)=dvv(i)+anorm(j,i)*dx(j)
       end do
       do j=1,mcnstr
          do k=1,nodim
             asl(i,j)=asl(i,j)+anorm(k,i)*anorm(k,j)*dtm3(k)
          end do
       end do
    end do

! ---

    naux=6*nodim+mcnstr
    allocate( aux(naux) )

    call DGELSS &
         (mcnstr,mcnstr,1,asl,mcnstr,dvv,mcnstr,fc,1.d-12,nrank,aux,naux,info)
    if ( info /= 0 ) stop "rattle','dgelss ended with info.ne.0"

    deallocate( aux )

! ---

    do i=1,mcnstr
       ylagr(i)=dvv(i)
    end do
    do j=1,nodim
       dvv(j)=0.0d0
       do i=1,mcnstr
          dvv(j)=dvv(j)-ylagr(i)*anorm(j,i)*dtm3(j)
       end do
    end do
    do j=1,nodim
       dx(j)=dx(j)+dvv(j)
    end do

    call put_1d_to_2d( Velocity, dx )

! ---

    do i=1,mcnstr
       fc(i)=0.0d0
       do j=1,nodim
          fc(i)=fc(i)+anorm(j,i)*dx(j)
       end do
    end do

    errx = sum( abs(fc) )
    if ( errx > tolx ) then
       call stop_program( " rattle| CONSTRAINT CAN NOT BE FULFILLED" )
    end if

    return
  END SUBROUTINE rattle


  SUBROUTINE put_1d_to_2d( tau, xin )
    implicit none
    real(8),intent(OUT) ::  tau(:,:)
    real(8),intent(IN)  :: xin(:)
    integer :: j,k,i
    k=0
    do i=1,size(tau,2)
       do j=1,size(tau,1)
          k=k+1
          tau(j,i) = xin(k)
       end do
    end do
  END SUBROUTINE put_1d_to_2d


  SUBROUTINE put_2d_to_1d( tau, xout )
    implicit none
    real(8),intent(IN)  :: tau(:,:)
    real(8),intent(OUT) :: xout(:)
    integer :: j,k,i
    k=0
    do i=1,size(tau,2)
       do j=1,size(tau,1)
          k=k+1
          xout(k) = tau(j,i)
       end do
    end do
  END SUBROUTINE put_2d_to_1d


  SUBROUTINE least_square(b,ldb,ndim,v)

    implicit none
    integer :: ldb,ndim, info
    real(8) :: b(:,:),v(:)
    integer :: nr,lw
    integer,parameter :: maxdis=20
    integer,parameter ::mrdiis=5
    real(8) :: diism(mrdiis+1,mrdiis+1)
    real(8) :: toleig=2.0d-16
    real(8) :: scr1(maxdis+1,maxdis+1),scr2(maxdis+1,maxdis+1)

    if ( ndim > maxdis+1 ) then
       write(*,*) "ndim=",ndim," maxdis+1=",maxdis+1
       call stop_program( "MDIM GREATER THAN MAXDIS+1" )
    end if

    lw=size( scr1 )
    call dgelss(ndim,ndim,1,b,ldb,v,ldb,scr1,toleig,nr,scr2,lw,info)

    if ( info /= 0 ) then
       call stop_program( "COULD NOT SOLVE DIIS EQUATION(least_square)" )
    end if

  END SUBROUTINE least_square


  SUBROUTINE calc_constraint( tau0, anorm, fc )

    implicit none
    real(8),intent(IN)  :: tau0(3,Natom)
    real(8),intent(OUT) :: anorm(:,:), fc(:)
    real(8) :: x1(3),x2(3),x3(3),dx(9)
    real(8),allocatable :: msk(:,:,:)
    integer :: kmax,k,n,i,j,iaa,ibb,icc
    integer :: lskptr(3,Natom)

    anorm(:,:)=0.0d0

    if ( ityp /= 1 .and. ityp /= 2 .and. ityp /= 3 ) then
       call stop_program( &
            "THIS TYPE OF CONSTRAINT HAS NOT BEEN IMPLIMENTED YET !!!" )
    end if

    select case( ityp )
    case( 1 ) ! constraints on distances

       kmax=6

       allocate( msk(nodim,mcnstr,kmax) ) ; msk=0.0d0

       do i=1,mcnstr

          x1(1:3)=tau0(1:3,ia(i))
          x2(1:3)=tau0(1:3,ib(i))
          call pbc_correction_distance( cval(i),x1,x2 )
          call constraint_func_distance( fc(i),fv(i),cval(i),x1,x2 )
          call diff_func_distance( dx,x1,x2 ) ! derivative of fc

          iaa=(ia(i)-1)*3
          ibb=(ib(i)-1)*3
          msk(:,:,:)    =0.0d0 ! mask for relevant dof
          msk(iaa+1,i,1)=1.0d0
          msk(iaa+2,i,2)=1.0d0
          msk(iaa+3,i,3)=1.0d0
          msk(ibb+1,i,4)=1.0d0
          msk(ibb+2,i,5)=1.0d0
          msk(ibb+3,i,6)=1.0d0
          do k=1,kmax
             do n=1,nodim
                anorm(n,i) = anorm(n,i) + dx(k) * msk(n,i,k)  
             end do
          end do

       end do ! i

       deallocate( msk )

    case( 2 ) ! constraints on angles

       kmax=9

       allocate( msk(nodim,mcnstr,kmax) ) ; msk=0.0d0

       do i=1,mcnstr

          x1(:) = tau0(:,ia(i))
          x2(:) = tau0(:,ib(i))
          x3(:) = tau0(:,ic(i))

          ! PBC correction may be necessary 
          ! call pbc_correction_angle( ... )
          call constraint_func_angle(fc(i),fv(i),cval(i),x1,x2,x3)
          call diff_func_angle(dx,x1,x2,x3)

          do k=1,kmax
             do n=1,nodim
                iaa=(ia(i)-1)*3
                ibb=(ib(i)-1)*3
                icc=(ic(i)-1)*3
                msk(:,:,:)    =0.0d0
                msk(iaa+1,i,1)=1.0d0
                msk(iaa+2,i,2)=1.0d0
                msk(iaa+3,i,3)=1.0d0
                msk(ibb+1,i,4)=1.0d0
                msk(ibb+2,i,5)=1.0d0
                msk(ibb+3,i,6)=1.0d0
                msk(icc+1,i,7)=1.0d0
                msk(icc+2,i,8)=1.0d0
                msk(icc+3,i,9)=1.0d0
                anorm(n,i)=anorm(n,i)+dx(k)*msk(n,i,k)
             end do
          end do
       end do
       deallocate( msk )

    case( 3 )

       lskptr(:,:)=0
       do i=1,Natom
          do k=1,3
             lskptr(k,i)=i*3-(3-k)
          end do
       end do

       do i=1,mcnstr
          call constraint_func_diff_Ncoordination &
               ( ia(i), cnpar(1,i), cnpar(2,i), tau0 &
                ,anorm(:,i), lskptr, fc(i), fv(i), cval(i) )
       end do

    end select

  END SUBROUTINE calc_constraint


  SUBROUTINE constraint_func_distance( fd, d, d0, x, y )
    implicit none
    real(8),intent(IN)  :: d0, x(3), y(3)
    real(8),intent(OUT) :: d, fd
    d  = sqrt( sum( (x-y)**2 ) )
    fd = d - d0
  END SUBROUTINE constraint_func_distance


  SUBROUTINE constraint_func_angle( ft, t, t0, x,y,z )
    implicit none
    real(8) :: ft, t, t0, x(3), y(3), z(3)
    real(8),parameter :: epsilon=1.d-12
    real(8) :: a(3), b(3), caa, cab, cbb, ac
    a(1) = x(1)-y(1)
    a(2) = x(2)-y(2)
    a(3) = x(3)-y(3)
    b(1) = y(1)-z(1)
    b(2) = y(2)-z(2)
    b(3) = y(3)-z(3)
    caa = sum( a*a )
    cab = sum( a*b )
    cbb = sum( b*b )
    if ( caa < epsilon .or. cbb < epsilon ) then
       ft=0.0d0
       t =t0
    else
       ac = cab/sqrt(caa*cbb)
       if ( ac < -1.0d0 ) ac=-1.0d0
       if ( ac >  1.0d0 ) ac= 1.0d0
       ft = ac+cos(t0)
       t = acos(-ac)
    end if
    return
  END SUBROUTINE constraint_func_angle


  SUBROUTINE diff_func_distance( dr,r1,r2 )
    implicit none
    real(8),intent(OUT) :: dr(6)
    real(8),intent(IN) :: r1(3),r2(3)
    real(8),parameter :: epsilon=1.d-12
    real(8) :: r,x,y,z
    x = r1(1) - r2(1)
    y = r1(2) - r2(2)
    z = r1(3) - r2(3)
    r = sqrt( x*x + y*y + z*z )
    if ( r < epsilon ) then
       dr(1) = 0.0d0
       dr(2) = 0.0d0
       dr(3) = 0.0d0
       dr(4) = 0.0d0
       dr(5) = 0.0d0
       dr(6) = 0.0d0
    else
       dr(1) =  x/r
       dr(2) =  y/r
       dr(3) =  z/r
       dr(4) = -x/r
       dr(5) = -y/r
       dr(6) = -z/r
    end if
  END SUBROUTINE diff_func_distance


  SUBROUTINE diff_func_angle( dt,x,y,z )
    implicit none
    real(8) :: dt(9), x(3), y(3), z(3)
    real(8),parameter :: epsilon=1.d-12
    real(8) :: caa, cab, cbb, ccc, a(3), b(3)
    a(1) = x(1)-y(1)
    a(2) = x(2)-y(2)
    a(3) = x(3)-y(3)
    b(1) = y(1)-z(1)
    b(2) = y(2)-z(2)
    b(3) = y(3)-z(3)
    caa = sum( a*a )
    cab = sum( a*b )
    cbb = sum( b*b )
    if ( caa < epsilon .or. cbb < epsilon ) then
       dt(:)=0.0d0
    else
       ccc = -1.0d0/sqrt(caa*cbb)
       dt(1) =  ccc*( cab/caa*a(1) - b(1) )
       dt(2) =  ccc*( cab/caa*a(2) - b(2) )
       dt(3) =  ccc*( cab/caa*a(3) - b(3) )
       dt(4) = -ccc*( cab/caa*a(1) - cab/cbb*b(1) + a(1) - b(1) )
       dt(5) = -ccc*( cab/caa*a(2) - cab/cbb*b(2) + a(2) - b(2) )
       dt(6) = -ccc*( cab/caa*a(3) - cab/cbb*b(3) + a(3) - b(3) )
       dt(7) = -ccc*( cab/cbb*b(1) - a(1) )
       dt(8) = -ccc*( cab/cbb*b(2) - a(2) )
       dt(9) = -ccc*( cab/cbb*b(3) - a(3) )
    end if
  END SUBROUTINE diff_func_angle


  SUBROUTINE pbc_correction_distance( cnstr_dist, x, y )
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
  END SUBROUTINE pbc_correction_distance


  SUBROUTINE constraint_func_diff_Ncoordination &
       ( ia_dum, kappa, rc, tscr, an, lsk, fci, fvi, cval_dum )

    implicit none
    integer,intent(IN)  :: ia_dum, lsk(:,:)
    real(8),intent(IN)  :: kappa, rc, cval_dum
    real(8),intent(IN)  :: tscr(:,:)
    real(8),intent(OUT) :: an(:),fci,fvi
    real(8) :: x0,y0,z0,dx,dy,dz,dd,ff,df,fff
    integer :: iat,i,k,l1,l2,l3
    real(8),save :: cval_save

    if ( ia_dum <= Natom ) then
       l1 = lsk(1,ia_dum)
       l2 = lsk(2,ia_dum)
       l3 = lsk(3,ia_dum)
    end if

    x0  = tscr(1,ia_dum)
    y0  = tscr(2,ia_dum)
    z0  = tscr(3,ia_dum)
    fvi = 0.0d0
    dd  = 0.0d0
    df  = 0.0d0
    ff  = 0.0d0
    an  = 0.0d0

    loop_iat: do iat=1,Natom

       do i=1,10
          if ( ki_atom(iat) == index(i) ) cycle loop_iat
       end do

       call pbc_distance( tscr(:,iat),x0,y0,z0, dx,dy,dz, dd )
       dd = sqrt(dd)

       if ( dd > 1.d-2 ) then

          df  = kappa*(dd-rc)
          fvi = fvi + 1.0d0/( exp(df) + 1.0d0 )
          ff  =-0.5d0*kappa/(cosh(df)+1.0d0)/dd

          k=lsk(1,iat)
          if ( k /= 0 ) an(k) = an(k) + ff*dx
          k=lsk(2,iat)
          if ( k /= 0 ) an(k) = an(k) + ff*dy
          k=lsk(3,iat)
          if ( k /= 0 ) an(k) = an(k) + ff*dz

          if ( l1 /= 0 ) an(l1) = an(l1) - ff*dx
          if ( l2 /= 0 ) an(l2) = an(l2) - ff*dy
          if ( l3 /= 0 ) an(l3) = an(l3) - ff*dz

       end if

    end do loop_iat

    if ( cval_dum == 999.9d0 ) then
       if ( itime == 1 ) cval_save=fvi
       fci = fvi - cval_save
    else
       fci = fvi - cval_dum
    end if

  END SUBROUTINE constraint_func_diff_Ncoordination


  SUBROUTINE pbc_distance(target,x0,y0,z0,dx,dy,dz,dd)
    implicit none
    real(8),intent(IN) :: target(3),x0,y0,z0
    real(8),intent(OUT) :: dd,dx,dy,dz
    integer :: i,fact1(27),fact2(27),fact3(27)
    real(8) :: x,y,z,dd_try
!           1 2 3 4 5 6 7 8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27      
    fact1=(/0,1,0,0,1,1,0,1,-1, 0, 0,-1, 1, 1,-1, 0, 0,-1,-1, 0, 1, 1,-1, 1,-1,-1,-1/)        
    fact2=(/0,0,1,0,1,0,1,1, 0,-1, 0, 1,-1, 0, 0,-1, 1,-1, 0,-1, 1,-1, 1,-1, 1,-1,-1/)        
    fact3=(/0,0,0,1,0,1,1,1, 0, 0,-1, 0, 0,-1, 1, 1,-1, 0,-1,-1,-1, 1, 1,-1,-1, 1,-1/)    

    dd=1.d100
    do i=1,27            
       x=target(1)-x0+fact1(i)*aa(1,1)+fact2(i)*aa(1,2)+fact3(i)*aa(1,3)
       y=target(2)-y0+fact1(i)*aa(2,1)+fact2(i)*aa(2,2)+fact3(i)*aa(2,3)
       z=target(3)-z0+fact1(i)*aa(3,1)+fact2(i)*aa(3,2)+fact3(i)*aa(3,3)
       dd_try = x*x + y*y + z*z
       if ( dd_try < dd ) then
          dx=x
          dy=y
          dz=z
          dd=dd_try
       end if
    end do

  END SUBROUTINE pbc_distance


  SUBROUTINE read_blue
    implicit none
    integer :: i,ispecies

    index(:)=0

    if ( myrank == 0 ) then

       open(888,file='blue_control.dat')
       read(888,*) ityp
       read(888,*) mcnstr

       allocate( ia(mcnstr)      ) ; ia=0
       allocate( ib(mcnstr)      ) ; ib=0
       allocate( ic(mcnstr)      ) ; ic=0
       allocate( cval(mcnstr)    ) ; cval=0.0d0
       allocate( cnpar(2,mcnstr) ) ; cnpar=0.0d0

       select case( ityp )
       case( 1 )

          read(888,*) (ia(i),i=1,mcnstr)
          read(888,*) (ib(i),i=1,mcnstr)
          read(888,*) (cval(i),i=1,mcnstr)

       case( 2 )

          read(888,*) (ia(i),i=1,mcnstr)
          read(888,*) (ib(i),i=1,mcnstr)
          read(888,*) (ic(i),i=1,mcnstr)
          read(888,*) (cval(i),i=1,mcnstr)

       case( 3 )

          read(888,*) ispecies
          read(888,*) (ia(i)     ,i=1,mcnstr)
          read(888,*) (cnpar(1,i),i=1,mcnstr)  !c_kappa
          read(888,*) (cnpar(2,i),i=1,mcnstr)  !c_rc
          read(888,*) (cval(i),i=1,mcnstr)
          read(888,*) (index(i),i=1,ispecies)
       end select

       close(888)

    end if

    call bcast_blue_data

  END SUBROUTINE read_blue


  SUBROUTINE bcast_blue_data 
    implicit none
    integer :: ierr

    call mpi_bcast(ityp  ,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(mcnstr,1,mpi_integer,0,mpi_comm_world,ierr)

    if ( .not.allocated(ia) ) then
       allocate( ia(mcnstr)      ) ; ia=0
       allocate( ib(mcnstr)      ) ; ib=0
       allocate( ic(mcnstr)      ) ; ic=0
       allocate( cval(mcnstr)    ) ; cval=0.0d0
       allocate( cnpar(2,mcnstr) ) ; cnpar=0.0d0
    end if

! ---

!   Bcast and allocation

    select case( ityp )
    case( 1 )

       call mpi_bcast(ia,mcnstr,mpi_integer,0,mpi_comm_world,ierr)
       call mpi_bcast(ib,mcnstr,mpi_integer,0,mpi_comm_world,ierr)
       call mpi_bcast(cval,mcnstr,mpi_real8,0,mpi_comm_world,ierr)

    case( 2 )

       call mpi_bcast(ia,mcnstr,mpi_integer,0,mpi_comm_world,ierr)
       call mpi_bcast(ib,mcnstr,mpi_integer,0,mpi_comm_world,ierr)
       call mpi_bcast(ic,mcnstr,mpi_integer,0,mpi_comm_world,ierr)
       call mpi_bcast(cval,mcnstr,mpi_real8,0,mpi_comm_world,ierr)

    case( 3 )

       call mpi_bcast(ia,mcnstr,mpi_integer,0,mpi_comm_world,ierr)
       call mpi_bcast(index,size(index),mpi_integer,0,mpi_comm_world,ierr)
       call mpi_bcast(cnpar,size(cnpar),mpi_real8,0,mpi_comm_world,ierr)
       call mpi_bcast(cval,mcnstr,mpi_real8,0,mpi_comm_world,ierr)

    end select

    allocate( xlagr(mcnstr)       ) ; xlagr=0.0d0
    allocate( ylagr(mcnstr)       ) ; ylagr=0.0d0
    allocate( fc(mcnstr)          ) ; fc=0.0d0
    allocate( fv(mcnstr)          ) ; fc=0.0d0
    allocate( ipvt(mcnstr)        ) ; ipvt=0
    nodim=3*Natom
    allocate( anorm(nodim,mcnstr) ) ; anorm=0.0d0

    return
  END SUBROUTINE bcast_blue_data


  SUBROUTINE dealloc_blue
    implicit none
    deallocate( ylagr  )
    deallocate( xlagr  )
    deallocate( ipvt   )
    deallocate( dtm3 )
    deallocate( dtm1 )
    deallocate( Rion0  )
    return
  END SUBROUTINE dealloc_blue


  SUBROUTINE write_blue_data( nfi, flag_io )
    implicit none
    integer,intent(IN) :: nfi
    logical,intent(IN) :: flag_io
    integer :: j
    logical,save :: init=.true.
    if ( flag_io ) then
       if ( init ) then
          open(889,file='CONSTRAINT',status="replace")
          init=.false.
       else
          open(889,file='CONSTRAINT',position="append")
       end if
       do j=1,mcnstr
          write(889,'(i7,2x,i4,5x,2(1pe20.10))') nfi,j,xlagr(j),fv(j)
       end do
       close(889)
    end if
    return
  END SUBROUTINE write_blue_data


END MODULE blue_moon_module



