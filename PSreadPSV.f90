MODULE PSreadPSV
  use VarPSMember
  implicit none

CONTAINS
  
!-------------------------------------------------------
  SUBROUTINE read_PSV( unit_ps,ielm,ddi,qqr,psi,phi,bet )
    implicit none
    integer,intent(IN) :: unit_ps,ielm
    character(30) :: file_name
    real(8) :: znuc
    integer :: nsmpl,ndlc,ndlc_1,ndata
    integer :: nl,l,i,j,m
    integer,allocatable :: nr(:)

    integer :: iost,iorb
    real(8) :: temp,h
    character(18) :: inbuf18
    character(80) :: inbuf

    integer :: ifchrg,verpot,ngauss
    real(8) :: zatom,Rc_in,a1,a2,a3,a4
    real(8) :: r1,r2,rend,dummy
    character(8) :: xc_pot

    real(8),allocatable :: rr(:),vl(:),cc(:),r(:)
    real(8),allocatable,intent(OUT) :: psi(:,:,:),phi(:,:,:),bet(:,:,:)
    real(8),allocatable,intent(OUT) :: ddi(:,:,:),qqr(:,:,:)
    real(8),allocatable :: a0(:),b0(:),c0(:),cdd_coef_0(:,:,:)

! Check pseudopotential type

    read(unit_ps,'(A)') inbuf

    if ( index(inbuf,'#PSV')==1 ) then

! Ver. 2.0 or more
!
       verpot=2
       read(unit_ps,*) xc_pot, zatom, ifchrg

! COMMENT LINE SKIP('#....')
!
       do i=1,10000
          read(unit_ps,'(A)') inbuf
          if ( index(inbuf,'#')==1 ) cycle
          backspace unit_ps
          exit
       end do
       if ( i>10000 ) stop "stop at psv"

       read(unit_ps,*) znuc

    else if ( index(inbuf,'#xPSV')==1 ) then
! xTAPP
!
       verpot=3
       ifchrg=0
       do i=1,10000
          read(unit_ps,'(A)') inbuf
          if ( index(inbuf,'#')==1 ) cycle
          backspace unit_ps
          exit
       end do
       if ( i>10000 ) stop "stop at psv"

       read(unit_ps,*) xc_pot
       read(unit_ps,*) znuc,zatom

    else

! Ver. 1.0 
!
       verpot=1
       ifchrg=0
       znuc=0
       write(6,'("WARNING! THE POTENTIAL TYPE IS OLD."/ &
                ,"XC-POTENTIAL TYPE: LDAPZ81 IS ASSUMED"/   &
                ,"We recommend rewrite or remake this potential data."/)')

       xc_pot='LDAPZ81'

! COMMENT LINE SKIP
!
       read(unit_ps,*)

    end if

    read(unit_ps,*) nl
    allocate( nr(nl) )
    read(unit_ps,*) nr(1:nl)
    read(unit_ps,*) nsmpl,Rc_in
    read(unit_ps,*) a1,a2,a3,a4

    temp=abs(a1+a3-1.0d0)
    if ( temp>1.0d-10 ) then
       write(6,*) ' LOCAL PARAMETER FITTING ERROR'
       write(6,*) 'A1,A2,A3,A4',a1,a2,a3,a4
       stop
    end if

! Read local potential
!
    read(unit_ps,*) ndlc

    allocate( rr(ndlc),vl(ndlc),cc(ndlc) )

    read(unit_ps,*,IOSTAT=iost) rr(1),vl(1),cc(1)
    backspace unit_ps

    read(unit_ps,*,IOSTAT=iost) rr(2),vl(2),cc(2)
    if ( rr(1)==rr(2) ) then
       iost=0
    else
       iost=1
       backspace unit_ps
       backspace unit_ps
    end if

    if ( iost==0 ) then
       do i=2,ndlc
          read(unit_ps,*) rr(i),vl(i),cc(i)
       end do
    else
       cc(:)=0.d0
       do i=2,ndlc
          read(unit_ps,*) rr(i),vl(i)
       end do
    end if

    if ( rr(1)<=0.d0 ) then
       write(6,*) ' ERROR IN LOACAL POT DATA rr(1)=',rr(1)
       stop 'stop at psv_format'
    end if

    h=exp( log(rr(ndlc)/rr(1))/(ndlc-1) )
    temp=h*rr(1)

    if ( abs(temp-rr(2))>1.d-12 ) then
       write(6,*) ' LOCAL POTENTIAL ERROR NDLC=',ndlc, 'H=',h
       write(6,*) ' rr(1),rr(2),r(3)=',rr(1),rr(2),rr(3)
       write(6,*) '          r(ndlc)=',rr(ndlc)
       write(6,*) ' temp=',temp,log(10.d0),log(exp(1.d0))
       stop
    end if

! Read nonlocal potential
!
    read(unit_ps,*)
    read(unit_ps,*) ndata,r1,rend,r2

    allocate( r(ndata) )

! MAKE R(I)
!
    if ( ndata/=nsmpl ) then
       write(6,*)' ERROR ON # OF DATA  NSMPL,NDATA=',nsmpl,ndata
       stop' MKQG0'
    end if
    h=exp( log(rend/r1)/(ndata-1) )
    r(1)=r1
    do i=2,ndata
       r(i)=r(i-1)*h
    end do
    if ( (abs(r(2)-r2)>1.d-6).or.(abs(rr(ndata)-rend)>1.d-8) ) then
       write(6,*)'ERROR  NDATA=',ndata
       write(6,'(" R(1),X1=",2e15.7," R(2),X2=",2e15.7/," R(NDATA),XEND=" &
                 ,2e15.7)') rr(1),r1,rr(2),r2,rr(ndata),rend
       stop
    end if
    temp=sqrt( sum( (rr(1:ndata)-r(1:ndata))**2 )/ndata )
    if ( temp>1.d-7 ) then
       write(*,*) "temp=",temp
       stop "Radial grid is inconsistent."
    end if

    m=maxval(nr(1:nl))
    allocate( psi(nsmpl,m,nl),phi(nsmpl,m,nl),bet(nsmpl,m,nl) )
    allocate( ddi(m,m,nl), qqr(m,m,nl) )

    do l=1,nl
       read(unit_ps,*)
       m=nr(l)
       read(unit_ps,'(3e20.12)') ( (ddi(i,j,l),i=1,j), j=1,m )
       do j=1,m
       do i=j+1,m
          ddi(i,j,l)=ddi(j,i,l)
       end do
       end do
       read(unit_ps,*) ( (qqr(i,j,l),i=1,j),j=1,m )
       do j=1,m
       do i=j+1,m
          qqr(i,j,l)=qqr(j,i,l)
       end do
       end do
       do j=1,m
          read(unit_ps,*)
          do i=1,nsmpl
             read(unit_ps,*) psi(i,j,l),phi(i,j,l),bet(i,j,l)
          end do
       end do
    end do

    select case( verpot )
    case default
       do j=1,10000
          read(unit_ps,'(A)') inbuf18
          if ( inbuf18=='### initial charge' ) then
             read(unit_ps,*) ngauss
             allocate( a0(ngauss),b0(ngauss),c0(ngauss) )
             do i=1,ngauss
                read(unit_ps,*) a0(i),b0(i),c0(i)
             end do
             exit
          end if
       end do
       if ( j>10000 ) stop "read_psv"
    case( 3 )
    end select

    i = max( ndlc,nsmpl,ndata ) + 1
    j = max( 1, sum(nr(1:nl)) )
    call ps_allocate(i,j)

    iorb=0
    do l=1,nl
       nrf(l,ielm)=nr(l)
       do j=1,nr(l)
          iorb=iorb+1
          lo(iorb,ielm)=l-1
          no(iorb,ielm)=j
          anorm(iorb,ielm)=abs( ddi(j,j,l) )
          inorm(iorb,ielm)=sign( 1.d0, ddi(j,j,l) )
          NRps(iorb,ielm)=nsmpl+1
          Rps(iorb,ielm)=Rc_in
       end do
    end do
    norb(ielm)=iorb
    nlf(ielm)=nl
    Mr(ielm)=max( ndlc,nsmpl,ndata )+1
    Zps(ielm)=znuc
    Zelement(ielm)=zatom

! r=0
    do i=1,ndlc
       vql(i+1,ielm)=vl(i)
       cdc(i+1,ielm)=cc(i)
       rad(i+1,ielm)=rr(i)
    end do
    rad(1,ielm)=0.d0
    vql(1,ielm)=vql(2,ielm)

    iorb=0
    do l=1,nl
       do j=1,nr(l)
          iorb=iorb+1
!          temp=sqrt( anorm(iorb,ielm) )
          do i=1,nsmpl
!             viod(i+1,iorb,ielm)=bet(i,j,l)*temp
             viod(i+1,iorb,ielm)=bet(i,j,l)
          end do
          viod(1,iorb,ielm)=0.d0
       end do
    end do

! dr/dx
    h=log( rr(ndlc)/rr(1) )/(ndlc-1)
    rab(1,ielm)=0.d0
    do i=1,ndlc
       rab(i+1,ielm)=rr(i)*h
    end do

    do i=1,ndlc
        rabr2(i,ielm)=rab(i,ielm)*(rr(i)**2)
    end do

    parloc(1,ielm)=a1
    parloc(2,ielm)=a2
    parloc(3,ielm)=a3
    parloc(4,ielm)=a4

!
! initial charge
!
    select case( verpot )
    case default
       temp=16.d0*atan(1.d0)
       cdd(:,ielm)=0.d0
       do l=1,ngauss
          do i=1,ndlc
             cdd(i+1,ielm)=cdd(i+1,ielm) &
                  +rr(i)**2*temp*(a0(l)+b0(l)*rr(i)**2)*exp(-c0(l)*rr(i)**2)
          end do
       end do
    case( 3 )
       file_name=trim(file_ps(ielm))//".ichr"
       open(unit_ps+1,file=file_name,status='old')
       read(unit_ps+1,*) ndlc_1
       if ( ndlc /= ndlc_1 ) stop "ndlc/=ndlc_1!!!"
       temp=16.d0*atan(1.d0)
       cdd(:,ielm)=0.d0
       do i=1,ndlc
          read(unit_ps+1,*) dummy,cdd(i+1,ielm)
          cdd(i+1,ielm)=temp*cdd(i+1,ielm)*rr(i)**2
       end do
       close(unit_ps+1)
    end select

    if ( max_ngauss == 0 ) then
       allocate( cdd_coef(3,ngauss,Nelement_PP) )
       cdd_coef=0.0d0
       max_ngauss=ngauss
    else if ( ngauss > max_ngauss ) then
       allocate( cdd_coef_0(3,max_ngauss,Nelement_PP) )
       cdd_coef_0(:,:,:)=cdd_coef(:,:,:)
       deallocate( cdd_coef )
       allocate( cdd_coef(3,ngauss,Nelement_PP) ) ; cdd_coef=0.0d0
       deallocate( cdd_coef_0 )
       cdd_coef(:,1:max_ngauss,:)=cdd_coef_0(:,:,:)
       max_ngauss=ngauss
    end if
    cdd_coef(1,1:ngauss,ielm)=a0(1:ngauss)
    cdd_coef(2,1:ngauss,ielm)=b0(1:ngauss)
    cdd_coef(3,1:ngauss,ielm)=c0(1:ngauss)

    write(*,*) "*** PSV format ***"
    write(*,*) zatom
    write(*,*) "# of radial mesh points =",Mr(ielm)
    write(*,*) "r(1),...,r(end)         =",rad(1,ielm),rad(Mr(ielm),ielm)
    write(*,*) "NRps                    =",NRps(1:norb(ielm),ielm)
    write(*,*) "# of orbitals           =",norb(ielm)
    write(*,*) "angular momentum        =",lo(1:norb(ielm),ielm)
    write(*,*) "cut off radius          =",Rps(1:norb(ielm),ielm)
    write(*,*) "Zps                     =",Zps(ielm)
    write(*,*) "xc_pot                  =   ",xc_pot
    write(*,*) "a1,a2,a3,a4             =",a1,a2
    write(*,*) "                         ",a3,a4
    write(*,*) "Reference               =",nr(1:nl)
    write(*,*) "anorm                   =",anorm(1:norb(ielm),ielm)
    write(*,*) "inorm                   =",inorm(1:norb(ielm),ielm)
    write(*,*) "sum(cdd)                =",sum( cdd(:,ielm)*rab(:,ielm) )
    write(*,*) "sum(cdc)                =",sum( cdc(:,ielm)*rab(:,ielm) )*temp

    if ( verpot /= 3 ) deallocate( c0,b0,a0 )
!    deallocate( qqr,ddi )
!    deallocate( bet,psi,phi )
    deallocate( r )
    deallocate( cc,vl,rr )
    deallocate( nr )

  END SUBROUTINE read_PSV

END MODULE PSreadPSV
