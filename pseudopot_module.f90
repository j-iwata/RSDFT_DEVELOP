MODULE pseudopot_module

  use ps_read_YB_module

  implicit none

  PRIVATE
  PUBLIC :: pselect,ippform,file_ps,inorm,NRps,norb,Npseudopot &
           ,Mr,lo,vql,cdd,cdc,rad,anorm,viod,Rps,Zps,parloc,rab &
           ,read_ppname_pseudopot,read_pseudopot &
           ,cdd_coef,read_param_pseudopot,read_param_oldformat_pseudopot &
           ,read_ppname_oldformat_pseudopot

  integer :: pselect,Npseudopot
  integer,allocatable :: ippform(:)
  character(30),allocatable :: file_ps(:)
  integer,allocatable :: inorm(:,:),NRps(:,:),norb(:),Mr(:),lo(:,:)
  real(8),allocatable :: vql(:,:),cdd(:,:),cdc(:,:),rad(:,:),parloc(:,:)
  real(8),allocatable :: anorm(:,:),viod(:,:,:),Rps(:,:),Zps(:),rab(:,:)
  real(8),allocatable :: cdd_coef(:,:,:)
  integer :: unit_ps,ielm,Nelement_PP
  integer :: max_psgrd=0, max_psorb=0, max_ngauss=0

CONTAINS


  SUBROUTINE read_ppname_pseudopot(MKI,rank,unit)
    implicit none
    integer,intent(IN) :: MKI,rank,unit
    integer :: i,j,ierr
    character(2) :: cbuf,ckey
    Nelement_PP=0
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=990) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:2) == "PP" ) Nelement_PP=Nelement_PP+1
       end do
990    continue
    end if
    call send_ppname_1(0)
    allocate( ippform(Nelement_PP),file_ps(Nelement_PP) )
    ippform(:)=0
    file_ps(:)=""
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:2) == "PP" ) then
             backspace(unit)
             do j=1,Nelement_PP
                read(unit,*) cbuf,ippform(j),file_ps(j)
             end do ! j
             exit
          end if
       end do
999    continue
       write(*,*) "Nelement_PP, Nelement=",Nelement_PP,MKI
       do i=1,Nelement_PP
          write(*,'(1x,"ippform, file_ps = ",i3,2x,a30,3x,3f10.5)') &
               ippform(i),file_ps(i)
       end do
       if ( Nelement_PP > MKI ) then
          write(*,'("WARNING(Nelement_PP>Nelement):First ",i1," potentials are used")') MKI
       end if
    end if
    call send_ppname_2(0)
    ierr=0
    if ( Nelement_PP < MKI ) ierr=-1
    do i=1,Nelement_PP
       if ( ippform(i) <= 0 .or. file_ps(i)=="" ) ierr=-1
    end do
    if ( ierr == -1 ) stop "stop@read_ppname_pseudopot"
    Nelement_PP=MKI
    Npseudopot =MKI
  END SUBROUTINE read_ppname_pseudopot

  SUBROUTINE read_ppname_oldformat_pseudopot(MKI,rank,unit)
    implicit none
    integer,intent(IN) :: MKI,rank,unit
    integer :: i
    Nelement_PP=MKI
    Npseudopot =MKI
    allocate( ippform(Nelement_PP),file_ps(Nelement_PP) )
    if ( rank == 0 ) then
       do i=1,Nelement_PP
          read(unit,*) ippform(i),file_ps(i)
       end do
       do i=1,Nelement_PP
          write(*,'(1x,"ippform, file_ps = ",i3,2x,a30,3x,3f10.5)') &
               ippform(i),file_ps(i)
       end do
    end if
    call send_ppname_2(0)
  END SUBROUTINE read_ppname_oldformat_pseudopot

  SUBROUTINE send_ppname_1(rank)
    implicit none
    integer,intent(IN) :: rank
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(Nelement_PP,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)    
  END SUBROUTINE send_ppname_1

  SUBROUTINE send_ppname_2(rank)
    implicit none
    integer,intent(IN) :: rank
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(file_ps,30*Nelement_PP,MPI_CHARACTER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(ippform,Nelement_PP,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_ppname_2


  SUBROUTINE read_param_pseudopot(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    character(7) :: cbuf,ckey
    integer :: i
    pselect=2
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:7) == "PSELECT" ) then
             backspace(unit)
             read(unit,*) cbuf,pselect
          end if
       end do
999    continue
       write(*,*) "pslect=",pselect
    end if
    call send_param_pseudopot(0)
  END SUBROUTINE read_param_pseudopot

  SUBROUTINE read_param_oldformat_pseudopot(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    if ( rank == 0 ) then
       read(unit,*) pselect
       write(*,*) "pselect=",pselect
    end if
    call send_param_pseudopot(0)
  END SUBROUTINE read_param_oldformat_pseudopot

  SUBROUTINE send_param_pseudopot(rank)
    implicit none
    integer,intent(IN) :: rank
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(pselect,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_param_pseudopot


  SUBROUTINE read_pseudopot(rank)
    implicit none
    integer,intent(IN) :: rank
    if ( rank == 0 ) then
       max_psgrd=0
       max_psorb=0
       do ielm=1,Nelement_PP
          unit_ps=33+ielm
          select case(ippform(ielm))
          case(1)
             call read_TM
          case(2)
             open(unit_ps,FILE=file_ps(ielm),STATUS='old')
             call read_PSV
             close(unit_ps)
          case(3)
             open(unit_ps,FILE=file_ps(ielm),STATUS='old')
             call ps_read_YB(unit_ps)
             close(unit_ps)
             call ps_allocate(ps_yb%nrr,ps_yb%norb)
             Mr(ielm)                 = ps_yb%nrr
             norb(ielm)               = ps_yb%norb
             Zps(ielm)                = ps_yb%znuc
             anorm(1:norb(ielm),ielm) = ps_yb%anorm(1:norb(ielm))
             inorm(1:norb(ielm),ielm) = ps_yb%inorm(1:norb(ielm))
             Rps(1:norb(ielm),ielm)   = ps_yb%Rc(1:norb(ielm))
             NRps(1:norb(ielm),ielm)  = ps_yb%NRc(1:norb(ielm))
             lo(1:norb(ielm),ielm)    = ps_yb%lo(1:norb(ielm))
             vql(1:Mr(ielm),ielm)     = ps_yb%vql(1:Mr(ielm))
             cdd(1:Mr(ielm),ielm)     = ps_yb%cdd(1:Mr(ielm))
             cdc(1:Mr(ielm),ielm)     = ps_yb%cdc(1:Mr(ielm))
             rad(1:Mr(ielm),ielm)     = ps_yb%rr(1:Mr(ielm))
             rab(1:Mr(ielm),ielm)     = ps_yb%rx(1:Mr(ielm))
             viod(1:Mr(ielm),1:norb(ielm),ielm) &
                                      = ps_yb%vps(1:Mr(ielm),1:norb(ielm))
          case default
             stop "ippform error"
          end select
       end do
    end if
    call send_pseudopot(rank)
  END SUBROUTINE read_pseudopot


  SUBROUTINE send_pseudopot(myrank)
    implicit none
    integer,intent(IN) :: myrank
    integer :: m,n,ierr
    include 'mpif.h'
    m=max_psgrd
    n=max_psorb
    call mpi_bcast(m,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Nelement_PP,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if ( myrank /= 0 ) then
       call ps_allocate(m,n)
    end if
    call mpi_bcast(Mr    ,Nelement_PP,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(norb  ,Nelement_PP,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Zps   ,Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(parloc,4*Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(anorm ,n*Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(inorm ,n*Nelement_PP,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Rps   ,n*Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(NRps  ,n*Nelement_PP,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(lo    ,n*Nelement_PP,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(vql   ,m*Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(cdd   ,m*Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(cdc   ,m*Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(rad   ,m*Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(rab   ,m*Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(viod  ,m*n*Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
!
    call mpi_bcast(max_ngauss,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if ( max_ngauss /= 0 ) then
       if ( myrank /= 0 ) then
          allocate( cdd_coef(3,max_ngauss,Nelement_PP) ) ; cdd_coef=0.0d0
       end if
       call mpi_bcast(cdd_coef,3*max_ngauss*Nelement_PP,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    end if
  END SUBROUTINE send_pseudopot


  SUBROUTINE read_TM
    implicit none
    stop
  END SUBROUTINE read_TM


  SUBROUTINE read_PSV
    implicit none
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
    real(8),allocatable :: psi(:,:,:),phi(:,:,:),bet(:,:,:)
    real(8),allocatable :: ddi(:,:,:),qqr(:,:,:)
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
       do j=1,nr(l)
          iorb=iorb+1
          lo(iorb,ielm)=l-1
          anorm(iorb,ielm)=abs( ddi(j,j,l) )
          inorm(iorb,ielm)=sign( 1.d0, ddi(j,j,l) )
          NRps(iorb,ielm)=nsmpl+1
          Rps(iorb,ielm)=Rc_in
       end do
    end do
    norb(ielm)=iorb
    Mr(ielm)=max( ndlc,nsmpl,ndata )+1
    Zps(ielm)=znuc

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
          temp=sqrt( anorm(iorb,ielm) )
          do i=1,nsmpl
             viod(i+1,iorb,ielm)=bet(i,j,l)*temp
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
    deallocate( qqr,ddi )
    deallocate( bet,psi,phi )
    deallocate( r )
    deallocate( cc,vl,rr )
    deallocate( nr )

  END SUBROUTINE read_PSV


  SUBROUTINE ps_allocate(n_grd,n_orb)
    implicit none
    integer,intent(IN) :: n_grd,n_orb
    integer :: mg,mo
    real(8),allocatable :: vql_tmp(:,:),cdd_tmp(:,:),rad_tmp(:,:)
    real(8),allocatable :: cdc_tmp(:,:),viod_tmp(:,:,:)
    real(8),allocatable :: anorm_tmp(:,:),Rps_tmp(:,:),rab_tmp(:,:)
    integer,allocatable :: inorm_tmp(:,:),lo_tmp(:,:),NRps_tmp(:,:)
    if ( max_psgrd==0 .or. max_psorb==0 ) then
       allocate( Mr(Nelement_PP)   ) ; Mr=0
       allocate( norb(Nelement_PP) ) ; norb=0
       allocate( Zps(Nelement_PP)  ) ; Zps=0.d0
       allocate( parloc(4,Nelement_PP)    ) ; parloc=0.d0
       allocate( anorm(n_orb,Nelement_PP) ) ; anorm=0.d0
       allocate( inorm(n_orb,Nelement_PP) ) ; inorm=0
       allocate( Rps(n_orb,Nelement_PP)   ) ; Rps=0.d0
       allocate( NRps(n_orb,Nelement_PP)  ) ; NRps=0
       allocate( lo(n_orb,Nelement_PP)    ) ; lo=0
       allocate( vql(n_grd,Nelement_PP)   ) ; vql=0.d0
       allocate( cdd(n_grd,Nelement_PP)   ) ; cdd=0.d0
       allocate( cdc(n_grd,Nelement_PP)   ) ; cdc=0.d0
       allocate( rad(n_grd,Nelement_PP)   ) ; rad=0.d0
       allocate( rab(n_grd,Nelement_PP)   ) ; rab=0.d0
       allocate( viod(n_grd,n_orb,Nelement_PP) ) ; viod=0.d0
       max_psgrd=n_grd
       max_psorb=n_orb
       return
    end if
    mg = max( max_psgrd, n_grd )
    mo = max( max_psorb, n_orb )
    if ( max_psgrd < mg ) then
       allocate( vql_tmp(mg,Nelement_PP) )
       allocate( cdd_tmp(mg,Nelement_PP) )
       allocate( rad_tmp(mg,Nelement_PP) )
       allocate( rab_tmp(mg,Nelement_PP) )
       allocate( cdc_tmp(mg,Nelement_PP) )
       vql_tmp(1:max_psgrd,1:Nelement_PP) = vql(1:max_psgrd,1:Nelement_PP)
       cdd_tmp(1:max_psgrd,1:Nelement_PP) = cdd(1:max_psgrd,1:Nelement_PP)
       rad_tmp(1:max_psgrd,1:Nelement_PP) = rad(1:max_psgrd,1:Nelement_PP)
       rab_tmp(1:max_psgrd,1:Nelement_PP) = rab(1:max_psgrd,1:Nelement_PP)
       cdc_tmp(1:max_psgrd,1:Nelement_PP) = cdc(1:max_psgrd,1:Nelement_PP)
       deallocate( cdc )
       deallocate( rab )
       deallocate( rad )
       deallocate( cdd )
       deallocate( vql )
       allocate( vql(mg,Nelement_PP) ) ; vql=0.d0
       allocate( cdd(mg,Nelement_PP) ) ; cdd=0.d0
       allocate( rad(mg,Nelement_PP) ) ; rad=0.d0
       allocate( rab(mg,Nelement_PP) ) ; rab=0.d0
       allocate( cdc(mg,Nelement_PP) ) ; cdc=0.d0
       vql(:,:)=vql_tmp(:,:)
       cdd(:,:)=cdd_tmp(:,:)
       rad(:,:)=rad_tmp(:,:)
       rab(:,:)=rab_tmp(:,:)
       cdc(:,:)=cdc_tmp(:,:)
       deallocate( cdc_tmp )
       deallocate( rab_tmp )
       deallocate( rad_tmp )
       deallocate( cdd_tmp )
       deallocate( vql_tmp )
       allocate( viod_tmp(mg,mo,Nelement_PP) )
       viod_tmp(1:max_psgrd,1:max_psorb,1:Nelement_PP) &
            = viod(1:max_psgrd,1:max_psorb,1:Nelement_PP)
       deallocate( viod )
       allocate( viod(mg,mo,Nelement_PP) ) ; viod=0.d0
       viod(:,:,:)=viod_tmp(:,:,:)
       deallocate( viod_tmp )
    end if
    if ( max_psorb < mo ) then
       allocate( anorm_tmp(mo,Nelement_PP) )
       allocate( inorm_tmp(mo,Nelement_PP) )
       allocate( lo_tmp(mo,Nelement_PP) )
       allocate( Rps_tmp(mo,Nelement_PP) )
       allocate( NRps_tmp(mo,Nelement_PP) )
       anorm_tmp(1:max_psorb,1:Nelement_PP) = anorm(1:max_psorb,1:Nelement_PP)
       inorm_tmp(1:max_psorb,1:Nelement_PP) = inorm(1:max_psorb,1:Nelement_PP)
       lo_tmp(1:max_psorb,1:Nelement_PP) = lo(1:max_psorb,1:Nelement_PP)
       Rps_tmp(1:max_psorb,1:Nelement_PP) = Rps(1:max_psorb,1:Nelement_PP)
       NRps_tmp(1:max_psorb,1:Nelement_PP) = NRps(1:max_psorb,1:Nelement_PP)
       deallocate( NRps )
       deallocate( Rps )
       deallocate( lo )
       deallocate( inorm )
       deallocate( anorm )
       allocate( anorm(mo,Nelement_PP) ) ; anorm=0.d0
       allocate( inorm(mo,Nelement_PP) ) ; inorm=0
       allocate( lo(mo,Nelement_PP)    ) ; lo=0
       allocate( Rps(mo,Nelement_PP)   ) ; Rps=0.d0
       allocate( NRps(mo,Nelement_PP)  ) ; NRps=0
       anorm(:,:) = anorm_tmp(:,:)
       inorm(:,:) = inorm_tmp(:,:)
       lo(:,:) = lo_tmp(:,:)
       Rps(:,:) = Rps_tmp(:,:)
       NRps(:,:) = NRps_tmp(:,:)
       deallocate( NRps_tmp )
       deallocate( Rps_tmp )
       deallocate( lo_tmp )
       deallocate( inorm_tmp )
       deallocate( anorm_tmp )
       if ( max_psgrd >= mg ) then
          allocate( viod_tmp(mg,mo,Nelement_PP) )
          viod_tmp(1:max_psgrd,1:max_psorb,1:Nelement_PP) &
               = viod(1:max_psgrd,1:max_psorb,1:Nelement_PP)
          deallocate( viod )
          allocate( viod(mg,mo,Nelement_PP) ) ; viod=0.d0
          viod(:,:,:)=viod_tmp(:,:,:)
          deallocate( viod_tmp )
       end if
    end if
    max_psgrd = mg
    max_psorb = mo
  END SUBROUTINE ps_allocate


END MODULE pseudopot_module
