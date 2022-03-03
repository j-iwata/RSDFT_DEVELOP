program cubegen_vrho
  implicit none

  integer,parameter :: u0=80, u1=1, u2=970, iu=2
  integer :: ML,ML1,ML2,ML3,zatom(9)
  integer :: i1,i2,i3,i,s,MI,MKI,n
  integer :: m1,m2,m3,n1,n2,n3
  integer,allocatable :: LL(:,:),Kion(:)
  real(8),allocatable :: rho(:,:,:,:),vloc(:,:,:,:)
  real(8),allocatable :: vh(:,:,:),vxc(:,:,:),rtmp(:)
  real(8),allocatable :: asi(:,:),rsi(:,:)
  real(8),parameter :: bohr=0.529177d0
  real(8) :: ax,aa(3,3),Va,dV
  character(64) :: cbuf
  logical :: flag_versioned=.false., flag_show_help
  integer :: ialgo_plot = 2  ! 1:cube 2:CHGCAR

! ---

  flag_show_help = .false.
  n = iargc()
  call getarg(1, cbuf)
  select case( trim(adjustl(cbuf)) )
  case( '-h', '--h', '-help', '--help' )
    flag_show_help = .true.
  case( '-f' )
    if ( n < 2 ) then
      flag_show_help = .true.
    else
      call getarg(2, cbuf)
      if ( trim(adjustl(cbuf)) == 'cub' .or. trim(adjustl(cbuf)) == 'cube' ) then
        ialgo_plot = 1
      else if ( trim(adjustl(cbuf)) == 'vasp' ) then
        ialgo_plot = 2
      else
        flag_show_help = .true.
      end if
    end if
  end select
  if ( flag_show_help ) then
    write(*,'(1x,"Usage:",/)')
    write(*,'(1x,"(help)")')
    write(*,'(2x,"$ cubegen_vrho -h",/)')
    write(*,'(1x,"(VASP format)")')
    write(*,'(2x,"$ cubegen_vrho")')
    write(*,'(2x,"$ cubegen_vrho -f vasp",/)')
    write(*,'(1x,"(GAUSSIAN cube format)")')
    write(*,'(2x,"$ cubegen_vrho -f cube")')
    stop
  end if

! ---

  open(u0,file='vrho.dat1',status='old',form='unformatted')

  read(u0) cbuf
  if ( cbuf(1:7) == "version" ) then
    write(*,*) "file format version: "//cbuf
    flag_versioned=.true.
  else
    write(*,*) "non-versioned (old) file format"
    flag_versioned=.false.
    rewind(u0)
 end if

  read(u0) ML,ML1,ML2,ML3
  write(*,*) "ML,ML1,ML2,ML3=",ML,ML1,ML2,ML3

  allocate( LL(3,ML) ) ; LL=0.0d0

  read(u0) LL

  if ( flag_versioned ) then

    read(u0) cbuf ; write(*,*) cbuf
    read(u0) aa
    read(u0) cbuf ; write(*,*) cbuf
    read(u0) MKI,MI
    allocate( asi(3,MI) ) ; asi=0.0d0
    allocate( rsi(3,MI) ) ; rsi=0.0d0
    allocate( Kion(MI)  ) ; Kion=0
    read(u0) asi
    read(u0) Kion
    read(u0) zatom(1:MKI)
    read(u0) cbuf ; write(*,*) cbuf
    read(u0)

  else

    read(u2,*) cbuf,ax
    read(u2,*) cbuf,aa(1:3,1)
    read(u2,*) cbuf,aa(1:3,2)
    read(u2,*) cbuf,aa(1:3,3)
    read(u2,*)
    read(u2,*) MKI,MI,zatom(1:MKI)

    allocate( asi(3,MI) ) ; asi=0.0d0
    allocate( rsi(3,MI) ) ; asi=0.0d0
    allocate( Kion(MI)  ) ; Kion=0

    do i=1,MI
      read(u2,*) Kion(i),asi(1:3,i)
    end do

    aa(:,:)=ax*aa(:,:)

  end if

  Va = aa(1,1)*aa(2,2)*aa(3,3)+aa(1,2)*aa(2,3)*aa(3,1) &
      +aa(1,3)*aa(2,1)*aa(3,2)-aa(1,3)*aa(2,2)*aa(3,1) &
      -aa(1,2)*aa(2,1)*aa(3,3)-aa(1,1)*aa(2,3)*aa(3,2)

! ---

!  do i=1,MI
!     asi(3,i) = asi(3,i) + 0.5d0
!     if ( asi(3,i) > 1.0d0 ) asi(3,i) = asi(3,i) - 1.0d0
!  end do

! ---

  rsi(:,:)=matmul( aa, asi )

! ---

  m1 = minval( LL(1,:) )
  m2 = minval( LL(2,:) )
  m3 = minval( LL(3,:) )
  n1 = maxval( LL(1,:) )
  n2 = maxval( LL(2,:) )
  n3 = maxval( LL(3,:) )

  allocate( rtmp(ML) ) ; rtmp=0.0d0
  allocate( rho(m1:n1,m2:n2,m3:n3,2)  ) ; rho=0.0d0
  allocate( vloc(m1:n1,m2:n2,m3:n3,2) ) ; vloc=0.0d0
  allocate( vh(m1:n1,m2:n2,m3:n3)     ) ; vh=0.0d0
  allocate( vxc(m1:n1,m2:n2,m3:n3)    ) ; vxc=0.0d0

  dV=abs(Va)/dble(ML)

  do s=1,2

    read(u0,END=10) rtmp
    do i=1,ML
      rho(LL(1,i),LL(2,i),LL(3,i),s)=rtmp(i)
    end do
    write(*,*) "for SPIN = ", s, " : sum rho  = ", sum( rho(:,:,:,s) )*dV

    call gen_plot( s, rho(:,:,:,s), "rho" )

    read(u0,END=10) rtmp
    do i=1,ML
      vloc(LL(1,i),LL(2,i),LL(3,i),s)=rtmp(i)
    end do

    call gen_plot( s, vloc(:,:,:,s), "vloc" )

    if ( s == 2 ) then
      rho(:,:,:,1) = rho(:,:,:,1) - rho(:,:,:,2)
      call gen_plot( 1, rho(:,:,:,1), "dspin" )
      write(*,*) "dspin(max,min)=",maxval(rho(:,:,:,1)),minval(rho(:,:,:,1))
      write(*,*) "sum  dspin  = ", sum( rho(:,:,:,1) )*dV
      write(*,*) "sum |dspin| = ", sum( abs(rho(:,:,:,1)) )*dV
    end if

    if ( s == 2 ) exit

    read(u0,END=10) rtmp
    do i=1,ML
      vh(LL(1,i),LL(2,i),LL(3,i))=rtmp(i)
    end do

    read(u0,END=10) rtmp
    do i=1,ML
      vxc(LL(1,i),LL(2,i),LL(3,i))=rtmp(i)
    end do

  end do ! s
10 continue

  close(u0)

contains

  subroutine gen_plot(s, f, fn0 )
    implicit none
    integer,intent(in) :: s
    real(8),intent(in) :: f(0:,0:,0:)
    character(*),intent(in) :: fn0
    select case( ialgo_plot )
    case(1)
      call gen_cube( s, f, fn0 ) 
    case(2)
      call gen_chgcar( s, f, fn0 )
    end select
  end subroutine gen_plot
  
  
  subroutine gen_cube( s, f, fn0 )
    implicit none
    integer,intent(in) :: s
    real(8),intent(in) :: f(0:,0:,0:)
    character(*),intent(in) :: fn0
    integer :: i,m3,i1,i2,i3,j3
    character(30) :: name, file_name
    real(8) :: r0(3), aa_del(3,3)
    character(1) :: cs
    !real(8),allocatable :: fsft(:,:,:)

    aa_del(:,1)=aa(:,1)/ML1
    aa_del(:,2)=aa(:,2)/ML2
    aa_del(:,3)=aa(:,3)/ML3

    name  = 'SYS1'
    r0(:) = 0.0d0

    write(cs,'(i1.1)') s

    file_name = fn0//"_"//cs//".cub"

!    allocate( fsft(0:ML1-1,0:ML2-1,0:ML3-1) ) ; fsft=0.0d0
!    m3=ML3/2
!    do i3=0,ML3-1
!       j3=mod(i3+m3,ML3)
!       do i2=0,ML2-1
!       do i1=0,ML1-1
!          fsft(i1,i2,j3) = f(i1,i2,i3)
!       end do
!       end do
!    end do

    open(iu,file=file_name)
    write(iu,100) name
    write(iu,100) name
    write(iu,110) MI,r0
    write(iu,110) ML1,aa_del(:,1)
    write(iu,110) ML2,aa_del(:,2)
    write(iu,110) ML3,aa_del(:,3)
    do i=1,MI
      write(iu,110) zatom(Kion(i)),real(zatom(Kion(i)),8),rsi(:,i)
    end do
    do i1=0,ML1-1
      do i2=0,ML2-1
        !write(iu,*) fsft(i1,i2,:)
        write(iu,*) f(i1,i2,:)
      end do
    end do
    close(iu)

100 format(a4)
110 format(i5,4f12.6)

!    deallocate( fsft )

  end subroutine gen_cube


  subroutine gen_chgcar( s, f, fn0 )
    implicit none
    integer,intent(in) :: s
    real(8),intent(in) :: f(0:,0:,0:)
    character(*),intent(in) :: fn0
    character(1) :: cs
    character(30) :: file_name
    integer :: i,j,k,m1,m2,m3
    real(8) :: c
    write(cs,'(i1.1)') s
    file_name = fn0//"_"//cs//".chgcar"
    open(iu,file=file_name)
    write(iu,*) file_name
    write(iu,'(1x,f15.8)') 1.0d0
    write(iu,'(1x,3f15.8)') aa(:,1)*bohr
    write(iu,'(1x,3f15.8)') aa(:,2)*bohr
    write(iu,'(1x,3f15.8)') aa(:,3)*bohr
    write(iu,'(1x,9(a6,2x))') (get_element_name(zatom(i)),i=1,MKI)
    write(iu,'(1x,9(i6,2x))') (count(Kion==i),i=1,MKI)
    write(iu,'("Direct")')
    do i = 1, MKI
      do j = 1, MI
        if ( Kion(j) == i ) write(iu,'(1x,3f15.8)') asi(:,j)
      end do
    end do
    write(iu,*)
    m1 = size(f,1)
    m2 = size(f,2)
    m3 = size(f,3)
    c  = Va
    write(*,*) "sum(f)*dV=",sum(f)*c/size(f)
    write(iu,*) m1,m2,m3
    write(iu,'((10g13.5))') (((f(i,j,k)*c,i=0,m1-1),j=0,m2-1),k=0,m3-1)
    close(iu)
  end subroutine gen_chgcar


  function get_element_name( z )

    implicit none
    character(2) :: get_element_name
    integer,intent(in) :: z
    character(2),save :: name(112)
  
    data name/ "H" ,"He", &
               "Li","Be","B" ,"C" ,"N" ,"O" ,"F" ,"Ne", &
               "Na","Mg","Al","Si","P" ,"S" ,"Cl","Ar", &
               "K" ,"Ca","Sc","Ti","V" ,"Cr","Mn","Fe","Co", &
               "Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr", &
               "Rb","Sr","Y" ,"Zr","Nb","Mo","Tc","Ru","Rh", &
               "Pd","Ag","Cd","In","Sn","Sb","Te","I" ,"Xe", &
               "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu", &
               "Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu", &
               "Hf","Ta","W" ,"Re","Os","Ir", &
               "Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn", &
               "Fr","Ra","Ac","Th","Pa","U" ,"Np","Pu","Am", &
               "Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf", &
               "Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn" /
  
    get_element_name = name(z)
  
  end function get_element_name
  
end program cubegen_vrho
