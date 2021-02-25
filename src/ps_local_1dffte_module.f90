module ps_local_1dffte_module

  use ggrid_module, only: construct_Ggrid,destruct_Ggrid,LLG,MGL,allgatherv_Ggrid
  use ffte_sub_module, only: npux,npuy,npuz,comm_fftx,comm_ffty,comm_fftz
  use pzfft3dv_test_module, only: zwork1_ffte, pzfft3dv_test
  use watch_module

  implicit none

  private
  public :: init_ps_local_1dffte
  public :: construct_ps_local_1dffte

  integer :: NGHT
  integer,allocatable :: LGHT(:,:),IGHT(:)
  complex(8),allocatable :: fg(:)

  logical :: flag_init_done = .false.

  integer :: ML_0,a1b,b1b,a2b,b2b,a3b,b3b,ab1,ab12
  integer :: ML,ML1,ML2,ML3

contains

  subroutine init_ps_local_1dffte( Ngrid, Igrid )
    implicit none
    integer,intent(in) :: Ngrid(0:3)
    integer,intent(in) :: Igrid(2,0:3)
    integer :: NGL,i1,i2,i3,n,i

    if ( flag_init_done ) return

    call write_border( 0, ' init_ps_local_1dffte(start)' )

    ML  = Ngrid(0)
    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)
    ML_0= Igrid(1,0)
    a1b = Igrid(1,1)
    b1b = Igrid(2,1)
    a2b = Igrid(1,2)
    b2b = Igrid(2,2)
    a3b = Igrid(1,3)
    b3b = Igrid(2,3)
    ab1 = b1b-a1b+1
    ab12= (b2b-a2b+1)*(b1b-a1b+1)

    call construct_Ggrid(0)

    NGL = size( LLG, 2 )

    n=0
    do i=1,NGL
      i1=mod( Ngrid(1)+LLG(1,i), Ngrid(1) )
      i2=mod( Ngrid(2)+LLG(2,i), Ngrid(2) )
      i3=mod( Ngrid(3)+LLG(3,i), Ngrid(3) )
      if ( a1b <= i1 .and. i1 <= b1b .and. &
           a2b <= i2 .and. i2 <= b2b .and. &
           a3b <= i3 .and. i3 <= b3b ) n=n+1
    end do
    NGHT=n

    allocate( LGHT(3,n) ); LGHT=0
    allocate( IGHT(n)   ); IGHT=0

    n=0
    do i=1,NGL
      i1=mod( Ngrid(1)+LLG(1,i), Ngrid(1) )
      i2=mod( Ngrid(2)+LLG(2,i), Ngrid(2) )
      i3=mod( Ngrid(3)+LLG(3,i), Ngrid(3) )
      if ( a1b <= i1 .and. i1 <= b1b .and. &
           a2b <= i2 .and. i2 <= b2b .and. &
           a3b <= i3 .and. i3 <= b3b ) then
        n=n+1
        LGHT(1,n)=i1
        LGHT(2,n)=i2
        LGHT(3,n)=i3
        IGHT(n)=i
      end if
    end do

    allocate( fg(NGL) ); fg=(0.0d0,0.0d0)

    call destruct_Ggrid

    flag_init_done = .true.

    call write_border( 0, ' init_ps_local_1dffte(end)' )

  end subroutine init_ps_local_1dffte

  subroutine construct_ps_local_1dffte( vqlg, SGK, Vion )
    implicit none
    real(8),intent(in) :: vqlg(:,:)
    complex(8),intent(in) :: SGK(:,:)
    real(8),intent(out) :: Vion(:)
    integer :: i,i1,i2,i3,ik,j,n,nelem,MG_0,MG_1
    real(8) :: ctt(0:3),ett(0:3)

    call write_border( 0, " construct_ps_local_1dffte(start)" )
    if ( .not.flag_init_done ) call stop_program('Call "init_ps_local_1dffte" first.')

    ctt(:)=0.d0
    ett(:)=0.d0

    call watch(ctt(0),ett(0))

    nelem = size( SGK, 2 )
    MG_0  = lbound( SGK, 1 )
    MG_1  = ubound( SGK, 1 )

!$OMP parallel private(i,j,ik)
!$OMP do
    do i=MG_0,MG_1
      j=MGL(i)
      fg(i)=vqlg(j,1)*SGK(i,1)
    end do
!$OMP end do
    do ik=2,nelem
!$OMP do
      do i=MG_0,MG_1
        j=MGL(i)
        fg(i)=fg(i)+vqlg(j,ik)*SGK(i,ik)
      end do
!$OMP end do
    end do
!$OMP end parallel

    call allgatherv_Ggrid(fg)

!$OMP parallel
!$OMP workshare
    zwork1_ffte(:,:,:)=(0.0d0,0.0d0)
!$OMP end workshare
!$OMP do
    do i=1,NGHT
      zwork1_ffte(LGHT(1,i),LGHT(2,i),LGHT(3,i)) = fg(IGHT(i))
    end do
!$OMP end do
!$OMP end parallel

    call watch(ctt(1),ett(1))

    call pzfft3dv_test( zwork1_ffte, 1 )
    !call pzfft3dv(zwork1_ffte,zwork2_ffte,ML1,ML2,ML3,comm_ffty,comm_fftz,npuy,npuz,1)

    call watch(ctt(2),ett(2))

!$OMP parallel do collapse(3) private(i)
    do i3=a3b,b3b
    do i2=a2b,b2b
    do i1=a1b,b1b
       i = 1 + i1-a1b + (i2-a2b)*ab1 + (i3-a3b)*ab12
       Vion(i)=real( zwork1_ffte(i1,i2,i3) )*ML
    end do
    end do
    end do
!$OMP end parallel do

    call watch(ctt(3),ett(3))

!    if ( disp_switch_parallel ) then
!      write(*,*) "time(construct_ps_loc_1dffte_1)",ctt(1)-ctt(0),ett(1)-ett(0)
!      write(*,*) "time(construct_ps_loc_1dffte_2)",ctt(2)-ctt(1),ett(2)-ett(1)
!      write(*,*) "time(construct_ps_loc_1dffte_3)",ctt(3)-ctt(2),ett(3)-ett(2)
!    end if

    call write_border( 0, " construct_ps_local_1dffte(end)" )

  end subroutine construct_ps_local_1dffte

end module ps_local_1dffte_module
