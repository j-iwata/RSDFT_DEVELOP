module hartree_sol_ffte_module

  use hartree_variables, only: E_hartree, Vh
  use bb_module, only: bb
  use ffte_sub_module, only: zwork1_ffte,zwork2_ffte,npuz,npuy,comm_fftx,comm_ffty,comm_fftz
  use rgrid_module, only: dV
  use ggrid_module, only: LLG, construct_ggrid, destruct_ggrid
  use parallel_module, only: comm_grid
  use watch_module
  use rsdft_mpi_module, only: rsdft_allreduce

  implicit none

  private
  public :: init_hartree_sol_ffte
  public :: calc_hartree_sol_ffte

  integer :: NGHT
  integer,allocatable :: LGHT(:,:)
  integer,allocatable :: IGHT(:,:)
  real(8),allocatable :: GGHT(:)

  integer :: a1b,b1b,a2b,b2b,a3b,b3b,ab1,ab12
  integer :: ML1,ML2,ML3

  logical :: flag_init_done = .false.

contains

  subroutine init_hartree_sol_ffte( Ngrid, Igrid )
    implicit none
    integer,intent(in) :: Ngrid(3)
    integer,intent(in) :: Igrid(2,3)
    integer :: n,i,i1,i2,i3,NLG
    real(8) :: pi4,g2

    if ( flag_init_done ) return

    call write_border( 0, ' init_hartree_sol_ffte(start)' )

    pi4 = 4.0d0*acos(-1.0d0)

    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

    a1b = Igrid(1,1)
    b1b = Igrid(2,1)
    a2b = Igrid(1,2)
    b2b = Igrid(2,2)
    a3b = Igrid(1,3)
    b3b = Igrid(2,3)
    ab1 = (b1b-a1b+1)
    ab12= (b1b-a1b+1)*(b2b-a2b+1)

    call construct_Ggrid(0)

    NLG = size( LLG, 2 )

    n=0
    do i=1,NLG
      i1=mod( Ngrid(1)+LLG(1,i), Ngrid(1) )
      i2=mod( Ngrid(2)+LLG(2,i), Ngrid(2) )
      i3=mod( Ngrid(3)+LLG(3,i), Ngrid(3) )
      if ( all(LLG(1:3,i)==0) ) cycle
      if ( a2b <= i2 .and. i2 <= b2b .and. a3b <= i3 .and. i3 <= b3b ) then
        n=n+1
      end if
    end do
    NGHT=n

    allocate( LGHT(3,NGHT) ); LGHT=0
    allocate( IGHT(3,NGHT) ); IGHT=0
    allocate( GGHT(NGHT) ); GGHT=0.0d0

    n=0
    do i=1,NLG
      i1=mod( Ngrid(1)+LLG(1,i), Ngrid(1) )
      i2=mod( Ngrid(2)+LLG(2,i), Ngrid(2) )
      i3=mod( Ngrid(3)+LLG(3,i), Ngrid(3) )
      if ( all(LLG(1:3,i)==0) ) cycle
      if ( a2b <= i2 .and. i2 <= b2b .and. a3b <= i3 .and. i3 <= b3b ) then
        n=n+1
        LGHT(1,n)=i1
        LGHT(2,n)=i2
        LGHT(3,n)=i3
        g2=( bb(1,1)*LLG(1,i)+bb(1,2)*LLG(2,i)+bb(1,3)*LLG(3,i) )**2 &
          +( bb(2,1)*LLG(1,i)+bb(2,2)*LLG(2,i)+bb(2,3)*LLG(3,i) )**2 &
          +( bb(3,1)*LLG(1,i)+bb(3,2)*LLG(2,i)+bb(3,3)*LLG(3,i) )**2
        GGHT(n)=pi4/g2
      end if
    end do

    call destruct_Ggrid

    flag_init_done = .true.

    call write_border( 0, ' init_hartree_sol_ffte(end)' )

  end subroutine init_hartree_sol_ffte

  subroutine calc_hartree_sol_ffte( rho )
    implicit none
    real(8),intent(in) :: rho(:,:)
    real(8) :: Eh0,ctt(0:5),ett(0:5)
    integer :: i,i1,i2,i3,ispin,nspin,ML_0
    complex(8),parameter :: z0=(0.0d0,0.0d0)
    include 'mpif.h'
    integer :: ierr

    call write_border( 1, " calc_hartree_sol_ffte(start)" )
    if ( .not.flag_init_done ) call stop_program('Call "init_hartree_sol_ffte" first.')

    nspin = size( rho, 2 )
    ML_0 = lbound( Vh, 1 )
    
    ctt(:)=0.d0
    ett(:)=0.d0

    call watch(ctt(0),ett(0))

    zwork1_ffte(:,:,:)=z0
!$OMP parallel do collapse(3) private(i)
    do i3=a3b,b3b
    do i2=a2b,b2b
    do i1=a1b,b1b
      i=1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
      zwork1_ffte(i1,i2,i3) = rho(i,1)
    end do
    end do
    end do
!$OMP end parallel do
    do ispin=2,nspin
!$OMP parallel do collapse(3) private(i)
      do i3=a3b,b3b
      do i2=a2b,b2b
      do i1=a1b,b1b
        i=1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
        zwork1_ffte(i1,i2,i3) = zwork1_ffte(i1,i2,i3) + rho(i,ispin)
      end do
      end do
      end do
!$OMP end parallel do
    end do !ispin

    !call rsdft_allreduce( zwork1_ffte, comm_fftx, out=zwork2_ffte )
    call MPI_Allreduce( zwork1_ffte, zwork2_ffte, size(zwork2_ffte), MPI_COMPLEX16, MPI_SUM, comm_fftx, ierr )

    call watch(ctt(1),ett(1))

    call pzfft3dv(zwork2_ffte,zwork1_ffte,ML1,ML2,ML3,comm_ffty,comm_fftz,npuy,npuz,-1)

    call watch(ctt(2),ett(2))

    zwork2_ffte(:,:,:)=z0
    do i=1,NGHT
      i1=LGHT(1,i)
      i2=LGHT(2,i)
      i3=LGHT(3,i)
      zwork2_ffte(i1,i2,i3) = zwork1_ffte(i1,i2,i3)*GGHT(i)
    end do

    call watch(ctt(3),ett(3))

    call pzfft3dv(zwork2_ffte,zwork1_ffte,ML1,ML2,ML3,comm_ffty,comm_fftz,npuy,npuz,1)

    call watch(ctt(4),ett(4))

!$OMP parallel do collapse(3) private(i)
    do i3=a3b,b3b
    do i2=a2b,b2b
    do i1=a1b,b1b
      i=ML_0+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
      Vh(i)=real( zwork1_ffte(i1,i2,i3) )
    end do
    end do
    end do
!$OMP end parallel do

    Eh0=0.0d0
    do ispin=1,nspin
!$OMP parallel do collapse(3) private(i) reduction(+:Eh0)
      do i3=a3b,b3b
      do i2=a2b,b2b
      do i1=a1b,b1b
        i=1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
        Eh0 = Eh0 + real( zwork1_ffte(i1,i2,i3) )*rho(i,ispin)
      end do
      end do
      end do
!$OMP end parallel do
    end do
    E_hartree = 0.5d0*Eh0*dV
    call rsdft_allreduce( E_hartree, comm_grid )

    call watch(ctt(5),ett(5))

!    call check_disp_switch( disp_switch, 0 )
!    if ( disp_switch ) then
!      write(*,*) "time(calc_hatree_sol_ffte_1)=",ctt(1)-ctt(0),ett(1)-ett(0)
!      write(*,*) "time(calc_hatree_sol_ffte_2)=",ctt(2)-ctt(1),ett(2)-ett(1)
!      write(*,*) "time(calc_hatree_sol_ffte_3)=",ctt(3)-ctt(2),ett(3)-ett(2)
!      write(*,*) "time(calc_hatree_sol_ffte_4)=",ctt(4)-ctt(3),ett(4)-ett(3)
!      write(*,*) "time(calc_hatree_sol_ffte_5)=",ctt(5)-ctt(4),ett(5)-ett(4)
!    end if

    call write_border( 1, " calc_hartree_sol_ffte(end)" )

  end subroutine calc_hartree_sol_ffte

end module hartree_sol_ffte_module
