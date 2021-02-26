module force_local_1dffte_module

  use rsdft_mpi_module, only: rsdft_allreduce, rsdft_allgatherv
  use rgrid_variables, only: Ngrid, Igrid, dV
  use ggrid_module, only: construct_Ggrid, destruct_Ggrid, get_Ggrid, LLG, MGL
  use bb_module, only: bb
  use pzfft3dv_test_module, only: zwork1_ffte,zwork2_ffte,pzfft3dv_test
  use density_module, only: rho
  use parallel_module, only: comm_grid
  use atom_module, only: ki_atom, aa_atom
  use ps_local_variables, only: vqlg
  use watch_module

  private
  public :: calc_force_local_1dffte

  logical :: flag_init_done = .false.

  integer :: MG,NGPS
  integer,allocatable :: LGPS(:,:),IGPS(:)
  integer :: MI_0,MI_1
  integer,allocatable :: icnta(:),idisa(:)
  complex(8),allocatable :: fg(:)
  integer,allocatable :: LLG_f(:,:)

contains

  subroutine calc_force_local_1dffte( force )
    implicit none
    real(8),intent(out) :: force(:,:)
    integer :: ispin,i,i1,i2,i3,ik,a,j,ierr,irank,N_MI,n
    integer :: ML1,ML2,ML3,ML,ML_0,MI
    integer :: a1b,b1b,a2b,b2b,a3b,b3b,ab1,ab12
    integer :: myrank, nprocs
    real(8) :: a1,a2,a3,pi2,Gr,Gx,Gy,Gz,Vcell
    real(8) :: ctt(0:9),ett(0:9)
    complex(8),allocatable :: zrho(:)
    complex(8),parameter :: z0=(0.0d0,0.0d0)
    real(8) :: zsum1,zsum2,zsum3,ztmp
    include 'mpif.h'

    call write_border( 1, ' calc_force_local_1dffte(start)' )

    force(:,:)=0.0d0

    ctt(:)=0.0d0
    ett(:)=0.0d0

    call MPI_Comm_size( MPI_COMM_WORLD, nprocs, ierr )
    call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr )

    MI    = size( force, 2 )
    ML    = Ngrid(0)
    ML1   = Ngrid(1)
    ML2   = Ngrid(2)
    ML3   = Ngrid(3)
    ML_0  = Igrid(1,0)
    pi2   = 2.0d0*acos(-1.0d0)
    Vcell = ML*dV
    a1b = Igrid(1,1)
    b1b = Igrid(2,1)
    a2b = Igrid(1,2)
    b2b = Igrid(2,2)
    a3b = Igrid(1,3)
    b3b = Igrid(2,3)
    ab1 = b1b-a1b+1
    ab12= (b1b-a1b+1)*(b2b-a2b+1)

    call watch(ctt(0),ett(0))

    if ( .not.flag_init_done ) then
      call construct_Ggrid(2)
      MG = size( LLG, 2 )
      n=0
      do i=1,MG
        i1=LLG(1,i)
        i2=LLG(2,i)
        i3=LLG(3,i)
        if ( a1b <= i1 .and. i1 <= b1b .and. &
             a2b <= i2 .and. i2 <= b2b .and. &
             a3b <= i3 .and. i3 <= b3b ) n=n+1
      end do
      NGPS=n
      allocate( LGPS(3,NGPS) ); LGPS=0
      allocate( IGPS(NGPS)   ); IGPS=0
      n=0
      do i=1,MG
        i1=LLG(1,i)
        i2=LLG(2,i)
        i3=LLG(3,i)
        if ( a1b <= i1 .and. i1 <= b1b .and. &
             a2b <= i2 .and. i2 <= b2b .and. &
             a3b <= i3 .and. i3 <= b3b       ) then
          n=n+1
          LGPS(1,n)=i1
          LGPS(2,n)=i2
          LGPS(3,n)=i3
          IGPS(n)=i
        end if
      end do
      call destruct_Ggrid
      allocate( icnta(0:nprocs-1) ); icnta=0
      allocate( idisa(0:nprocs-1) ); idisa=0
      N_MI = MI/nprocs
      icnta(0:nprocs-1) = N_MI
      n = MI - N_MI*nprocs
      if ( n>0 ) then
        do irank=0,n-1
          icnta(irank)=icnta(irank)+1
        end do
      end if
      do irank=0,nprocs-1
        idisa(irank) = sum( icnta(0:irank) ) - icnta(irank)
      end do
      MI_0 = idisa(myrank)+1
      MI_1 = idisa(myrank)+icnta(myrank)
      idisa(:)=idisa(:)*3
      icnta(:)=icnta(:)*3
      allocate( fg(MG) ) ; fg=z0
      flag_init_done = .true.
    end if

    if ( .not.allocated(LLG_f) ) then
      allocate( LLG_f(3,MG) ) ; LLG_f=0
      call get_Ggrid(0,LLG_f)
    end if

    call watch(ctt(1),ett(1))

!$OMP parallel private(i)
!$OMP do collapse(3)
    do i3=a3b,b3b
    do i2=a2b,b2b
    do i1=a1b,b1b
      i=ML_0+i1-a1b+(i2-a2b)*ab1+(i3-a3b)*ab12
      zwork1_ffte(i1,i2,i3) = rho(i,1)
    end do
    end do
    end do
!$OMP end do
    do ispin=2,size(rho,2)
!$OMP do collapse(3)
      do i3=a3b,b3b
      do i2=a2b,b2b
      do i1=a1b,b1b
        i=ML_0+i1-a1b+(i2-a2b)*ab1+(i3-a3b)*ab12
        zwork1_ffte(i1,i2,i3) = zwork1_ffte(i1,i2,i3) + rho(i,ispin)
      end do
      end do
      end do
!$OMP end do
    end do !ispin
!$OMP end parallel

    call watch(ctt(2),ett(2))

    call pzfft3dv_test( zwork1_ffte, -1 )

    call watch(ctt(4),ett(4))

    fg(:)=z0
!$OMP parallel do
    do i=1,NGPS
      fg(IGPS(i))=conjg( zwork1_ffte(LGPS(1,i),LGPS(2,i),LGPS(3,i)) )
    end do
!$OMP end parallel do

    call watch(ctt(5),ett(5))

    call rsdft_allreduce( fg, comm_grid )

    call watch(ctt(6),ett(6))

!    if ( .not.allocated(LLG_f) ) then
!      allocate( LLG_f(3,MG) ) ; LLG_f=0
!      call get_Ggrid(0,LLG_f)
!    end if

    call watch(ctt(7),ett(7))

    do a=MI_0,MI_1

      ik=ki_atom(a)
      a1=pi2*aa_atom(1,a)
      a2=pi2*aa_atom(2,a)
      a3=pi2*aa_atom(3,a)

      zsum1=z0
      zsum2=z0
      zsum3=z0
!$OMP parallel do reduction(+:zsum1,zsum2,zsum3) private( Gx,Gy,Gz,Gr,j,ztmp )
      do i=1,MG
        Gx=bb(1,1)*LLG_f(1,i)+bb(1,2)*LLG_f(2,i)+bb(1,3)*LLG_f(3,i)
        Gy=bb(2,1)*LLG_f(1,i)+bb(2,2)*LLG_f(2,i)+bb(2,3)*LLG_f(3,i)
        Gz=bb(3,1)*LLG_f(1,i)+bb(3,2)*LLG_f(2,i)+bb(3,3)*LLG_f(3,i)
        Gr=a1*LLG_f(1,i)+a2*LLG_f(2,i)+a3*LLG_f(3,i)
        j=MGL(i)
        ztmp=-vqlg(j,ik)*dcmplx(sin(Gr),cos(Gr))*fg(i)
        zsum1=zsum1+Gx*ztmp
        zsum2=zsum2+Gy*ztmp
        zsum3=zsum3+Gz*ztmp
      end do
!$OMP end parallel do

      force(1,a) = -zsum1*dV
      force(2,a) = -zsum2*dV
      force(3,a) = -zsum3*dV

    end do ! a
!!$OMP end parallel do

    call watch(ctt(8),ett(8))

    if ( nprocs <= MI ) then
      call rsdft_allgatherv( force(:,MI_0:MI_1), force, icnta, idisa, MPI_COMM_WORLD )
    else
      call rsdft_allreduce( force )
    end if

    call watch(ctt(9),ett(9))

!    if ( disp_switch_parallel ) then
!      write(*,*) "time(force_local_1dffte_1)=",ctt(1)-ctt(0),ett(1)-ett(0)
!      write(*,*) "time(force_local_1dffte_2)=",ctt(2)-ctt(1),ett(2)-ett(1)
!      write(*,*) "time(force_local_1dffte_3)=",ctt(3)-ctt(2),ett(3)-ett(2)
!      write(*,*) "time(force_local_1dffte_4)=",ctt(4)-ctt(3),ett(4)-ett(3)
!      write(*,*) "time(force_local_1dffte_5)=",ctt(5)-ctt(4),ett(5)-ett(4)
!      write(*,*) "time(force_local_1dffte_6)=",ctt(6)-ctt(5),ett(6)-ett(5)
!      write(*,*) "time(force_local_1dffte_7)=",ctt(7)-ctt(6),ett(7)-ett(6)
!      write(*,*) "time(force_local_1dffte_8)=",ctt(8)-ctt(7),ett(8)-ett(7)
!      write(*,*) "time(force_local_1dffte_9)=",ctt(9)-ctt(8),ett(9)-ett(8)
!    end if

    call write_border( 1, ' calc_force_local_1dffte(end)' )

  end subroutine calc_force_local_1dffte

end module force_local_1dffte_module
