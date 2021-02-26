module force_local_sol_module

  use parallel_module
  use bb_module, only: bb
  use rgrid_module, only: Ngrid, Igrid, dV
  use ggrid_module, only: NGgrid, MGL, LLG, get_ggrid &
                         ,construct_ggrid, destruct_ggrid
  use density_module, only: rho
  use atom_module, only: aa_atom, ki_atom
  use ps_local_variables, only: vqlg
  use ffte_sub_module
  use watch_module
  use fft_module
  use rsdft_mpi_module
  use force_local_1dffte_module
  use force_local_fftw_module

  implicit none

  private
  public :: calc_force_local_sol

  logical :: first_time2 = .true.
  integer :: NGPS
  integer,allocatable :: LGPS(:,:),IGPS(:)
  integer :: MI_0,MI_1
  integer,allocatable :: icnta(:),idisa(:)
  complex(8),allocatable :: fg(:)
  integer,allocatable :: LLG_f(:,:)

contains

  subroutine calc_force_local_sol( force )
    implicit none
    real(8),intent(out) :: force(:,:)
    call write_border( 1, ' calc_force_ps_local_sol(start)' )
    select case( iswitch_fft )
    case( 'FFTE', 'FFTE1' )
      call calc_force_local_ffte( force )
    case( 'FFTE2' )
      call calc_force_local_1dffte( force )
    case( 'FFTW', 'FFTW1' )
      call calc_force_local_fftw( force )
    case default
      call calc_force_local_fft0( force )
    end select
    call write_border( 1, ' calc_force_local_sol(end)' )
  end subroutine calc_force_local_sol


  subroutine calc_force_local_fft0( force )
    implicit none
    real(8),intent(out) :: force(:,:)
    integer :: a,i,j,ik,ispin,irank,i1,i2,i3,n,l1,l2,l3,m1,m2,m3,j1,j2,j3
    integer :: MI,MI_0,MI_1,ML1,ML2,ML3,ML_0,ML_1,ML,MG,N_MI,ierr
    integer,allocatable :: icnt(:),idis(:)
    real(8) :: a1,a2,a3,pi2,Gr,Gx,Gy,Gz,Vcell
    real(8),allocatable :: work(:)
    complex(8),allocatable :: zrho3(:,:,:),zrho(:),zwork(:,:,:)
    complex(8),parameter :: z0=(0.0d0,0.0d0)
    complex(8) :: zsum1,zsum2,zsum3,ztmp

    call write_border( 1, ' calc_force_local_fft0(start)' )

    force(:,:)=0.0d0

    MI    = size( force, 2 )
    MG    = NGgrid(0)
    ML    = Ngrid(0)
    ML1   = Ngrid(1)
    ML2   = Ngrid(2)
    ML3   = Ngrid(3)
    ML_0  = Igrid(1,0)
    ML_1  = Igrid(2,0)
    pi2   = 2.d0*acos(-1.d0)
    Vcell = ML*dV

    allocate( icnt(0:nprocs-1) )
    allocate( idis(0:nprocs-1) )
    N_MI = MI/nprocs
    icnt(0:nprocs-1) = N_MI
    n = MI - N_MI*nprocs
    if ( n>0 ) then
      do irank=0,n-1
        icnt(irank)=icnt(irank)+1
      end do
    end if
    do irank=0,nprocs-1
      idis(irank) = sum( icnt(0:irank) ) - icnt(irank)
    end do
    MI_0 = idis(myrank)+1
    MI_1 = idis(myrank)+icnt(myrank)
    deallocate( idis, icnt )

    call init_fft

    allocate( work(ML_0:ML_1) ) ; work=0.0d0
    do ispin=1,size(rho,2)
      work(:) = work(:) + rho(:,ispin)
    end do
    call d1_to_z3_fft( work, zrho3 )
    deallocate( work )

    call forward_fft( zrho3, zwork )

    if ( allocated(zwork) ) deallocate( zwork )

    call finalize_fft

    call construct_Ggrid(0)

    allocate( zrho(MG) ) ; zrho=z0

    do i=1,NGgrid(0)
      i1=mod(ML1+LLG(1,i),ML1)
      i2=mod(ML2+LLG(2,i),ML2)
      i3=mod(ML3+LLG(3,i),ML3)
      zrho(i) = conjg( zrho3(i1,i2,i3) )
    end do

    if ( allocated(zrho3) ) deallocate( zrho3 )

    do a=MI_0,MI_1

      ik=ki_atom(a)
      a1=pi2*aa_atom(1,a)
      a2=pi2*aa_atom(2,a)
      a3=pi2*aa_atom(3,a)

      zsum1=z0
      zsum2=z0
      zsum3=z0
      do i=1,NGgrid(0)
        Gx=bb(1,1)*LLG(1,i)+bb(1,2)*LLG(2,i)+bb(1,3)*LLG(3,i)
        Gy=bb(2,1)*LLG(1,i)+bb(2,2)*LLG(2,i)+bb(2,3)*LLG(3,i)
        Gz=bb(3,1)*LLG(1,i)+bb(3,2)*LLG(2,i)+bb(3,3)*LLG(3,i)
        Gr=a1*LLG(1,i)+a2*LLG(2,i)+a3*LLG(3,i)
        j=MGL(i)
        ztmp=-vqlg(j,ik)*dcmplx(sin(Gr),cos(Gr))*zrho(i)
        zsum1=zsum1+Gx*ztmp
        zsum2=zsum2+Gy*ztmp
        zsum3=zsum3+Gz*ztmp
      end do
      force(1,a) = -zsum1*Vcell
      force(2,a) = -zsum2*Vcell
      force(3,a) = -zsum3*Vcell

    end do ! a

    call destruct_Ggrid

    deallocate( zrho )

    allocate( work(3*MI) )
    n=0
    do a=1,MI
      do i=1,3
        n=n+1
        work(n)=force(i,a)
      end do
    end do
    call mpi_allreduce(work,force,n,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    deallocate( work )

    call write_border( 1, ' calc_force_local_fft0(end)' )

  end subroutine calc_force_local_fft0


  subroutine calc_force_local_ffte( force )
    implicit none
    real(8),intent(out) :: force(:,:)
    integer :: ispin,i,i1,i2,i3,ik,a,j,ierr,irank,N_MI,n
    integer :: ML1,ML2,ML3,ML,MG,ML_0,MI
    integer :: a1b,b1b,a2b,b2b,a3b,b3b,ab1,ab12
    real(8) :: a1,a2,a3,pi2,Gr,Gx,Gy,Gz,Vcell
    real(8) :: ctt(0:9),ett(0:9)
    complex(8),allocatable :: zrho(:)
    complex(8),parameter :: z0=(0.0d0,0.0d0)
    real(8) :: zsum1,zsum2,zsum3,ztmp
    include 'mpif.h'

    call write_border( 1, ' calc_force_local_ffte(start)' )

    force(:,:)=0.0d0

    ctt(:)=0.0d0
    ett(:)=0.0d0

    MI    = size( force, 2 )
    MG    = NGgrid(0)
    ML    = Ngrid(0)
    ML1   = Ngrid(1)
    ML2   = Ngrid(2)
    ML3   = Ngrid(3)
    ML_0  = Igrid(1,0)
    pi2   = 2.d0*acos(-1.d0)
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

    if ( first_time2 ) then
      call construct_Ggrid(2)
      n=0
      do i=1,NGgrid(0)
        i1=LLG(1,i)
        i2=LLG(2,i)
        i3=LLG(3,i)
        if ( a1b <= i1 .and. i1 <= b1b .and. &
             a2b <= i2 .and. i2 <= b2b .and. &
             a3b <= i3 .and. i3 <= b3b       ) then
          n=n+1
        end if
      end do
      allocate( LGPS(3,n) ) ; LGPS=0
      allocate( IGPS(n)   ) ; IGPS=0
      n=0
      do i=1,NGgrid(0)
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
      NGPS=n
      call destruct_Ggrid
      allocate( icnta(0:nprocs-1) ) ; icnta=0
      allocate( idisa(0:nprocs-1) ) ; idisa=0
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
      first_time2=.false.
    end if

    call watch(ctt(1),ett(1))

!$OMP parallel private(i)
!$OMP workshare
    zwork1_ffte(:,:,:)=z0
!$OMP end workshare
    do ispin=1,size(rho,2)
!$OMP do collapse(3)
      do i3=a3b,b3b
      do i2=a2b,b2b
      do i1=a1b,b1b
        i=ML_0+i1-a1b+(i2-a2b)*ab1+(i3-a3b)*ab12
        zwork1_ffte(i1,i2,i3)=zwork1_ffte(i1,i2,i3)+rho(i,ispin)
      end do
      end do
      end do
!$OMP end do
    end do
!$OMP end parallel

    call watch(ctt(2),ett(2))

    call mpi_allreduce(zwork1_ffte,zwork2_ffte,ML1*(b2b-a2b+1)*(b3b-a3b+1) &
         ,RSDFT_MPI_COMPLEX16,mpi_sum,comm_fftx,ierr)

    call watch(ctt(3),ett(3))

    call pzfft3dv(zwork2_ffte,zwork1_ffte,ML1,ML2,ML3,comm_ffty,comm_fftz,npuy,npuz,-1)

    call watch(ctt(4),ett(4))

    fg(:)=z0
!$OMP parallel do
    do i=1,NGPS
      fg(IGPS(i))=conjg( zwork1_ffte(LGPS(1,i),LGPS(2,i),LGPS(3,i)) )
    end do
!$OMP end parallel do

    call watch(ctt(5),ett(5))

    call rsdft_allreduce_sum( fg, comm_grid )

    call watch(ctt(6),ett(6))

    if ( .not.allocated(LLG_f) ) then
      allocate( LLG_f(3,NGgrid(0)) ) ; LLG_f=0
      call get_Ggrid(0,LLG_f)
    end if

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
      do i=1,NGgrid(0)
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
      call rsdft_allreduce_sum( force, MPI_COMM_WORLD )
    end if

    call watch(ctt(9),ett(9))

!    if ( disp_switch_parallel ) then
!      write(*,*) "time(force_local_ffte_1)=",ctt(1)-ctt(0),ett(1)-ett(0)
!      write(*,*) "time(force_local_ffte_2)=",ctt(2)-ctt(1),ett(2)-ett(1)
!      write(*,*) "time(force_local_ffte_3)=",ctt(3)-ctt(2),ett(3)-ett(2)
!      write(*,*) "time(force_local_ffte_4)=",ctt(4)-ctt(3),ett(4)-ett(3)
!      write(*,*) "time(force_local_ffte_5)=",ctt(5)-ctt(4),ett(5)-ett(4)
!      write(*,*) "time(force_local_ffte_6)=",ctt(6)-ctt(5),ett(6)-ett(5)
!      write(*,*) "time(force_local_ffte_7)=",ctt(7)-ctt(6),ett(7)-ett(6)
!      write(*,*) "time(force_local_ffte_8)=",ctt(8)-ctt(7),ett(8)-ett(7)
!      write(*,*) "time(force_local_ffte_9)=",ctt(9)-ctt(8),ett(9)-ett(8)
!    end if

    call write_border( 1, ' calc_force_local_ffte(end)' )

  end subroutine calc_force_local_ffte

end module force_local_sol_module
