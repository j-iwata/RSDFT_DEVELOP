MODULE force_local_fftw_module

  use parallel_module
  use bb_module, only: bb
  use rgrid_module, only: Ngrid, Igrid, dV
  use ggrid_module, only: NGgrid, MGL, LLG, get_ggrid &
                         ,construct_ggrid, destruct_ggrid
  use density_module, only: rho
  use atom_module, only: aa_atom, ki_atom
  use ps_local_variables, only: vqlg
  use fftw_module
  use watch_module
  use fft_module
  use,intrinsic :: iso_c_binding

  implicit none

  PRIVATE
  PUBLIC :: calc_force_local_fftw

  logical :: first_time2 = .true.
  integer :: NGPS
  integer,allocatable :: LGPS(:,:),IGPS(:)
  integer :: MI_0,MI_1
  integer,allocatable :: icnta(:),idisa(:)
  complex(8),allocatable :: fg(:)
  integer,allocatable :: LLG_f(:,:)
  complex(8),allocatable :: zwork3(:,:,:)
  real(8),allocatable :: work(:)

CONTAINS


  SUBROUTINE calc_force_local_fftw( MI, force )
    implicit none
    integer,intent(IN) :: MI
    real(8),intent(OUT) :: force(3,MI)
    integer :: ispin,i,i1,i2,i3,ik,a,j,ierr,irank,N_MI,n
    integer :: ML1,ML2,ML3,ML,MG,ML_0,j1,j2,j3
    integer :: a1b,b1b,a2b,b2b,a3b,b3b,ab1,ab12
    real(8) :: a1,a2,a3,pi2,Gr,Gx,Gy,Gz,Vcell
    real(8) :: ctt(0:9),ett(0:9)
    complex(8),allocatable :: zrho(:) !, zrho3(:,:,:)
    complex(8),parameter :: z0=(0.d0,0.d0)
    real(8) :: zsum1,zsum2,zsum3,ztmp
    include 'mpif.h'
#ifdef _FFTW_
    include 'fftw3-mpi.f03'
#endif

    force(:,:)=0.d0

#ifdef _FFTW_

    ctt(:)=0.0d0
    ett(:)=0.0d0

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
             LGPS(1,n)=i1+1
             LGPS(2,n)=i2+1
             LGPS(3,n)=i3+1-ML3_c0
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
       allocate( zwork3(0:ML1-1,0:ML2-1,0:ML3-1) ) ; zwork3=z0
       allocate( work(size(rho,1)) ) ; work=0.0d0
       first_time2=.false.
    end if

    call watch(ctt(1),ett(1))

    work(:)=rho(:,1)
    do n=2,size(rho,2)
       work(:)=work(:)+rho(:,n)
    end do
    call d1_to_z3_fft( work, zwork3 )

    do i3=1,N_ML3_c
       do i2=1,ML2_c
       do i1=1,ML1_c
          j1=i1-1
          j2=i2-1
          j3=i3-1+ML3_c0
          zwork3_ptr0(i1,i2,i3) = zwork3(j1,j2,j3)
       end do
       end do
    end do

    call watch(ctt(2),ett(2))

    call fftw_mpi_execute_dft( plan_forward, zwork3_ptr0, zwork3_ptr1 )

    call watch(ctt(3),ett(3))

    fg(:)=z0
!$OMP parallel do
    do i=1,NGPS
       fg(IGPS(i))=conjg( zwork3_ptr1(LGPS(1,i),LGPS(2,i),LGPS(3,i)) )
    end do
!$OMP end parallel do

    call watch(ctt(4),ett(4))

    call mpi_allreduce(MPI_IN_PLACE,fg,MG,MPI_COMPLEX16 &
         ,MPI_SUM,comm_grid,ierr)

    call watch(ctt(5),ett(5))

    if ( .not.allocated(LLG_f) ) then
       allocate( LLG_f(3,NGgrid(0)) ) ; LLG_f=0
       call get_Ggrid(0,LLG_f)
    end if
!    call construct_Ggrid(0)

    call watch(ctt(6),ett(6))

    do a=MI_0,MI_1
!!$OMP parallel do private( ik,a1,a2,a3,zsum1,zsum2,zsum3,i,Gr,Gx,Gy,Gz,j,ztmp )
!    do a=1,MI

       ik=ki_atom(a)
       a1=pi2*aa_atom(1,a)
       a2=pi2*aa_atom(2,a)
       a3=pi2*aa_atom(3,a)

       zsum1=z0
       zsum2=z0
       zsum3=z0
!$OMP parallel do reduction(+:zsum1,zsum2,zsum3) private( Gx,Gy,Gz,j,ztmp )
       do i=1,NGgrid(0)
!       do i=MG_0,MG_1
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

    call watch(ctt(7),ett(7))

!    call destruct_Ggrid

    if ( nprocs <= MI ) then
       call mpi_allgatherv(force(1,MI_0),icnta(myrank),MPI_REAL8,force &
                          ,icnta,idisa,MPI_REAL8,MPI_COMM_WORLD,ierr)
    else
       call mpi_allreduce(MPI_IN_PLACE,force,size(force),mpi_real8 &
                         ,mpi_sum,mpi_comm_world,ierr)
    end if

    call watch(ctt(8),ett(8))

    if ( disp_switch_parallel ) then
       write(*,*) "time(force_local_fftw1)=",ctt(1)-ctt(0),ett(1)-ett(0)
       write(*,*) "time(force_local_fftw2)=",ctt(2)-ctt(1),ett(2)-ett(1)
       write(*,*) "time(force_local_fftw3)=",ctt(3)-ctt(2),ett(3)-ett(2)
       write(*,*) "time(force_local_fftw4)=",ctt(4)-ctt(3),ett(4)-ett(3)
       write(*,*) "time(force_local_fftw5)=",ctt(5)-ctt(4),ett(5)-ett(4)
       write(*,*) "time(force_local_fftw6)=",ctt(6)-ctt(5),ett(6)-ett(5)
       write(*,*) "time(force_local_fftw7)=",ctt(7)-ctt(6),ett(7)-ett(6)
       write(*,*) "time(force_local_fftw8)=",ctt(8)-ctt(7),ett(8)-ett(7)
    end if

#endif

  END SUBROUTINE calc_force_local_fftw


END MODULE force_local_fftw_module
