module construct_vion_vh_floc_module

  use rgrid_variables, only: Ngrid, dV, Igrid
  use ffte_sub_module, only: zwork1_ffte, zwork2_ffte &
                           , comm_fftx, comm_ffty, comm_fftz, npuy, npuz
  use atom_module, only: ki_atom, aa_atom
  use bb_module, only: bb
  use ggrid_module, only: construct_Ggrid, MGL
  use ps_local_variables, only: vqlg
  use parallel_module, only: RSDFT_MPI_COMPLEX16, comm_grid, construct_id_ir_parallel, myrank
  use rsdft_mpi_module, only: rsdft_allreduce_sum, rsdft_allgatherv

  implicit none

  private
  public :: construct_vion_vh_floc

  logical :: first_time=.true.
  integer,allocatable :: LLG(:,:)
  integer,allocatable :: LGHT(:,:), IGHT(:)
  real(8),allocatable :: GGHT(:)
  integer :: a1b,b1b,a2b,b2b,a3b,b3b,ab1,ab12
  complex(8),allocatable :: fg(:), GGVL(:)
  integer,allocatable :: icnta(:), idisa(:)
  integer :: MI_0, MI_1

contains

  subroutine construct_vion_vh_floc( rho, Vion, Vh, force_local, Eh )
    implicit none
    real(8),intent(in)  :: rho(:,:)
    real(8),intent(out) :: Vion(:), Vh(:), force_local(:,:)
    real(8),intent(out) :: Eh
    real(8) :: a1,a2,a3,pi2,Gx,Gy,Gz,Gr
    real(8),allocatable :: work(:,:)
    complex(8),parameter :: z0=(0.0d0,0.0d0), zi=(0.0d0,1.0d0)
    complex(8) :: zsum1,zsum2,zsum3,ztmpv,ztmpf
    integer :: i1,i2,i3,is,i,ii,j,ik,a,ierr
    integer :: Nspin, Natom
    include 'mpif.h'

    Nspin = size( rho, 2 )
    Natom = size( force_local, 2 )
    pi2   = 2.0d0*acos(-1.0d0)

    if ( first_time ) then
       call initial_preparation
       call construct_id_ir_parallel( idisa, icnta, Natom, MPI_COMM_WORLD, MI_0, MI_1 )
       idisa=idisa*3
       icnta=icnta*3
       first_time = .false.
    end if

! ---

    zwork1_ffte(:,:,:)=z0
    do is=1,Nspin
!$OMP parallel do collapse(3) private(i)
       do i3=a3b,b3b
       do i2=a2b,b2b
       do i1=a1b,b1b
          i=1+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
          zwork1_ffte(i1,i2,i3)=zwork1_ffte(i1,i2,i3)+rho(i,is)
       end do
       end do
       end do
!$OMP end parallel do
    end do

    call MPI_Allreduce(zwork1_ffte,zwork2_ffte,size(zwork2_ffte) &
         ,RSDFT_MPI_COMPLEX16,MPI_SUM,comm_fftx,ierr)

    call pzfft3dv(zwork2_ffte,zwork1_ffte,Ngrid(1),Ngrid(2),Ngrid(3) &
         ,comm_ffty,comm_fftz,npuy,npuz,-1)

    fg(:)=z0
!$OMP parallel do
    do i=1,size(IGHT)
       i1=LGHT(1,i)
       if ( i1 < a1b .or. b1b < i1 ) cycle
       i2=LGHT(2,i)
       i3=LGHT(3,i)
       ii=IGHT(i)
       fg(ii)=conjg( zwork1_ffte(i1,i2,i3) )
    end do
!$OMP end parallel do
    call rsdft_allreduce_sum( fg, comm_grid )

! ---

    GGVL(:)=z0

    do a=MI_0,MI_1

       ik=ki_atom(a)
       a1=pi2*aa_atom(1,a)
       a2=pi2*aa_atom(2,a)
       a3=pi2*aa_atom(3,a)

       zsum1=z0
       zsum2=z0
       zsum3=z0
       do i=1,size(MGL)
          i1=LLG(1,i)
          i2=LLG(2,i)
          i3=LLG(3,i)
          Gx=bb(1,1)*i1+bb(1,2)*i2+bb(1,3)*i3
          Gy=bb(2,1)*i1+bb(2,2)*i2+bb(2,3)*i3
          Gz=bb(3,1)*i1+bb(3,2)*i2+bb(3,3)*i3
          Gr=a1*i1+a2*i2+a3*i3
          j=MGL(i)
          ztmpv = vqlg(j,ik)*dcmplx( cos(Gr),-sin(Gr) )
          GGVL(i) = GGVL(i) + ztmpv
          ztmpf =-ztmpv*fg(i)*zi
          zsum1 = zsum1 + Gx*ztmpf
          zsum2 = zsum2 + Gy*ztmpf
          zsum3 = zsum3 + Gz*ztmpf
       end do ! i

       force_local(1,a) = -zsum1*dV
       force_local(2,a) = -zsum2*dV
       force_local(3,a) = -zsum3*dV

    end do ! a

    if ( size(icnta) <= Natom ) then
       call rsdft_allgatherv &
            ( force_local(:,MI_0:MI_1), force_local, icnta, idisa, MPI_COMM_WORLD )
    else
       call rsdft_allreduce_sum( force_local, MPI_COMM_WORLD )
    end if

    call rsdft_allreduce_sum( GGVL, MPI_COMM_WORLD )

! ---

    Eh=0.0d0
    zwork2_ffte(:,:,:)=z0
    do i=1,size(IGHT)
       ii=IGHT(i)
       i1=LGHT(1,i)
       i2=LGHT(2,i)
       i3=LGHT(3,i)
       ztmpv=zwork1_ffte(i1,i2,i3)
       ztmpf=ztmpv*GGHT(i)
       Eh=Eh+conjg(ztmpv)*ztmpf
       zwork2_ffte(i1,i2,i3)=ztmpf+GGVL(ii)*Ngrid(0)
    end do

    call rsdft_allreduce_sum( Eh, comm_ffty )
    call rsdft_allreduce_sum( Eh, comm_fftz )

    Eh=-Eh*0.5d0*dV/Ngrid(0)

    call pzfft3dv(zwork2_ffte,zwork1_ffte,Ngrid(1),Ngrid(2),Ngrid(3) &
         ,comm_ffty,comm_fftz,npuy,npuz,1)

!$OMP parallel do collapse(3) private(i)
    do i3=a3b,b3b
    do i2=a2b,b2b
    do i1=a1b,b1b
       i = 1 + i1-a1b + (i2-a2b)*ab1 + (i3-a3b)*ab12
       Vion(i)=real( zwork1_ffte(i1,i2,i3) )
    end do
    end do
    end do
!$OMP end parallel do

    Vh(:)=0.0d0

  end subroutine construct_vion_vh_floc


  subroutine initial_preparation
    implicit none
    integer :: n,i,i1,i2,i3,j1,j2,j3
    real(8) :: pi4, g2

    pi4 = 4.0d0*acos(-1.0d0)

    a1b = Igrid(1,1)
    b1b = Igrid(2,1)
    a2b = Igrid(1,2)
    b2b = Igrid(2,2)
    a3b = Igrid(1,3)
    b3b = Igrid(2,3)
    ab1 = b1b-a1b+1
    ab12= ab1*(b2b-a2b+1)

    call construct_Ggrid( 0, LLG )

    n=0
    do i=1,size(LLG,2)
       i1=mod( Ngrid(1)+LLG(1,i), Ngrid(1) )
       i2=mod( Ngrid(2)+LLG(2,i), Ngrid(2) )
       i3=mod( Ngrid(3)+LLG(3,i), Ngrid(3) )
       if ( a2b <= i2 .and. i2 <= b2b .and. a3b <= i3 .and. i3 <= b3b ) then
          n=n+1
       end if
    end do

    allocate( LGHT(3,n) ); LGHT=0
    allocate( IGHT(n)   ); IGHT=0
    allocate( GGHT(n)   ); GGHT=0.0d0

    n=0
    do i=1,size(LLG,2)
       i1=LLG(1,i)
       i2=LLG(2,i)
       i3=LLG(3,i)
       j1=mod( Ngrid(1)+i1, Ngrid(1) )
       j2=mod( Ngrid(2)+i2, Ngrid(2) )
       j3=mod( Ngrid(3)+i3, Ngrid(3) )
       if ( a2b<=j2.and.j2<=b2b .and. a3b<=j3.and.j3<=b3b ) then
          n=n+1
          LGHT(1,n)=j1
          LGHT(2,n)=j2
          LGHT(3,n)=j3
          IGHT(n)=i
          if ( i1/=0 .or. i2/=0 .or. i3/=0 ) then
             g2=( bb(1,1)*i1+bb(1,2)*i2+bb(1,3)*i3 )**2 &
               +( bb(2,1)*i1+bb(2,2)*i2+bb(2,3)*i3 )**2 &
               +( bb(3,1)*i1+bb(3,2)*i2+bb(3,3)*i3 )**2
             GGHT(n)=pi4/g2
          end if
       end if
    end do

    allocate( fg(size(LLG,2)) ); fg=(0.0d0,0.0d0)
    allocate( GGVL(size(LLG,2)) ); GGVL=(0.0d0,0.0d0)

  end subroutine initial_preparation


end module construct_vion_vh_floc_module
