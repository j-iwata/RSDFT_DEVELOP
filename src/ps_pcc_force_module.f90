MODULE ps_pcc_force_module

  use atom_module
  use bb_module
  use ggrid_module
  use rgrid_module
  use parallel_module
  use xc_module, only: Vxc
  use ps_pcc_module, only: cdcg
  use fft_module
  use rsdft_mpi_module

  implicit none

  PRIVATE
  PUBLIC :: calc_ps_pcc_force

CONTAINS


  SUBROUTINE calc_ps_pcc_force( MI, force )
    implicit none
    integer,intent(IN) :: MI
    real(8),intent(OUT) :: force(3,MI)
    integer :: a,ik,i,j,ierr,i1,i2,i3,j1,j2,j3,irank,n,n1,n2
    integer :: MG,ML,ML1,ML2,ML3,N_MI,MI_0,MI_1,MSP_0,MSP_1
    integer,allocatable :: icnt(:),idis(:)
    real(8) :: pi2,a1,a2,a3,Gx,Gy,Gz,Gr,Vcell
    real(8),allocatable :: w1(:)
    complex(8) :: zsum1,zsum2,zsum3,ztmp
    complex(8),allocatable :: z1(:,:,:),zvxc(:),zvxc3(:,:,:)
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    include 'mpif.h'

    force(:,:) = 0.0d0

    pi2 = 2.0d0*acos(-1.0d0)
    MG  = NGgrid(0)
    ML  = Ngrid(0)
    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

    Vcell = ML*dV

! ---

    allocate( icnt(0:nprocs-1), idis(0:nprocs-1) )
    N_MI = Natom/nprocs
    icnt(0:nprocs-1) = N_MI
    n = Natom - N_MI*nprocs
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

! ---

    call init_fft

    n     = size(Vxc,1)
    n1    = Igrid(1,0)
    n2    = Igrid(2,0)
    MSP_0 = id_spin(myrank_s) + 1
    MSP_1 = id_spin(myrank_s) + ir_spin(myrank_s)

    allocate( w1(n1:n2) ) ; w1=0.0d0

    do j=MSP_0,MSP_1
       w1(n1:n2) = w1(n1:n2) + Vxc(n1:n2,j)
    end do
    call rsdft_allreduce_sum( w1(n1:n2), comm_spin )

    call d1_to_z3_fft( w1, zvxc3 )

    deallocate( w1 )

    call forward_fft( zvxc3, z1 )

    if ( allocated(z1) ) deallocate( z1 )

    call finalize_fft

! ---

    call construct_Ggrid(0)

    allocate( zvxc(MG) ) ; zvxc=zero

    do i=1,MG
       i1=mod(ML1+LLG(1,i),ML1)
       i2=mod(ML2+LLG(2,i),ML2)
       i3=mod(ML3+LLG(3,i),ML3)
       zvxc(i) = conjg( zvxc3(i1,i2,i3) )
    end do

    if ( allocated(zvxc3) ) deallocate(zvxc3)

    do a=MI_0,MI_1

       ik=ki_atom(a)
       a1=pi2*aa_atom(1,a)
       a2=pi2*aa_atom(2,a)
       a3=pi2*aa_atom(3,a)

       zsum1=zero
       zsum2=zero
       zsum3=zero
       do i=1,MG
          Gx=bb(1,1)*LLG(1,i)+bb(1,2)*LLG(2,i)+bb(1,3)*LLG(3,i)
          Gy=bb(2,1)*LLG(1,i)+bb(2,2)*LLG(2,i)+bb(2,3)*LLG(3,i)
          Gz=bb(3,1)*LLG(1,i)+bb(3,2)*LLG(2,i)+bb(3,3)*LLG(3,i)
          Gr=a1*LLG(1,i)+a2*LLG(2,i)+a3*LLG(3,i)
          j=MGL(i)
          ztmp=cdcg(j,ik)*dcmplx(sin(Gr),cos(Gr))*zvxc(i)
          zsum1=zsum1+Gx*ztmp
          zsum2=zsum2+Gy*ztmp
          zsum3=zsum3+Gz*ztmp
       end do
       force(1,a) = zsum1*Vcell
       force(2,a) = zsum2*Vcell
       force(3,a) = zsum3*Vcell

    end do ! a

    call destruct_Ggrid

    deallocate( zvxc )

    call rsdft_allreduce_sum( force, MPI_COMM_WORLD )

  END SUBROUTINE calc_ps_pcc_force


END MODULE ps_pcc_force_module
