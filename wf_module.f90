MODULE wf_module

  use rgrid_module, only: dV
  use array_bound_module, only: ML_0,ML_1,MB,MB_0,MB_1 &
                               ,MBZ,MBZ_0,MBZ_1,MSP,MSP_0,MSP_1

  implicit none

  PRIVATE
  PUBLIC :: unk,esp,res,init_wf,test_on_wf,gather_wf

#ifdef _DRSDFT_
  real(8),parameter :: zero=0.d0
  real(8),allocatable :: unk(:,:,:,:)
#else
  complex(8),parameter :: zero=(0.d0,0.d0)
  complex(8),allocatable :: unk(:,:,:,:)
#endif

  real(8),allocatable :: esp(:,:,:)
  real(8),allocatable :: res(:,:,:)

CONTAINS

  SUBROUTINE init_wf
    integer :: s,k,n,i
    integer,allocatable :: ir(:)
    real(8) :: u1,u2

    if ( allocated(res) ) deallocate(res)
    if ( allocated(esp) ) deallocate(esp)
    if ( allocated(unk) ) deallocate(unk)

    allocate( unk(ML_0:ML_1,MB,MBZ_0:MBZ_1,MSP_0:MSP_1) ) ; unk=zero
    allocate( esp(MB,MBZ,MSP) ) ; esp=0.d0
    allocate( res(MB,MBZ,MSP) ) ; res=0.d0

    call random_seed( size=n )
    allocate( ir(n) )
    ir(:)=MB_0+ML_0
    call random_seed( put=ir )
    deallocate( ir )
    do s=MSP_0,MSP_1
       do k=MBZ_0,MBZ_1
          do n=MB_0,MB_1
             do i=ML_0,ML_1
                call random_number(u1)
                call random_number(u2)
                unk(i,n,k,s)=dcmplx(u1,u2)
             end do
          end do
       end do
    end do
  END SUBROUTINE init_wf

  SUBROUTINE test_on_wf(disp_switch)
    logical,intent(IN) :: disp_switch
    include 'mpif.h'
    integer :: ierr,s,k,m,n
    complex(8),allocatable :: uu(:,:),vv(:,:)
    allocate( uu(MB,MB) )
    allocate( vv(MB,MB) )
    do s=1,MSP
    do k=1,MBZ
       uu=(0.d0,0.d0)
       do n=1,MB
       do m=1,n
#ifdef _DRSDFT_
          uu(m,n)=sum( unk(:,m,k,s)*unk(:,n,k,s) )*dV
#else
          uu(m,n)=sum( conjg(unk(:,m,k,s))*unk(:,n,k,s) )*dV
#endif
       end do
       end do
       call mpi_allreduce(uu,vv,MB*MB,MPI_COMPLEX16,MPI_SUM &
            ,MPI_COMM_WORLD,ierr)
       do n=1,MB
       do m=1,n
          if ( disp_switch ) then
          write(*,'(1x,i2,i5,2i7,2g25.16)') &
               s,k,m,n,real(vv(m,n)),aimag(vv(m,n))
          end if
       end do
       end do
    end do
    end do
    deallocate( vv )
    deallocate( uu )
  END SUBROUTINE test_on_wf


  SUBROUTINE gather_wf
    use parallel_module
    integer :: k,s,ierr
    ir_band(:)=ir_band(:)*(ML_1-ML_0+1)
    id_band(:)=id_band(:)*(ML_1-ML_0+1)
    do s=MSP_0,MSP_1
    do k=MBZ_0,MBZ_1
       call mpi_allgatherv( unk(ML_0,MB_0,k,s),ir_band(myrank_b),MPI_REAL8 &
            ,unk(ML_0,1,k,s),ir_band,id_band,MPI_REAL8,comm_band,ierr )
    end do
    end do
    ir_band(:)=ir_band(:)/(ML_1-ML_0+1)
    id_band(:)=id_band(:)/(ML_1-ML_0+1)
  END SUBROUTINE gather_wf

END MODULE wf_module
