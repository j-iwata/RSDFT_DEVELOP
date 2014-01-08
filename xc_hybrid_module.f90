MODULE xc_hybrid_module

  implicit none

  PRIVATE
  PUBLIC :: read_xc_hybrid, update_xc_hybrid &
           ,Nsweep_hybrid,npart,reduce_num,omega,R_hf,alpha_hf &
           ,iflag_hf,iflag_pbe0,iflag_hse,iflag_lcwpbe &
           ,VFunk,unk_hf &
           ,tr_switch,MB_fock,MBZ_fock,occ_fock,occ_factor,kbb_fock &
           ,icount_sweep_hybrid

  integer :: npart,reduce_num
  integer :: Nsweep_hybrid
  real(8) :: R_hf,omega

  integer :: icount_sweep_hybrid=0

  integer :: iflag_hf=0
  integer :: iflag_pbe0=0
  integer :: iflag_hse=0
  integer :: iflag_lcwpbe=0

  complex(8),allocatable :: VFunk(:,:,:,:)
  complex(8),allocatable :: unk_hf(:,:,:,:)

  integer :: tr_switch
  integer :: MB_fock,MBZ_fock
  real(8),parameter :: alpha_hf = 0.25d0
  real(8),allocatable :: occ_focK(:,:,:),kbb_fock(:,:)
  real(8) :: occ_factor

CONTAINS

  SUBROUTINE read_xc_hybrid(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    character(2) :: cbuf,ckey
    integer :: i
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey == "HF" ) then
             backspace(unit)
             read(unit,*) cbuf,omega,Nsweep_hybrid,npart,reduce_num
          end if
       end do
999    continue
       R_hf=omega
       write(*,*) "----- Parameters for Hybrid_XC -----"
       write(*,*) "Partition number of mpi_allgather",npart
       write(*,*) "Partition number of mpi_allreduce",reduce_num
       write(*,*) "Nsweep_hybrid =",Nsweep_hybrid
       write(*,*) "R_hf, omega   =",R_hf,omega
    end if
    call send_xc_hybrid
  END SUBROUTINE read_xc_hybrid

  SUBROUTINE send_xc_hybrid
    implicit none
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(npart        ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(reduce_num   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Nsweep_hybrid,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(R_hf      ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(omega     ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_xc_hybrid

  SUBROUTINE update_xc_hybrid(mb_0,mb_1,mk_0,ms_0,f)
    implicit none
    integer,intent(IN) :: mb_0,mb_1,mk_0,ms_0
    complex(8),intent(IN)  :: f(:,:,:,:)
    integer :: ml,mb,mk,ms,n,k,s
    ml=size(f,1)
    mb=size(f,2)
    mk=size(f,3)
    ms=size(f,4)
    unk_hf(:,:,:,:)=(0.0d0,0.0d0)
    do s=1,ms
       do k=1,mk
          do n=mb_0,mb_1
             unk_hf(:,n,mk_0-1+k,ms_0-1+s)=f(:,n,k,s)
          end do
       end do
    end do
  END SUBROUTINE update_xc_hybrid

END MODULE xc_hybrid_module
