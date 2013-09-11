MODULE localpot_module

  use rgrid_module
  use array_bound_module

  implicit none

  PRIVATE
  PUBLIC :: Vloc,init_localpot,op_localpot,read_localpot

  real(8),allocatable :: Vloc(:,:)

CONTAINS

  SUBROUTINE init_localpot
    allocate( Vloc(ML_0:ML_1,MSP_0:MSP_1) )
    Vloc=0.d0
  END SUBROUTINE init_localpot

  SUBROUTINE op_localpot(s,mm,nn,f,vf)
    integer,intent(IN) :: s,mm,nn
#ifdef _DRSDFT_
    real(8),intent(IN) :: f(mm,nn)
    real(8),intent(INOUT) :: vf(mm,nn)
#else
    complex(8),intent(IN) :: f(mm,nn)
    complex(8),intent(INOUT) :: vf(mm,nn)
#endif
    integer :: n,i
!$OMP parallel
    do n=1,nn
!$OMP do
       do i=1,mm
          vf(i,n) = vf(i,n) + Vloc(ML_0-1+i,s)*f(i,n)
       end do
!$OMP end do
    end do
!$OMP end parallel
  END SUBROUTINE op_localpot

  SUBROUTINE read_localpot(file_name,rank)
    character(*),intent(IN) :: file_name
    integer,intent(IN) :: rank
    integer,parameter :: unit=80
    integer :: ML_tmp(0:3),ierr,i1,i2,i3,i,s
    integer,allocatable :: LL_tmp(:,:)
    real(8),allocatable :: rtmp(:),rtmp3(:,:,:)
    include 'mpif.h'
    if ( rank == 0 ) then
       open(unit,file=file_name,form='unformatted')
       read(unit) ML_tmp(0:3)
       write(*,*) "ML_TMP=",ML_tmp
    end if
    call mpi_bcast(ML_tmp,4,mpi_integer,0,mpi_comm_world,ierr)
    if ( ML_tmp(0) /= Ngrid(0) ) stop "at read_localpot"
    allocate( LL_tmp(3,ML_tmp(0)) )
    if ( rank == 0 ) then
       read(unit) LL_tmp(:,:)
    end if
    call mpi_bcast(LL_tmp,3*ML_tmp(0),mpi_integer,0,mpi_comm_world,ierr)
    allocate( rtmp3(0:ML_tmp(1)-1,0:ML_tmp(2)-1,0:ML_tmp(3)-1) )
    allocate( rtmp(ML_tmp(0)) )
    do s=1,MSP
       if ( rank == 0 ) then
          read(unit) rtmp(:) ! rho
          read(unit) rtmp(:) ! Vloc
       end if
       call mpi_bcast(rtmp,ML_tmp(0),mpi_real8,0,mpi_comm_world,ierr)
       do i=1,ML_tmp(0)
          rtmp3( LL_tmp(1,i),LL_tmp(2,i),LL_tmp(3,i) ) = rtmp(i)
       end do
       if ( MSP_0 <= s .and. s <= MSP_1 ) then
          i=ML_0-1
          do i3=Igrid(1,3),Igrid(2,3)
          do i2=Igrid(1,2),Igrid(2,2)
          do i1=Igrid(1,1),Igrid(2,1)
             i=i+1
             Vloc(i,s) = rtmp3(i1,i2,i3)
          end do
          end do
          end do
       end if
    end do ! s
    if ( rank == 0 ) then
       close(unit)
    end if
    deallocate( rtmp   )
    deallocate( rtmp3  )
    deallocate( LL_tmp )
  END SUBROUTINE read_localpot

END MODULE localpot_module
