MODULE sqerr_g_module

  use parallel_module, only: comm_grid, comm_spin
  use bb_module, only: bb
  use fft_module, only: init_fft,finalize_fft,forward_fft,d1_to_z3_fft

  implicit none

  PRIVATE
  PUBLIC :: calc_sqerr_g
  PUBLIC :: init_sqerr_g

  real(8),allocatable :: fold(:,:,:)

CONTAINS


  SUBROUTINE init_sqerr_g( m, n, nmax, f1_in, f2_in )
    implicit none
    integer,intent(IN) :: m,n,nmax
    real(8),intent(IN) :: f1_in(m,n),f2_in(m,n)
    integer :: n0,ierr
    integer,parameter :: ndat=2
    include 'mpif.h'
    call write_border( 1, " init_sqerr_g(start)" )
    if ( allocated(fold) ) deallocate( fold )
    allocate( fold(m,nmax,ndat) ) ; fold=0.0d0
    n0=size( f1_in(:,:)  )
    call mpi_allgather( f1_in,n0,MPI_REAL8,fold(1,1,1),n0,MPI_REAL8,comm_spin,ierr )
    call mpi_allgather( f2_in,n0,MPI_REAL8,fold(1,1,2),n0,MPI_REAL8,comm_spin,ierr )
    call write_border( 1, " init_sqerr_g(end)" )
  END SUBROUTINE init_sqerr_g


  SUBROUTINE calc_sqerr_g( m, n, nmax, ndat, f_in )
    implicit none
    integer,intent(IN)  :: m,n,nmax,ndat
    real(8),intent(IN)  :: f_in(m,n,ndat)
    integer :: ierr,i,n0,j
    real(8),allocatable :: f(:,:,:),ftmp(:)
    complex(8),allocatable :: zw(:,:,:),zw0(:,:,:),zw1(:,:,:)
    include 'mpif.h'

    call write_border( 1, " calc_sqerr_g(start)" )

    allocate( f(m,nmax,ndat) ) ; f=0.0d0

    n0=size( f_in(:,:,1) )
    do i=1,ndat
       call mpi_allgather( f_in(1,1,i),n0,MPI_REAL8,f(1,1,i),n0,MPI_REAL8,comm_spin,ierr )
    end do

    call init_fft

    allocate( ftmp(m) ) ; ftmp=0.0d0

    do i=1,ndat
       do j=1,nmax
!          call d1_to_z3_fft( fold(:,j,i), zw0 )
!          call forward_fft( zw0, zw )
!          call d1_to_z3_fft( f(:,j,i), zw1 )!
!          call forward_fft( zw1, zw )
!          zw1=zw1-zw0
          ftmp(:)=f(:,j,i)-fold(:,j,i)
          call d1_to_z3_fft( ftmp, zw1 )
          call forward_fft( zw1, zw )
          call analyze_residue( zw1 )
       end do
    end do

    deallocate( ftmp )

    call finalize_fft

    fold=f

    deallocate( f )

    call write_border( 1, " calc_sqerr_g(end)" )

    return
  END SUBROUTINE calc_sqerr_g


  SUBROUTINE analyze_residue( zw )
    implicit none
    complex(8),intent(IN) :: zw(0:,0:,0:)
    integer :: ML1,ML2,ML3,i1,i2,i3,j1,j2,j3,k1,k2,k3
    real(8) :: dif_tot,dif_min,dif_max,dif,gmin,gmax
    real(8) :: gg,gx,gy,gz,gg_min,gg_max,dif_gg_min,dif_gg_max
    logical :: disp

    ML1=size(zw,1)
    ML2=size(zw,2)
    ML3=size(zw,3)

    dif_tot=0.0d0
    dif_min= 1.d100
    dif_max=-1.d100

    dif_gg_min=0.0d0
    dif_gg_max=0.0d0
    gg_min= 1.d100
    gg_max=-1.d100

    do i3=0,ML3-1
    do i2=0,ML2-1
    do i1=0,ML1-1
       j1 = i1-ML1 ; if ( j1 < -(ML1-1)/2 ) j1=i1
       j2 = i2-ML2 ; if ( j2 < -(ML2-1)/2 ) j2=i2
       j3 = i3-ML3 ; if ( j3 < -(ML3-1)/2 ) j3=i3
       k1 = mod(j1+ML1,ML1)
       k2 = mod(j2+ML2,ML2)
       k3 = mod(j3+ML3,ML3)
       if ( k1/=i1 .or. k2/=i2 .or. k3/=i3 ) then
          write(*,*) k1,k2,k3,i1,i2,i3
          call stop_program_f( "stop@sqerr_g_module" )
       end if
       gx=j1*bb(1,1)+j2*bb(1,2)+j3*bb(1,3)
       gy=j1*bb(2,1)+j2*bb(2,2)+j3*bb(2,3)
       gz=j1*bb(3,1)+j2*bb(3,2)+j3*bb(3,3)
       gg=gx*gx+gy*gy+gz*gz
       dif=abs(zw(i1,i2,i3))**2
       dif_tot=dif_tot+dif
       if ( j1==0 .and. j2==0 .and. j3==0 ) cycle
       if ( dif < dif_min ) then
          dif_min=dif
          gmin=gg
       end if
       if ( dif > dif_max ) then
          dif_max=dif
          gmax=gg
       end if
       if ( gg < gg_min ) then
          gg_min=gg
          dif_gg_min=dif
       end if
       if ( gg > gg_max ) then
          gg_max=gg
          dif_gg_max=dif
       end if
    end do
    end do
    end do
    call check_disp_switch( disp, 0 )
    if ( disp ) then
       write(*,*)
       write(*,'(1x,"   max_diff=",es14.5," at gg=",f10.5)') dif_max, gmax
       write(*,'(1x,"diff@gg_min=",es14.5," at gg=",f10.5)') dif_gg_min,gg_min
       write(*,'(1x,"   min_diff=",es14.5," at gg=",f10.5)') dif_min, gmin
       write(*,'(1x,"diff@gg_max=",es14.5," at gg=",f10.5)') dif_gg_max,gg_max
       write(*,'(1x," diff_total=",es14.5)') dif_tot
    end if
  END SUBROUTINE analyze_residue


END MODULE sqerr_g_module
