MODULE kinetic_fft_module

  use parallel_module
  use rgrid_module
  use bb_module
  use fft_module

  implicit none

  PRIVATE
  PUBLIC :: op_kinetic_fft

CONTAINS


  SUBROUTINE op_kinetic_fft(k,tpsi,htpsi,n1,n2,ib1,ib2)
    implicit none
    integer,intent(IN) :: k,n1,n2,ib1,ib2
#ifdef _DRSDFT_
    real(8),intent(IN) :: tpsi(n1:n2,ib1:ib2)
    real(8),intent(INOUT) :: htpsi(n1:n2,ib1:ib2)
    real(8),parameter :: zero=0.0d0
    integer,parameter :: TYPE_MAIN=MPI_REAL8
    real(8),allocatable :: work(:)
#else
    complex(8),intent(IN) :: tpsi(n1:n2,ib1:ib2)
    complex(8),intent(INOUT) :: htpsi(n1:n2,ib1:ib2)
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    integer,parameter :: TYPE_MAIN=MPI_COMPLEX16
    complex(8),allocatable :: work(:)
#endif

    real(8) :: Gx,Gy,Gz,GG
    integer :: ierr,i,ib,irank,i1,i2,i3,j1,j2,j3
    integer :: ML,ML1,ML2,ML3
    complex(8),allocatable :: zwork0(:,:,:),zwork1(:,:,:)

    ML  = Ngrid(0)
    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

    allocate( work(ML) ) ; work=zero
    allocate( zwork0(0:ML1-1,0:ML2-1,0:ML3-1) ) ; zwork0=(0.0d0,0.0d0)
    allocate( zwork1(0:ML1-1,0:ML2-1,0:ML3-1) ) ; zwork1=(0.0d0,0.0d0)

    call init_fft

    do ib=ib1,ib2

       call mpi_allgatherv(tpsi(n1,ib),n2-n1+1,TYPE_MAIN &
            ,work,ir_grid,id_grid,TYPE_MAIN,comm_grid,ierr)

       i=0
       irank=-1
       do i3=1,node_partition(3)
       do i2=1,node_partition(2)
       do i1=1,node_partition(1)
          irank=irank+1
          do j3=pinfo_grid(5,irank),pinfo_grid(5,irank)+pinfo_grid(6,irank)-1
          do j2=pinfo_grid(3,irank),pinfo_grid(3,irank)+pinfo_grid(4,irank)-1
          do j1=pinfo_grid(1,irank),pinfo_grid(1,irank)+pinfo_grid(2,irank)-1
             i=i+1
             zwork0(j1,j2,j3)=work(i)
          end do
          end do
          end do
       end do
       end do
       end do

       call forward_fft( zwork0, zwork1 )

       zwork1(:,:,:)=(0.0d0,0.0d0)
       do i3=-ML3/2,ML3/2
       do i2=-ML2/2,ML2/2
       do i1=-ML1/2,ML1/2
!       do i3=-ML3/2,(ML3-1)/2
!       do i2=-ML2/2,(ML2-1)/2
!       do i1=-ML1/2,(ML1-1)/2
!       do i3=-(ML3-1)/2,ML3/2
!       do i2=-(ML2-1)/2,ML2/2
!       do i1=-(ML1-1)/2,ML1/2
!       do i3=-(ML3-1)/2,(ML3-1)/2
!       do i2=-(ML2-1)/2,(ML2-1)/2
!       do i1=-(ML1-1)/2,(ML1-1)/2
          Gx = bb(1,1)*i1 + bb(1,2)*i2 + bb(1,3)*i3
          Gy = bb(2,1)*i1 + bb(2,2)*i2 + bb(2,3)*i3
          Gz = bb(3,1)*i1 + bb(3,2)*i2 + bb(3,3)*i3
          GG = Gx*Gx + Gy*Gy + Gz*Gz
          j1 = mod(i1+ML1,ML1)
          j2 = mod(i2+ML2,ML2)
          j3 = mod(i3+ML3,ML3)
          zwork1(j1,j2,j3) = 0.5d0*GG*zwork0(j1,j2,j3)
       end do ! i1
       end do ! i2
       end do ! i3

       call backward_fft( zwork1, zwork0 )

       i=n1-1
       do i3=Igrid(1,3),Igrid(2,3)
       do i2=Igrid(1,2),Igrid(2,2)
       do i1=Igrid(1,1),Igrid(2,1)
          i=i+1
          htpsi(i,ib) = htpsi(i,ib) + zwork1(i1,i2,i3)
       end do
       end do
       end do

    end do ! ib

    call finalize_fft

    deallocate( zwork0 )
    deallocate( zwork1 )
    deallocate( work )

  END SUBROUTINE op_kinetic_fft


END MODULE kinetic_fft_module
