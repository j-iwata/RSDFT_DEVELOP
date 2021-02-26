module ffte_sub_module

  implicit none

  private
  public :: init_ffte_sub, free_ffte_sub, comm_fftx, comm_ffty, comm_fftz &
           ,npux, npuy, npuz, zwork1_ffte, zwork2_ffte

  integer :: comm_fftx, comm_ffty, comm_fftz
  integer :: npux, npuy, npuz

  complex(8),allocatable :: zwork1_ffte(:,:,:)
  complex(8),allocatable :: zwork2_ffte(:,:,:)

contains

  subroutine init_ffte_sub(ig,ng,np,comm)
    implicit none
    integer,intent(in) :: ig(2,3),ng(3),np(3),comm
    integer :: ix,iy,iz,icolor,ierr,nprocs,i,j,n
    integer :: myrnk,irnk,icolor_x,icolor_y,icolor_z
    complex(8) :: z1(1),z2(1)
    complex(8) :: z0=(0.0d0,0.0d0)
    include 'mpif.h'

    call write_border( 0, ' init_ffte_sub(start)' )

!------------------------------- parameter check

    ierr=0
    ix=mod(ng(1),np(2))
    iy=mod(ng(2),np(2))
    if ( ix /= 0 .or. iy /= 0 ) then
       write(*,'(1x,"Both NX and NY must be divisible by NPUY",6i5)') ng,np
       ierr=1
    end if
    iy=mod(ng(2),np(3))
    iz=mod(ng(3),np(3))
    if ( iy /= 0 .or. iz /= 0 ) then
       write(*,'(1x,"Both NY and NZ must be divisible by NPUZ",6i5)') ng,np
       ierr=1
    end if
    if ( ierr /= 0 ) call stop_program( "stop@init_ffte_sub(1)" )

!---

    ierr=0
    do j=1,3
       n=ng(j)
       do i=1,ng(j)
          if ( mod(n,2) == 0 ) then
             n=n/2
          else if ( mod(n,3) == 0 ) then
             n=n/3
          else if ( mod(n,5) == 0 ) then
             n=n/5
          else
             if ( n /= 1 ) then
                ierr=ierr+1
                exit
             end if
          end if
       end do ! i
    end do ! j
    if ( ierr /= 0 ) then
       write(*,*) "ng(1:3)=",(ng(i),i=1,3)
       call stop_program( "ng should be composit numbers of 2,3,5" )
    end if

!------------------------------- parameter check (end)

    call MPI_Comm_rank( comm, myrnk, ierr )

    icolor_x=0
    icolor_y=0
    icolor_z=0
    irnk=-1
    do iz = 1, np(3)
    do iy = 1, np(2)
    do ix = 1, np(1)
      irnk = irnk + 1
      if ( irnk == myrnk ) then
        icolor_x = 1 + (iy-1) + np(2)*(iz-1)
        icolor_y = 1 + (iz-1) + np(3)*(ix-1)
        icolor_z = 1 + (ix-1) + np(1)*(iy-1)
      end if
    end do
    end do
    end do

    call MPI_Comm_split( comm, icolor_x, 0, comm_fftx, ierr )
    call MPI_Comm_split( comm, icolor_y, 0, comm_ffty, ierr )
    call MPI_Comm_split( comm, icolor_z, 0, comm_fftz, ierr )
    call MPI_Comm_size(comm_fftx, npux, ierr)
    call MPI_Comm_size(comm_ffty, npuy, ierr)
    call MPI_Comm_size(comm_fftz, npuz, ierr)

    call pzfft3dv(z1,z2,ng(1),ng(2),ng(3),comm_ffty,comm_fftz,npuy,npuz,0)

    allocate( zwork1_ffte(0:ng(1)-1,ig(1,2):ig(2,2),ig(1,3):ig(2,3)) ) ; zwork1_ffte=z0
    allocate( zwork2_ffte(0:ng(1)-1,ig(1,2):ig(2,2),ig(1,3):ig(2,3)) ) ; zwork2_ffte=z0

    call write_border( 0, ' init_ffte_sub(end)' )

    return



    ix = ig(1,1)/( ng(1)/np(1) )
    iy = ig(1,2)/( ng(2)/np(2) )
    iz = ig(1,3)/( ng(3)/np(3) )
    nprocs=np(1)*np(2)*np(3)
    icolor=iy+iz*np(2)
    call mpi_comm_split(comm, icolor, 0, comm_fftx, ierr)
    icolor=iz+ix*nprocs
    call mpi_comm_split(comm, icolor, 0, comm_ffty, ierr)
    icolor=iy+ix*nprocs
    call mpi_comm_split(comm, icolor, 0, comm_fftz, ierr)
    call mpi_comm_size(comm_fftx, npux, ierr)
    call mpi_comm_size(comm_ffty, npuy, ierr)
    call mpi_comm_size(comm_fftz, npuz, ierr)
    call pzfft3dv(z1,z2,ng(1),ng(2),ng(3),comm_ffty,comm_fftz,npuy,npuz,0)
    iy=ig(1,2)+ng(2)/np(2)-1
    iz=ig(1,3)+ng(3)/np(3)-1
    allocate( zwork1_ffte(0:ng(1)-1,ig(1,2):iy,ig(1,3):iz) ) ; zwork1_ffte=z0
    allocate( zwork2_ffte(0:ng(1)-1,ig(1,2):iy,ig(1,3):iz) ) ; zwork2_ffte=z0

  end subroutine init_ffte_sub

  subroutine free_ffte_sub
    implicit none
    integer :: ierr
    include 'mpif.h'
    call MPI_Comm_free(comm_fftz,ierr)
    call MPI_Comm_free(comm_ffty,ierr)
    call MPI_Comm_free(comm_fftx,ierr)
    deallocate( zwork2_ffte )
    deallocate( zwork1_ffte )
  end subroutine free_ffte_sub

end module ffte_sub_module
