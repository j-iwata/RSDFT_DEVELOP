MODULE band_analysis_module

  use wf_module
  use parallel_module
  use rgrid_module
  use aa_module
  use array_bound_module

  implicit none

  PRIVATE
  PUBLIC :: band_analysis

CONTAINS

  SUBROUTINE band_analysis
    implicit none

    integer,allocatable :: indx(:,:,:),ichk(:)
    integer :: nk,ML3,m,n,k,nrnk,i,j,i1,i2,i3,j3,ierr,icount
    real(8),allocatable :: omax(:,:,:),f(:,:),f1(:,:)
    real(8) :: pi2,c3,h3,kz,kr,c1,c2,e1,e2
    complex(8),parameter ::z0=(0.d0,0.d0)
    complex(8),allocatable :: ovlp0(:),ovlp(:)

    pi2=2.d0*acos(-1.d0)
    c3=1.d0/Ngrid(3)
    h3=Hgrid(3)
    ML3=Ngrid(3)
    nrnk=6
    nk=19
    kz=pi2/aa(3,3)

    allocate( ovlp0(0:nk-1)   ) ; ovlp0=z0
    allocate( ovlp(0:nk-1)    ) ; ovlp=z0
    allocate( indx(nrnk,0:nk-1,MB) ) ; indx=0
    allocate( omax(nrnk,0:nk-1,MB) ) ; omax=0.d0
    allocate( ichk(MB)        ) ; ichk=0
    allocate( f(0:Ngrid(3)-1,MB) ) ; f=0.d0
    allocate( f1(0:Ngrid(3)-1,MB) ) ; f1=0.d0

    icount=0

    do n=1 , 20 !884,MB
       do m=n,n
          ovlp0(:)=z0
          i=Igrid(1,0)-1
          do i3=Igrid(1,3),Igrid(2,3)
          do i2=Igrid(1,2),Igrid(2,2)
          do i1=Igrid(1,1),Igrid(2,1)
             i=i+1
             do k=0,nk/2
                kr=k*kz*i3*h3
                ovlp0(k)=ovlp0(k) &
                     +dcmplx(cos(kr),sin(kr))*unk(i,n,1,1)
             end do
!             f(i3,n)=f(i3,n)+unk(i,n,1,1)**2
          end do
          end do
          end do
!          call mpi_allreduce(f(:,n),f1(:,n),Ngrid(3),mpi_real8 &
!               ,mpi_sum,comm_grid,ierr)
          call mpi_allreduce(ovlp0,ovlp,nk,MPI_COMPLEX16,MPI_SUM,comm_grid,ierr)
          do k=0,nk/2
             c1=abs(ovlp(k))**2
             if ( myrank == 0 ) write(*,'(1x,2i6,3g20.10,f15.10)') n,k,c1,ovlp(k),esp(n,1,1)
!             if ( c1 > 1.d0 .and. myrank == 0 ) &
!                  write(*,'(1x,3i6,2f15.10)') n,m,k,c1,esp(m,1,1)
!             do i=1,nrnk
!                if ( omax(i,k,n) < c1 ) then
!                   do j=nrnk,i+1,-1
!                      omax(j,k,n)=omax(j-1,k,n)
!                      indx(j,k,n)=indx(j-1,k,n)
!                   end do
!                   omax(i,k,n)=c1
!                   indx(i,k,n)=m
!                   exit
!                end if
!             end do ! i
          end do ! k
       end do ! m

!       i1=indx(1,1,n)
!       i2=indx(2,1,n)
!       e1=esp(i1,1,1)
!       e2=esp(i2,1,1)
!       c1=omax(1,1,n)
!       c2=omax(2,1,n)

!       if ( abs(e2-e1) > 1.d-2 ) cycle 

!       if ( ichk(i1) /= 0 .or. ichk(i2) /= 0 .or. ichk(n) /= 0 ) cycle
!       ichk(i1)=ichk(i1)+1
!       ichk(i2)=ichk(i2)+1
!       ichk(n )=ichk(n )+1
!       icount=icount+1

       if ( myrank == 0 ) then
!          write(*,*) c2-c1,e2-e1,icount,count(ichk/=0)
!          write(*,'(1x,4i6,2f15.10)') n,0,1,n,omax(1,0,n),esp(n,1,1)
!          do k=1,nk/2
!             do i=1,nrnk
!                m=indx(i,k,n)
!                if ( m /= 0 ) then
!                   write(*,'(1x,4i6,2f15.10)') n,k,i,m,omax(i,k,n),esp(m,1,1)
!                else
!                   write(*,'(1x,4i6,2f15.10)') n,k,i,m,omax(i,k,n)
!                end if
!             end do
!          end do
!          write(*,*)
!          write(*,*)

!          do i3=0,Ngrid(3)-1
!             write(100,*) i3,i3*Hgrid(3),f1(i3,n)
!          end do
!          write(100,*)
!          write(100,*)

       end if

    end do ! n

    if ( myrank == 0 ) then
    do m=1,MB
       do n=1,MB
       do i=1,nrnk
          if ( indx(i,1,n) == m ) write(*,*) m,n,omax(i,1,n) 
       end do
       end do
    end do
    end if

    deallocate( f1,f )
    deallocate( ichk )
    deallocate( omax )
    deallocate( indx )
    deallocate( ovlp,ovlp0 )

    write(*,*) "end bandanalysis",myrank
    return

  END SUBROUTINE band_analysis

END MODULE band_analysis_module
