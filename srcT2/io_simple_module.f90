MODULE io_simple_module

  use aa_module, only: aa
  use rgrid_module, only: Ngrid,Igrid
  use parallel_module, only: comm_grid,myrank
  use rsdft_mpi_module

  implicit none

  PRIVATE
  PUBLIC :: write_io_simple

CONTAINS


  SUBROUTINE write_io_simple( f, file_name_in )

    implicit none
    real(8),intent(IN) :: f(:)
    character(*),optional,intent(IN) :: file_name_in
    real(8),allocatable :: w(:,:,:)
    integer :: i,i1,i2,i3
    integer,parameter :: u=17
    character(40) :: file_name

    call write_border( 0, " io_simple(start)" )

    allocate( w(0:Ngrid(1)-1,0:Ngrid(2)-1,0:Ngrid(3)-1) ) ; w=0.0d0

    i=0
    do i3=Igrid(1,3),Igrid(2,3)
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
       i=i+1
       w(i1,i2,i3)=f(i)
    end do
    end do
    end do

    call rsdft_allreduce_sum( w, comm_grid )

    file_name="io_simple.dat"
    if ( present(file_name_in) ) file_name=file_name_in

    if ( myrank == 0 ) then
       open(u,file=file_name,form='unformatted')
       write(u) Ngrid(1:3)
       write(u) aa(:,:)
       write(u) w(:,:,:)
       close(u)
    end if

    deallocate( w )

    call write_border( 0, " io_simple(end)" )

  END SUBROUTINE write_io_simple


END MODULE io_simple_module
