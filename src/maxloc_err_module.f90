MODULE maxloc_err_module

  use grid_module, only: convert_1d_to_3d_grid

  implicit none

  PRIVATE
  PUBLIC :: show_maxloc_err

  logical,PUBLIC :: flag_maxloc_err=.false.

CONTAINS


  SUBROUTINE show_maxloc_err( f, g, n )

    implicit none
    real(8),intent(IN) :: f(:), g(:)
    integer,intent(IN) :: n
    integer :: m,i,j
    integer,allocatable :: loc(:),loc3(:,:)
    logical,allocatable :: msk(:)
    logical :: disp
    real(8),allocatable :: err(:)

    if ( .not.flag_maxloc_err ) return

    allocate( msk(size(f)) ) ; msk=.true.
    allocate( loc(n)    ) ; loc=0
    allocate( loc3(3,n) ) ; loc3=0
    allocate( err(n)    ) ; err=0.0d0

    do j=1,n
       loc(j) = maxloc( abs(f-g), 1, msk )
       err(j) = abs( f(loc(j)) - g(loc(j)) )
       msk(loc(j))=.false.
    end do

    do j=1,n
       call convert_1d_to_3d_grid( loc(j), loc3(:,j) )
    end do

    call check_disp_switch( disp, 0 )
    if ( disp ) then
       write(*,*) "--- grid points of maximum error ---"
       do j=1,n
          write(*,'(1x,i1," (",3i4,")",i8,es14.5)') j,loc3(:,j),loc(j),err(j)
       end do
    end if

    deallocate( err  )
    deallocate( loc3 )
    deallocate( loc  )
    deallocate( msk  )

  END SUBROUTINE show_maxloc_err


END MODULE maxloc_err_module
