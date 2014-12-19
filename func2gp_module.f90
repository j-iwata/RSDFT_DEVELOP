MODULE func2gp_module

  use esm_rgrid_module
  use rgrid_module
  use rgrid_mol_module, LL_mol => LL
  use parallel_module
  
  implicit none

  PRIVATE
  PUBLIC :: func2gp_c_esm,func2gp_r_esm,func2gp_r,func2gp_r_mol

CONTAINS


  SUBROUTINE func2gp_c_esm(unit,n1,n2,f)
    implicit none
    integer,intent(IN) :: unit,n1,n2
    complex(8),intent(IN) :: f(n1:n2)
    integer :: i,i1,i2,i3,a1,a2,a3,b1,b2,b3,ierr,m
    real(8),allocatable :: w(:,:,:,:)
    real(8),parameter :: zero=(0.d0,0.d0)

    a1 = -Ngrid(1)/2
    b1 =  Ngrid(1)/2+1
    a2 = -Ngrid(2)/2
    b2 =  Ngrid(2)/2+1
    a3 = -Ngrid(3)/2
    b3 =  Ngrid(3)/2+1

    allocate( w(a1:b1,a2:b2,a3:b3,2) )
    w=zero

    do i=n1,n2
       i1=LL_ESM(1,i)
       i2=LL_ESM(2,i)
       i3=LL_ESM(3,i)
       w(i1,i2,i3,1) = abs( f(i) )**2
    end do

    m=size( w(:,:,:,1) )
    call mpi_allreduce(w(:,:,:,1),w(:,:,:,2),m,mpi_real8,mpi_sum,comm_grid,ierr)

    if ( myrank == 0 ) then
       rewind unit
       do i=a1,b1
          write(unit,'(1x,4g20.10)') i*Hgrid(1),w(i,0,0,2),w(0,i,0,2),w(0,0,i,2)
       end do
    end if

    deallocate( w )

  END SUBROUTINE func2gp_c_esm

  SUBROUTINE func2gp_r_esm(unit,n1,n2,f)
    implicit none
    integer,intent(IN) :: unit,n1,n2
    real(8),intent(IN) :: f(n1:n2)
    integer :: i,i1,i2,i3,a1,a2,a3,b1,b2,b3,ierr,m
    real(8),allocatable :: w(:,:,:,:)
    real(8),parameter :: zero=(0.d0,0.d0)

    a1 = -Ngrid(1)/2
    b1 =  Ngrid(1)/2+1
    a2 = -Ngrid(2)/2
    b2 =  Ngrid(2)/2+1
    a3 = -Ngrid(3)/2
    b3 =  Ngrid(3)/2+1

    allocate( w(a1:b1,a2:b2,a3:b3,2) )
    w=zero

    do i=n1,n2
       i1=LL_ESM(1,i)
       i2=LL_ESM(2,i)
       i3=LL_ESM(3,i)
       w(i1,i2,i3,1) = f(i) 
    end do

    m=size( w(:,:,:,1) )
    call mpi_allreduce(w(:,:,:,1),w(:,:,:,2),m,mpi_real8,mpi_sum,comm_grid,ierr)

    if ( myrank == 0 ) then
       rewind unit
       do i=a1,b1
          write(unit,'(1x,4g20.10)') i*Hgrid(1),w(i,0,0,2),w(0,i,0,2),w(0,0,i,2)
       end do
    end if

    deallocate( w )

  END SUBROUTINE func2gp_r_esm


  SUBROUTINE func2gp_r(unit,n1,n2,f)
    implicit none
    integer,intent(IN) :: unit,n1,n2
    real(8),intent(IN) :: f(n1:n2)
    integer :: i,i1,i2,i3,a1,a2,a3,b1,b2,b3,ierr,m
    real(8),allocatable :: w(:,:,:,:)
    real(8),parameter :: zero=0.0d0

    a1 = 0
    b1 = Ngrid(1)-1
    a2 = 0
    b2 = Ngrid(2)-1
    a3 = 0
    b3 = Ngrid(3)-1

    allocate( w(a1:b1,a2:b2,a3:b3,2) )
    w=zero

    i=n1-1
    do i3=Igrid(1,3),Igrid(2,3)
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
       i=i+1
       w(i1,i2,i3,1) = f(i) 
    end do
    end do
    end do

    m=size( w(:,:,:,1) )
    call mpi_allreduce(w(:,:,:,1),w(:,:,:,2),m,mpi_real8,mpi_sum,comm_grid,ierr)

    if ( myrank == 0 ) then
       rewind unit
       do i=a1,b1
          write(unit,'(1x,4g20.10)') i*Hgrid(1),w(i,0,0,2),sum(w(i,:,:,2))
       end do
       write(unit,*)
       write(unit,*)
       do i=a2,b2
          write(unit,'(1x,4g20.10)') i*Hgrid(2),w(0,i,0,2),sum(w(:,i,:,2))
       end do
       write(unit,*)
       write(unit,*)
       do i=a3,b3
          write(unit,'(1x,4g20.10)') i*Hgrid(3),w(0,0,i,2),sum(w(:,:,i,2))
       end do
    end if

    deallocate( w )

  END SUBROUTINE func2gp_r


  SUBROUTINE func2gp_r_mol(unit,n1,n2,f)
    implicit none
    integer,intent(IN) :: unit,n1,n2
    real(8),intent(IN) :: f(n1:n2)
    integer :: i,i1,i2,i3,a1,a2,a3,b1,b2,b3,ierr,m
    real(8),allocatable :: w(:,:,:,:)
    real(8),parameter :: zero=0.0d0

    a1 =-Ngrid(1)
    b1 = Ngrid(1)
    a2 =-Ngrid(2)
    b2 = Ngrid(2)
    a3 =-Ngrid(3)
    b3 = Ngrid(3)

    allocate( w(a1:b1,a2:b2,a3:b3,2) )
    w=zero

    do i=n1,n2
       w(LL_mol(1,i),LL_mol(2,i),LL_mol(3,i),1) = f(i) 
    end do

    m=size( w(:,:,:,1) )
    call mpi_allreduce(w(:,:,:,1),w(:,:,:,2),m,mpi_real8,mpi_sum,comm_grid,ierr)

    if ( myrank == 0 ) then
       rewind unit
       do i=a1,b1
          write(unit,'(1x,4g20.10)') i*Hgrid(1),w(i,0,0,2),sum(w(i,:,:,2))
       end do
       write(unit,*)
       write(unit,*)
       do i=a2,b2
          write(unit,'(1x,4g20.10)') i*Hgrid(2),w(0,i,0,2),sum(w(:,i,:,2))
       end do
       write(unit,*)
       write(unit,*)
       do i=a3,b3
          write(unit,'(1x,4g20.10)') i*Hgrid(3),w(0,0,i,2),sum(w(:,:,i,2))
       end do
    end if

    deallocate( w )

  END SUBROUTINE func2gp_r_mol


END MODULE func2gp_module
