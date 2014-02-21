MODULE rgrid_mol_module

  implicit none

  PRIVATE
  PUBLIC :: LL,KK,Hsize &
           ,read_rgrid_mol &
           ,init_rgrid_mol,mesh_div_1 &
           ,construct_rgrid_mol,destruct_rgrid_mol &
           ,construct_boundary_rgrid_mol,destruct_boundary_rgrid_mol &
           ,map_g2p_rgrid_mol

  integer :: Box_Shape
  real(8) :: Hsize,Rsize,Zsize
  integer,allocatable :: LL(:,:),KK(:,:)

  real(8),parameter :: eps = 1.d-10

CONTAINS


  SUBROUTINE read_rgrid_mol(rank,unit)
    integer,intent(IN)  :: rank,unit
    if ( rank == 0 ) then
       read(unit,*) Box_Shape
       select case(Box_Shape)
       case(1)
          read(unit,*) Hsize,Rsize
       case(2)
          read(unit,*) Hsize,Rsize,Zsize
       end select
       write(*,*) "Box_Shape=",Box_Shape
       write(*,*) "Hsize=",Hsize
       write(*,*) "Rsize=",Rsize
       write(*,*) "Zsize=",Zsize
    end if
    call send_rgrid_mol(rank)
  END SUBROUTINE read_rgrid_mol

  SUBROUTINE send_rgrid_mol(rank)
    integer,intent(IN)  :: rank
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(Box_Shape,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Hsize,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Rsize,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Zsize,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_rgrid_mol


  SUBROUTINE init_rgrid_mol(Ngrid,Hgrid,aa,bb,disp_switch)
    integer,intent(OUT) :: Ngrid(3)
    real(8),intent(OUT) :: Hgrid(3),aa(3,3),bb(3,3)
    logical,optional,intent(IN)  :: disp_switch
    call init_mesh_rsmol( Ngrid(1),Ngrid(2),Ngrid(3) )
    Hgrid(1)=Hsize
    Hgrid(2)=Hsize
    Hgrid(3)=Hsize
    aa(:,:)=0.d0
    bb(:,:)=0.d0
    aa(1,1)=Hsize*Ngrid(1)
    aa(2,2)=Hsize*Ngrid(2)
    aa(3,3)=Hsize*Ngrid(3)
    bb(1,1)=2.d0*acos(-1.d0)/aa(1,1)
    bb(2,2)=2.d0*acos(-1.d0)/aa(2,2)
    bb(3,3)=2.d0*acos(-1.d0)/aa(3,3)
    if ( present(disp_switch) ) then
       if ( disp_switch ) then
       write(*,*) "Ngrid(1:3)=",Ngrid(1:3)
       write(*,*) "Hgrid(1:3)=",Hgrid(1:3)
       write(*,*) "aa"
       write(*,*) aa(1:3,1)
       write(*,*) aa(1:3,2)
       write(*,*) aa(1:3,3)
       write(*,*) "bb"
       write(*,*) bb(1:3,1)
       write(*,*) bb(1:3,2)
       write(*,*) bb(1:3,3)
       end if
    end if
  END SUBROUTINE init_rgrid_mol

  SUBROUTINE init_mesh_rsmol(mx,my,mz)
    integer,intent(OUT) :: mx,my,mz
    select case(Box_Shape)
    case(1)
       mx = nint(Rsize/Hsize)
       if ( mx*Hsize < Rsize ) mx=mx+1
       my = mx
       mz = mx
    case(2)
       mx = nint(Rsize/Hsize)
       if ( mx*Hsize < Rsize ) mx=mx+1
       my = mx
       mz = nint(Zsize/Hsize)
       if ( mz*Hsize < Zsize ) mz=mz+1
    end select
    mx = 2*mx+1
    my = 2*my+1
    mz = 2*mz+1
  END SUBROUTINE init_mesh_rsmol


  SUBROUTINE mesh_div_1(np_grid,np,Ngrid,pinfo_grid,disp_switch)
    implicit none
    integer,intent(IN)  :: np_grid,np(3),Ngrid(3)
    logical,intent(IN)  :: disp_switch
    integer,intent(OUT) :: pinfo_grid(8,0:np_grid-1)
    integer,parameter :: max_loop=600
    integer :: iloop,ix,iy,iz,i1,i2,i3,m,n,mx0,my0,mz0
    integer :: mx,my,mz,i,j
    integer,allocatable :: itmp(:),mxp(:),myp(:),mzp(:)
    real(8) :: r2,Rc2,z,ave,var

    allocate( itmp(0:np_grid-1) ) ; itmp=0
    allocate( mxp(np(1)),myp(np(2)),mzp(np(3)) )

    mxp(:) = Ngrid(1)/np(1)
    myp(:) = Ngrid(2)/np(2)
    mzp(:) = Ngrid(3)/np(3)
    mxp(1) = mxp(1)+Ngrid(1)-sum(mxp)
    myp(1) = myp(1)+Ngrid(2)-sum(myp)
    mzp(1) = mzp(1)+Ngrid(3)-sum(mzp)

    Rc2 = Rsize**2

    mx = (Ngrid(1)-1)/2
    my = (Ngrid(2)-1)/2
    mz = (Ngrid(3)-1)/2

    do iloop=1,max_loop

       itmp(:)=0

       n=-1

       mz0=-mz
       do i3=1,np(3)
       my0=-my
       do i2=1,np(2)
       mx0=-mx
       do i1=1,np(1)

          n=n+1

          select case(Box_Shape)
          case(1)
             do iz=mz0,mz0+mzp(i3)-1
             do iy=my0,my0+myp(i2)-1
             do ix=mx0,mx0+mxp(i1)-1
                r2=Hsize*Hsize*(ix*ix+iy*iy+iz*iz)
                if( r2 <= Rc2+eps )then
                   itmp(n)=itmp(n)+1
                end if
             end do
             end do
             end do
          case(2)
             do iz=mz0,mz0+mzp(i3)-1
                z=iz*Hsize
                if ( abs(z) > Zsize+eps ) cycle
                do iy=my0,my0+myp(i2)-1
                do ix=mx0,mx0+mxp(i1)-1
                   r2=Hsize*Hsize*(ix*ix+iy*iy)
                   if ( r2 <= Rc2+eps ) then
                      itmp(n)=itmp(n)+1
                   end if
                end do
                end do
             end do
          end select

       mx0=mx0+mxp(i1)
       end do
       my0=my0+myp(i2)
       end do
       mz0=mz0+mzp(i3)
       end do

       if ( iloop == max_loop ) exit

       i=0 ; j=0
       do n=0,np_grid-1
          if ( itmp(i) < itmp(n) ) i=n
          if ( itmp(j) > itmp(n) ) j=n
       end do

       n=-1
       do i3=1,np(3)
       do i2=1,np(2)
       do i1=1,np(1)
          n=n+1
          select case( mod(iloop,3) )
          case(0)
             if ( n == i ) mxp(i1)=mxp(i1)-1
             if ( n == j ) mxp(i1)=mxp(i1)+1
          case(1)
             if ( n == i ) myp(i2)=myp(i2)-1
             if ( n == j ) myp(i2)=myp(i2)+1
          case(2)
             if ( n == i ) mzp(i3)=mzp(i3)-1
             if ( n == j ) mzp(i3)=mzp(i3)+1
          end select
       end do
       end do
       end do

    end do ! iloop

    n=-1
    do i3=1,np(3)
    do i2=1,np(2)
    do i1=1,np(1)
       n=n+1
       pinfo_grid(1,n) = sum( mxp(1:i1) ) - mxp(i1) - mx - 1
       pinfo_grid(2,n) = mxp(i1)
       pinfo_grid(3,n) = sum( myp(1:i2) ) - myp(i2) - my - 1
       pinfo_grid(4,n) = myp(i2)
       pinfo_grid(5,n) = sum( mzp(1:i3) ) - mzp(i3) - mz - 1
       pinfo_grid(6,n) = mzp(i3)
       pinfo_grid(7,n) = sum( itmp(0:n) ) - itmp(n)
       pinfo_grid(8,n) = itmp(n)
    end do
    end do
    end do

    if ( disp_switch ) then
       write(*,'(1x,a5,2x,3a4,2x,a8)') "rank","mxp","myp","mzp","#"
       n=-1
       do i3=1,np(3)
       do i2=1,np(2)
       do i1=1,np(1)
          n=n+1
!          write(*,'(1x,i5,2x,3i4,2x,i8)') n,mxp(i1),myp(i2),mzp(i3),itmp(n)
          pinfo_grid(2:6:2,n)=pinfo_grid(2:6:2,n)+pinfo_grid(1:5:2,n)
          pinfo_grid(1:5:2,n)=pinfo_grid(1:5:2,n)+1
          write(*,'(1x,i5,2x,6i4,2x,2i8)') n,pinfo_grid(1:8,n)
          pinfo_grid(1:5:2,n)=pinfo_grid(1:5:2,n)-1
          pinfo_grid(2:6:2,n)=pinfo_grid(2:6:2,n)-pinfo_grid(1:5:2,n)
       end do
       end do
       end do
    end if

    deallocate( mzp,myp,mxp )
    deallocate( itmp )

    return
  END SUBROUTINE mesh_div_1


  SUBROUTINE construct_rgrid_mol
    use rgrid_module, only: Igrid
    use array_bound_module, only: ML_0,ML_1
    integer :: i,ix,iy,iz
    real(8) :: r2,Rc2,z,x,y
    if ( allocated(LL) ) deallocate(LL)
    allocate( LL(3,ML_0:ML_1) ) ; LL=0
    Rc2=Rsize*Rsize
    select case(Box_Shape)
    case(1)
       i=ML_0-1
       do iz=Igrid(1,3),Igrid(2,3)
       do iy=Igrid(1,2),Igrid(2,2)
       do ix=Igrid(1,1),Igrid(2,1)
          r2=Hsize*Hsize*(ix*ix+iy*iy+iz*iz)
          if( r2 <= Rc2+eps )then
             i=i+1
             LL(1,i)=ix
             LL(2,i)=iy
             LL(3,i)=iz
          end if
       end do
       end do
       end do
    case(2)
       i=ML_0-1
       do iz=Igrid(1,3),Igrid(2,3)
          z=iz*Hsize
          if ( abs(z) > Zsize+eps ) cycle
          do iy=Igrid(1,2),Igrid(2,2)
          do ix=Igrid(1,1),Igrid(2,1)
             r2=Hsize*Hsize*(ix*ix+iy*iy)
             if ( r2 <= Rc2+eps ) then
                i=i+1
                LL(1,i)=ix
                LL(2,i)=iy
                LL(3,i)=iz
             end if
          end do
          end do
       end do
    end select
  END SUBROUTINE construct_rgrid_mol

  SUBROUTINE destruct_rgrid_mol
    if ( allocated(LL) ) deallocate(LL)
  END SUBROUTINE destruct_rgrid_mol


  SUBROUTINE construct_boundary_rgrid_mol(Md)
    use array_bound_module, only: ML_0,ML_1
    integer,intent(IN) :: Md
    integer :: i,ix,iy,iz,m,m1,m2,m3,m4,m5,m6
    integer,allocatable :: icheck(:,:,:)
    real(8) :: z,r2,Rc2,H2
    m1 = minval( LL(1,:) ) - Md
    m2 = maxval( LL(1,:) ) + Md
    m3 = minval( LL(2,:) ) - Md
    m4 = maxval( LL(2,:) ) + Md
    m5 = minval( LL(3,:) ) - Md
    m6 = maxval( LL(3,:) ) + Md
    allocate( icheck(m1:m2,m3:m4,m5:m6) ) ; icheck=0
    Rc2 = Rsize*Rsize
    H2  = Hsize*Hsize
    do m=-Md,Md
       if ( m == 0 ) cycle
       select case(Box_Shape)
       case(1)
          do i=ML_0,ML_1
             ix = LL(1,i)
             iy = LL(2,i)
             iz = LL(3,i)
             r2 = H2*( (ix+m)*(ix+m) + iy*iy + iz*iz )
             if ( r2 > Rc2+eps ) icheck(ix+m,iy,iz)=icheck(ix+m,iy,iz)+1
             r2 = H2*( ix*ix + (iy+m)*(iy+m) + iz*iz )
             if ( r2 > Rc2+eps ) icheck(ix,iy+m,iz)=icheck(ix,iy+m,iz)+1
             r2 = H2*( ix*ix + iy*iy + (iz+m)*(iz+m) )
             if ( r2 > Rc2+eps ) icheck(ix,iy,iz+m)=icheck(ix,iy,iz+m)+1
          end do
       case(2)
          do i=ML_0,ML_1
             ix = LL(1,i)
             iy = LL(2,i)
             iz = LL(3,i)
             z  = iz*Hsize
             r2 = H2*( (ix+m)*(ix+m) + iy*iy )
             if ( r2 > Rc2+eps .and. abs(z) <= Zsize+eps ) icheck(ix+m,iy,iz)=icheck(ix+m,iy,iz)+1
             r2 = H2*( ix*ix + (iy+m)*(iy+m) )
             if ( r2 > Rc2+eps .and. abs(z) <= Zsize+eps ) icheck(ix,iy+m,iz)=icheck(ix,iy+m,iz)+1
             z  = (iz+m)*Hsize
             r2 = H2*( ix*ix + iy*iy )
             if ( r2 <= Rc2+eps .and. abs(z) > Zsize+eps ) icheck(ix,iy,iz+m)=icheck(ix,iy,iz+m)+1
          end do
       end select
    end do

    if ( allocated(KK) ) deallocate(KK)

    m = count( icheck(m1:m2,m3:m4,m5:m6) > 0 )
    allocate( KK(3,m) )
    KK=0

    i=0
    do iz=m5,m6
    do iy=m3,m4
    do ix=m1,m2
       if ( icheck(ix,iy,iz) > 0 ) then
          i=i+1
          KK(1,i)=ix
          KK(2,i)=iy
          KK(3,i)=iz
       end if
    end do
    end do
    end do

  END SUBROUTINE construct_boundary_rgrid_mol

  SUBROUTINE destruct_boundary_rgrid_mol
    if ( allocated(KK) ) deallocate(KK)
  END SUBROUTINE destruct_boundary_rgrid_mol


  SUBROUTINE map_g2p_rgrid_mol(a1,b1,a2,b2,a3,b3,map_g2p,np,pinfo_grid)
    integer,intent(IN)  :: a1,b1,a2,b2,a3,b3,np,pinfo_grid(8,0:np-1)
    integer,intent(OUT) :: map_g2p(a1:b1,a2:b2,a3:b3)
    integer :: n,i1,i2,i3
    real(8) :: r2,Rc2,z
    Rc2=Rsize*Rsize
    map_g2p(:,:,:)=-1
    do n=0,np-1
       do i3=pinfo_grid(5,n)+1,pinfo_grid(5,n)+pinfo_grid(6,n)
       do i2=pinfo_grid(3,n)+1,pinfo_grid(3,n)+pinfo_grid(4,n)
       do i1=pinfo_grid(1,n)+1,pinfo_grid(1,n)+pinfo_grid(2,n)
          select case(Box_Shape)
          case(1)
             r2=Hsize*Hsize*(i1*i1+i2*i2+i3*i3)
             if ( r2 <= Rc2+eps ) then
                map_g2p(i1,i2,i3)=n
             end if
          case(2)
             r2=Hsize*Hsize*(i1*i1+i2*i2)
             z=abs(i3*Hsize)
             if ( r2 <= Rc2+eps .and. z <= Zsize+eps ) then
                map_g2p(i1,i2,i3)=n
             end if
          end select
       end do
       end do
       end do
    end do
  END SUBROUTINE map_g2p_rgrid_mol


END MODULE rgrid_mol_module
