!--------1---------2---------3---------4---------5---------6---------7--

      SUBROUTINE init_fdparallel_mol
      use global_variables
      implicit none

!      integer,allocatable :: Lxyz(:,:,:)
      integer :: ix,iy,iz !,j,k,n,m,a,il,ik,ierr,n1,n2,iloop
      integer :: i,i1,i2,i3,n,j !,k,m,n
      real(8) :: H,r2,Rc2,z !,r,x,y,r0
      integer :: a1,a2,a3,b1,b2,b3
      integer,allocatable :: map_grid_2_pinfo(:,:,:,:)
      integer,allocatable :: neighbor_info(:,:,:)
      real(8) :: ctime0,ctime1,etime0,etime1
      real(8) :: mem,memax
      real(8),parameter :: eps=1.d-10
!      integer :: Mx,My,Mz
      integer :: ix0,ix1,iy0,iy1,iz0,iz1,jx,jy,jz
!      integer :: ntmp,i1,i2,i3,ilp,ikp
!      integer :: Mx00,My00,Mz00,Mx11,My11,Mz11
!      integer :: mxi,myi,mzi
!      integer,allocatable :: Mxp(:),Myp(:),Mzp(:)
      integer,allocatable :: itmp(:),jtmp(:)
!      integer,allocatable :: ip(:,:,:)
!      integer,allocatable :: itmp2(:,:)
!      integer :: sp,fp,ns,jrank
      integer :: irank
!      integer :: ireqs(18),ireqr(18),nreq,ireq(100)
!      integer :: itags,itagr,istatus(mpi_status_size)
!      integer :: icheck_1,icheck_2,icheck_3,icheck_4

      call watch(ctime0,etime0)

      mem   = 0.d0
      memax = 0.d0
      Rc2   = Rsize**2
      H     = H1

      a1=MLI(1,1)-Md ; a2=MLI(1,2)-Md ; a3=MLI(1,3)-Md
      b1=MLI(2,1)+Md ; b2=MLI(2,2)+Md ; b3=MLI(2,3)+Md

      allocate( map_grid_2_pinfo(a1:b1,a2:b2,a3:b3,2) ) ; map_grid_2_pinfo(:,:,:,:)=-1
      mem=mem+bsintg*(b1-a1+1)*(b2-a2+1)*(b3-a3+1)*2 ; memax=max(mem,memax)

      do n=0,np_grid-1
         a1=pinfo_grid_2d(1,n) ; b1=a1+pinfo_grid_2d(2,n)-1
         a2=pinfo_grid_2d(3,n) ; b2=a2+pinfo_grid_2d(4,n)-1
         a3=pinfo_grid_2d(5,n) ; b3=a3+pinfo_grid_2d(6,n)-1
         do i3=a3,b3
         do i2=a2,b2
         do i1=a1,b1
            select case(SYStype)
            case(1)
               r2=H*H*(i1*i1+i2*i2+i3*i3)
               if ( r2<Rc2+eps ) then
                  map_grid_2_pinfo(i1,i2,i3,1)=n
               end if
            case(2)
               r2=H*H*(i1*i1+i2*i2)
               z=abs(i3*H)
               if ( r2<Rc2+eps .and. z<Zsize+eps ) then
                  map_grid_2_pinfo(i1,i2,i3,1)=n
               end if
            end select
         end do
         end do
         end do
      end do

!
! --- FD Boundary conditions ( eqdiv ) ---
!

!- allocate -----------------------------
      call gv_alloc("n_neighbor")
!----------------------------------------
!- allocate ----------------------------------------------------
      allocate( neighbor_info(8,6,100) )
      mem=mem+bsintg*size(neighbor_info) ; memax=max(mem,memax)
      allocate( itmp(0:np_grid-1) ) ; itmp=-1
      allocate( jtmp(0:np_grid-1) ) ; jtmp=0
      mem=mem+bsintg*np_grid*2 ; memax=max(mem,memax)
!---------------------------------------------------------------


      N_neighbor(1:6)=0
      neighbor_info(1:8,1:6,1:100)=0
      neighbor_info(1:2,:,:)=MPI_PROC_NULL


      iy0=    pinfo_grid_2d(3,myrank_g)
      iy1=iy0+pinfo_grid_2d(4,myrank_g)-1
      iz0=    pinfo_grid_2d(5,myrank_g)
      iz1=iz0+pinfo_grid_2d(6,myrank_g)-1

      ix0=pinfo_grid_2d(1,myrank_g)-Md
      ix1=pinfo_grid_2d(1,myrank_g)-1
      itmp(:)=-1
      jtmp(:)=0
      do ix=ix0,ix1
         if ( ix<-ML1 .or. ML1<ix ) cycle
         do iz=iz0,iz1
         do iy=iy0,iy1
            irank=map_grid_2_pinfo(ix,iy,iz,1)
            if ( 0<=irank .and. irank<=np_grid-1 ) then
               jtmp(irank)=jtmp(irank)+1
            end if
         end do
         end do
      end do
      n=0
      do i=0,np_grid-1
         if ( jtmp(i)>0 ) then
            n=n+1
            itmp(n)=i
         end if
      end do
      N_neighbor(1)=n
      do i=1,n
         irank=itmp(i)
         ix=max( ix0, pinfo_grid_2d(1,irank) )
         jx=min( ix1, pinfo_grid_2d(1,irank)+pinfo_grid_2d(2,irank)-1 )
         iy=max( iy0, pinfo_grid_2d(3,irank) )
         jy=min( iy1, pinfo_grid_2d(3,irank)+pinfo_grid_2d(4,irank)-1 )
         iz=max( iz0, pinfo_grid_2d(5,irank) )
         jz=min( iz1, pinfo_grid_2d(5,irank)+pinfo_grid_2d(6,irank)-1 )
         neighbor_info( 1, 1, i )=myrank_g
         neighbor_info( 2, 1, i )=irank
         neighbor_info( 3, 1, i )=ix
         neighbor_info( 4, 1, i )=jx
         neighbor_info( 5, 1, i )=iy
         neighbor_info( 6, 1, i )=jy
         neighbor_info( 7, 1, i )=iz
         neighbor_info( 8, 1, i )=jz
      end do

      ix0=pinfo_grid_2d(1,myrank_g)+pinfo_grid_2d(2,myrank_g)
      ix1=ix0+Md-1
      itmp(:)=-1
      jtmp(:)=0
      do ix=ix0,ix1
         if ( ix<-ML1 .or. ML1<ix ) cycle
         do iz=iz0,iz1
         do iy=iy0,iy1
            irank=map_grid_2_pinfo(ix,iy,iz,1)
            if ( irank>=0 ) then
               jtmp(irank)=jtmp(irank)+1
            end if
         end do
         end do
      end do
      n=0
      do i=0,np_grid-1
         if ( jtmp(i)>0 ) then
            n=n+1
            itmp(n)=i
         end if
      end do
      N_neighbor(2)=n
      do i=1,n
         irank=itmp(i)
         ix=max( ix0, pinfo_grid_2d(1,irank) )
         jx=min( ix1, pinfo_grid_2d(1,irank)+pinfo_grid_2d(2,irank)-1 )
         iy=max( iy0, pinfo_grid_2d(3,irank) )
         jy=min( iy1, pinfo_grid_2d(3,irank)+pinfo_grid_2d(4,irank)-1 )
         iz=max( iz0, pinfo_grid_2d(5,irank) )
         jz=min( iz1, pinfo_grid_2d(5,irank)+pinfo_grid_2d(6,irank)-1 )
         neighbor_info( 1, 2, i )=myrank_g
         neighbor_info( 2, 2, i )=irank
         neighbor_info( 3, 2, i )=ix
         neighbor_info( 4, 2, i )=jx
         neighbor_info( 5, 2, i )=iy
         neighbor_info( 6, 2, i )=jy
         neighbor_info( 7, 2, i )=iz
         neighbor_info( 8, 2, i )=jz
      end do

      ix0=    pinfo_grid_2d(1,myrank_g)
      ix1=ix0+pinfo_grid_2d(2,myrank_g)-1
      iz0=    pinfo_grid_2d(5,myrank_g)
      iz1=iz0+pinfo_grid_2d(6,myrank_g)-1

      iy0=pinfo_grid_2d(3,myrank_g)-Md
      iy1=pinfo_grid_2d(3,myrank_g)-1
      itmp(:)=-1
      jtmp(:)=0
      do iy=iy1,iy0,-1
         if ( iy<-ML2 .or. ML2<iy ) cycle
         do iz=iz0,iz1
         do ix=ix0,ix1
            irank=map_grid_2_pinfo(ix,iy,iz,1)
            if ( irank>=0 ) then
               jtmp(irank)=jtmp(irank)+1
            end if
         end do
         end do
      end do
      n=0
      do i=0,np_grid-1
         if ( jtmp(i)>0 ) then
            n=n+1
            itmp(n)=i
         end if
      end do
      N_neighbor(3)=n
      do i=1,n
         irank=itmp(i)
         ix=max( ix0, pinfo_grid_2d(1,irank) )
         jx=min( ix1, pinfo_grid_2d(1,irank)+pinfo_grid_2d(2,irank)-1 )
         iy=max( iy0, pinfo_grid_2d(3,irank) )
         jy=min( iy1, pinfo_grid_2d(3,irank)+pinfo_grid_2d(4,irank)-1 )
         iz=max( iz0, pinfo_grid_2d(5,irank) )
         jz=min( iz1, pinfo_grid_2d(5,irank)+pinfo_grid_2d(6,irank)-1 )
         neighbor_info( 1, 3, i )=myrank_g
         neighbor_info( 2, 3, i )=irank
         neighbor_info( 3, 3, i )=ix
         neighbor_info( 4, 3, i )=jx
         neighbor_info( 5, 3, i )=iy
         neighbor_info( 6, 3, i )=jy
         neighbor_info( 7, 3, i )=iz
         neighbor_info( 8, 3, i )=jz
      end do

      iy0=pinfo_grid_2d(3,myrank_g)+pinfo_grid_2d(4,myrank_g)
      iy1=iy0+Md-1
      itmp(:)=-1
      jtmp(:)=0
      do iy=iy0,iy1
         if ( iy<-ML2 .or. ML2<iy ) cycle
         do iz=iz0,iz1
         do ix=ix0,ix1
            irank=map_grid_2_pinfo(ix,iy,iz,1)
            if ( irank>=0 ) then
               jtmp(irank)=jtmp(irank)+1
            end if
         end do
         end do
      end do
      n=0
      do i=0,np_grid-1
         if ( jtmp(i)>0 ) then
            n=n+1
            itmp(n)=i
         end if
      end do
      N_neighbor(4)=n
      do i=1,n
         irank=itmp(i)
         ix=max( ix0, pinfo_grid_2d(1,irank) )
         jx=min( ix1, pinfo_grid_2d(1,irank)+pinfo_grid_2d(2,irank)-1 )
         iy=max( iy0, pinfo_grid_2d(3,irank) )
         jy=min( iy1, pinfo_grid_2d(3,irank)+pinfo_grid_2d(4,irank)-1 )
         iz=max( iz0, pinfo_grid_2d(5,irank) )
         jz=min( iz1, pinfo_grid_2d(5,irank)+pinfo_grid_2d(6,irank)-1 )
         neighbor_info( 1, 4, i )=myrank_g
         neighbor_info( 2, 4, i )=irank
         neighbor_info( 3, 4, i )=ix
         neighbor_info( 4, 4, i )=jx
         neighbor_info( 5, 4, i )=iy
         neighbor_info( 6, 4, i )=jy
         neighbor_info( 7, 4, i )=iz
         neighbor_info( 8, 4, i )=jz
      end do


      ix0=    pinfo_grid_2d(1,myrank_g)
      ix1=ix0+pinfo_grid_2d(2,myrank_g)-1
      iy0=    pinfo_grid_2d(3,myrank_g)
      iy1=iy0+pinfo_grid_2d(4,myrank_g)-1

      iz0=pinfo_grid_2d(5,myrank_g)-Md
      iz1=pinfo_grid_2d(5,myrank_g)-1
      itmp(:)=-1
      jtmp(:)=0
      do iz=iz1,iz0,-1
         if ( iz<-ML3 .or. ML3<iz ) cycle
         do iy=iy0,iy1
         do ix=ix0,ix1
            irank=map_grid_2_pinfo(ix,iy,iz,1)
            if ( irank>=0 ) then
               jtmp(irank)=jtmp(irank)+1
            end if
         end do
         end do
      end do
      n=0
      do i=0,np_grid-1
         if ( jtmp(i)>0 ) then
            n=n+1
            itmp(n)=i
         end if
      end do
      N_neighbor(5)=n
      do i=1,n
         irank=itmp(i)
         ix=max( ix0, pinfo_grid_2d(1,irank) )
         jx=min( ix1, pinfo_grid_2d(1,irank)+pinfo_grid_2d(2,irank)-1 )
         iy=max( iy0, pinfo_grid_2d(3,irank) )
         jy=min( iy1, pinfo_grid_2d(3,irank)+pinfo_grid_2d(4,irank)-1 )
         iz=max( iz0, pinfo_grid_2d(5,irank) )
         jz=min( iz1, pinfo_grid_2d(5,irank)+pinfo_grid_2d(6,irank)-1 )
         neighbor_info( 1, 5, i )=myrank_g
         neighbor_info( 2, 5, i )=irank
         neighbor_info( 3, 5, i )=ix
         neighbor_info( 4, 5, i )=jx
         neighbor_info( 5, 5, i )=iy
         neighbor_info( 6, 5, i )=jy
         neighbor_info( 7, 5, i )=iz
         neighbor_info( 8, 5, i )=jz
      end do

      iz0=pinfo_grid_2d(5,myrank_g)+pinfo_grid_2d(6,myrank_g)
      iz1=iz0+Md-1
      itmp(:)=-1
      jtmp(:)=0
      do iz=iz0,iz1
         if ( iz<-ML3 .or. ML3<iz ) cycle
         do iy=iy0,iy1
         do ix=ix0,ix1
            irank=map_grid_2_pinfo(ix,iy,iz,1)
            if ( irank>=0 ) then
               jtmp(irank)=jtmp(irank)+1
            end if
         end do
         end do
      end do
      n=0
      do i=0,np_grid-1
         if ( jtmp(i)>0 ) then
            n=n+1
            itmp(n)=i
         end if
      end do
      N_neighbor(6)=n
      do i=1,n
         irank=itmp(i)
         ix=max( ix0, pinfo_grid_2d(1,irank) )
         jx=min( ix1, pinfo_grid_2d(1,irank)+pinfo_grid_2d(2,irank)-1 )
         iy=max( iy0, pinfo_grid_2d(3,irank) )
         jy=min( iy1, pinfo_grid_2d(3,irank)+pinfo_grid_2d(4,irank)-1 )
         iz=max( iz0, pinfo_grid_2d(5,irank) )
         jz=min( iz1, pinfo_grid_2d(5,irank)+pinfo_grid_2d(6,irank)-1 )
         neighbor_info( 1, 6, i )=myrank_g
         neighbor_info( 2, 6, i )=irank
         neighbor_info( 3, 6, i )=ix
         neighbor_info( 4, 6, i )=jx
         neighbor_info( 5, 6, i )=iy
         neighbor_info( 6, 6, i )=jy
         neighbor_info( 7, 6, i )=iz
         neighbor_info( 8, 6, i )=jz
      end do

!- deallocate ------------------------------------------
      deallocate( jtmp,itmp ) ; mem=mem-bsintg*np_grid*2
!-------------------------------------------------------
!- allocate -----------------------
      call gv_alloc("fdinfo_send")
!----------------------------------

      fdinfo_send(7,:,:)=MPI_PROC_NULL
      fdinfo_recv(8,:,:)=MPI_PROC_NULL

      do i=1,6
         do n=1,n_neighbor(i)
            ix=neighbor_info( 3, i, n )
            jx=neighbor_info( 4, i, n )
            iy=neighbor_info( 5, i, n )
            jy=neighbor_info( 6, i, n )
            iz=neighbor_info( 7, i, n )
            jz=neighbor_info( 8, i, n )
            fdinfo_recv(1,n,i)=ix
            fdinfo_recv(2,n,i)=jx
            fdinfo_recv(3,n,i)=iy
            fdinfo_recv(4,n,i)=jy
            fdinfo_recv(5,n,i)=iz
            fdinfo_recv(6,n,i)=jz
            fdinfo_recv(7,n,i)=neighbor_info( 1, i, n ) !=myrank_g
            fdinfo_recv(8,n,i)=neighbor_info( 2, i, n ) !=irank
            fdinfo_recv(9,n,i)=(jx-ix+1)*(jy-iy+1)*(jz-iz+1)
            select case(i)
            case(1,2)
               fdinfo_recv(10,n,i)=(jy-iy+1)*(jz-iz+1)
            case(3,4)
               fdinfo_recv(10,n,i)=(jx-ix+1)*(jz-iz+1)
            case(5,6)
               fdinfo_recv(10,n,i)=(jx-ix+1)*(jy-iy+1)
            end select
         end do
      end do

!- deallocate ---------------------------------------------------------
      mem=mem-bsintg*size(neighbor_info) ; deallocate( neighbor_info )
!----------------------------------------------------------------------
!- deallocate ---------------------------------------------------------------
      mem=mem-bsintg*size(map_grid_2_pinfo) ; deallocate( map_grid_2_pinfo )
!----------------------------------------------------------------------------

      call watch(ctime1,etime1)
      if (DISP_SWITCH) then
         write(*,*) "TIME(INIT_FDPARALLEL_MOL)=",ctime1-ctime0,etime1-etime0
         write(*,*) "MEM(MB)=",mem,memax*B2MB
      end if

      return

      END SUBROUTINE init_fdparallel_mol
