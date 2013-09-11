!--------1---------2---------3---------4---------5---------6---------7--
!
! Grid generation and numbering
!
      SUBROUTINE init_mesh
      use global_variables
      implicit none

      integer,allocatable :: Lxyz(:,:,:)
      integer :: ix,iy,iz,i !,j,k,n,m,a,il,ik,ierr,n1,n2,iloop
      real(8) :: H,r,x,y,z,r0,r2,Rc2
      real(8) :: ctime0,ctime1,etime0,etime1
      real(8) :: mem,memax
      real(8),parameter :: eps=1.d-10
      integer :: Mx,My,Mz
!      integer :: ix0,ix1,iy0,iy1,iz0,iz1,jx,jy,jz
!      integer :: ntmp,i1,i2,i3,ilp,ikp
!      integer :: Mx00,My00,Mz00,Mx11,My11,Mz11
!      integer :: mxi,myi,mzi
!      integer,allocatable :: Mxp(:),Myp(:),Mzp(:),itmp(:),jtmp(:)
!      integer,allocatable :: ip(:,:,:)
!      integer,allocatable :: itmp2(:,:)
!      integer :: sp,fp,jrank,irank,ns
!      integer :: ireqs(18),ireqr(18),nreq,ireq(100)
!      integer :: itags,itagr,istatus(mpi_status_size)
      integer :: icheck_1 !,icheck_2,icheck_3,icheck_4

      call watch(ctime0,etime0)
      if (DISP_SWITCH) then
         write(*,'(a60," init_mesh")') repeat("-",60)
      end if

      mem   = 0.d0
      memax = 0.d0
      H     = H1

!
! --- Mesh, Volume element ---
!

      H2=H1
      H3=H1

      dV=Jacobian*H1*H2*H3

      zdV=dV

!
! --- cut off ---
!

      qf=Pi/H1
      Ecut=qf*qf

      if ( DISP_SWITCH ) then
         write(*,*) "qf, Ecut=",qf,Ecut
      end if

!
! --- Box shape ---
!

      select case( SYStype )
      case(1)
         if (DISP_SWITCH) write(*,*) "Spherical Box"
         Mx=int(Rsize/H) ; if ( Mx*H<Rsize ) Mx=Mx+1
         My=Mx ; Mz=Mx
         if(DISP_SWITCH)write(*,*) "Rsize=",Rsize,Mx*H
      case(2)
         if (DISP_SWITCH) write(*,*) "Cylindrical Box"
         Mz=int(Zsize/H) ; if ( Mz*H<Zsize ) Mz=Mz+1
         Mx=int(Rsize/H) ; if ( Mx*H<Rsize ) Mx=Mx+1
         My=Mx
      case default
         goto 900
      end select

      icheck_1 = (2*Mx+1)*(2*My+1)*(2*Mz+1)

      if (DISP_SWITCH) then
         write(*,*) "Mx, My, Mz  =",Mx ,My ,Mz
         write(*,'(1x,"2Mx+1,2My+1,2Mz+1 =",3i6,i12)') 2*Mx+1,2*My+1,2*Mz+1
         write(*,'(1x,"(2Mx+1)*(2My+1)*(2Mz+1) =",i12)') icheck_1
      end if

!
! --- Mx+Md ---
!

      ML1=Mx ; ML2=My ; ML3=Mz

      Mx=Mx+Md ; My=My+Md ; Mz=Mz+Md

      if (DISP_SWITCH) then
         write(*,*) "( replace Mx=>Mx+Md, My=>My+Md, Mz=>Mz+Md )"
         write(*,*) "Mx My Mz      =", Mx,My,Mz
         write(*,*) "ML1,ML2,ML3   =",ML1,ML2,ML3
         write(*,*) "# of grid points =",(2*Mx+1)*(2*My+1)*(2*Mz+1)
      end if

!
! --- symmetry ---
!

      ML_irreducible=0
      MK_irreducible=0
      if ( isymmetry==1 ) then
         call init_mesh_sym
         return
      end if

!
! --- Division of area to inner(=2) and outer(=1) one ---
!

!- allocate ------------------------------------------------------------
      allocate( Lxyz(-Mx:Mx,-My:My,-Mz:Mz) ) ; Lxyz=0
      mem=mem+bsintg*(2*Mx+1)*(2*My+1)*(2*Mz+1) ; memax=max(mem,memax)
!-----------------------------------------------------------------------

      ML=0
      MK=0
      Rc2=Rsize**2

      do iz=-Mz+Md,Mz-Md
      do iy=-My+Md,My-Md
      do ix=-Mx+Md,Mx-Md
         select case(SYStype)
         case(1)
            r2=H*H*(ix*ix+iy*iy+iz*iz)
            if ( r2<Rc2+eps ) then
               do i=-Md,Md
                  Lxyz(ix+i,iy,iz)=-1
                  Lxyz(ix,iy+i,iz)=-1
                  Lxyz(ix,iy,iz+i)=-1
               end do
            end if
         case(2)
            r2=H*H*(ix*ix+iy*iy)
            z=abs(iz*H)
            if ( r2<Rc2+eps .and. z<Zsize+eps ) then
               do i=-Md,Md
                  Lxyz(ix+i,iy,iz)=-1
                  Lxyz(ix,iy+i,iz)=-1
                  Lxyz(ix,iy,iz+i)=-1
               end do
            end if
         end select
      end do
      end do
      end do
      do iz=-Mz+Md,Mz-Md
      do iy=-My+Md,My-Md
      do ix=-Mx+Md,Mx-Md
         select case(SYStype)
         case(1)
            r2=H*H*(ix*ix+iy*iy+iz*iz)
            if ( r2<Rc2+eps ) then
               Lxyz(ix,iy,iz)=-2
            end if
         case(2)
            r2=H*H*(ix*ix+iy*iy)
            z=abs(iz*H)
            if ( r2<Rc2+eps .and. z<Zsize+eps ) then
               Lxyz(ix,iy,iz)=-2
            end if
         end select
      end do
      end do
      end do

      ML=count(Lxyz==-2)
      MK=count(Lxyz==-1)

      if(DISP_SWITCH) then
         write(*,*) "ML=",ML
         write(*,*) "MK=",MK
      end if

!- deallocate ---------------------------------------------------------
      deallocate( Lxyz ) ; mem=mem-bsintg*(2*Mx+1)*(2*My+1)*(2*Mz+1)
!----------------------------------------------------------------------

      call watch(ctime1,etime1)
      if (DISP_SWITCH) then
         write(*,*) "TIME(INIT_MESH)=",ctime1-ctime0,etime1-etime0
         write(*,*) "MEM(MB)=",mem,memax*B2MB
      end if

      return

 900  call stop_program1("init_mesh",1)

      END SUBROUTINE init_mesh
