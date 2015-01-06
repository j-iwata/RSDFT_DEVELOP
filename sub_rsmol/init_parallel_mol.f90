!--------1---------2---------3---------4---------5---------6---------7--
!
! Grid generation and numbering
!
      SUBROUTINE init_parallel_mol
      use global_variables
      implicit none

      integer :: a1,a2,a3,b1,b2,b3
      integer :: i,i1,i2,i3,j,k,m,n
      real(8) :: H,r2,Rc2,z
      real(8) :: ctime0,ctime1,etime0,etime1
      real(8) :: mem,memax
      real(8),parameter :: eps=1.d-10
      integer :: Mx,My,Mz,ierr
      integer,allocatable :: Mxp(:),Myp(:),Mzp(:)
      integer,allocatable :: itmp(:),jtmp(:),LLL1(:,:,:),LLL0(:,:,:)
      integer,allocatable :: ip(:,:,:)
      integer :: icheck_1,icheck_2,icheck_3,icheck_4

      if ( isymmetry==1 ) then
         call init_parallel_mol_sym
         return
      end if

      call watch(ctime0,etime0)

      mem   = 0.d0
      memax = 0.d0

      Mx = ML1
      My = ML2
      Mz = ML3

      H = H1

      Rc2 = Rsize**2

      icheck_1=(2*Mx+1)*(2*My+1)*(2*Mz+1)

!
! --- processor grid ---
!
      allocate( LLp(3,0:np_grid-1), LLLp(np1,np2,np3) )
      mem=mem+bsintg*(3*np_grid+np_grid) ; memax=max(mem,memax)

      n=-1
      do i3=1,np3
      do i2=1,np2
      do i1=1,np1
         n=n+1
         LLp(1,n)=i1
         LLp(2,n)=i2
         LLp(3,n)=i3
         LLLp(i1,i2,i3)=n
      end do
      end do
      end do

!
! --- make real-space grid ---
!

!- allocate ----------------------------------------------
      allocate( Mxp(np1),Myp(np2),Mzp(np3) )
      allocate( ip(3,2,0:np_grid-1) ) ; ip=0
      mem=mem+bsintg*(np1+np2+np3)+bsintg*(6*np_grid)
      memax=max(mem,memax)
!---------------------------------------------------------

      Mxp(1:np1)=(2*Mx+1)/np1
      Myp(1:np2)=(2*My+1)/np2
      Mzp(1:np3)=(2*Mz+1)/np3

      m=(2*Mx+1)-sum(Mxp)
      if(m/=0)then
         do i=1,abs(m)
            if( mod(i,2)==1 )then
               Mxp(np1)=Mxp(np1)+sign(1,m)
            else
               Mxp(1)=Mxp(1)+sign(1,m)
            end if
         end do
      end if

      m=(2*My+1)-sum(Myp)
      if(m/=0)then
         do i=1,abs(m)
            if( mod(i,2)==1 )then
               Myp(np2)=Myp(np2)+sign(1,m)
            else
               Myp(1)=Myp(1)+sign(1,m)
            end if
         end do
      end if

      m=(2*Mz+1)-sum(Mzp)
      if(m/=0)then
         do i=1,abs(m)
            if( mod(i,2)==1 )then
               Mzp(np3)=Mzp(np3)+sign(1,m)
            else
               Mzp(1)=Mzp(1)+sign(1,m)
            end if
         end do
      end if

      if (DISP_SWITCH) then
         write(*,'(1x,"(before)")')
         write(*,'(1x,"Mxp=",10i6)') Mxp(1:np1)
         write(*,'(1x,"Myp=",10i6)') Myp(1:np2)
         write(*,'(1x,"Mzp=",10i6)') Mzp(1:np3)
         write(*,'(1x,"sum(Mxp),sum(Myp),sum(Mzp)=",3i6)') sum(Mxp),sum(Myp),sum(Mzp)
      end if

! --- adjust the block size of each node (for equal division) ----

      select case(iswitch_eqdiv)
      case default

         call mesh_div_1(icheck_2)
         n=-1
         do i3=1,np3
         do i2=1,np2
         do i1=1,np1
            n=n+1
            pinfo_grid_2d(1,n)=-Mx+sum( Mxp(1:i1) )-Mxp(i1)
            pinfo_grid_2d(2,n)= Mxp(i1)
            pinfo_grid_2d(3,n)=-My+sum( Myp(1:i2) )-Myp(i2)
            pinfo_grid_2d(4,n)= Myp(i2)
            pinfo_grid_2d(5,n)=-Mz+sum( Mzp(1:i3) )-Mzp(i3)
            pinfo_grid_2d(6,n)= Mzp(i3)
         end do
         end do
         end do

         if (DISP_SWITCH) then
            write(*,'(1x,"(after)")')
            write(*,'(1x,"Mxp=",10i6)') Mxp(1:np1)
            write(*,'(1x,"Myp=",10i6)') Myp(1:np2)
            write(*,'(1x,"Mzp=",10i6)') Mzp(1:np3)
            write(*,'(1x,"sum(Mxp),sum(Myp),sum(Mzp)=",3i6)') sum(Mxp),sum(Myp),sum(Mzp)
         end if
         m=sum(Mxp)*sum(Myp)*sum(Mzp)
         if ( m/=icheck_1 ) goto 900

      case(1)

         call mesh_div_2
         n=-1
         do i3=1,np3
         do i2=1,np2
         do i1=1,np1
            n=n+1
            pinfo_grid_2d(1,n)=ip(1,1,n)
            pinfo_grid_2d(2,n)=ip(1,2,n)-ip(1,1,n)+1
            pinfo_grid_2d(3,n)=ip(2,1,n)
            pinfo_grid_2d(4,n)=ip(2,2,n)-ip(2,1,n)+1
            pinfo_grid_2d(5,n)=ip(3,1,n)
            pinfo_grid_2d(6,n)=ip(3,2,n)-ip(3,1,n)+1
         end do
         end do
         end do

      end select

!- deallocate ---------------------------------------------------
      deallocate( ip ) ; mem=mem-bsintg*6*np_grid
      deallocate( Mxp,Myp,Mzp ) ; mem=mem-bsintg*(np1+np2+np3)
      deallocate( LLLp,LLp ) ; mem=mem-bsintg*(np_grid+3*np_grid)
!-----------------------------------------------------------------
!- allocate -----------------------------------------------------------
      allocate( itmp(0:np_grid-1) ) ; itmp=0
      allocate( jtmp(0:np_grid-1) ) ; jtmp=0
      mem=mem+bsintg*np_grid*2 ; memax=max(mem,memax)
!      Mx=ML1+Md ; My=ML2+Md ; Mz=ML3+Md
!      allocate( LLL1(-Mx:Mx,-My:My,-Mz:Mz) ) ; LLL1=0
!      mem=mem+bsintg*(2*Mx+1)*(2*My+1)*(2*Mz+1) ; memax=max(mem,memax)
!----------------------------------------------------------------------

      i=0
      do n=0,np_grid-1
         a1=pinfo_grid_2d(1,n) ; b1=a1+pinfo_grid_2d(2,n)-1
         a2=pinfo_grid_2d(3,n) ; b2=a2+pinfo_grid_2d(4,n)-1
         a3=pinfo_grid_2d(5,n) ; b3=a3+pinfo_grid_2d(6,n)-1
!-allocate-
         allocate( LLL0(a1-Md:b1+Md,a2-Md:b2+Md,a3-Md:b3+Md) )
         mem=mem+bsintg*size(LLL0) ; memax=max(mem,memax)
         LLL0(:,:,:)=0
         select case(SYStype)
         case(1)
            m=0
            do i3=a3,b3
            do i2=a2,b2
            do i1=a1,b1
               r2=H*H*(i1*i1+i2*i2+i3*i3)
               if ( r2<Rc2+eps ) then
                  i=i+1
                  m=m+1
                  do j=-Md,Md
                     LLL0(i1+j,i2,i3)=1
                     LLL0(i1,i2+j,i3)=1
                     LLL0(i1,i2,i3+j)=1
                  end do
               end if
            end do
            end do
            end do
            do i3=a3-Md,b3+Md
            do i2=a2-Md,b2+Md
            do i1=a1-Md,b1+Md
               r2=H*H*(i1*i1+i2*i2+i3*i3)
               if ( r2<Rc2+eps ) then
                  LLL0(i1,i2,i3)=-1
               end if
            end do
            end do
            end do
         case(2)
            m=0
            do i3=a3,b3
            do i2=a2,b2
            do i1=a1,b1
               r2=H*H*(i1*i1+i2*i2)
               z=abs(i3*H)
               if ( r2<Rc2+eps .and. z<Zsize+eps ) then
                  i=i+1
                  m=m+1
                  do j=-Md,Md
                     LLL0(i1+j,i2,i3)=1
                     LLL0(i1,i2+j,i3)=1
                     LLL0(i1,i2,i3+j)=1
                  end do
               end if
            end do
            end do
            end do
            do i3=a3-Md,b3+Md
            do i2=a2-Md,b2+Md
            do i1=a1-Md,b1+Md
               r2=H*H*(i1*i1+i2*i2)
               z=abs(i3*H)
               if ( r2<Rc2+eps .and. z<Zsize+eps ) then
                  LLL0(i1,i2,i3)=-1
               end if
            end do
            end do
            end do
         end select
         pinfo_grid_2d(7,n)=i-m
         pinfo_grid_2d(8,n)=m
         if ( DISP_SWITCH ) then
            write(*,*) n,count(LLL0/=0),count(LLL0<0),count(LLL0>0)
         end if
         jtmp(n)=count(LLL0(:,:,:)>0)
         mem=mem-bsintg*size(LLL0) ; deallocate( LLL0 )
      end do

      do n=0,np_grid-1
         itmp(n)=sum(jtmp(0:n))-jtmp(n)
      end do

!- allocate ------------------
      call gv_alloc("idisp2")
!-----------------------------

      do n=0,nprocs-1
         idisp2(n)=itmp(id_class(n,0))
         ircnt2(n)=jtmp(id_class(n,0))
      end do
      if ( DISP_SWITCH ) then
         write(*,'(1x,a6,2x,2a8)') "nprocs","idisp2","ircnt2"
         do n=0,nprocs-1
            write(*,'(1x,i6,2x,2i8)') n,idisp2(n),ircnt2(n)
         end do
      end if

!- deallocate ---------------------------------------------------------
      mem=mem-bsintg*(size(jtmp)+size(itmp)) ; deallocate( jtmp,itmp )
!----------------------------------------------------------------------

      call watch(ctime1,etime1)
      if (DISP_SWITCH) then
         write(*,*) "TIME(INIT_PARALLEL_MOL)=",ctime1-ctime0,etime1-etime0
         write(*,*) "MEM(MB)=",mem,memax*B2MB
      end if

      return

 900  call stop_program1("init_parallel_mol",1)


      CONTAINS


      SUBROUTINE mesh_div_1(icheck_2)
      implicit none
      integer,intent(INOUT) :: icheck_2
      integer,parameter :: max_loop=600
      integer :: iloop,ix,iy,iz,i1,i2,i3,m,n,Mx0,My0,Mz0
      integer,allocatable :: itmp(:)
      real(8),parameter :: eps=1.d-10
      real(8) :: r2,Rc2,z

      if ( np_grid<=1 ) return

      allocate( itmp(0:np_grid-1) ) ; itmp=0
      mem=mem+bsintg*np_grid ; memax=max(mem,memax)

      Rc2=Rsize**2

      do iloop=1,max_loop

         itmp(:)=0

         n=-1

         Mz0=-Mz
         do i3=1,np3
         My0=-My
         do i2=1,np2
         Mx0=-Mx
         do i1=1,np1
            n=n+1
            select case(SYStype)
            case(1)
               do iz=Mz0,Mz0+Mzp(i3)-1
               do iy=My0,My0+Myp(i2)-1
               do ix=Mx0,Mx0+Mxp(i1)-1
                  r2=H*H*(ix*ix+iy*iy+iz*iz)
                  if( r2<=Rc2+eps )then
                     itmp(n)=itmp(n)+1
                  end if
               end do
               end do
               end do
            case(2)
               do iz=Mz0,Mz0+Mzp(i3)-1
                  z=iz*H
                  if ( abs(z)>Zsize+eps ) cycle
                  do iy=My0,My0+Myp(i2)-1
                  do ix=Mx0,Mx0+Mxp(i1)-1
                     r2=H*H*(ix*ix+iy*iy)
                     if ( r2<=Rc2+eps ) then
                        itmp(n)=itmp(n)+1
                     end if
                  end do
                  end do
               end do
            end select
            Mx0=Mx0+Mxp(i1)
         end do
            My0=My0+Myp(i2)
         end do
            Mz0=Mz0+Mzp(i3)
         end do

         i=0 ; j=0
         do n=0,np_grid-1
            if( itmp(i)<itmp(n) ) i=n
            if( itmp(j)>itmp(n) ) j=n
         end do

         if( mod(iloop,3)==1 )then
            Mzp(LLp(3,i))=Mzp(LLp(3,i))-1
            Mzp(LLp(3,j))=Mzp(LLp(3,j))+1
         else if( mod(iloop,3)==2 )then
            Myp(LLp(2,i))=Myp(LLp(2,i))-1
            Myp(LLp(2,j))=Myp(LLp(2,j))+1
         else if( mod(iloop,3)==0 )then
            Mxp(LLp(1,i))=Mxp(LLp(1,i))-1
            Mxp(LLp(1,j))=Mxp(LLp(1,j))+1
         end if

      end do ! iloop

      if(DISP_SWITCH)then
         write(*,*) "iloop=",iloop
         write(*,'(5("(",i4,")",i6,2x))') (n,itmp(n),n=0,np_grid-1)
         write(*,*) "tot. ave.",sum(itmp),sum(itmp)/real(np_grid)
      end if

      icheck_2=sum(itmp)

      deallocate( itmp ) ; mem=mem-bsintg*np_grid

      return
      END SUBROUTINE mesh_div_1


      SUBROUTINE mesh_div_2
      implicit none
      integer,parameter :: max_loop=100
      integer :: i1,i2,i3,i,j,iloop,iloc(1),m1,m2,m3
      integer :: ix,jx,iy,jy,iz,jz,i_dif
      integer :: min_dif,max_val,icount
      integer,allocatable :: LLLL(:,:,:),itmp(:),jtmp(:),mtmp(:,:)
      real(8),parameter :: eps=1.d-10
      real(8) :: Rc2,r2,z

      n   = max( np1, np2, np3 )
      Rc2 = Rsize**2

      allocate( itmp(n) )   ; itmp=0
      allocate( jtmp(n) )   ; jtmp=0
      allocate( mtmp(n,3) ) ; mtmp=0
      mem=mem+bsintg*n*2+bsintg*n*3 ; memax=max(mem,memax)

      allocate( LLLL(-Mx:Mx,-My:My,-Mz:Mz) ) ; LLLL=0
      mem=mem+bsintg*(2*Mx+1)*(2*My+1)*(2*Mz+1) ; memax=max(mem,memax)

      select case(SYStype)
      case(1)
         do i3=-Mz,Mz
         do i2=-My,My
         do i1=-Mx,Mx
            r2=H*H*(i1*i1+i2*i2+i3*i3)
            if ( r2<=Rc2+eps ) LLLL(i1,i2,i3)=1
         end do
         end do
         end do
      case(2)
         do i3=-Mz,Mz
            z=H*i3
            if ( abs(z)>Zsize+eps ) cycle
            do i2=-My,My
            do i1=-Mx,Mx
               r2=H*H*(i1*i1+i2*i2)
               if ( r2<=Rc2+eps ) LLLL(i1,i2,i3)=1
            end do
            end do
         end do
      end select

      if (DISP_SWITCH) then
         write(*,*) "count(LLLL==1)=",count(LLLL==1)
      end if

      min_dif = sum(LLLL)
      max_val = min_dif
      icount  = 0
      jtmp(:) = 0

      loop_x : do iloop=1,max_loop
         itmp(:)=0
         do m1=1,np1
            ix=-Mx+sum(Mxp(1:m1))-Mxp(m1)
            jx=ix+Mxp(m1)-1
            do i1=ix,jx
               itmp(m1)=itmp(m1)+sum(LLLL(i1,:,:))
            end do
         end do
         iloc(:)=maxloc( itmp(1:np1) ) ; i=iloc(1)
         iloc(:)=minloc( itmp(1:np1) ) ; j=iloc(1)
         i_dif=itmp(i)-itmp(j)
         if ( i_dif<min_dif ) then
            min_dif=i_dif
            max_val=maxval(itmp(1:np1))
            mtmp(1:np1,1)=Mxp(1:np1)
            jtmp(1:np1)=itmp(1:np1)
            icount=0
         else if ( i_dif==min_dif ) then
            if ( maxval(itmp(1:np1))==max_val ) then
               icount=icount+1
            else if ( maxval(itmp(1:np1))<max_val ) then
               max_val=maxval(itmp(1:np1))
               mtmp(1:np1,1)=Mxp(1:np1)
               jtmp(1:np1)=itmp(1:np1)
               icount=0
            end if
         end if
         if ( icount>9 ) exit
         Mxp(i)=Mxp(i)-1
         Mxp(j)=Mxp(j)+1
      end do loop_x

      i=maxval(jtmp(1:np1))
      j=minval(jtmp(1:np1))
      if (DISP_SWITCH) then
         write(*,'(1x,"iloop =",i4,2x,"max min dif =",3i6)') iloop,i,j,i-j
         write(*,'(1x,"rank_x",16i8)') (i,i=0,np1-1)
         write(*,'(1x,"      ",16i8)') jtmp(1:np1)
         write(*,'(1x,"      ",16i8)') mtmp(1:np1,1)
         write(*,*) "total :",sum(jtmp(1:np1)),sum(mtmp(1:np1,1))
      end if
      do i3=1,np3
      do i2=1,np2
      do i1=1,np1
         n=LLLp(i1,i2,i3)
         ip(1,1,n)=-Mx+sum(mtmp(1:i1,1))-mtmp(i1,1)
         ip(1,2,n)=ip(1,1,n)+mtmp(i1,1)-1
      end do
      end do
      end do

      do m1=1,np1

         min_dif = sum(LLLL)
         max_val = min_dif
         icount  = 0
         jtmp(:) = 0

         loop_y : do iloop=1,max_loop
            itmp(:)=0
            do m2=1,np2
               n=LLLp(m1,m2,1)
               ix=ip(1,1,n)
               jx=ip(1,2,n)
               iy=-My+sum(Myp(1:m2))-Myp(m2)
               jy=iy+Myp(m2)-1
               do i2=iy,jy
                  itmp(m2)=itmp(m2)+sum(LLLL(ix:jx,i2,:))
               end do
            end do
            iloc(:)=maxloc(itmp(1:np2)) ; i=iloc(1)
            iloc(:)=minloc(itmp(1:np2)) ; j=iloc(1)
            i_dif=itmp(i)-itmp(j)
            if ( i_dif<min_dif ) then
               min_dif=i_dif
               max_val=maxval(itmp(1:np2))
               mtmp(1:np2,2)=Myp(1:np2)
               jtmp(1:np2)=itmp(1:np2)
               icount=0
            else if ( i_dif==min_dif ) then
               if ( maxval(itmp(1:np2))==max_val ) then
                  icount=icount+1
               else if ( maxval(itmp(1:np2))<max_val ) then
                  max_val=maxval(itmp(1:np2))
                  mtmp(1:np2,2)=Myp(1:np2)
                  jtmp(1:np2)=itmp(1:np2)
                  icount=0
               end if
            end if
            if ( icount>9 ) exit
            Myp(i)=Myp(i)-1
            Myp(j)=Myp(j)+1
         end do loop_y

         i=maxval(jtmp(1:np2)) ; j=minval(jtmp(1:np2))
         if (DISP_SWITCH) then
            write(*,'(1x,"iloop =",i4,2x,"max min dif =",3i6)') iloop,i,j,i-j
            write(*,'(1x,"rank_y",16i8)') (i,i=0,np2-1)
            write(*,'(1x,"      ",16i8)') jtmp(1:np2)
            write(*,'(1x,"      ",16i8)') mtmp(1:np2,2)
            write(*,*) "total :",sum(jtmp(1:np2)),sum(mtmp(1:np2,2))
         end if

         do i3=1,np3
         do i2=1,np2
            n=LLLp(m1,i2,i3)
            ip(2,1,n)=-My+sum(mtmp(1:i2,2))-mtmp(i2,2)
            ip(2,2,n)=ip(2,1,n)+mtmp(i2,2)-1
         end do
         end do

      end do ! m1


      do m1=1,np1
      do m2=1,np2

         min_dif = sum(LLLL)
         max_val = min_dif
         icount  = 0
         jtmp(:) = 0

         loop_z : do iloop=1,max_loop
            itmp(:)=0
            do m3=1,np3
               n=LLLp(m1,m2,m3)
               ix=ip(1,1,n)
               jx=ip(1,2,n)
               iy=ip(2,1,n)
               jy=ip(2,2,n)
               iz=-Mz+sum(Mzp(1:m3))-Mzp(m3)
               jz=iz+Mzp(m3)-1
               do i3=iz,jz
                  itmp(m3)=itmp(m3)+sum(LLLL(ix:jx,iy:jy,i3))
               end do
            end do
            iloc(:)=maxloc(itmp(1:np3)) ; i=iloc(1)
            iloc(:)=minloc(itmp(1:np3)) ; j=iloc(1)
            i_dif=itmp(i)-itmp(j)
            if ( i_dif<min_dif ) then
               min_dif=i_dif
               max_val=maxval(itmp(1:np3))
               mtmp(1:np3,3)=Mzp(1:np3)
               jtmp(1:np3)=itmp(1:np3)
               icount=0
            else if ( i_dif==min_dif ) then
               if ( maxval(itmp(1:np3))==max_val ) then
                  icount=icount+1
               else if ( maxval(itmp(1:np3))<max_val ) then
                  max_val=maxval(itmp(1:np3))
                  mtmp(1:np3,3)=Mzp(1:np3)
                  jtmp(1:np3)=itmp(1:np3)
                  icount=0
               end if
            end if
            if ( icount>9 ) exit
            Mzp(i)=Mzp(i)-1
            Mzp(j)=Mzp(j)+1
         end do loop_z

         i=maxval(jtmp(1:np3)) ; j=minval(jtmp(1:np3))
         if (DISP_SWITCH) then
            write(*,'(1x,"iloop =",i4,2x,"max min dif =",3i6)') iloop,i,j,i-j
            write(*,'(1x,"rank_z",16i8)') (i,i=0,np3-1)
            write(*,'(1x,"      ",16i8)') jtmp(1:np3)
            write(*,'(1x,"      ",16i8)') mtmp(1:np3,3)
            write(*,*) "total :",sum(jtmp(1:np3)),sum(mtmp(1:np3,3))
         end if

         do i3=1,np3
            n=LLLp(m1,m2,i3)
            ip(3,1,n)=-Mz+sum(mtmp(1:i3,3))-mtmp(i3,3)
            ip(3,2,n)=ip(3,1,n)+mtmp(i3,3)-1
         end do

      end do ! m2
      end do ! m1

      mem=mem-bsintg*size(LLLL) ; deallocate( LLLL )
      mem=mem-bsintg*size(mtmp) ; deallocate( mtmp )
      mem=mem-bsintg*size(jtmp) ; deallocate( jtmp )
      mem=mem-bsintg*size(itmp) ; deallocate( itmp )

      return
      END SUBROUTINE mesh_div_2


      END SUBROUTINE init_parallel_mol
