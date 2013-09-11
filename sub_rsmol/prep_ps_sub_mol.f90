!--------1---------2---------3---------4---------5---------6---------7--
!
      SUBROUTINE Make_MinimalBox_Mol(Rc,mm1,mm2,mm3,mm4)
      use global_variables
      implicit none
      real(8),intent(IN)  :: Rc
      integer,intent(OUT) :: mm1,mm2,mm3,mm4
      real(8) :: a1,a2,c1,c2,c3
      real(8) :: x,y,z,r2,R0,R0R0,Rmax
      real(8) :: mem,memax,ctime0,etime0,ctime1,etime1
      integer :: m,n,m1,m2,m3,mm,mm0,isize,i1,i2,i3,mmmax,j1,j2,j3
      integer :: mt1,mt2,mt3
      integer :: n1a,n1b,n2a,n2b,n3a,n3b
      integer,allocatable :: jtmp(:,:,:),jtmp0(:,:,:)

      call bwatch(ctime0,etime0)

      if (DISP_SWITCH) then
         write(*,*) "--- SUBROUTINE 'Make_MinimalBox' ---"
         write(*,*) "--- (in 'prep_ps_sub.f')         ---" 
      end if

      mem   = 0.d0
      memax = 0.d0

      R0 = Rc + max(H1,H2,H3) + 1.d-8
      if (DISP_SWITCH) write(*,*) "R_search=",R0,R0*R0

      R0R0 = R0*R0

      mm1=0
      mm2=0
      mm3=0
      mm4=0

!
! --- 1st way ---
!

      mm1=nint(R0/H1)+1
      mm2=nint(R0/H2)+1
      mm3=nint(R0/H3)+1

      if (DISP_SWITCH) then
         write(*,*) "mm1,mm2,mm3 (before)=",mm1,mm2,mm3
      end if

!
! --- 2nd way ---
!

      mm0=(2*mm1+1)*(2*mm2+1)*(2*mm3+1)
      allocate( jtmp(-mm1:mm1,-mm2:mm2,-mm3:mm3) ) ; jtmp=0
      mem=mem+bsintg*mm0 ; memax=max(mem,memax)

      mmmax=1000
      Rmax=0.d0
      m1=0
      m2=0
      m3=0
      m=1
      n=1

      jtmp(0,0,0)=1

      do mm=1,mmmax

         do i3=-m3,m3
            c3=i3*H3
         do i2=-m2,m2
            c2=i2*H2
            do i1=-m1-1,-m1-1
               c1=i1*H1
               do j3=0,Ndense-1
                  z=c3+j3*H3d
               do j2=0,Ndense-1
                  y=c2+j2*H2d
               do j1=0,Ndense-1
                  x=c1+j1*H1d
                  r2=x*x+y*y+z*z
                  if ( r2<=R0R0 ) then
                     m=m+1
                     Rmax=max(Rmax,r2)
                     jtmp(i1,i2,i3)=jtmp(i1,i2,i3)+1
                  end if
               end do
               end do
               end do
            end do
            do i1=m1+1,m1+1
               c1=i1*H1
               do j3=0,Ndense-1
                  z=c3+j3*H3d
               do j2=0,Ndense-1
                  y=c2+j2*H2d
               do j1=0,Ndense-1
                  x=c1+j1*H1d
                  r2=x*x+y*y+z*z
                  if ( r2<=R0R0 ) then
                     m=m+1
                     Rmax=max(Rmax,r2)
                     jtmp(i1,i2,i3)=jtmp(i1,i2,i3)+1
                  end if
               end do
               end do
               end do
            end do
         end do
         end do
         m1=m1+1

         do i3=-m3,m3
            c3=i3*H3
         do i1=-m1,m1
            c1=i1*H1
            do i2=-m2-1,-m2-1
               c2=i2*H2
               do j3=0,Ndense-1
                  z=c3+j3*H3d
               do j2=0,Ndense-1
                  y=c2+j2*H2d
               do j1=0,Ndense-1
                  x=c1+j1*H1d
                  r2=x*x+y*y+z*z
                  if ( r2<=R0R0 ) then
                     m=m+1
                     Rmax=max(Rmax,r2)
                     jtmp(i1,i2,i3)=jtmp(i1,i2,i3)+1
                  end if
               end do
               end do
               end do
            end do
            do i2=m2+1,m2+1
               c2=i2*H2
               do j3=0,Ndense-1
                  z=c3+j3*H3d
               do j2=0,Ndense-1
                  y=c2+j2*H2d
               do j1=0,Ndense-1
                  x=c1+j1*H1d
                  r2=x*x+y*y+z*z
                  if ( r2<=R0R0 ) then
                     m=m+1
                     Rmax=max(Rmax,r2)
                     jtmp(i1,i2,i3)=jtmp(i1,i2,i3)+1
                  end if
               end do
               end do
               end do
            end do
         end do
         end do
         m2=m2+1

         do i2=-m2,m2
            c2=i2*H2
         do i1=-m1,m1
            c1=i1*H1
            do i3=-m3-1,-m3-1
               c3=i3*H3
               do j3=0,Ndense-1
                  z=c3+j3*H3d
               do j2=0,Ndense-1
                  y=c2+j2*H2d
               do j1=0,Ndense-1
                  x=c1+j1*H1d
                  r2=x*x+y*y+z*z
                  if ( r2<=R0R0 ) then
                     m=m+1
                     Rmax=max(Rmax,r2)
                     jtmp(i1,i2,i3)=jtmp(i1,i2,i3)+1
                  end if
               end do
               end do
               end do
            end do
            do i3=m3+1,m3+1
               c3=i3*H3
               do j3=0,Ndense-1
                  z=c3+j3*H3d
               do j2=0,Ndense-1
                  y=c2+j2*H2d
               do j1=0,Ndense-1
                  x=c1+j1*H1d
                  r2=x*x+y*y+z*z
                  if ( r2<=R0R0 ) then
                     m=m+1
                     Rmax=max(Rmax,r2)
                     jtmp(i1,i2,i3)=jtmp(i1,i2,i3)+1
                  end if
               end do
               end do
               end do
            end do
         end do
         end do
         m3=m3+1

         if (DISP_SWITCH) then
            write(*,'(1x,i4,2i12,3i6,f10.5)') mm,m,n,m1,m2,m3,Rmax
         end if

         if ( m==n ) exit
         n=m

         if ( m1==mm1 .or. m2==mm2 .or. m3==mm3 ) then

            allocate( jtmp0(-mm1:mm1,-mm2:mm2,-mm3:mm3) )
            mem=mem+bsintg*size(jtmp0) ; memax=max(mem,memax)

            jtmp0(:,:,:)=jtmp(:,:,:)

            mt1=mm1 ; if ( m1==mm1 ) mt1=mt1+2
            mt2=mm2 ; if ( m2==mm2 ) mt2=mt2+2
            mt3=mm3 ; if ( m3==mm3 ) mt3=mt3+2

            mem=mem-bsintg*size(jtmp) ; deallocate( jtmp )

            allocate( jtmp(-mt1:mt1,-mt2:mt2,-mt3:mt3) ) ; jtmp=0
            mem=mem+bsintg*size(jtmp) ; memax=max(mem,memax)

            jtmp(-mm1:mm1,-mm2:mm2,-mm3:mm3) = jtmp0(-mm1:mm1,-mm2:mm2,-mm3:mm3)

            mm1=mt1
            mm2=mt2
            mm3=mt3

            mem=mem-bsintg*size(jtmp0) ; deallocate( jtmp0 )

         end if

      end do

      isize = count(jtmp>0)
      Rmax  = sqrt(Rmax)

      if (DISP_SWITCH) then
         write(*,*) "Rmax=",Rmax,isize
      end if

!- allocate ---------------------------------------------
      if ( allocated(map_grid_ion) ) then
         Max_mem=Max_mem-bsintg*size(map_grid_ion)
         deallocate( map_grid_ion )
      end if
      if ( allocated(mcube_grid_ion) ) then
         Max_mem=Max_mem-bsintg*size(mcube_grid_ion)
         deallocate( mcube_grid_ion )
      end if
      allocate( map_grid_ion(3,isize)  ) ; map_grid_ion=0
      allocate( mcube_grid_ion(2,3) ) ; mcube_grid_ion=0
      Max_mem=Max_mem+bsintg*3*isize
      Max_mem=Max_mem+bsintg*2*3
!--------------------------------------------------------

      n=0
      do i3=-mm3,mm3
      do i2=-mm2,mm2
      do i1=-mm1,mm1
         if ( jtmp(i1,i2,i3)>0 ) then
            n=n+1
            map_grid_ion(1,n)=i1
            map_grid_ion(2,n)=i2
            map_grid_ion(3,n)=i3
         end if
      end do
      end do
      end do

      if( n/=isize ) then
         write(*,*) "i/=isize!!!",n,isize
         goto 900
      end if

      M_grid_ion = n

      mcube_grid_ion(1,1)=minval( map_grid_ion(1,1:M_grid_ion) )
      mcube_grid_ion(1,2)=minval( map_grid_ion(2,1:M_grid_ion) )
      mcube_grid_ion(1,3)=minval( map_grid_ion(3,1:M_grid_ion) )
      mcube_grid_ion(2,1)=maxval( map_grid_ion(1,1:M_grid_ion) )
      mcube_grid_ion(2,2)=maxval( map_grid_ion(2,1:M_grid_ion) )
      mcube_grid_ion(2,3)=maxval( map_grid_ion(3,1:M_grid_ion) )

      mem=mem-bsintg*size(jtmp) ; deallocate( jtmp )

      call bwatch(ctime1,etime1)

      if (DISP_SWITCH) then
         write(*,*) "TIME(MAKE_MINIMAL_BOX_MOL)",ctime1-ctime0,etime1-etime0
         write(*,*) "mcube_grid_ion(1,1:3)=",mcube_grid_ion(1,1:3)
         write(*,*) "mcube_grid_ion(2,1:3)=",mcube_grid_ion(2,1:3)
         write(*,*) "mem(MB)",mem,memax*B2MB
      end if

      return

 900  call stop_program

      END SUBROUTINE Make_MinimalBox_Mol
