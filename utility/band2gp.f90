
      PROGRAM bans2gp
      implicit none
      integer :: mb,mbz,mb1,mb2,mb3,mb4,ib
      integer :: k,n,idummy,iloc(2,2),iv,iv2,ic,ik,m,md,i
      real(8),parameter :: HT=27.2116d0
      real(8),allocatable :: e(:,:),kbb(:)
      real(8) :: eg,ev,ec,tmp1,tmp2,kbb_tmp
      real(8) :: dk,minv(3),mass(3),cfd(0:6)
      real(8) :: aa(3,3),bb(3,3),ax,va,ak(3),aak,b(3)
      character(32) :: form

!------kuchida 2010-10-29
      integer :: nbk
      real(8),allocatable :: kpd(:) 


      aa(:,:) = 0.d0
      bb(:,:) = 0.d0
      ev      = 0.d0

!------kuchida 2010-10-29
      read(*,*) nbk
      allocate( kpd(nbk) )
      Do i=1,nbk
       read(*,*) kpd(i)
      End Do



      read(*,*) mb,mbz,mb1,mb2,mb3,mb4
      write(*,*) "MB ,MBZ=",mb ,mbz
      write(*,*) "MB1,MB2=",mb1,mb2,"(",mb2-mb1+1,")"
      write(*,*) "MBv,MBc=",mb3,mb4

      if ( mb1<1.or.mb<mb1.or.mb2<1.or.mb<mb2.or.mb2<mb1 ) then
         mb1=1
         mb2=mb
         write(*,*) "MB1,MB2=",mb1,mb2,"(",mb2-mb1+1,")"
      end if

      allocate( e(mb,-mbz:mbz) ) ; e=0.d0
      allocate( kbb(-mbz:mbz)  ) ; kbb=0.d0

      do k=0,mbz-1
         do n=1,mb
            read(*,*) kbb_tmp,idummy,e(n,k)
         end do
         kbb(k)=kbb_tmp
         write(*,*) k,e(1,k),e(mb,k)
      end do

      do k=-1,-mbz+1
         do n=1,mb
            e(n,k)=e(n,-k)
         end do
         kbb(k)=kbb(-k)
      end do

      if ( 0<mb3 .and. mb3<mb ) then
         ev=maxval( e(1:mb3,0:mbz-1) )
         ec=minval( e(mb4:mb,0:mbz-1) )
      else
         ev=0.d0
         ec=0.d0
      end if
      eg=ec-ev

!c      iloc(:,1)=maxloc( e(1:iv,1:mbz) )
!c      iloc(:,2)=minloc( e(ic:mb,1:mbz) ) ; iloc(1,2)=iloc(1,2)+ic-1
!c      aa(:,:)=ax*aa(:,:)
!c      va=aa(1,1)*aa(2,2)*aa(3,3)+aa(1,2)*aa(2,3)*aa(3,1)
!c     &  +aa(1,3)*aa(2,1)*aa(3,2)-aa(1,3)*aa(2,2)*aa(3,1)
!c     &  -aa(1,2)*aa(2,1)*aa(3,3)-aa(1,1)*aa(2,3)*aa(3,2)
!c      bb(1,1) = aa(2,2)*aa(3,3) - aa(3,2)*aa(2,3) ! x-component
!c      bb(2,1) = aa(3,2)*aa(1,3) - aa(1,2)*aa(3,3) ! y
!c      bb(3,1) = aa(1,2)*aa(2,3) - aa(2,2)*aa(1,3) ! z
!c      bb(1,2) = aa(2,3)*aa(3,1) - aa(3,3)*aa(2,1)
!c      bb(2,2) = aa(3,3)*aa(1,1) - aa(1,3)*aa(3,1)
!c      bb(3,2) = aa(1,3)*aa(2,1) - aa(2,3)*aa(1,1)
!c      bb(1,3) = aa(2,1)*aa(3,2) - aa(3,1)*aa(2,2)
!c      bb(2,3) = aa(3,1)*aa(1,2) - aa(1,1)*aa(3,2)
!c      bb(3,3) = aa(1,1)*aa(2,2) - aa(2,1)*aa(1,2)
!c      bb(:,:) = 8.d0*atan(1.d0)/va*bb(:,:)
!c      ak(1:3) = b(1)*bb(1:3,1) + b(2)*bb(1:3,2) + b(3)*bb(1:3,3)
!c      aak=sqrt(sum(ak(1:3)*ak(1:3)))
!c      dk=kbb(2)-kbb(1)
!c      kbb(:)=kbb(:)/dk
!c      dk=aak/real(mbz,8)
!c      kbb(:)=kbb(:)*dk
!c      write(*,*) "dk,aak=",dk,aak
!c      do i=1,3
!c         write(*,*) aa(:,i)/ax
!c      end do
!c      do i=1,3
!c         write(*,*) ax/(atan(1.d0)*8.d0)*bb(:,i)
!c      end do
!c      do n=1,mb
!c         write(*,'(1x,i3,2x,3f15.8)') n,e(n,0:2)
!c      end do

      write(*,*) "Eg(eV)=",eg*HT
      write(*,*) "Ev(eV)=",ev*HT,mb3
      write(*,*) "Ec(eV)=",ec*HT,mb4

!
! --- effective mass ---
!
!c      cfd(:)=0.d0
!cc      md=1 ; cfd(0)=-2.d0  ; cfd(1)=1.d0
!c      md=2 ; cfd(0)=-2.5d0 ; cfd(1)=4.d0/3.d0 ; cfd(2)=-1.d0/12.d0
!c      do k=-3,mbz-md
!c         minv(1)=cfd(0)*e(iv2,k)
!c         minv(2)=cfd(0)*e(iv ,k)
!c         minv(3)=cfd(0)*e(ic ,k)
!c         do m=1,md
!c            minv(1) = minv(1) + cfd(m)*( e(iv2,k+m) + e(iv2,k-m) )
!c            minv(2) = minv(2) + cfd(m)*( e(iv ,k+m) + e(iv ,k-m) )
!c            minv(3) = minv(3) + cfd(m)*( e(ic ,k+m) + e(ic ,k-m) )
!c         end do
!c         minv(1:3)=minv(1:3)/dk**2
!c         mass(1:3)=1.d0/minv(1:3)
!cc         write(*,'(1x,i4,2x,4f16.8)') k,minv(1:2),mass(1:2)
!c         write(*,'(1x,i4,2x,5f10.5)') k,e(iv,k)*HT,e(ic,k)*HT,mass(1:3)
!c      end do

!
! --- for gnuplot ---
!
      if ( mb2-mb1+1>50 ) then
         write(*,*) "WARNING : check write format"
      end if

!---------------
         rewind 51
         do k=0,mbz-1
            write(51,'(1x,f12.7,50f15.7)') kbb(k)/kbb(mbz-1),e(mb1:mb2,k)*HT
         end do

!---------------

      rewind 53
      Do i=1,nbk-1
        Write(53,*) kpd(i)/kbb(mbz-1), -10
        Write(53,*) kpd(i)/kbb(mbz-1),   0
        Write(53,*) kpd(i)/kbb(mbz-1),  30
        Write(53,*) 
      End Do
      Write(53,*) 

      rewind 52
         do ib=mb1,mb2
            do k=0,mbz-1
               write(52,'(1x,f12.7,f15.7)') kbb(k)/kbb(mbz-1),e(ib,k)*HT
            end do
            write(52,*) ' '  
         end do

      deallocate( kpd )
      stop
      END PROGRAM bans2gp
