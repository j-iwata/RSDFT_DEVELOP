PROGRAM wfiomap

  implicit none

  integer :: u1=5, u2=10
  integer :: ML1,ML2,ML3
  integer :: MBN,MBZ,MSP,i1,i2,i3,i4,i5,i6,i,irank,jrank,n
  integer :: np_old(6),np_new(6),nprocs_old,nprocs_new
  integer,allocatable :: ntmp(:,:),Igrid(:,:,:),Jgrid(:,:,:)
  integer,allocatable :: itmp(:),iomap(:,:),nmap(:)
  integer :: a1b,b1b,a2b,b2b,a3b,b3b,a1a,b1a,a2a,b2a,a3a,b3a
  integer :: a4b,b4b,a5b,b5b,a6b,b6b,a4a,b4a,a5a,b5a,a6a,b6a

  read(u1,*) ML1,ML2,ML3
  read(u1,*) MBN,MBZ,MSP
  write(*,*) "ML1,ML2,ML3",ML1,ML2,ML3
  write(*,*) "MBN,MBZ,MSP",MBN,MBZ,MSP

  np_old=1
  np_new=1
  read(u1,*) np_old(1:6)
  read(u1,*) np_new(1:6)
  write(*,*) "np_old",np_old
  write(*,*) "np_new",np_new

  nprocs_old=1
  nprocs_new=1
  do i=1,6
     nprocs_old=nprocs_old*np_old(i)
     nprocs_new=nprocs_new*np_new(i)
  end do
  write(*,*) "nprocs_old",nprocs_old
  write(*,*) "nprocs_new",nprocs_new

  i1=maxval( np_old )
  i2=maxval( np_new )
  i3=max( i1,i2 )
  allocate( ntmp(0:i3-1,6) )
!
  ntmp=0
  do i3=0,ML3-1
     i=mod(i3,np_old(3))
     ntmp(i,3)=ntmp(i,3)+1
  end do
  do i2=0,ML2-1
     i=mod(i2,np_old(2))
     ntmp(i,2)=ntmp(i,2)+1
  end do
  do i1=0,ML1-1
     i=mod(i1,np_old(1))
     ntmp(i,1)=ntmp(i,1)+1
  end do
!
  do i4=1,MBN
     i=mod(i4-1,np_old(4))
     ntmp(i,4)=ntmp(i,4)+1
  end do
  do i5=1,MBZ
     i=mod(i5-1,np_old(5))
     ntmp(i,5)=ntmp(i,5)+1
  end do
  do i6=1,MSP
     i=mod(i6-1,np_old(6))
     ntmp(i,6)=ntmp(i,6)+1
  end do
!
!  do i=1,6
!     write(*,'(1x,10i6)') ntmp(0:np_old(i)-1,i)
!  end do
!
  allocate( Igrid(2,6,0:nprocs_old-1) ) ; Igrid=0

  irank=-1
  do i6=0,np_old(6)-1
  do i5=0,np_old(5)-1
  do i4=0,np_old(4)-1
     do i3=0,np_old(3)-1
     do i2=0,np_old(2)-1
     do i1=0,np_old(1)-1
        irank=irank+1
        Igrid(1,1,irank) = sum( ntmp(0:i1,1) ) - ntmp(i1,1)
        Igrid(2,1,irank) = Igrid(1,1,irank) + ntmp(i1,1) - 1
        Igrid(1,2,irank) = sum( ntmp(0:i2,2) ) - ntmp(i2,2)
        Igrid(2,2,irank) = Igrid(1,2,irank) + ntmp(i2,2) - 1
        Igrid(1,3,irank) = sum( ntmp(0:i3,3) ) - ntmp(i3,3)
        Igrid(2,3,irank) = Igrid(1,3,irank) + ntmp(i3,3) - 1
        Igrid(1,4,irank) = sum( ntmp(0:i4,4) ) - ntmp(i4,4) + 1
        Igrid(2,4,irank) = Igrid(1,4,irank) + ntmp(i4,4) - 1
        Igrid(1,5,irank) = sum( ntmp(0:i5,5) ) - ntmp(i5,5) + 1
        Igrid(2,5,irank) = Igrid(1,5,irank) + ntmp(i5,5) - 1
        Igrid(1,6,irank) = sum( ntmp(0:i6,6) ) - ntmp(i6,6) + 1
        Igrid(2,6,irank) = Igrid(1,6,irank) + ntmp(i6,6) - 1
     end do
     end do
     end do
  end do
  end do
  end do

!  do irank=0,nprocs_old-1
!     write(*,'(1x,i6,2x,6(2i6,1x))') irank,Igrid(:,:,irank)
!  end do

! ---

  ntmp=0
  do i3=0,ML3-1
     i=mod(i3,np_new(3))
     ntmp(i,3)=ntmp(i,3)+1
  end do
  do i2=0,ML2-1
     i=mod(i2,np_new(2))
     ntmp(i,2)=ntmp(i,2)+1
  end do
  do i1=0,ML1-1
     i=mod(i1,np_new(1))
     ntmp(i,1)=ntmp(i,1)+1
  end do
  do i4=1,MBN
     i=mod(i4-1,np_new(4))
     ntmp(i,4)=ntmp(i,4)+1
  end do
  do i5=1,MBZ
     i=mod(i5-1,np_new(5))
     ntmp(i,5)=ntmp(i,5)+1
  end do
  do i6=1,MSP
     i=mod(i6-1,np_new(6))
     ntmp(i,6)=ntmp(i,6)+1
  end do

!  do i=1,6
!     write(*,'(1x,10i6)') ntmp(0:np_new(i)-1,i)
!  end do

!---

  allocate( Jgrid(2,6,0:nprocs_new-1) ) ; Jgrid=0

  irank=-1
  do i6=0,np_new(6)-1
  do i5=0,np_new(5)-1
  do i4=0,np_new(4)-1
     do i3=0,np_new(3)-1
     do i2=0,np_new(2)-1
     do i1=0,np_new(1)-1
        irank=irank+1
        Jgrid(1,1,irank) = sum( ntmp(0:i1,1) ) - ntmp(i1,1)
        Jgrid(2,1,irank) = Jgrid(1,1,irank) + ntmp(i1,1) - 1
        Jgrid(1,2,irank) = sum( ntmp(0:i2,2) ) - ntmp(i2,2)
        Jgrid(2,2,irank) = Jgrid(1,2,irank) + ntmp(i2,2) - 1
        Jgrid(1,3,irank) = sum( ntmp(0:i3,3) ) - ntmp(i3,3)
        Jgrid(2,3,irank) = Jgrid(1,3,irank) + ntmp(i3,3) - 1
        Jgrid(1,4,irank) = sum( ntmp(0:i4,4) ) - ntmp(i4,4) + 1
        Jgrid(2,4,irank) = Jgrid(1,4,irank) + ntmp(i4,4) - 1
        Jgrid(1,5,irank) = sum( ntmp(0:i5,5) ) - ntmp(i5,5) + 1
        Jgrid(2,5,irank) = Jgrid(1,5,irank) + ntmp(i5,5) - 1
        Jgrid(1,6,irank) = sum( ntmp(0:i6,6) ) - ntmp(i6,6) + 1
        Jgrid(2,6,irank) = Jgrid(1,6,irank) + ntmp(i6,6) - 1
     end do
     end do
     end do
  end do
  end do
  end do

!  do irank=0,nprocs_new-1
!     write(*,'(1x,i6,2x,6(2i6,1x))') irank,Jgrid(:,:,irank)
!  end do

!---

  allocate( itmp(nprocs_old)             ) ; itmp=0
  allocate( iomap(nprocs_old,nprocs_new) ) ; iomap=0
  allocate( nmap(nprocs_new)             ) ; nmap=0

  open(u2,file="input_wfiodir")
  write(u2,*) nprocs_new ! # of new processes

  do jrank=0,nprocs_new-1

     a1b = Jgrid(1,1,jrank)
     b1b = Jgrid(2,1,jrank)
     a2b = Jgrid(1,2,jrank)
     b2b = Jgrid(2,2,jrank)
     a3b = Jgrid(1,3,jrank)
     b3b = Jgrid(2,3,jrank)
     a4b = Jgrid(1,4,jrank)
     b4b = Jgrid(2,4,jrank)
     a5b = Jgrid(1,5,jrank)
     b5b = Jgrid(2,5,jrank)
     a6b = Jgrid(1,6,jrank)
     b6b = Jgrid(2,6,jrank)

     write(*,'(1x,"new rank: ",i6,1x,12i6)') jrank,a1b,b1b,a2b,b2b,a3b,b3b &
          ,a4b,b4b,a5b,b5b,a6b,b6b

     n=0
     do irank=0,nprocs_old-1

        a1a = Igrid(1,1,irank)
        b1a = Igrid(2,1,irank)
        a2a = Igrid(1,2,irank)
        b2a = Igrid(2,2,irank)
        a3a = Igrid(1,3,irank)
        b3a = Igrid(2,3,irank)
        a4a = Igrid(1,4,irank)
        b4a = Igrid(2,4,irank)
        a5a = Igrid(1,5,irank)
        b5a = Igrid(2,5,irank)
        a6a = Igrid(1,6,irank)
        b6a = Igrid(2,6,irank)

        if ( (a1b <= b1a .and. b1a <= b1b .or. a1b <= a1a .and. a1a <= b1b) &
        .and.(a2b <= b2a .and. b2a <= b2b .or. a2b <= a2a .and. a2a <= b2b) &
        .and.(a3b <= b3a .and. b3a <= b3b .or. a3b <= a3a .and. a3a <= b3b) &
        .and.(a4b <= b4a .and. b4a <= b4b .or. a4b <= a4a .and. a4a <= b4b) &
        .and.(a5b <= b5a .and. b5a <= b5b .or. a5b <= a5a .and. a5a <= b5b) &
        .and.(a6b <= b6a .and. b6a <= b6b .or. a6b <= a6a .and. a6a <= b6b) &
           ) then
           write(*,'(1x,"    rank: ",i6,1x,12i6)') irank &
                ,a1a,b1a,a2a,b2a,a3a,b3a,a4a,b4a,a5a,b5a,a6a,b6a
           n=n+1
           itmp(n)=irank
           iomap(n,jrank+1) = irank
        end if

     end do ! irnak(old)

     nmap(jrank+1) = n

     write(u2,*) n ! # of related old processes
     do i=1,n
        write(u2,*) itmp(i)
     end do

  end do ! jrank(new)

  do jrank=0,nprocs_new-1
     write(u2,'(1x,"new_rank: ",i6,1x,12i6)') jrank,Jgrid(:,:,jrank)
     do i=1,nmap(jrank+1)
        irank=iomap(i,jrank+1)
        write(u2,'(1x,"    rank: ",i6,1x,12i6)') irank,Igrid(:,:,irank)
     end do
  end do
        
  close(u2)

  deallocate( nmap  )
  deallocate( iomap )
  deallocate( itmp  )

!---

  deallocate( Jgrid )
  deallocate( Igrid )
  deallocate( ntmp  )

END PROGRAM wfiomap
