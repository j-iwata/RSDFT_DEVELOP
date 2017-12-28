PROGRAM ufld2gp

  implicit none
  integer,parameter :: u1=1, u2=10, u3=11
  integer,parameter :: maxloop=1000000000
  integer :: loop,ik,k,n,i,j,s,ne,nk,nb,ns,kmax,nmax,nmin,m
  integer :: nkp(0:9),nbp(0:9),nsp(0:9),kmaxp(0:9),nmaxp(0:9),nminp(0:9)
  integer :: kminp(0:9),k0,gmax(0:9),ig,ngp(0:9),ng,jg,ig1,ig2
  integer :: loop_read, n_loop_read
  integer :: unit_conv=1, unit_uf=30
  character(1) :: c1
  character(5) :: cbuf
  character(12) :: file_name
  character(128) :: cbuf_long
  real(8),parameter :: aB_=0.529177d0, HT_=27.2116d0
  real(8) :: aB, HT
  real(8),allocatable :: kxyz_pc(:,:),kxyz_sc(:,:)
  real(8),allocatable :: esp(:,:,:,:), weight(:,:,:,:)
  real(8) :: emax,emin,e_mergin,eta,de,pi,f(2),e,x,dummy(4),d
  real(8) :: z1,z2,ee(2),ww(2)
  logical :: flag

! ---

  ne = 500

  eta = 0.10d0/27.2116d0

  pi = acos(-1.0d0)

! ---

!  write(*,*) "eta(eV), ne="
!  read(*,*) eta, ne
!  eta = eta/HT_

! ---

  nsp(:)=0
  nbp(:)=0
  nkp(:)=0
  ngp(:)=0
  kminp(:)=1000000000
  kmaxp(:)=0
  nmaxp(:)=0
  nminp(:)=1000000000
  gmax(:)=0

  do loop_read=0,9

     if ( loop_read == 0 ) then
        file_name="band_ufld_z"
     else
        write(c1,'(i1)') loop_read
        file_name = "band_ufld_z"//c1
     end if

     inquire( FILE=file_name, EXIST=flag )
     if ( .not.flag ) exit

     write(*,*) loop_read, file_name

     open(u1,file=file_name,status="old")

     do loop=1,maxloop

        read(u1,*,END=1) cbuf,k

        if ( cbuf == "iktrj" ) then

           kminp(loop_read)=min(k,kminp(loop_read))
           kmaxp(loop_read)=max(k,kmaxp(loop_read))
           nkp(loop_read)=nkp(loop_read)+1
           if ( nbp(loop_read) /= 0 ) then
              nminp(loop_read)=min(nminp(loop_read),nbp(loop_read))
           end if
           nbp(loop_read)=0

        else

           backspace(u1)
           read(u1,*) i,j,ig
           gmax(loop_read)=max(gmax(loop_read),ig)

           nbp(loop_read)=nbp(loop_read)+1
           nmaxp(loop_read)=max(nmaxp(loop_read),nbp(loop_read))

        end if

     end do ! loop
1    continue

     ngp(loop_read)=gmax(loop_read)+1

     rewind u1
     read(u1,*)
     read(u1,'(a)') cbuf_long
     read(cbuf_long,*,END=2) i,j,dummy(1:4)
     nsp(loop_read)=2
2    continue
     if ( nsp(loop_read) == 0 ) nsp(loop_read)=1

     if ( loop_read > 0 .and. kminp(loop_read) == 1 ) then
        kminp(loop_read)=kminp(loop_read)+kmaxp(loop_read-1)
        kmaxp(loop_read)=kmaxp(loop_read)+kmaxp(loop_read-1)
     end if

     write(*,'(a50,"loop_read=",i1)') repeat("-",50),loop_read
     write(*,*) file_name, loop_read
     write(*,*) "ns, nk, nb=",nsp(loop_read),nkp(loop_read),nbp(loop_read)
     write(*,*) "kmin,max=",kminp(loop_read),kmaxp(loop_read)
     write(*,*) "nmin=",nminp(loop_read)/ngp(loop_read),mod(nminp(loop_read),ngp(loop_read))
     write(*,*) "nmax=",nmaxp(loop_read)/ngp(loop_read),mod(nmaxp(loop_read),ngp(loop_read))
     write(*,*) "gmax, ng=",gmax(loop_read),ngp(loop_read)

     nminp(loop_read)=nminp(loop_read)/ngp(loop_read)
     nmaxp(loop_read)=nmaxp(loop_read)/ngp(loop_read)
     nbp(loop_read)=nbp(loop_read)/ngp(loop_read)

     close(u1)

  end do ! loop_read

  n_loop_read=loop_read-1

  nk=sum(nkp)
  ns=maxval(nsp)
  nmin=minval(nmaxp(0:n_loop_read))
  nmax=maxval(nmaxp)
  ng=maxval(ngp)

  write(*,*)
  write(*,*) "ns, nk=",ns,nk
  write(*,*) "nmin,nmax=",nmin,nmax
  write(*,*) "ng=",ng
  write(*,*)

! ---

  allocate( kxyz_pc(3,0:nk)       ) ; kxyz_pc=0.0d0
  allocate( kxyz_sc(3,nk  )       ) ; kxyz_sc=0.0d0
  allocate( esp(nmax,ng,nk,ns)    ) ; esp=0.0d0
  allocate( weight(nmax,ng,nk,ns) ) ; weight=0.0d0

! ---

  do loop_read=0,n_loop_read

     if ( loop_read == 0 ) then
        file_name="band_ufld_z"
     else
        write(c1,'(i1)') loop_read
        file_name = "band_ufld_z"//c1
     end if

     k0=0
!     if ( loop_read > 0 ) k0=kmaxp(loop_read-1)

     open(u1,file=file_name,status="old")

     do ik=1,nkp(loop_read)

        read(u1,*) cbuf,k,kxyz_pc(1:3,k+k0),kxyz_sc(1:3,k+k0)

        k=k+k0

!        if ( k == kminp(loop_read) .or. k == kmaxp(loop_read) ) then
        if ( .true. ) then
           write(*,'(1x,3i4,2(3f15.8,2x))') &
                loop_read, k, k0, kxyz_pc(:,k), kxyz_sc(:,k)
        end if

        do ig=1,ngp(loop_read)
        do n=1,nmaxp(loop_read)
           read(u1,*) i,j,jg,(esp(n,ig,k,s),weight(n,ig,k,s),s=1,nsp(loop_read))
           read(u1,*,END=9) cbuf
           backspace(u1)
           if ( cbuf == "iktrj" ) exit
        end do
        end do

     end do
9    continue

     close(u1)

  end do ! loop_read

! ---

  write(*,*)
  do k=1,nk
     write(*,'(1x,i4,2(3f15.8,2x))') k,kxyz_pc(:,k),kxyz_sc(:,k)
  end do
  write(*,*)
  do loop=0,n_loop_read
     do k=2,nk
        d=sum((kxyz_pc(:,k)-kxyz_pc(:,k-1))**2)
        if ( d < 1.d-8 ) then
           write(*,'(1x,i4,2(3f15.8,2x))') k,kxyz_pc(:,k),kxyz_pc(:,k-1)
           do ik=k,nk-1
              kxyz_pc(:,ik)=kxyz_pc(:,ik+1)
              kxyz_sc(:,ik)=kxyz_sc(:,ik+1)
              esp(:,:,ik,:)=esp(:,:,ik+1,:)
              weight(:,:,ik,:)=weight(:,:,ik+1,:)
           end do
           nk=nk-1
           exit
        end if
     end do ! k
  end do ! loop
  write(*,*)
  do k=1,nk
     write(*,'(1x,i4,2(3f15.8,2x))') k,kxyz_pc(:,k),kxyz_sc(:,k)
  end do
  write(*,*)

! ---

!  nk=sum(nkp)
  write(*,*) "nk(new)=",nk

! ---

  do loop=1,maxloop

!  write(*,*) "input z range in lattice coordinate [0-1] = ? ?"
!  read(*,*) z1,z2
!  if ( z1 < 0.0d0 .or. z2 < 0.0d0 ) exit

  write(*,*) "input two grid indices in z direction [0 - ngz-1] = ? ?"
  write(*,*) "( negative values: quit )"
  read(*,*) ig1,ig2
  if ( ig1 < 0 .or. ig2 < 0 ) exit
  ig1=ig1+1
  ig2=ig2+1


!  ig1 = nint( z1*ng ) ; if ( ig1 == 0 ) ig1=1
!  ig2 = nint( z2*ng )

  write(*,*) "ig1=",ig1,dble(ig1)/ng,z1
  write(*,*) "ig2=",ig2,dble(ig2)/ng,z2

  do k=1,nk

     write(unit_uf,'(1x,"iktrj=",i6,2x,2(3f10.6,2x))') k,kxyz_pc(1:3,k),kxyz_sc(1:3,k)

     do n=1,nmin

        do s=1,ns
           ee(s) = esp(n,ig1,k,s)
           ww(s) = sum( weight(n,ig1:ig2,k,s) )
        end do
        write(unit_uf,'(1x,2i6,2(f20.15,2x,g20.10,2x))') k,n,( ee(s),ww(s),s=1,ns )

     end do

  end do

  unit_uf=unit_uf+1

  end do ! loop

! ---

  deallocate( weight  )
  deallocate( esp     )
  deallocate( kxyz_sc )
  deallocate( kxyz_pc )

END PROGRAM ufld2gp
