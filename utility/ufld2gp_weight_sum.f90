PROGRAM ufld2gp_weight_sum

  implicit none
  integer,parameter :: u1=1, u2=10, u3=11
  integer,parameter :: maxloop=1000000
  integer :: loop,ik,k,n,i,j,s,ne,nk,nb,ns,kmax,nmax,nmin
  integer :: nkp(0:9),nbp(0:9),nsp(0:9),kmaxp(0:9),nmaxp(0:9),nminp(0:9)
  integer :: kminp(0:9),k0
  integer :: loop_read, n_loop_read
  integer :: unit_conv=1
  character(1) :: c1
  character(5) :: cbuf
  character(10) :: file_name
  character(128) :: cbuf_long
  real(8),parameter :: aB_=0.529177d0, HT_=27.2116d0
  real(8) :: aB, HT
  real(8),allocatable :: kxyz_pc(:,:),kxyz_sc(:,:)
  real(8),allocatable :: esp(:,:,:), weight(:,:,:)
  real(8) :: emax,emin,e_mergin,eta,de,pi,f(2),e,x,dummy(4),d
  real(8) :: weight_tmp(2)
  logical :: flag

! ---

  ne = 500

  eta = 0.10d0/27.2116d0

  pi = acos(-1.0d0)

! ---

  write(*,*) "eta(eV), ne="
  read(*,*) eta, ne
  eta = eta/HT_

! ---

  nsp(:)=0
  nbp(:)=0
  nkp(:)=0
  kminp(:)=1000000000
  kmaxp(:)=0
  nmaxp(:)=0
  nminp(:)=1000000000

  do loop_read=0,9

     if ( loop_read == 0 ) then
        file_name="band_ufld"
     else
        write(c1,'(i1)') loop_read
        file_name = "band_ufld"//c1
     end if

     inquire( FILE=file_name, EXIST=flag )
     if ( .not.flag ) exit

     write(*,*) loop_read,file_name

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
           nbp(loop_read)=nbp(loop_read)+1
           nmaxp(loop_read)=max(nmaxp(loop_read),nbp(loop_read))
        end if
     end do ! loop
1    continue

     rewind u1
     read(u1,*)
     read(u1,'(a)') cbuf_long
     read(cbuf_long,*,END=2) i,j,dummy(1:4)
     nsp(loop_read)=2
2    continue
     if ( nsp(loop_read) == 0 ) nsp(loop_read)=1

!     if ( loop_read > 0 .and. kminp(loop_read) == 1 ) then
!        kminp(loop_read)=kminp(loop_read)+kmaxp(loop_read-1)
!        kmaxp(loop_read)=kmaxp(loop_read)+kmaxp(loop_read-1)
!     end if

     write(*,'(a50,"loop_read=",i1)') repeat("-",50),loop_read
     write(*,*) file_name, loop_read
     write(*,*) "ns, nk, nb=",nsp(loop_read),nkp(loop_read),nbp(loop_read)
     write(*,*) "kmin,max=",kminp(loop_read),kmaxp(loop_read)
     write(*,*) "nmin,nmax=",nminp(loop_read),nmaxp(loop_read)

     close(u1)

  end do ! loop_read

  n_loop_read=loop_read-1

!  nk=sum(nkp)
  nk=maxval(nkp) ; if ( any(nkp(1:n_loop_read)/=nk) ) stop "stop(1)"
  ns=maxval(nsp)
  nmin=minval(nmaxp(0:n_loop_read))
  nmax=maxval(nmaxp)

  write(*,*) "ns, nk=",ns,nk
  write(*,*) "nmin,nmax=",nmin,nmax

! ---

  allocate( kxyz_pc(3,0:nk)    ) ; kxyz_pc=0.0d0
  allocate( kxyz_sc(3,nk)      ) ; kxyz_sc=0.0d0
  allocate( esp(nmax,nk,ns)    ) ; esp=0.0d0
  allocate( weight(nmax,nk,ns) ) ; weight=0.0d0

! ---

  do loop_read=0,n_loop_read

     if ( loop_read == 0 ) then
        file_name="band_ufld"
     else
        write(c1,'(i1)') loop_read
        file_name = "band_ufld"//c1
     end if

     k0=0
!     if ( loop_read > 0 ) k0=kmaxp(loop_read-1)

     open(u1,file=file_name,status="old")

     do ik=1,nkp(loop_read)

        read(u1,*) cbuf,k,kxyz_pc(1:3,k+k0),kxyz_sc(1:3,k+k0)

        k=k+k0

        if ( k == kminp(loop_read) .or. k == kmaxp(loop_read) ) then
           write(*,'(1x,2i4,2(3f15.8,2x))') &
                loop_read, k, kxyz_pc(:,k), kxyz_sc(:,k)
        end if

        do n=1,nmaxp(loop_read)
           read(u1,*) i,j,(esp(n,k,s),weight_tmp(s),s=1,nsp(loop_read))
           do s=1,nsp(loop_read)
              weight(n,k,s)=weight(n,k,s)+weight_tmp(s)
           end do
           read(u1,*,END=9) cbuf
           backspace(u1)
           if ( cbuf == "iktrj" ) exit
        end do

     end do
9    continue

!     if ( loop_read > 0 ) then
!        d=sum((kxyz_pc(:,kminp(loop_read))-kxyz_pc(:,kmaxp(loop_read-1)))**2)
!        if ( d < 1.d-8 ) then
!           write(*,*) "omit the first k points"
!           do k=kminp(loop_read),kmaxp(loop_read)-1
!              kxyz_pc(:,k)=kxyz_pc(:,k+1)
!              esp(:,k,:)=esp(:,k+1,:)
!              weight(:,k,:)=weight(:,k+1,:)
!           end do
!           kminp(loop_read)=kminp(loop_read)+1
!           nkp(loop_read)=nkp(loop_read)-1
!        end if
!     end if

     close(u1)

  end do ! loop_read

!  nk=sum(nkp)
!  write(*,*) "nk(new)=",nk

! ---

  emin = minval( esp(1:nmin,:,:) )
  emax = maxval( esp(1:nmin,:,:) )

  write(*,*) "emin, emax=",emin,emax

  e_mergin = (emax-emin)*0.1d0

  emax = emax + e_mergin
  emin = emin - e_mergin

  write(*,*) "emin, emax, e_mergin=",emin,emax,e_mergin

  de = (emax-emin)/ne

  write(*,*) "ne, de=",ne,de

  write(*,*) "eta=",eta

! ---

  write(*,*) "Chose the unit [ 0: Hartree a.u.,  1: eV & angstrome^-1 ]"
  read(*,*) unit_conv

  if ( unit_conv == 0 ) then
     aB=1.0d0
     HT=1.0d0
  else
     aB=aB_
     HT=HT_
  end if

! ---

  rewind u2
  rewind u3

  kxyz_pc(:,0) = kxyz_pc(:,1)
  x=0.0d0

  do k=1,nk

     x = x + sqrt( sum( (kxyz_pc(:,k)-kxyz_pc(:,k-1))**2 ) )/aB

     do i=1,ne

        e = emin + (i-1)*de

        do s=1,ns
           f(s)=0.0d0
           do n=1,nmin
              f(s) = f(s) + eta/( (e-esp(n,k,s))**2 + eta**2 )*weight(n,k,s)
           end do
           f(s)=f(s)/pi
        end do

        write(u2,'(1x,4f25.15)') x,e*HT,(f(s),s=1,ns)
        write(u3,'(1x,4f25.15)') x,e*HT,(log(f(s)),s=1,ns)

     end do ! i

     write(u2,*)
     write(u3,*)

  end do ! k

! ---

  deallocate( weight  )
  deallocate( esp     )
  deallocate( kxyz_sc )
  deallocate( kxyz_pc )

END PROGRAM ufld2gp_weight_sum
