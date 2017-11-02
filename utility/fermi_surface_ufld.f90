PROGRAM fermi_surface_ufld

  !use libtetrabz, only: libtetrabz_dos
  !use libtetrabz_mpi, only: libtetrabz_mpi_dos

  implicit none

  integer,parameter :: u1=1, u2=2, u5=5, u99=99, u10=10
  integer,parameter :: max_loop=10000000
  integer :: ltetra
  integer :: k,i,j,k1,k2,k3,i1,i2,i3,j1,j2,j3,n,n1,n2,ne,s,idummy
  integer :: nk, n_irwedge, n_whole, indx_range(2,3)
  integer :: nge(3),ngw(3),ngf(3),ngg(3),nfac,nbnd
  integer :: nband,nspin,nkpnt,q1,q2,q3
  integer :: myrank,ie,Nextend(3)
  integer,allocatable :: kbz(:,:),ktmp(:,:,:),kgrd(:,:,:)
  integer,allocatable :: itmp(:,:),kxyz(:,:,:,:)
  real(8),allocatable :: eval(:,:,:),occp(:,:,:)
  real(8),allocatable :: wbz(:),kpnt(:,:),weight(:)
  real(8),allocatable :: eig(:,:,:,:),wght(:,:,:,:,:),e0(:)
  real(8),allocatable :: www(:,:,:,:),ufnld(:,:,:,:)
  real(8),allocatable :: ufnld_ext(:,:,:)
  real(8),parameter :: HT=27.2116d0
  real(8) :: emin,emax,de,e1,e2,e,enk,eta,emargin
  real(8) :: bb(3,3), ktry(3), err, etime(0:1),et(0:6)
  real(8) :: aa_pc(3,3),ax_pc,bb_pc(3,3)
  real(8) :: rdummy(3),dk(3)
  character(80) :: cbuf
  character(16) :: file_bz, file_eigv
  character(8) :: cbuf3(3)
!  include 'mpif.h'
  complex(8),allocatable :: zw0(:,:,:), zw1(:,:,:), zw2(:,:,:)
  complex(8),allocatable :: zwork(:,:,:), zwfft(:,:,:)
  complex(8),parameter :: zi=(0.0d0,1.0d0), z0=(0.0d0,0.0d0)
  real(8) :: kr
  real(8),parameter :: pi=3.141592653589793d0
  real(8),allocatable :: ee0(:,:,:), ww0(:,:,:)
  logical :: flag_exist

!  call MPI_INIT(i)
!  call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, i )
  myrank=0
!  etime(0)=mpi_wtime()

  file_bz   = "bz_info"
!  file_eigv = "eigenvalues"
  file_eigv = "band_ufld"

!------------------------------

  call read_aa_pc( ax_pc, aa_pc )

!------------------------------

  open(u1,file=file_bz,status='old')
  read(u1,*)
  read(u1,*) bb(1:3,1)
  read(u1,*) bb(1:3,2)
  read(u1,*) bb(1:3,3)
  read(u1,'(a)') cbuf
  call get_num_from_chr( cbuf, nk )
  read(u1,'(a)') cbuf
  call get_num_from_chr( cbuf, n_irwedge )
  read(u1,'(a)') cbuf
  call get_num_from_chr( cbuf, n_whole )

  read(u1,*)
  read(u1,*)

  allocate( kbz(3,n_irwedge) ) ; kbz=0
  allocate( wbz(n_irwedge)   ) ; wbz=0.0d0

  do i=1,n_irwedge
     read(u1,*) k, kbz(1:3,i), wbz(i)
  end do

  do i=1,3
     read(u1,'(a)') cbuf
     call get_nums_from_chr( cbuf, indx_range(:,i) )
  end do

  allocate( ktmp(indx_range(1,1):indx_range(2,1) &
                ,indx_range(1,2):indx_range(2,2) &
                ,indx_range(1,3):indx_range(2,3)) ) ; ktmp=0

  read(u1,*)

  do i=1,n_whole
     read(u1,*) k,k1,k2,k3,ktmp(k1,k2,k3)
  end do

  close(u1)

! ---

  k=max( size(ktmp,1), size(ktmp,2), size(ktmp,3) )
  allocate( itmp(k,3) ) ; itmp=0 

  loop_3 : do k3=indx_range(1,3),indx_range(2,3)
           do k2=indx_range(1,2),indx_range(2,2)
              if ( any(ktmp(:,k2,k3)/=0) ) exit loop_3
           end do
           end do loop_3
  nge(1) = count( ktmp(:,k2,k3) /= 0 )
  i=0
  do k1=indx_range(1,1),indx_range(2,1)
     if ( ktmp(k1,k2,k3) /= 0 ) then
        i=i+1
        itmp(i,1)=k1
     end if
  end do

  loop_1 : do k1=indx_range(1,1),indx_range(2,1)
           do k3=indx_range(1,3),indx_range(2,3)
              if ( any(ktmp(k1,:,k3)/=0) ) exit loop_1
           end do
           end do loop_1
  nge(2) = count( ktmp(k1,:,k3) /= 0 )
  i=0
  do k2=indx_range(1,2),indx_range(2,2)
     if ( ktmp(k1,k2,k3) /= 0 ) then
        i=i+1
        itmp(i,2)=k2
     end if
  end do

  loop_2 : do k2=indx_range(1,2),indx_range(2,2)
           do k1=indx_range(1,1),indx_range(2,1)
              if ( any(ktmp(k1,k2,:)/=0) ) exit loop_2
           end do
           end do loop_2
  nge(3) = count( ktmp(k1,k2,:) /= 0 )
  i=0
  do k3=indx_range(1,3),indx_range(2,3)
     if ( ktmp(k1,k2,k3) /= 0 ) then
        i=i+1
        itmp(i,3)=k3
     end if
  end do

  write(*,*) "size(ktmp),nge",size(ktmp,1),nge(1)
  write(*,*) "size(ktmp),nge",size(ktmp,2),nge(2)
  write(*,*) "size(ktmp),nge",size(ktmp,3),nge(3)

! ---

  allocate( kgrd(nge(1),nge(2),nge(3))   ) ; kgrd=0
  allocate( kxyz(3,nge(1),nge(2),nge(3)) ) ; kxyz=0

  do i3=1,nge(3)
  do i2=1,nge(2)
  do i1=1,nge(1)
     kgrd(i1,i2,i3) = ktmp(itmp(i1,1),itmp(i2,2),itmp(i3,3))
     kxyz(1,i1,i2,i3) = itmp(i1,1)
     kxyz(2,i1,i2,i3) = itmp(i2,2)
     kxyz(3,i1,i2,i3) = itmp(i3,3)
  end do
  end do
  end do

  write(*,'(1x,"itmp(:,1):",16i3)') itmp(1:nge(1),1)
  write(*,'(1x,"itmp(:,2):",16i3)') itmp(1:nge(2),2)
  write(*,'(1x,"itmp(:,3):",16i3)') itmp(1:nge(3),3)

! ---

  open(u2,file=file_eigv)

  k1=0
  k2=0
  n1=0
  do i=1,max_loop
     read(u2,*,END=9) cbuf
     if ( cbuf == "iktrj=" ) then
        backspace(u2)
        read(u2,*,END=9) cbuf,k
        k1=max(k1,k)
        k2=k2+1
     else
        read(u2,*,END=9) j,n
        n1=max(n1,n)
     end if
  end do
9 continue

  nband = n1
  nkpnt = k1
  write(*,*) "nband, nkpnt =", nband, nkpnt

  allocate( eig(nband,nge(1),nge(2),nge(3)) ) ; eig=0.0d0
  allocate( www(nband,nge(1),nge(2),nge(3)) ) ; www=0.0d0

  aa_pc(:,:)=ax_pc*aa_pc(:,:)

  call calc_bb( aa_pc, bb_pc )

  rewind u2  ! file_eigv

  do i=1,nkpnt
     read(u2,*,END=8) cbuf3(1:2), ktry(1:3)
     rdummy=matmul(transpose(aa_pc),ktry)/(2.0d0*pi)*nge
     write(*,'(1x,i4,3f10.5)') i,rdummy
     i1 = nint( rdummy(1) )
     i2 = nint( rdummy(2) )
     i3 = nint( rdummy(3) )
     j1 = mod(i1+nge(1),nge(1))+1
     j2 = mod(i2+nge(2),nge(2))+1
     j3 = mod(i3+nge(3),nge(3))+1
!     write(*,'(1x,3(3i3,2x))') i1,i2,i3,j1,j2,j3,kxyz(:,j1,j2,j3)
     do j=1,nband
        read(u2,*) rdummy(1:2), eig(j,j1,j2,j3), www(j,j1,j2,j3)
     end do
     k1 = mod(-i1+nge(1),nge(1))+1
     k2 = mod(-i2+nge(2),nge(2))+1
     k3 = mod(-i3+nge(3),nge(3))+1
     do j=1,nband
        eig(j,k1,k2,k3)=eig(j,j1,j2,j3)
        www(j,k1,k2,k3)=www(j,j1,j2,j3)
     end do
  end do
8 continue

  close(u2)

  rewind 20
  i3=0
  do i2=-(nge(2)-1)/2,nge(2)/2
  do i1=-(nge(1)-1)/2,nge(1)/2
     j1=mod(i1+nge(1),nge(1))+1
     j2=mod(i2+nge(2),nge(2))+1
     j3=mod(i3+nge(3),nge(3))+1
     write(20,'(1x,2i4,3g20.10)') i1,i2,eig(1,j1,j2,j3),www(1,j1,j2,j3) &
          ,eig(1,j1,j2,j3)*www(1,j1,j2,j3)
  end do
  write(20,*)
  end do
  rewind 21
  i3=0
  i2=0
  do i1=-(nge(1)-1)/2,nge(1)/2
     j1=mod(i1+nge(1),nge(1))+1
     j2=mod(i2+nge(2),nge(2))+1
     j3=mod(i3+nge(3),nge(3))+1
     write(21,'(1x,i4,2g20.10)') i1,eig(1,j1,j2,j3),www(1,j1,j2,j3)
  end do

!------------------------------

  emin = minval( eig )
  emax = maxval( eig )

  if ( myrank == 0 ) then
     write(*,*) "emin     (HT,eV)=",emin,emin*HT
     write(*,*) "emax     (HT,eV)=",emax,emax*HT
     write(*,*) "emax-emin(HT,eV)=",emax-emin,(emax-emin)*HT
  end if

!------------------------------ FFT

  nk=nge(1)
  write(*,*) "nk=",nk

  write(*,*) "nfac=?"
  read(*,*) nfac
!  nfac=8
  write(*,*) "nfac=",nfac
  ngf(1:3) = nge(1:3)*nfac

  ne=300
  write(*,*) "eta(eV)=?"
  read(*,*) eta
  write(*,*) "e1,e2(HT)=?"
  read(*,*) e1,e2
  write(*,*) "e1(ht,eV)=",e1,e1*HT
  write(*,*) "e2(ht,eV)=",e2,e2*HT
  de=(e2-e1)/ne
  write(*,*) "de(ht,eV)=",de,de*HT
  emargin=1.0d0/HT

  allocate( ee0(-nge(1):nge(1),-nge(2):nge(2),-nge(3):nge(3)) )
  ee0=0.0d0
  allocate( ww0(-nge(1):nge(1),-nge(2):nge(2),-nge(3):nge(3)) )
  ww0=0.0d0

  allocate( zw0(-nge(1):nge(1),-nge(2):nge(2),-nge(3):nge(3)) )
  zw0=(0.0d0,0.0d0)
  allocate( zw1(-ngf(1):ngf(1),-ngf(2):ngf(2),-ngf(3):ngf(3)) )
  zw1=(0.0d0,0.0d0)
  allocate( zw2(-ngf(1):ngf(1),-ngf(2):ngf(2),-ngf(3):ngf(3)) )
  zw2=(0.0d0,0.0d0)

  allocate( ufnld(-ngf(1):ngf(1),-ngf(2):ngf(2),-ngf(3):ngf(3),ne) )
  ufnld=0.0d0

  allocate( zwork(ngf(1),ngf(2),ngf(3)) ) ; zwork=(0.0d0,0.0d0)
  allocate( zwfft(ngf(1),ngf(2),ngf(3)) ) ; zwfft=(0.0d0,0.0d0)

  do n=1,nband

     if ( maxval(eig(n,:,:,:)) < e1-emargin .or. &
          minval(eig(n,:,:,:)) > e2+emargin ) cycle

!     et(0)=mpi_wtime()

     do j3=-(nge(3)-1)/2,nge(3)/2
     do j2=-(nge(2)-1)/2,nge(2)/2
     do j1=-(nge(1)-1)/2,nge(1)/2
        i1=mod(j1+nge(1),nge(1))+1
        i2=mod(j2+nge(2),nge(2))+1
        i3=mod(j3+nge(3),nge(3))+1
        ee0(j1,j2,j3) = eig(n,i1,i2,i3)
        ww0(j1,j2,j3) = www(n,i1,i2,i3)
     end do
     end do
     end do

!     et(1)=mpi_wtime()

     call check_file( nfac, n, flag_exist )

     if ( flag_exist ) then

        call read_file( nfac, n, zw1, zw2 )

     else

        call dft_forward( nk, nfac, nge, ngf, ee0, zw0, zw1 )
        call dft_forward( nk, nfac, nge, ngf, ww0, zw0, zw2 )

        call write_file( nfac, n, zw1, zw2 )

     end if

!     et(2)=mpi_wtime()

     do ie=1,ne
        e=e1+(ie-1)*de
        do i3=-ngf(3)/2,ngf(3)/2
           if ( i3 /= 0 ) cycle
        do i2=-ngf(2)/2,ngf(2)/2
        do i1=-ngf(1)/2,ngf(1)/2
           enk=real( zw1(i1,i2,i3) )
           ufnld(i1,i2,i3,ie) = ufnld(i1,i2,i3,ie) &
                + real(zw2(i1,i2,i3)) * eta/( (e-enk)**2 + eta**2 )/pi
        end do
        end do
        end do
     end do ! ie

!     et(3)=mpi_wtime()

     write(*,'(1x,i5,2f15.10,2x,3f8.3)') &
          n, minval(eig(n,:,:,:)),maxval(eig(n,:,:,:)) &
          ,et(1)-et(0),et(2)-et(1),et(3)-et(2)

  end do ! n

! ---

  write(*,*) "Do you need zone extention ? ( Nextend(1:3)=?,?,? )"
  read(*,*) Nextend(1:3)
!  Nextend(1:3)=(/ 2, 2, 1 /)

  ngg(1:3) = ngf(1:3)*Nextend(1:3)
  write(*,'(1x,6i4)') ngf,ngg

  allocate( ufnld_ext(-ngg(1):ngg(1),-ngg(2):ngg(2),-ngg(3):ngg(3)) )
  ufnld_ext=0.0d0

  do i3=-ngf(3)/2,ngf(3)/2
  do i2=-ngf(2)/2,ngf(2)/2
  do i1=-ngf(1)/2,ngf(1)/2
     ufnld_ext(i1,i2,i3) = sum( ufnld(i1,i2,i3,:) )
  end do
  end do
  end do
  call periodic_image( ngf(1)/2,ngf(2)/2,ngf(3)/2,ngg(1),ngg(2),ngg(3),ufnld_ext )

  dk(1) = sqrt( sum( bb_pc(:,1)**2 ) )/ngf(1)
  dk(2) = sqrt( sum( bb_pc(:,2)**2 ) )/ngf(2)

  write(*,*) '|b1|',sqrt(sum(bb_pc(:,1)**2))
  write(*,*) '|b2|',sqrt(sum(bb_pc(:,2)**2))

  !i1=25
  !i2=25
  !write(*,*) "i1,i2=? ?"
  !read(*,*) i1,i2
  !rdummy(1) = i1*bb_pc(1,1)/ngf(1)+i2*bb_pc(1,2)/ngf(2)
  !rdummy(2) = i1*bb_pc(2,1)/ngf(1)+i2*bb_pc(2,2)/ngf(2)
  !write(*,*) sqrt(sum(rdummy(1:2)**2))

  rewind 10
  do i2=-ngf(2)/2,ngf(2)/2
  do i1=-ngf(1)/2,ngf(1)/2
     write(10,'(1x,2i5,2g20.10)') i1, i2, sum( ufnld(i1,i2,0,:) ) &
                                  ,log(abs(sum(ufnld(i1,i2,0,:))))
!     write(10,'(1x,2f10.5,2g20.10)') i1*dk(1), i2*dk(2) &
!                                 , sum( ufnld(i1,i2,0,:) ) &
!                                  ,log(abs(sum(ufnld(i1,i2,0,:))))
  end do
  write(10,*)
  end do

  rewind 12
  do i2=-ngg(2)/2,ngg(2)/2
  do i1=-ngg(1)/2,ngg(1)/2
     write(12,'(1x,2i5,2g20.10)') i1, i2, ufnld_ext(i1,i2,0) &
                 ,log(abs(ufnld_ext(i1,i2,0)))
!     write(12,'(1x,2f10.5,2g20.10)') i1*dk(1), i2*dk(2), ufnld_ext(i1,i2,0) &
!                 ,log(abs(ufnld_ext(i1,i2,0)))
  end do
  write(12,*)
  end do

!  rewind 11
!  i3=0
!  do i2=-ngf(2)/2,ngf(2)/2
!  do i1=-ngf(1)/2,ngf(1)/2
!     write(11,'(1x,2i4,3g20.10)') i1,i2,real(zw1(i1,i2,i3)) &
!                                      ,real(zw2(i1,i2,i3)) &
!                                      ,real(zw1(i1,i2,i3))*real(zw2(i1,i2,i3))
!  end do
!  write(11,*)
!  end do

  deallocate( zw0, zw1 )
  deallocate( zwfft )
  deallocate( zwork )

!------------------------------

!  etime(1)=mpi_wtime()
!  if ( myrank == 0 ) write(*,*) "TIME=",etime(1)-etime(0)

!  call MPI_FINALIZE( i )

CONTAINS


  SUBROUTINE get_num_from_chr( cbuf, n )
    implicit none
    character(*),intent(IN) :: cbuf
    integer,intent(OUT) :: n
    integer :: i
    do i=1,len_trim(cbuf)
       if ( cbuf(i:i) == ":" ) exit
    end do
    read(cbuf(i+1:),*) n
  END SUBROUTINE get_num_from_chr


  SUBROUTINE get_nums_from_chr( cbuf, n )
    implicit none
    character(*),intent(IN) :: cbuf
    integer,intent(OUT) :: n(:)
    integer :: i
    do i=1,len_trim(cbuf)
       if ( cbuf(i:i) == ":" ) exit
    end do
    read(cbuf(i+1:),*) n(:)
  END SUBROUTINE get_nums_from_chr


  SUBROUTINE read_aa_pc( ax, aa )
    implicit none
    real(8),intent(OUT) :: ax, aa(3,3)
    character(10) :: cbuf
    rewind 1
10  read(1,*,end=99) cbuf
    if ( cbuf == "BANDUF" ) then
       read(1,*) ax
       read(1,*) aa(1:3,1)
       read(1,*) aa(1:3,2)
       read(1,*) aa(1:3,3)
       return
    end if
    goto 10
99  stop "BANDUF is not found"
  END SUBROUTINE read_aa_pc


  SUBROUTINE dft_forward( nk, nfac, nge, ngf, ff0, zw0, zw1 )
    implicit none
    integer,intent(IN) :: nk, nfac, nge(3),ngf(3)
    real(8),intent(IN) :: ff0(-nge(1):nge(1),-nge(2):nge(2),-nge(3):nge(3))
    complex(8),intent(OUT) :: zw0(-nge(1):nge(1),-nge(2):nge(2),-nge(3):nge(3))
    complex(8),intent(OUT) :: zw1(-ngf(1):ngf(1),-ngf(2):ngf(2),-ngf(3):ngf(3))
    complex(8),parameter :: z0=(0.0d0,0.0d0),zi=(0.0d0,1.0d0)
    integer :: j1,j2,j3,k1,k2,k3
    integer,allocatable :: itmp3(:,:,:)
    real(8) :: kR
    real(8),parameter :: pi=3.141592653589793d0
    zw0=z0
    do j3=-(nge(3)-1)/2,nge(3)/2
    do j2=-(nge(2)-1)/2,nge(2)/2
    do j1=-(nge(1)-1)/2,nge(1)/2
       do k3=-(nge(3)-1)/2,nge(3)/2
       do k2=-(nge(2)-1)/2,nge(2)/2
       do k1=-(nge(1)-1)/2,nge(1)/2
          kR=2.0d0*pi*( k1*j1 + k2*j2 + k3*j3 )/dble(nk)
          zw0(j1,j2,j3)=zw0(j1,j2,j3)+ff0(k1,k2,k3)*exp(zi*kR)
       end do
       end do
       end do
    end do
    end do
    end do
    zw0(:,:,:)=zw0(:,:,:)/dble( nge(1)*nge(2)*nge(3) )

    do j3=-(nge(3)-1)/2,nge(3)/2
    do j2=-(nge(2)-1)/2,nge(2)/2
    do j1=-(nge(1)-1)/2,nge(1)/2
       zw0(-j1,-j2,-j3) = conjg( zw0(j1,j2,j3) )
    end do
    end do
    end do
    do j3=-(nge(3)-1)/2,nge(3)/2
    do j2=-(nge(2)-1)/2,nge(2)/2
       zw0(-nge(1)/2,j2,j3) = zw0(nge(1)/2,j2,j3)
    end do
    end do
    do j3=-(nge(3)-1)/2,nge(3)/2
    do j1=-(nge(1)-1)/2,nge(1)/2
       zw0(j1,-nge(2)/2,j3) = zw0(j1,nge(2)/2,j3)
    end do
    end do
    do j3=-(nge(3)-1)/2,nge(3)/2
       zw0(-nge(1)/2,-nge(2)/2,j3) = zw0(nge(1)/2,nge(2)/2,j3)
    end do

    allocate( itmp3(-nge(1):nge(1),-nge(2):nge(2),-nge(3):nge(3)) )
    itmp3=0
    do j3=-(nge(3)-1)/2,nge(3)/2
    do j2=-(nge(2)-1)/2,nge(2)/2
    do j1=-(nge(1)-1)/2,nge(1)/2
       itmp3(j1,j2,j3)=1
    end do
    end do
    end do

    do j1=-(nge(1)-1)/2,nge(1)/2
       itmp3(j1, nge(2)/2,0)=2
       itmp3(j1,-nge(2)/2,0)=2
    end do
    do j2=-(nge(2)-1)/2,nge(2)/2
       itmp3( nge(1)/2,j2,0)=2
       itmp3(-nge(1)/2,j2,0)=2
    end do
    itmp3( nge(1)/2, nge(2)/2,0)=4
    itmp3(-nge(1)/2, nge(2)/2,0)=4
    itmp3( nge(1)/2,-nge(2)/2,0)=4
    itmp3(-nge(1)/2,-nge(2)/2,0)=4

    zw1=z0
    do k3=-ngf(3)/2,ngf(3)/2
    do k2=-ngf(2)/2,ngf(2)/2
    do k1=-ngf(1)/2,ngf(1)/2
       do j3=-nge(3)/2,nge(3)/2
       do j2=-nge(2)/2,nge(2)/2
       do j1=-nge(1)/2,nge(1)/2
          kR = 2.0d0*pi*( k1*j1 + k2*j2 + k3*j3 )/dble(nk*nfac)
          zw1(k1,k2,k3)=zw1(k1,k2,k3)+zw0(j1,j2,j3)*exp(-zi*kR) &
               /dble(itmp3(j1,j2,j3))
       end do
       end do
       end do
    end do
    end do
    end do
    deallocate( itmp3 )
  END SUBROUTINE dft_forward


  SUBROUTINE periodic_image( m1,m2,m3,n1,n2,n3,f )
    implicit none
    integer,intent(IN) :: m1,m2,m3,n1,n2,n3
    real(8),intent(INOUT) :: f(-n1:n1,-n2:n2,-n3:n3)
    integer :: i1,i2,i3,j1,j2,j3
    real(8),allocatable :: g(:,:,:)
    allocate( g(0:2*m1-1,0:2*m2-1,0:2*m3-1) ) ; g=0.0d0
    do i3=-m3+1,m3
    do i2=-m2+1,m2
    do i1=-m1+1,m1
       j1=mod(i1+2*m1,2*m1)
       j2=mod(i2+2*m2,2*m2)
       j3=mod(i3+2*m3,2*m3)
       g(j1,j2,j3) = f(i1,i2,i3)
    end do
    end do
    end do
    do i3=-n3,n3
       j3=i3 ; call sub1( 2*m3, j3 )
       do i2=-n2,n2
          j2=i2 ; call sub1( 2*m2, j2 )
          do i1=-n1,n1
             j1=i1 ; call sub1( 2*m1, j1 )
             !write(*,'(1x,6i4)') i1,i2,i3,j1,j2,j3
             f(i1,i2,i3) = g(j1,j2,j3)
          end do
       end do
    end do
    deallocate( g )
  END SUBROUTINE periodic_image

  SUBROUTINE sub1( m, i )
    implicit none
    integer,intent(IN) :: m
    integer,intent(INOUT) :: i
    integer :: j
    j=i
1   if ( j < -m ) then
       j=j+m
    else if ( m < j ) then
       j=j-m
    else
       i=mod(j+m,m)
       return
    end if
    goto 1
  END SUBROUTINE sub1


  SUBROUTINE check_file( nfac, n, file_exist )
    implicit none
    integer,intent(IN) :: nfac, n
    logical,intent(OUT) :: file_exist
    integer,parameter :: u=80
    character(8),parameter :: file_name="file_chk"
    integer,parameter :: ndata_max=10000
    integer :: nfac_tmp, n_tmp, i
    logical :: flag
    file_exist=.false.
    inquire( file=file_name, EXIST=flag )
    if ( flag ) then
       open(u,file=file_name,status="old")
       read(u,*) nfac_tmp
       if ( nfac_tmp == nfac ) then
          do i=1,ndata_max
             read(u,*,END=10) n_tmp
             if ( n_tmp == n ) then
                file_exist = .true.
                exit
             end if
          end do
10        continue
       end if
       close(u)
    end if
  END SUBROUTINE check_file


  SUBROUTINE read_file( nfac, n, z1, z2 )
    implicit none
    integer,intent(IN) :: nfac, n
    complex(8),intent(INOUT) :: z1(:,:,:),z2(:,:,:)
    integer,parameter :: u=81
    character(8),parameter :: file_name="file_dat"
    integer,parameter :: ndata_max=10000
    integer :: nfac_tmp, n_tmp, i
    open(u,file=file_name,status="old",form="unformatted")
    read(u) nfac_tmp
    do i=1,ndata_max
       read(u,END=10) n_tmp
       if ( n_tmp == n ) then
          read(u) z1
          read(u) z2
          exit
       end if
       read(u)
       read(u)
    end do
10  continue
    close(u)
  END SUBROUTINE read_file


  SUBROUTINE write_file( nfac, n, z1, z2 )
    implicit none
    integer,intent(IN) :: nfac, n
    complex(8),intent(IN) :: z1(:,:,:),z2(:,:,:)
    integer,parameter :: u1=82, u2=83
    character(8),parameter :: file_n1="file_chk"
    character(8),parameter :: file_n2="file_dat"
    integer,parameter :: ndata_max=10000
    integer :: nfac_tmp, i
    logical :: flag

    inquire( FILE=file_n1, EXIST=flag )

    if ( flag ) then
       open(u1,file=file_n1,status="old")
       read(u1,*) nfac_tmp
       if ( nfac_tmp /= nfac ) flag=.false.
       close(u1)
    end if

    if ( flag ) then

       open(u1,file=file_n1,status="old",position="append")
       write(u1,*) n
       close(u1)

       open(u2,file=file_n2,status="old",form="unformatted",position="append")
       write(u2) n
       write(u2) z1
       write(u2) z2
       close(u2)

    else

       open(u1,file=file_n1)
       write(u1,*) nfac
       write(u1,*) n
       close(u1)

       open(u2,file=file_n2,form="unformatted")
       write(u2) nfac
       write(u2) n
       write(u2) z1
       write(u2) z2
       close(u2)

    end if
  END SUBROUTINE write_file


  SUBROUTINE calc_bb(aa,bb_out)
    implicit none
    real(8),intent(IN)  :: aa(3,3)
    real(8),intent(OUT) :: bb_out(3,3)
    real(8) :: Vaa,PI2
    PI2 = 2.d0*acos(-1.d0)
    Vaa = aa(1,1)*aa(2,2)*aa(3,3)+aa(1,2)*aa(2,3)*aa(3,1) &
         +aa(1,3)*aa(2,1)*aa(3,2)-aa(1,3)*aa(2,2)*aa(3,1) &
         -aa(1,2)*aa(2,1)*aa(3,3)-aa(1,1)*aa(2,3)*aa(3,2)
    bb_out(1,1) = aa(2,2)*aa(3,3) - aa(3,2)*aa(2,3)
    bb_out(2,1) = aa(3,2)*aa(1,3) - aa(1,2)*aa(3,3)
    bb_out(3,1) = aa(1,2)*aa(2,3) - aa(2,2)*aa(1,3)
    bb_out(1,2) = aa(2,3)*aa(3,1) - aa(3,3)*aa(2,1)
    bb_out(2,2) = aa(3,3)*aa(1,1) - aa(1,3)*aa(3,1)
    bb_out(3,2) = aa(1,3)*aa(2,1) - aa(2,3)*aa(1,1)
    bb_out(1,3) = aa(2,1)*aa(3,2) - aa(3,1)*aa(2,2)
    bb_out(2,3) = aa(3,1)*aa(1,2) - aa(1,1)*aa(3,2)
    bb_out(3,3) = aa(1,1)*aa(2,2) - aa(2,1)*aa(1,2)
    bb_out(:,:)=bb_out(:,:)*PI2/Vaa
  END SUBROUTINE calc_bb


END PROGRAM fermi_surface_ufld
