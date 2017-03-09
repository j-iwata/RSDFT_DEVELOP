PROGRAM ldos2gp

  !use libtetrabz, only: libtetrabz_dos
  use libtetrabz_mpi, only: libtetrabz_mpi_dos

  implicit none

  integer,parameter :: u1=1, u2=2, u3=3, u5=5, u99=99, u10=10
  integer,parameter :: max_loop=10000000
  integer :: ltetra
  integer :: k,i,k1,k2,k3,i1,i2,i3,n,ne,s,j
  integer :: nk, n_irwedge, n_whole, indx_range(2,3)
  integer :: nge(3),ngw(3)
  integer :: nband,nspin,nkpnt
  integer :: myrank
  integer,allocatable :: kbz(:,:),ktmp(:,:,:),kgrd(:,:,:)
  integer,allocatable :: itmp(:,:)
  real(8),allocatable :: eval(:,:,:),occp(:,:,:)
  real(8),allocatable :: wbz(:),kpnt(:,:),weight(:)
  real(8),allocatable :: eig(:,:,:,:),wght(:,:,:,:,:),e0(:)
  real(8),parameter :: HT=27.2116d0
  real(8) :: emin,emax,de,e1,e2
  real(8) :: bb(3,3), ktry(3), err, etime(0:1)
  character(80) :: cbuf
  character(16) :: file_bz, file_eigv
  character(8) :: cbuf3(3)
  include 'mpif.h'
  integer :: ngrid(0:3),mb,mb1,mb2,mk,ms,mmbz,jtmp(14:21)
  integer :: io_ctrl, oc, type_wf, icheck_format
  integer,allocatable :: LL(:,:)
  real(8),allocatable :: occ(:,:,:),kbb(:,:)
  real(8) :: aa(3,3), dummy33(3,3), w1,w2,w3
  real(8),allocatable :: dtmp(:)
  real(8),allocatable :: rho_x(:,:,:),rho_y(:,:,:),rho_z(:,:,:)
  complex(8),allocatable :: ztmp(:), wf(:,:,:)
  complex(8),parameter :: z0=(0.0d0,0.0d0)
  real(8) :: hgrid(3),dV,Vaa,dVYZ,dVZX,dVXY

  call MPI_INIT(i)
  call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, i )

  etime(0)=mpi_wtime()

  file_bz   = "bz_info"
  file_eigv = "eigenvalues"

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

  read(u1,'(a)') cbuf
  call get_nums_from_chr( cbuf, indx_range(:,1) )
  read(u1,'(a)') cbuf
  call get_nums_from_chr( cbuf, indx_range(:,2) )
  read(u1,'(a)') cbuf
  call get_nums_from_chr( cbuf, indx_range(:,3) )

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

! ---

  if ( myrank == 0 ) then
     write(*,'(1x,"indx_range:",3(2i4,2x))') indx_range(:,:)
     write(*,'(1x,"nge       :",3i4)') nge
  end if

! ---

  allocate( kgrd(nge(1),nge(2),nge(3)) ) ; kgrd=0

  do i3=1,nge(3)
  do i2=1,nge(2)
  do i1=1,nge(1)
     kgrd(i1,i2,i3) = ktmp(itmp(i1,1),itmp(i2,2),itmp(i3,3))
  end do
  end do
  end do

! ---

  if ( myrank == 0 ) then

  open(u3,file="wf.dat1",status="old",form="unformatted")

  read(u3) icheck_format
  if ( icheck_format > 0 ) stop "old format of wf.dat1 is not supported!"

  read(u3) ngrid(0:3)
  read(u3) mb, mb1, mb2
  read(u3) mk, ms, mmbz
  read(u3) io_ctrl, oc, type_wf

  allocate( LL(3,ngrid(0)) ) ; LL=0
  allocate( occ(mb,mk,ms)  ) ; occ=0.0d0
  allocate( kbb(3,mk)      ) ; kbb=0.0d0

  read(u3) LL
  read(u3) occ
  read(u3) aa, dummy33, kbb

  deallocate( kbb )
  deallocate( occ )

  if ( type_wf == 0 ) then
     allocate( ztmp(ngrid(0)) ) ; ztmp=z0
  else
     allocate( dtmp(ngrid(0)) ) ; dtmp=0.0d0
  end if

  allocate( wf(0:ngrid(1)-1,0:ngrid(2)-1,0:ngrid(3)-1) ) ; wf=z0

  allocate( rho_x(0:ngrid(1)-1,mb1:mb2,mk) ) ; rho_x=0.0d0
  allocate( rho_y(0:ngrid(2)-1,mb1:mb2,mk) ) ; rho_y=0.0d0
  allocate( rho_z(0:ngrid(3)-1,mb1:mb2,mk) ) ; rho_z=0.0d0
  read(u3)
  read(u3)
  read(u3)
  read(u3)
  read(u3)

  do s=1,ms
  do k=1,mk
  do n=mb1,mb2

     if ( type_wf == 0 ) then
        read(u3) ztmp
        do i=1,ngrid(0)
           wf( LL(1,i),LL(2,i),LL(3,i) ) = ztmp(i)
        end do
     else
        read(u3) dtmp
        do i=1,ngrid(0)
           wf( LL(1,i),LL(2,i),LL(3,i) ) = dtmp(i)
        end do
     end if

     do i3=0,ngrid(3)-1
     do i2=0,ngrid(2)-1
     do i1=0,ngrid(1)-1
        rho_x(i1,n,k) = rho_x(i1,n,k) + abs( wf(i1,i2,i3) )**2
        rho_y(i2,n,k) = rho_y(i2,n,k) + abs( wf(i1,i2,i3) )**2
        rho_z(i3,n,k) = rho_z(i3,n,k) + abs( wf(i1,i2,i3) )**2
     end do
     end do
     end do

  end do
  end do
  end do

  deallocate( wf )
  if ( allocated(ztmp) ) deallocate( ztmp )
  if ( allocated(dtmp) ) deallocate( dtmp )
  deallocate( LL )

  close(u3)

  end if

! ---

  call MPI_BCAST( ngrid, 4, MPI_INTEGER, 0, MPI_COMM_WORLD, i )
  call MPI_BCAST( mb1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, i )
  call MPI_BCAST( mb2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, i )
  call MPI_BCAST( mk , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, i )

  if ( myrank /= 0 ) then
     allocate( rho_x(0:ngrid(1)-1,mb1:mb2,mk) ) ; rho_x=0.0d0
     allocate( rho_y(0:ngrid(2)-1,mb1:mb2,mk) ) ; rho_y=0.0d0
     allocate( rho_z(0:ngrid(3)-1,mb1:mb2,mk) ) ; rho_z=0.0d0
  end if

  call MPI_BCAST( rho_x, size(rho_x), MPI_REAL8, 0, MPI_COMM_WORLD, i )
  call MPI_BCAST( rho_y, size(rho_y), MPI_REAL8, 0, MPI_COMM_WORLD, i )
  call MPI_BCAST( rho_z, size(rho_z), MPI_REAL8, 0, MPI_COMM_WORLD, i )

  call MPI_BCAST( aa, 9, MPI_REAL8, 0, MPI_COMM_WORLD, i )

! ---

  Vaa = aa(1,1)*aa(2,2)*aa(3,3)+aa(1,2)*aa(2,3)*aa(3,1) &
       +aa(1,3)*aa(2,1)*aa(3,2)-aa(1,3)*aa(2,2)*aa(3,1) &
       -aa(1,2)*aa(2,1)*aa(3,3)-aa(1,1)*aa(2,3)*aa(3,2)

  dV = Vaa/ngrid(0)

  do i=1,3
     hgrid(i) = sqrt( sum(aa(:,i)**2) )/ngrid(i)
  end do

  dVYZ = dV/hgrid(1)
  dVZX = dV/hgrid(2)
  dVXY = dV/hgrid(3)

! ---

  open(u2,file=file_eigv)

  nkpnt=0

  do i=1,max_loop
     read(u2,*,END=10) cbuf
     if ( cbuf == "Nband" ) then
        backspace(u2)
        read(u2,*) cbuf3, nband, nspin, k
        nkpnt=nkpnt+1
     end if
  end do
10 continue

  allocate( eval(nband,nkpnt,nspin) ) ; eval=0.0d0
  allocate( kpnt(3,nkpnt) ) ; kpnt=0.0d0
  allocate( weight(nkpnt) ) ; weight=0.0d0

  rewind u2

  do k=1,nkpnt
     read(u2,*)
     read(u2,*) kpnt(1:3,k),weight(k)
     do n=1,nband
        read(u2,*) i,( eval(n,k,s), s=1,nspin )
     end do
  end do

  close(u2)

! ---

  allocate( eig(nband,nge(1),nge(2),nge(3)) ) ; eig=0.0d0

  do i3=1,nge(3)
  do i2=1,nge(2)
  do i1=1,nge(1)

     i = kgrd(i1,i2,i3)

     ktry(1:3) = dble(kbz(1:3,i))/dble(nk)

     do k=1,nkpnt
        err = sum( (ktry-kpnt(:,k))**2 )
        if ( err < 1.d-10 ) then
           do n=1,nband
              eig(n,i1,i2,i3) = eval(n,k,1)
           end do
           exit
        end if
     end do ! k
     if ( .not.(err < 1.d-10) ) stop "error(1)"

  end do ! i1
  end do ! i2
  end do ! i3

! ---

  write(*,*) "count(eig/=0.0d0):",count(eig/=0.0d0),size(eig)

!------------------------------

  emin = minval( eval )
  emax = maxval( eval )

  if ( myrank == 0 ) then
     write(*,*) "emin     (HT,eV)=",emin,emin*HT
     write(*,*) "emax     (HT,eV)=",emax,emax*HT
     write(*,*) "emax-emin(HT,eV)=",emax-emin,(emax-emin)*HT
  end if

!------------------------------

  if ( myrank == 0 ) then
     write(*,*) "# of energy grid: ne="
     read(u5,*) ne
     write(*,*) ne
  end if
  call MPI_BCAST( ne, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, i )

  if ( myrank == 0 ) then
     write(*,*) "energy range: e1,e2= ( 0,0: default values are set )"
     read(u5,*) e1,e2
  end if
  call MPI_BCAST( e1, 1, MPI_REAL8, 0, MPI_COMM_WORLD, i )
  call MPI_BCAST( e2, 1, MPI_REAL8, 0, MPI_COMM_WORLD, i )
  if ( e1 == 0.0d0 ) e1=emin
  if ( e2 == 0.0d0 ) e2=emax
  if ( myrank == 0 ) write(*,*) e1,e2

  if ( myrank == 0 ) then
     write(*,*) "ltetra [1:linear tetrahedron, 2:opt-linear tetrahedron]"
     read(u5,*) ltetra
     write(*,*) ltetra
  end if
  call MPI_BCAST( ltetra, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, i )

! ---

  ngw(1:3) = nge(1:3)

  allocate( wght(ne,nband,ngw(1),ngw(2),ngw(3)) ) ; wght=0.0d0
  allocate( e0(ne) ) ; e0=0.0d0


  de=(emax-emin)/ne
  do i=1,ne
     e0(i) = emin + (i-1)*de
  end do

  !call libtetrabz_dos( ltetra, bb, nband, nge, eig, ngw, wght, ne, e0 )
  call libtetrabz_mpi_dos &
       ( ltetra, MPI_COMM_WORLD, bb, nband, nge, eig, ngw, wght, ne, e0 )

  if ( myrank == 0 ) then

     rewind 10
     do i=1,ne
        write(10,*) e0(i), sum( wght(i,:,:,:,:) )
     end do

     rewind 11
     rewind 12
     rewind 13
     do i=1,ne

        do j=0,ngrid(1)-1
           w1=0.0d0
           do i3=1,nge(3)
           do i2=1,nge(2)
           do i1=1,nge(1)
              k=kgrd(i1,i2,i3)
              do n=mb1,mb2
                 w1=w1+wght(i,n,i1,i2,i3)*rho_x(j,n,k)*dVYZ
              end do
           end do
           end do
           end do
           write(11,*) e0(i), j*hgrid(1), w1
        end do

        do j=0,ngrid(2)-1
           w2=0.0d0
           do i3=1,nge(3)
           do i2=1,nge(2)
           do i1=1,nge(1)
              k=kgrd(i1,i2,i3)
              do n=mb1,mb2
                 w2=w2+wght(i,n,i1,i2,i3)*rho_y(j,n,k)*dVZX
              end do
           end do
           end do
           end do
           write(12,*) e0(i), j*hgrid(2), w2
        end do

        do j=0,ngrid(3)-1
           w3=0.0d0
           do i3=1,nge(3)
           do i2=1,nge(2)
           do i1=1,nge(1)
              k=kgrd(i1,i2,i3)
              do n=mb1,mb2
                 w3=w3+wght(i,n,i1,i2,i3)*rho_z(j,n,k)*dVXY
              end do
           end do
           end do
           end do
           write(13,*) e0(i), j*hgrid(3), w3
        end do

     end do ! i

  end if

  etime(1)=mpi_wtime()
  if ( myrank == 0 ) write(*,*) "TIME=",etime(1)-etime(0)

  call MPI_FINALIZE( i )

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


END PROGRAM ldos2gp
