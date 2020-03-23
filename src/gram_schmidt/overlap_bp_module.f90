module overlap_bp_module

  use sl_tools_module, only: slinfo

  implicit none

  private
  public :: init_overlap_bp
  public :: calc_overlap_bp

  interface calc_overlap_bp
     module procedure d_calc_overlap_bp
  end interface

  integer :: comm_g, nprocs_g, myrank_g
  integer :: comm_b, nprocs_b, myrank_b
  real(8) :: dV
  integer,allocatable :: ir_bp(:), id_bp(:)
  integer :: nband

  logical :: has_init_done = .false.

contains

  subroutine init_overlap_bp( ng, nb, alpha, comm1, comm2 )
    implicit none
    integer,intent(in) :: ng, nb
    real(8),intent(in) :: alpha
    integer,intent(in) :: comm1, comm2
    integer :: ierr
    include 'mpif.h'
    if ( has_init_done ) return
!    call write_border( 1, " init_overlap_bp(start)" )
    dV=alpha
    comm_g=comm1
    comm_b=comm2
    call MPI_Comm_size( comm_g, nprocs_g, ierr )
    call MPI_Comm_rank( comm_g, myrank_g, ierr )
    call MPI_Comm_size( comm_b, nprocs_b, ierr )
    call MPI_Comm_rank( comm_b, myrank_b, ierr )
    allocate( ir_bp(0:nprocs_b-1) ); ir_bp=0
    allocate( id_bp(0:nprocs_b-1) ); id_bp=0
    call MPI_Allreduce(nb,nband,1,MPI_INTEGER,MPI_SUM,comm_b,ierr)
    call load_div( nband, nprocs_b, ir_bp, id_bp )
    has_init_done = .true.
!    call write_border( 1, " init_overlap_bp(end)" )
  end subroutine init_overlap_bp


  subroutine d_calc_overlap_bp( Ssub, u, sl, UPLO )

    implicit none
    real(8),intent(out)     :: Ssub(:,:)
    real(8),intent(inout)   :: u(:,:)
    type(slinfo),intent(in) :: sl
    character(1),intent(in) :: UPLO
    real(8),parameter :: zero=0.0d0
    real(8),allocatable :: b(:,:),w(:,:)
    integer,allocatable :: irank_s(:), irank_r(:)
    integer :: ir,ic,i,j,i0,j0,istep,nstep,is,js,i1,j1,k
    integer :: m,n,mb0_0,mb1_0,nb0,nb1,np,np1,ip,jp,inp
    integer :: itag,ierr,ireq(2),comm1,comm2,nprocs
    integer :: is_min,is_max,ir_min,ir_max,me1,me2
    include 'mpif.h'
    integer :: istatus(MPI_STATUS_SIZE,2),myrank
    real(8),allocatable :: Sgp(:,:),Stmp(:,:),S(:,:)
    real(8),allocatable :: buff(:)
    integer,allocatable :: icheck_s(:,:,:,:),icheck_r(:,:,:,:),irank(:,:)
    integer,allocatable :: ir_gp(:),id_gp(:)

!    call write_border( 1, "d_calc_overlap_bp(start)" )

    m     = size( u, 1 )    ! grid size ( local )
    n     = size( u, 2 )    ! band size ( local )
    np    = nprocs_b        ! # of processes for band-parallel
    np1   = nprocs_g        ! # of processes for grid-parallel
    comm1 = comm_g
    comm2 = comm_b
    nstep = np/2

    allocate( w(m,n) ); w=zero
    allocate( b(m,n) ); b=zero

! ---

    allocate( irank_s(nstep) ); irank_s=0
    allocate( irank_r(nstep) ); irank_r=0

    if ( UPLO == "L" .or. UPLO == "l" ) then
       do istep=1,nstep
          irank_s(istep) = mod( myrank_b-istep+np, np )
          irank_r(istep) = mod( myrank_b+istep   , np )
       end do
    else if ( UPLO == "U" .or. UPLO == "u" ) then
       do istep=1,nstep
          irank_r(istep) = mod( myrank_b-istep+np, np )
          irank_s(istep) = mod( myrank_b+istep   , np )
       end do
    end if
!    write(*,'(1x,"s",i2,2x,10i2)') myrank_b,irank_s
!    write(*,'(1x,"r",i2,2x,10i2)') myrank_b,irank_r

! ---

    nb0 = id_bp(myrank_b) + 1
    nb1 = nb0 + ir_bp(myrank_b) - 1

    allocate( ir_gp(0:nprocs_g-1) ); ir_gp=0
    allocate( id_gp(0:nprocs_g-1) ); id_gp=0

    call load_div( nband, nprocs_g, ir_gp, id_gp )

    mb0_0 = id_gp(myrank_g) + 1
    mb1_0 = mb0_0 + ir_gp(myrank_g) - 1

    call MPI_Comm_size( MPI_COMM_WORLD, nprocs, ierr )
    call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr )
!    write(*,*) "mb0_0,mb0_1",mb0_0,mb1_0,myrank
!    write(*,*) "size(Ssub)",size(Ssub,1),size(Ssub,2)
!    call stop_program("")

    allocate( Stmp(n,n) ); Stmp=0.0d0
    allocate( Sgp(mb0_0:mb1_0,nb0:nb1) ); Sgp=0.0d0

! ---

    itag=0

    ic=id_bp(myrank_b)+1
    ir=ic

    b(:,1:n) = u(:,1:n)

    do istep=0,nstep

       w(:,:) = b(:,:)  ! --> same size as u(:,:) the wave function

       if ( istep < nstep ) then
          call MPI_Irecv( b,size(b),MPI_REAL8,irank_r(1),itag,comm2,ireq(1),ierr )
          call MPI_Isend( w,size(w),MPI_REAL8,irank_s(1),itag,comm2,ireq(2),ierr )
       end if

       if ( istep == 0 ) then
          call DSYRK( UPLO, 'T', n, m, dV, u, m, zero, Stmp, size(Stmp,1) )
       else
          call DGEMM( 'T', 'N', n, n, m, dV, w, m, u, m, zero, Stmp, size(Stmp,1) )
       end if

       call MPI_Allreduce(MPI_IN_PLACE,Stmp,size(Stmp),MPI_REAL8,MPI_SUM,comm1,ierr) 

       do j=1,n
       do i=1,n
          i0=ir +i-1
          j0=nb0+j-1
          if ( mb0_0 <= i0 .and. i0 <= mb1_0 ) Sgp(i0,j0)=Stmp(i,j)
       end do
       end do

       if ( UPLO == "L" .or. UPLO == "l" ) then
          ir=ir+n ; if ( ir > nband ) ir=ir-nband
       else if ( UPLO == "U" .or. UPLO == "u" ) then
          ir=ir-n ; if ( ir < 1 ) ir=ir+nband
       end if

       if ( istep < nstep ) call mpi_waitall( 2, ireq, istatus, ierr )

    end do ! istep

    deallocate( Stmp )
    deallocate( irank_r )
    deallocate( irank_s )
    deallocate( b )
    deallocate( w )

! ---

    ir=0
    is=0

!    do j=nb0  ,nb1
!    do i=mb0_0,mb1_0
!       if ( Sgp(i,j) == 0.0d0 ) then
!          if ( i <= j ) ir=ir+1
!       else
!          if ( i > j ) is=is+1
!       end if
!    end do
!    end do

    ir_min=nband
    ir_max=1
    is_min=nband
    is_max=1
    do i=mb0_0,mb1_0
       if ( i < nb0 ) then
          if ( any(Sgp(i,:)==0.0d0) ) then
             ir=ir+1
             ir_min=min(ir_min,i)
             ir_max=max(ir_max,i)
          end if
       end if
       if ( i > nb1 ) then
          if ( any(Sgp(i,:)/=0.0d0) ) then
             is=is+1
             is_min=min(is_min,i)
             is_max=max(is_max,i)
          end if
       end if
    end do

    if ( ir == 0 ) then
       ir_min=0
       ir_max=0
    end if
    if ( is == 0 ) then
       is_min=0
       is_max=0
    end if

!    write(*,'(1x,i4,3x,"(ir)",3i6,3x,"(is)",3i6)') &
!         myrank,ir,ir_min,ir_max,is,is_min,is_max

    allocate( icheck_s(0:np1-1,0:np-1,4,0:nprocs-1) ); icheck_s=0
    allocate( icheck_r(0:np1-1,0:np-1,4,0:nprocs-1) ); icheck_r=0

    if ( is_min > 1 ) then
       do jp=0,nprocs_b-1
          j0=id_bp(jp)+1
          j1=id_bp(jp)+ir_bp(jp)-1
       do ip=0,nprocs_g-1
          i0=id_gp(ip)+1
          i1=id_gp(ip)+ir_gp(ip)-1
          do is=is_min,is_max
             if ( j0 <= is .and. is <= j1 ) then
                do js=nb0,nb1
                   if ( i0 <= js .and. js <= i1 ) then
!                      write(*,'(1x,2i4,2x,3i4,2x,3i4,2x,i4)') &
!                           is,js,ip,i0,i1,jp,j0,j1,myrank
                      if ( icheck_s(ip,jp,1,myrank) == 0 ) then
                         icheck_s(ip,jp,1,myrank)=is
                         icheck_s(ip,jp,3,myrank)=js
                      end if
                      icheck_s(ip,jp,2,myrank)=is
                      icheck_s(ip,jp,4,myrank)=js
                   end if
                end do
             end if
          end do
       end do
       end do
    end if

    call MPI_Allreduce( icheck_s, icheck_r, size(icheck_r), MPI_INTEGER &
         , MPI_SUM, MPI_COMM_WORLD, ierr )

!    do jp=0,np-1
!       do ip=0,np1-1
!          if ( icheck_s(ip,jp,1,myrank) /= 0 ) then
!             write(*,'(1x,"s",9i4)') ip,jp,icheck_s(ip,jp,:,myrank),myrank
!          end if
!       end do
!    end do

!    ip=myrank_g
!    jp=myrank_b
!    do i=0,nprocs-1
!       if ( icheck_r(ip,jp,1,i) /= 0 ) then
!          write(*,'(1x,"r",9i4)') ip,jp,icheck_r(ip,jp,:,i),i,myrank
!       end if
!    end do

!    write(*,*) count(Sgp/=0.0d0),myrank
!    call stop_program("")

    me1=myrank_g
    me2=myrank_b

    allocate( irank(0:np1-1,0:np-1) ); irank=0
    irank(me1,me2)=myrank
    call MPI_Allreduce(MPI_IN_PLACE,irank,size(irank),MPI_INTEGER,MPI_SUM &
         ,MPI_COMM_WORLD,ierr)

    allocate( buff(size(Sgp)) ); buff=0.0d0

    do jp=0,np-1
    do ip=0,np1-1
       do inp=0,nprocs-1
          k=0
          if ( inp==myrank .and. icheck_s(ip,jp,1,inp)/=0 ) then
             i0=icheck_s(ip,jp,1,inp)
             i1=icheck_s(ip,jp,2,inp)
             j0=icheck_s(ip,jp,3,inp)
             j1=icheck_s(ip,jp,4,inp)
             k=0
             do j=j0,j1
             do i=i0,i1
                k=k+1
                buff(k)=Sgp(i,j)
             end do
             end do
             call MPI_Isend(buff,k,MPI_REAL8,irank(ip,jp),1,MPI_COMM_WORLD,ireq(1),ierr)
!             write(*,*) "isend",ip,jp,myrank,irank(ip,jp)
          end if
          if ( ip==me1 .and. jp==me2 .and. icheck_r(ip,jp,1,inp)/=0 ) then
             i0=icheck_r(ip,jp,1,inp)
             i1=icheck_r(ip,jp,2,inp)
             j0=icheck_r(ip,jp,3,inp)
             j1=icheck_r(ip,jp,4,inp)
             k=(i1-i0+1)*(j1-j0+1)
             call MPI_Irecv(buff,k,MPI_REAL8,inp,1,MPI_COMM_WORLD,ireq(1),ierr)
!             write(*,*) "irecv",ip,jp,inp,myrank
          end if
          if ( k > 0 ) call MPI_Wait(ireq(1),istatus,ierr)
          if ( ip==me1 .and. jp==me2 .and. icheck_r(ip,jp,1,inp)/=0 ) then
             i0=icheck_r(ip,jp,1,inp)
             i1=icheck_r(ip,jp,2,inp)
             j0=icheck_r(ip,jp,3,inp)
             j1=icheck_r(ip,jp,4,inp)
!             write(*,'("aaa",16i4)') j0,j1,i0,i1,myrank,count(Sgp(j0:j1,i0:i1)/=0.0d0),count(b(j0:j1,i0:i1)/=0.0d0),count(b/=0.0d0),mb0_0,mb1_0,nb0,nb1
!             do i=i0,i1
!             do j=j0,j1
!                k=1+(i-i0)+(i1-i0+1)*(j-j0)
!                Sgp(j,i)=buff(k)
!             end do
!             end do
             k=0
             do j=j0,j1
             do i=i0,i1
                k=k+1
                Sgp(j,i)=buff(k)
             end do
             end do
          end if
       end do
    end do
    end do

!    write(*,*) count(Sgp/=0.0d0),myrank
!    call stop_program("")

    deallocate( buff )
    deallocate( irank )
    deallocate( icheck_r )
    deallocate( icheck_s )

! ---

#ifdef test

    allocate( S(nband,nband) ) ; S=zero
    call mochikae( UPLO,comm1,comm2,nband,mb0_0,mb1_0,nb0,nb1,Sgp,S )
    call d_distribute_matrix( sl, S, Ssub )
    deallocate( S )

#else

!    allocate( S(nband,nband) ) ; S=zero
!    call mochikae( UPLO,comm1,comm2,nband,mb0_0,mb1_0,nb0,nb1,Sgp,S )
!    deallocate( S )

    !call d_distribute_matrix_2( sl, nband,nband, mb0_0,mb1_0,nb0,nb1, Sgp, Ssub )
    call d_distribute_matrix_3( sl, nband,nband, mb0_0,mb1_0,nb0,nb1, Sgp, Ssub )

#endif

    deallocate( Sgp )

!    call write_border( 1, "d_calc_overlap_bp(end)" )

  end subroutine d_calc_overlap_bp


  subroutine load_div( m, np, ir_out, id_out )
    implicit none
    integer,intent(in) :: m, np
    integer,intent(out) :: ir_out(0:np-1),id_out(0:np-1)
    integer :: i,ip
    ir_out=0.0d0
    do i=1,m
       ip=mod(i-1,np)
       ir_out(ip)=ir_out(ip)+1
    end do
    do ip=0,np-1
       id_out(ip)=sum(ir_out(0:ip))-ir_out(ip)
    end do
  end subroutine load_div


  SUBROUTINE d_distribute_matrix( sl, H, Hsub )
    implicit none
    type(slinfo),intent(IN) :: sl
    real(8),intent(IN) :: H(:,:)
    real(8),intent(OUT) :: Hsub(:,:)
    integer :: m,n,i0,j0,i1,j1,ip,jp,is,js,i,j,ii,jj
!    call write_border( 1, "d_distribute_matrix(start)" )
    m = size( H, 1 )
    n = size( H, 2 )
    jp=-1
    js= 0
    do j0=1,n,sl%nbsize
       j1=min(j0+sl%nbsize-1,n)
       jj=j1-j0+1
       jp=mod(jp+1,sl%npcol) ; if ( jp /= sl%mycol ) cycle
       ip=-1
       is= 0
       do i0=1,m,sl%mbsize
          i1=min(i0+sl%mbsize-1,m)
          ii=i1-i0+1
          ip=mod(ip+1,sl%nprow) ; if ( ip /= sl%myrow ) cycle
          do j=1,jj
          do i=1,ii
             Hsub(is+i,js+j) = H(i0+i-1,j0+j-1)
          end do
          end do
          is = is + ii
       end do
       js = js + jj
    end do
!    call write_border( 1, "d_distribute_matrix(end)" )
  END SUBROUTINE d_distribute_matrix


  SUBROUTINE d_distribute_matrix_2( sl, m,n, mb0,mb1,nb0,nb1, H, Hsub )

    implicit none
    type(slinfo),intent(IN) :: sl
    integer,intent(IN) :: m,n,mb0,mb1,nb0,nb1
    real(8),intent(IN) :: H(mb0:mb1,nb0:nb1)
    real(8),intent(OUT) :: Hsub(:,:)
    integer :: i,i0,i1,ii,ip,is,ib,j,j0,j1,jj,jp,js,jb,ierr,myrank
    integer,allocatable :: is_(:,:),js_(:,:)
    real(8),allocatable :: work(:,:),buff(:,:)
    include 'mpif.h'

!    call write_border( 1, "d_distribute_matrix_2(start)" )

    Hsub(:,:)=0.0d0

    allocate( work(size(Hsub,1),size(Hsub,2)) ) ; work=0.0d0
    allocate( buff(size(Hsub,1),size(Hsub,2)) ) ; buff=0.0d0
    allocate( is_(0:sl%nprow-1,0:sl%npcol-1) ) ; is_=0
    allocate( js_(0:sl%nprow-1,0:sl%npcol-1) ) ; js_=0

    jp=-1
    js_=0
    call MPI_COMM_RANK(mpi_comm_world,myrank,ierr)

    do j0=1,n,sl%nbsize

       j1=min(j0+sl%nbsize-1,n)
       jj=j1-j0+1
       jp=mod(jp+1,sl%npcol)

       ip=-1
       is_(:,jp)=0

       do i0=1,m,sl%mbsize

          i1=min(i0+sl%mbsize-1,m)
          ii=i1-i0+1
          ip=mod(ip+1,sl%nprow)

          is = is_(ip,jp)
          js = js_(ip,jp)

          work(:,:)=0.0d0
          do j=1,jj
          do i=1,ii
             ib=i0+i-1
             jb=j0+j-1
             if ( mb0 <= ib .and. ib <= mb1 .and. nb0 <= jb .and. jb <= nb1 ) then
                work(is+i,js+j) = work(is+i,js+j) + H(ib,jb)
             end if
             if ( ib /= jb .and. mb0 <= jb .and. jb <= mb1 .and. nb0 <= ib .and. ib <= nb1 ) then
                work(is+i,js+j) = work(is+i,js+j) + H(jb,ib)
             end if
          end do
          end do

!          do i=0,size(sl%map_1to2,2)-1
!             if ( ip == sl%map_1to2(1,i) .and. jp == sl%map_1to2(2,i) ) exit
!          end do
          i = sl%map_2to1(ip,jp)
          call MPI_Reduce( work, buff, size(buff), MPI_REAL8, MPI_SUM, i, MPI_COMM_WORLD, ierr )
          if ( ip == sl%myrow .and. jp == sl%mycol ) then
             Hsub = Hsub + buff
          end if

          is_(ip,jp) = is_(ip,jp) + ii

       end do

       js_(:,jp) = js_(:,jp) + jj

    end do

    deallocate( js_,is_ )
    deallocate( buff )
    deallocate( work )

!    call write_border( 1, "d_distribute_matrix_2(end)" )

  END SUBROUTINE d_distribute_matrix_2


  SUBROUTINE d_distribute_matrix_3( sl, m,n, mb0,mb1,nb0,nb1, H, Hsub )
    implicit none
    type(slinfo),intent(IN) :: sl
    integer,intent(IN) :: m,n,mb0,mb1,nb0,nb1
    real(8),intent(IN) :: H(mb0:mb1,nb0:nb1)
    real(8),intent(OUT) :: Hsub(:,:)
    integer :: ictxt,ia,ja
!    call write_border( 1, "d_distribute_matrix_3(start)" )
    call pdgemr2d(n,n,H,1,1,sl%descb,Hsub,1,1,sl%desca,sl%icontxt_b)
!    call pdgemr2d(n,n,Hsub,1,1,sl%desca,Hsub,1,1,sl%desca,ictxt)
!    call pdgemr2d(n,n,H,1,1,sl%descb,H,1,1,sl%descb,sl%icontxt_b)
!    call write_border( 1, "d_distribute_matrix_3(end)" )
  END SUBROUTINE d_distribute_matrix_3


  SUBROUTINE mochikae( UPLO,comm1,comm2,nband,mb0_0,mb1_0,nb0,nb1,Sgp,S )
    implicit none
    character(1),intent(IN) :: UPLO
    integer,intent(IN) :: comm1,comm2,nband,mb0_0,mb1_0,nb0,nb1
    real(8),intent(INOUT) :: Sgp(mb0_0:mb1_0,nb0:nb1)
    real(8),intent(INOUT) :: S(nband,nband)
    include 'mpif.h'
    integer :: ierr,ic,ir
!    call write_border( 1, "mochikae(start)" )
    S(mb0_0:mb1_0,nb0:nb1) = Sgp(mb0_0:mb1_0,nb0:nb1)
    call MPI_ALLREDUCE( MPI_IN_PLACE,S,size(S),MPI_REAL8,MPI_SUM,comm1,ierr )
    if ( UPLO == "L" .or. UPLO == "l" ) then
       do ic=1,nband
       do ir=ic+1,nband
          if ( S(ic,ir) /= 0.0d0 ) then
             S(ir,ic)=S(ic,ir)
             S(ic,ir)=0.0d0
          end if
       end do
       end do
    else if ( UPLO == "U" .or. UPLO == "u" ) then
       do ic=1,nband
       do ir=1,ic-1
          if ( S(ic,ir) /= 0.0d0 ) then
             S(ir,ic)=S(ic,ir)
             S(ic,ir)=0.0d0
          end if
       end do
       end do
    end if
    call MPI_ALLREDUCE(MPI_IN_PLACE,S,size(S),MPI_REAL8,MPI_SUM,comm2,ierr)
    Sgp(mb0_0:mb1_0,nb0:nb1) = S(mb0_0:mb1_0,nb0:nb1)
!    call write_border( 1, "mochikae(end)" )
  END SUBROUTINE mochikae


end module overlap_bp_module
