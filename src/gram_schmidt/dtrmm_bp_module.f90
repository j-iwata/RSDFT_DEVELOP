module dtrmm_bp_module

  use sl_tools_module, only: slinfo

  implicit none

  private
  public :: init_dtrmm_bp
  public :: dtrmm_bp

  integer :: comm_g, nprocs_g, myrank_g
  integer :: comm_b, nprocs_b, myrank_b
  integer,allocatable :: ir_bp(:), id_bp(:)
  logical :: has_init_done = .false.

contains

  subroutine init_dtrmm_bp( ng, nb, comm1, comm2 )
    implicit none
    integer,intent(in) :: ng, nb
    integer,intent(in) :: comm1, comm2
    integer :: nband,ierr
    include 'mpif.h'
    if ( has_init_done ) return
!    call write_border( 1, " init_dtrmm_bp(start)" )
    comm_g=comm1
    comm_b=comm2
    call MPI_Comm_size( comm_g, nprocs_g, ierr )
    call MPI_Comm_rank( comm_g, myrank_g, ierr )
    call MPI_Comm_size( comm_b, nprocs_b, ierr )
    call MPI_Comm_rank( comm_b, myrank_b, ierr )
    call MPI_Allreduce(nb,nband,1,MPI_INTEGER,MPI_SUM,comm_b,ierr)
    call load_div( nband, nprocs_b, ir_out=ir_bp, id_out=id_bp )
    has_init_done = .true.
!    call write_border( 1, " init_dtrmm_bp(end)" )
  end subroutine init_dtrmm_bp


  subroutine dtrmm_bp( UPLO, sl, u, Ssub )

    implicit none
    character(1),intent(in) :: UPLO
    type(slinfo),intent(in) :: sl
    real(8),intent(inout)   :: u(:,:)
    real(8),intent(in)      :: Ssub(:,:)
    real(8),allocatable :: S(:,:)
    real(8),parameter :: zero=0.0d0, one=1.0d0
    real(8),allocatable :: b(:,:),w(:,:),ubak(:,:),w1(:,:)
    real(8),allocatable :: Sgp(:,:),Stmp(:,:)
    integer :: m,n,np,nn,istep,comm2,nband,mrnk2,irnk2,nrnk2,mrnk1
    integer :: nstep,ir,nr,ic,i,j,ir0,ir1,nr0,icc
    integer,allocatable :: irank_s(:), irank_r(:)
    integer :: ierr,ireq(2),itag
    include 'mpif.h'
    integer :: istatus(MPI_STATUS_SIZE,2)
    integer :: mb0_0,mb1_0,nb0,nb1
    real(8) :: c

!    call write_border( 1, "dtrmm_bp(start)" )

    m  = size( u, 1 )
    n  = size( u, 2 )
    np = nprocs_b

!    if ( nprocs_g == 1 .and. nprocs_b == 1 ) then
!       call DTRMM( 'R', 'U', 'N', 'N', m, n, one, Ssub, n, u, m )
!       call write_border( 1, "dtrmm_bp(return)" )
!       return
!    end if

    nn = maxval( ir_bp )
    comm2 = comm_b
    mrnk2 = myrank_b
    mrnk1 = myrank_g
    nstep = np - (np-1)/2

!    if ( mrnk1 == 0 .and. mrnk2 == 0 ) write(*,*) "nstep=",nstep

    nband = sum( ir_bp )
    ir0 = id_bp( myrank_b ) + 1
    ir1 = ir0 + ir_bp( myrank_b ) - 1
    nr0 = ir1 - ir0 + 1

    allocate( ubak(m,n) ); ubak=0.0d0
    allocate( w(m,nn)   ); w=zero
    allocate( b(m,nn)   ); b=zero
    allocate( w1(m,nn)  ); w1=zero

! ---

    allocate( irank_s(0:nstep) ) ; irank_s=0
    allocate( irank_r(0:nstep) ) ; irank_r=0

    do istep=0,nstep
       irank_r(istep) = myrank_b-1
       irank_s(istep) = myrank_b+1
       if ( irank_s(istep) >= np ) irank_s(istep)=MPI_PROC_NULL
       if ( istep >= mrnk2       ) irank_r(istep)=MPI_PROC_NULL
       if ( istep >  mrnk2       ) irank_s(istep)=MPI_PROC_NULL
       if ( istep > (np-1)/2-1   ) irank_s(istep)=MPI_PROC_NULL
       if ( istep > (np-1)/2-1   ) irank_r(istep)=MPI_PROC_NULL
    end do
    do istep=0,nstep
       if ( istep > mrnk2 ) then
          if ( istep >= mrnk2+2 ) irank_s(istep)=mrnk2+1
          if ( istep >= mrnk2+1 .and. mrnk2 > 0 ) irank_r(istep)=mrnk2-1
       end if
    end do
    irank_s(nstep)=MPI_PROC_NULL
    irank_r(nstep)=MPI_PROC_NULL
    if ( mrnk2 <= (np-1)/2-1  ) irank_s(nstep) = np -(np-1)/2 + mrnk2
    if ( mrnk2 >= np-(np-1)/2 ) irank_r(nstep) = mrnk2 - np + (np-1)/2
!    if ( v(1)%pinfo%me == 0 ) then
!       write(*,'(1x,"s",i2,2x,10i2)') v(2)%pinfo%me,irank_s(0:nstep)
!       write(*,'(1x,"r",i2,2x,10i2)') v(2)%pinfo%me,irank_r(0:nstep)
!    end if

! ---

    call load_div( nband, nprocs_g, myrank_g, mb0_0, mb1_0 )

    nb0 = id_bp(myrank_b) + 1
    nb1 = nb0 + ir_bp(myrank_b) - 1

    allocate( Sgp(mb0_0:mb1_0,nb0:nb1) ); Sgp=0.0d0
    allocate( Stmp(n,n) ); Stmp=0.0d0

!    call scatter_matrix( S, mb0_0,mb1_0,nb0,nb1, Sgp )

    call d_gather_matrix_3( sl, Ssub, nband,mb0_0,mb1_0,nb0,nb1, Sgp )

!    write(*,*) "sum(S**2)",sum(S**2)
!    c=sum(Sgp**2)
!    !c=sum(Ssub**2)
!    call MPI_ALLREDUCE(MPI_IN_PLACE,c,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
!    write(*,*) "c=",c
!    call stop_program("")

! ---

    itag=0

    ic=id_bp(myrank_b)+1
    ir=ic

    ubak(:,:) = u(:,:)

    w(:,1:n) = u(:,1:n)


    do istep=0,nstep

       if ( istep == 0 ) then

          call MPI_Irecv( b,size(b),MPI_REAL8,irank_r(istep),itag,comm2,ireq(1),ierr )
          call MPI_Isend( w,size(w),MPI_REAL8,irank_s(istep),itag,comm2,ireq(2),ierr )

!          call DTRMM( 'R', 'U', 'N', 'N', m, n, one, S(ir,ic), nband, u, m )
          call get_stmp( Sgp,mb0_0,mb1_0,nb0,nb1,ir,ir+n-1,ic,ic+n-1,comm_g,Stmp )
          call DTRMM( 'R', 'U', 'N', 'N', m, n, one, Stmp, n, u, m )

          ir=ir-n
          if ( ir <= 0 ) ir=nband+ir

          call MPI_Waitall( 2, ireq, istatus, ierr )
          w=b

       else if ( istep < nstep ) then

          if ( istep <= mrnk2 ) then

             call MPI_Irecv( b,size(b),MPI_REAL8,irank_r(istep),itag,comm2,ireq(1),ierr )
             call MPI_Isend( w,size(w),MPI_REAL8,irank_s(istep),itag,comm2,ireq(2),ierr )

!             call DGEMM( 'N', 'N', m, n, n, one, w, m, S(ir,ic), nband, one, u, m )
             call get_stmp( Sgp,mb0_0,mb1_0,nb0,nb1,ir,ir+n-1,ic,ic+n-1,comm_g,Stmp )
             call DGEMM( 'N', 'N', m, n, n, one, w, m, Stmp, n, one, u, m )

             ir=ir-n
             if ( ir <= 0 ) ir=nband+ir

             call MPI_Waitall( 2, ireq, istatus, ierr )
             w=b

          else if ( istep > mrnk2 ) then

             call MPI_Irecv( b,size(b),MPI_REAL8,irank_r(istep),itag,comm2,ireq(1),ierr )
             call MPI_Isend( w,size(w),MPI_REAL8,irank_s(istep),itag,comm2,ireq(2),ierr )

!             do j=1,n
!             do i=1,n
!                Stmp(i,j)=S(ic+i-1,ir+j-1)
!             end do
!             end do
             call get_stmp( Sgp,mb0_0,mb1_0,nb0,nb1,ir,ir+n-1,ic,ic+n-1,comm_g,Stmp )
             call DGEMM( 'N', 'N', m, n, n, one, ubak, m, Stmp, n, zero, w1, m )

             ir=ir-n
             if ( ir <= 0 ) ir=nband+ir

             call MPI_Waitall( 2, ireq, istatus, ierr )
             if ( mrnk2 == 0 ) then
                w=w1
             else
                w=w1+b
             end if

          end if

       else if ( istep == nstep ) then

          call MPI_Irecv( b,size(b),MPI_REAL8,irank_r(istep),itag,comm2,ireq(1),ierr )
          call MPI_Isend( w,size(w),MPI_REAL8,irank_s(istep),itag,comm2,ireq(2),ierr )
          call MPI_Waitall( 2, ireq, istatus, ierr )

          if ( mrnk2 >= np-(np-1)/2 ) u = u + b

       end if

    end do ! istep

! ---

    deallocate( Stmp )
    deallocate( Sgp  )

! ---

    deallocate( irank_r )
    deallocate( irank_s )

    deallocate( w1   )
    deallocate( b    )
    deallocate( w    )
    deallocate( ubak )

!    call write_border( 1, "dtrmm_bp(end)" )

  end subroutine dtrmm_bp


  SUBROUTINE dtrmm_bp_0( UPLO, u, S )

    implicit none
    character(1),intent(IN) :: UPLO
    real(8),intent(INOUT)   :: u(:,:)
    real(8),intent(IN)      :: S(:,:)
    real(8),parameter :: zero=0.0d0, one=1.0d0
    real(8),allocatable :: b(:,:),w(:,:),ubak(:,:),w1(:,:)
    real(8),allocatable :: Sgp(:,:),Stmp(:,:)
    integer :: m,n,np,nn,istep,comm2,nband,mrnk2,irnk2,nrnk2,mrnk1
    integer :: nstep,ir,nr,ic,i,j,ir0,ir1,nr0,icc
    integer,allocatable :: irank_s(:), irank_r(:)
    integer :: ierr,ireq(2),itag
    include 'mpif.h'
    integer :: istatus(MPI_STATUS_SIZE,2)
    integer :: mb0_0,mb1_0,nb0,nb1

    call write_border( 1, "dtrmm_bp_0(start)" )
#ifdef TEST
    m  = size( u, 1 )
    n  = size( u, 2 )
    np = v(2)%pinfo%np

    if ( np == 1 ) call DTRMM( 'R', 'U', 'N', 'N', m, n, one, S, n, u, m )

    nn = maxval( v(2)%pinfo%ir )
    comm2 = v(2)%pinfo%comm
    mrnk2 = v(2)%pinfo%me
    mrnk1 = v(1)%pinfo%me
    nstep = np - (np-1)/2

    if ( mrnk1 == 0 .and. mrnk2 == 0 ) write(*,*) "nstep=",nstep

    nband = sum( v(2)%pinfo%ir )
    ir0 = v(2)%pinfo%id( v(2)%pinfo%me ) + 1
    ir1 = ir0 + v(2)%pinfo%ir( v(2)%pinfo%me ) - 1
    nr0 = ir1 - ir0 + 1

    allocate( ubak(m,n) ) ; ubak=0.0d0
    allocate( w(m,nn)   ) ; w=zero
    allocate( b(m,nn)   ) ; b=zero
    allocate( w1(m,nn)  ) ; w1=zero

! ---

    allocate( irank_s(0:nstep) ) ; irank_s=0
    allocate( irank_r(0:nstep) ) ; irank_r=0

    do istep=0,nstep
       irank_r(istep) = v(2)%pinfo%me-1
       irank_s(istep) = v(2)%pinfo%me+1
       if ( irank_s(istep) >= np ) irank_s(istep)=MPI_PROC_NULL
       if ( istep >= mrnk2       ) irank_r(istep)=MPI_PROC_NULL
       if ( istep >  mrnk2       ) irank_s(istep)=MPI_PROC_NULL
       if ( istep > (np-1)/2-1   ) irank_s(istep)=MPI_PROC_NULL
       if ( istep > (np-1)/2-1   ) irank_r(istep)=MPI_PROC_NULL
    end do
    do istep=0,nstep
       if ( istep > mrnk2 ) then
          if ( istep >= mrnk2+2 ) irank_s(istep)=mrnk2+1
          if ( istep >= mrnk2+1 .and. mrnk2 > 0 ) irank_r(istep)=mrnk2-1
       end if
    end do
    irank_s(nstep)=MPI_PROC_NULL
    irank_r(nstep)=MPI_PROC_NULL
    if ( mrnk2 <= (np-1)/2-1  ) irank_s(nstep) = np -(np-1)/2 + mrnk2
    if ( mrnk2 >= np-(np-1)/2 ) irank_r(nstep) = mrnk2 - np + (np-1)/2
!    if ( v(1)%pinfo%me == 0 ) then
!       write(*,'(1x,"s",i2,2x,10i2)') v(2)%pinfo%me,irank_s(0:nstep)
!       write(*,'(1x,"r",i2,2x,10i2)') v(2)%pinfo%me,irank_r(0:nstep)
!    end if

! ---

    call load_div( nband, v(1)%pinfo%np, v(1)%pinfo%me, mb0_0, mb1_0 )

    nb0 = v(2)%pinfo%id(v(2)%pinfo%me)+1
    nb1 = nb0+v(2)%pinfo%ir(v(2)%pinfo%me)-1

    allocate( Sgp(mb0_0:mb1_0,nb0:nb1) ) ; Sgp=0.0d0
    allocate( Stmp(n,n) ) ; Stmp=0.0d0

    call scatter_matrix( S, mb0_0,mb1_0,nb0,nb1, Sgp )

! ---

    itag=0

    ic=v(2)%pinfo%id(v(2)%pinfo%me)+1
    ir=ic

    ubak(:,:) = u(:,:)

    w(:,1:n) = u(:,1:n)


    do istep=0,nstep

       if ( istep == 0 ) then

          call mpi_irecv( b,size(b),MPI_REAL8,irank_r(istep),itag,comm2,ireq(1),ierr )
          call mpi_isend( w,size(w),MPI_REAL8,irank_s(istep),itag,comm2,ireq(2),ierr )

          call DTRMM( 'R', 'U', 'N', 'N', m, n, one, S(ir,ic), nband, u, m )
!          call get_stmp( Sgp,mb0_0,mb1_0,nb0,nb1,ir,ir+n-1,ic,ic+n-1,v(1)%pinfo%comm,Stmp )
!          call DTRMM( 'R', 'U', 'N', 'N', m, n, one, Stmp, n, u, m )

          ir=ir-n
          if ( ir <= 0 ) ir=nband+ir

          call mpi_waitall( 2, ireq, istatus, ierr )
          w=b

       else if ( istep < nstep ) then

          if ( istep <= mrnk2 ) then

             call mpi_irecv( b,size(b),MPI_REAL8,irank_r(istep),itag,comm2,ireq(1),ierr )
             call mpi_isend( w,size(w),MPI_REAL8,irank_s(istep),itag,comm2,ireq(2),ierr )

             call DGEMM( 'N', 'N', m, n, n, one, w, m, S(ir,ic), nband, one, u, m )
             !call get_stmp( Sgp,mb0_0,mb1_0,nb0,nb1,ir,ir+n-1,ic,ic+n-1,v(1)%pinfo%comm,Stmp )
             !call DGEMM( 'N', 'N', m, n, n, one, w, m, Stmp, n, one, u, m )

             ir=ir-n
             if ( ir <= 0 ) ir=nband+ir

             call mpi_waitall( 2, ireq, istatus, ierr )
             w=b

          else if ( istep > mrnk2 ) then

             call mpi_irecv( b,size(b),MPI_REAL8,irank_r(istep),itag,comm2,ireq(1),ierr )
             call mpi_isend( w,size(w),MPI_REAL8,irank_s(istep),itag,comm2,ireq(2),ierr )

             do j=1,n
             do i=1,n
                Stmp(i,j)=S(ic+i-1,ir+j-1)
             end do
             end do
             call get_stmp( Sgp,mb0_0,mb1_0,nb0,nb1,ir,ir+n-1,ic,ic+n-1,v(1)%pinfo%comm,Stmp )
             call DGEMM( 'N', 'N', m, n, n, one, ubak, m, Stmp, n, zero, w1, m )

             ir=ir-n
             if ( ir <= 0 ) ir=nband+ir

             call mpi_waitall( 2, ireq, istatus, ierr )
             if ( mrnk2 == 0 ) then
                w=w1
             else
                w=w1+b
             end if

          end if

       else if ( istep == nstep ) then

          call mpi_irecv( b,size(b),MPI_REAL8,irank_r(istep),itag,comm2,ireq(1),ierr )
          call mpi_isend( w,size(w),MPI_REAL8,irank_s(istep),itag,comm2,ireq(2),ierr )
          call mpi_waitall( 2, ireq, istatus, ierr )

          if ( mrnk2 >= np-(np-1)/2 ) u = u + b

       end if

    end do ! istep

! ---

    deallocate( Stmp )
    deallocate( Sgp  )

! ---

    deallocate( irank_r )
    deallocate( irank_s )

    deallocate( w1   )
    deallocate( b    )
    deallocate( w    )
    deallocate( ubak )
#endif
    call write_border( 1, "dtrmm_bp_0(end)" )

  END SUBROUTINE dtrmm_bp_0


  subroutine load_div( m, np, me, m0, m1, ir_out, id_out )
    implicit none
    integer,intent(in) :: m, np
    integer,optional,intent(in)  :: me
    integer,optional,intent(out) :: m0, m1
    integer,allocatable,optional,intent(out) :: ir_out(:), id_out(:)
    integer,allocatable :: ir(:)
    integer :: i,ip
    allocate( ir(0:np-1) ) ; ir=0
    do i=1,m
       ip=mod(i-1,np)
       ir(ip)=ir(ip)+1
    end do
    if ( present(me) .and. present(m0) .and. present(m1) ) then
       m0 = sum( ir(0:me) ) - ir(me) + 1
       m1 = m0 + ir(me) - 1
    end if
    if ( present(ir_out) .and. present(id_out) ) then
       if ( allocated(ir_out) ) deallocate(ir_out)
       if ( allocated(id_out) ) deallocate(id_out)
       allocate( ir_out(0:np-1) ); ir_out=0
       allocate( id_out(0:np-1) ); id_out=0
       ir_out=ir
       do ip=0,np-1
          id_out(ip) = sum(ir_out(0:ip)) - ir_out(ip)
       end do
    end if
    deallocate( ir )
  end subroutine load_div


  SUBROUTINE scatter_matrix( S, mb0,mb1,nb0,nb1, Sgp )
    implicit none
    integer :: mb0,mb1,nb0,nb1
    real(8),intent(IN)  :: S(:,:)
    real(8),intent(OUT) :: Sgp(mb0:mb1,nb0:nb1)
    integer :: i,j,n,is,ie,js,je,nblk
    n=size(S,1)
    nblk=nb1-nb0+1
    Sgp(:,:)=0.0d0
    do js=1,n,nblk
       je=js+nblk-1
       do is=1,n,nblk
          ie=is+nblk-1
          if ( is == js ) then
             do j=js,je
             do i=js,j
                if ( mb0 <= i .and. i <= mb1 .and. nb0 <= j .and. j <= nb1 ) then
                   Sgp(i,j)=S(i,j)
                end if
             end do ! i
             end do ! j
          else if ( is < js ) then
             do j=js,je
             do i=is,ie
                if ( mb0 <= i .and. i <= mb1 .and. nb0 <= j .and. j <= nb1 ) then
                   Sgp(i,j)=S(i,j)
                end if
             end do ! i
             end do ! j
          else if ( is > js ) then
             do j=js,je
             do i=is,ie
                if ( mb0 <= i .and. i <= mb1 .and. nb0 <= j .and. j <= nb1 ) then
                   Sgp(i,j)=S(js+i-is,is+j-js)
                end if
             end do ! i
             end do ! j
          end if
       end do
    end do
  END SUBROUTINE scatter_matrix


  SUBROUTINE get_Stmp( Sgp, mb0,mb1,nb0,nb1, ir0,ir1,ic0,ic1,comm, Stmp )
    implicit none
    integer,intent(IN)  :: mb0,mb1,nb0,nb1,ir0,ir1,ic0,ic1,comm
    real(8),intent(IN)  :: Sgp(mb0:mb1,nb0:nb1)
    real(8),intent(OUT) :: Stmp(ir0:ir1,ic0:ic1)
    integer :: ir,ic
    include 'mpif.h'
    Stmp=0.0d0
    do ic=nb0,nb1
    do ir=mb0,mb1
       if ( ir0 <= ir .and. ir <= ir1 .and. ic0 <= ic .and. ic <= ic1 ) then
          Stmp(ir,ic)=Sgp(ir,ic)
       end if
    end do
    end do
    call MPI_ALLREDUCE( MPI_IN_PLACE, Stmp, size(Stmp), MPI_REAL8, MPI_SUM, comm, ir )
  END SUBROUTINE get_Stmp


  SUBROUTINE d_gather_matrix_2( sl, Hsub, n,mb0,mb1,nb0,nb1, H )
    implicit none
    type(slinfo) :: sl
    integer,intent(IN) :: n,mb0,mb1,nb0,nb1
    real(8),intent(IN) :: Hsub(:,:)
    real(8),intent(OUT) :: H(mb0:mb1,nb0:nb1)
    real(8),allocatable :: Htmp(:,:)
    integer :: i0,j0,i1,j1,ip,jp,is,js,i,j,ii,jj,ib,jb,ierr
    integer :: myrank
    include 'mpif.h'

    call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )

    H(:,:) = 0.0d0

    allocate( Htmp(sl%mbsize,sl%nbsize) ) ; Htmp=0.0d0

    jp=-1
    js= 0
    do j0=1,n,sl%nbsize
       j1=min(j0+sl%nbsize-1,n)
       jj=j1-j0+1
       jp=mod(jp+1,sl%npcol)
       ip=-1
       is= 0
       do i0=1,n,sl%mbsize
          i1=min(i0+sl%mbsize-1,n)
          ii=i1-i0+1
          ip=mod(ip+1,sl%nprow)
          if ( ip == sl%myrow .and. jp == sl%mycol ) then
             do j=1,jj
             do i=1,ii
                Htmp(i,j) = Hsub(is+i,js+j)
             end do
             end do
          end if
          call MPI_BCAST(Htmp,size(Htmp),MPI_REAL8,sl%map_2to1(ip,jp),MPI_COMM_WORLD,ierr)
          do j=1,jj
          do i=1,ii
             ib=i0+i-1
             jb=j0+j-1
             if ( mb0<=ib.and.ib<=mb1 .and. nb0<=jb.and.jb<=nb1 ) then
                H(ib,jb) = Htmp(i,j)
             end if
          end do
          end do
          if ( ip == sl%myrow ) is = is + ii
       end do
       if ( jp == sl%mycol ) js = js + jj
    end do

    deallocate( Htmp )

! ---

  END SUBROUTINE d_gather_matrix_2


  SUBROUTINE d_gather_matrix_3( sl, Hsub, n,mb0,mb1,nb0,nb1, H )
    implicit none
    type(slinfo) :: sl
    integer,intent(IN) :: n,mb0,mb1,nb0,nb1
    real(8),intent(IN) :: Hsub(:,:)
    real(8),intent(OUT) :: H(mb0:mb1,nb0:nb1)
    real(8),allocatable :: Htmp(:,:)
    integer :: i0,j0,i1,j1,ip,jp,is,js,i,j,ii,jj,ib,jb,ierr
    integer :: myrank,nblk,iblk,jblk,ia,ja,i00,j00
    include 'mpif.h'

    call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )

    H(:,:) = 0.0d0

    allocate( Htmp(sl%mbsize,sl%nbsize) ) ; Htmp=0.0d0

    nblk = nb1 - nb0 + 1

    jp=-1
    js= 0
    do j0=1,n,sl%nbsize
       j1=min(j0+sl%nbsize-1,n)
       jj=j1-j0+1
       jp=mod(jp+1,sl%npcol)
       ip=-1
       is= 0
       do i0=1,n,sl%mbsize
          i1=min(i0+sl%mbsize-1,n)
          ii=i1-i0+1
          ip=mod(ip+1,sl%nprow)
          if ( ip == sl%myrow .and. jp == sl%mycol ) then
             do j=1,jj
             do i=1,ii
                Htmp(i,j) = Hsub(is+i,js+j)
             end do
             end do
          end if
          call MPI_BCAST(Htmp,size(Htmp),MPI_REAL8,sl%map_2to1(ip,jp),MPI_COMM_WORLD,ierr)
          do j=1,jj
          do i=1,ii
             ib=i0+i-1
             jb=j0+j-1
             if ( ib > jb ) cycle
             if ( mb0<=ib.and.ib<=mb1 .and. nb0<=jb.and.jb<=nb1 ) then
                H(ib,jb) = Htmp(i,j)
             end if
             if ( ib == jb ) cycle
             iblk=(jb-1)/nblk+1
             jblk=(ib-1)/nblk+1
             i00=(iblk-1)*nblk+1
             j00=(jblk-1)*nblk+1
             ia=jb-i00+1
             ja=ib-j00+1
             ib=ja+i00-1
             jb=ia+j00-1
             if ( mb0<=ib.and.ib<=mb1 .and. nb0<=jb.and.jb<=nb1 ) then
                H(ib,jb) = Htmp(i,j)
             end if
          end do
          end do
          if ( ip == sl%myrow ) is = is + ii
       end do
       if ( jp == sl%mycol ) js = js + jj
    end do

    deallocate( Htmp )

! ---

  END SUBROUTINE d_gather_matrix_3


END MODULE dtrmm_bp_module
