module z_ps_nloc2_op_module

  implicit none

  private
  public :: z_op_ps_nloc2
  public :: z_op_ps_nloc2_hp
  public :: z_reset_init_flag_ps_nloc2_op

  integer,allocatable :: ompMJJ1(:,:), ompMJJ2(:,:)
  integer,allocatable :: omplns1(:,:), omplns2(:,:)
  integer :: ompnzlma1,ompnzlma2,ompmyrank,ompnprocs,ompn1,ompn2

  integer,allocatable :: n_i2j(:)
  integer,allocatable :: i2j(:,:,:)

!$omp threadprivate( ompmyrank,ompnprocs,ompnzlma1,ompnzlma2,ompn1,ompn2 )

  logical :: init_done = .false.

contains

  subroutine z_reset_init_flag_ps_nloc2_op
    implicit none
    init_done = .false.
  end subroutine z_reset_init_flag_ps_nloc2_op


  subroutine z_init_op_ps_nloc2_hp
    !$ use omp_lib
    use rgrid_module, only: Igrid
    use ps_nloc2_variables, only: nzlma, lma_nsend, sendmap, JJP, MJJ
    use parallel_module, only: np_grid
    implicit none
    integer :: i,j,k,ompblock,ompblock0,mm

    call write_border( 0, " z_init_op_ps_nloc2_hp(start)" )

    mm = Igrid(2,0) - Igrid(1,0) + 1

!$omp parallel private( i,j,k,ompblock,ompblock0 )

    ompnprocs = 1
!$  ompnprocs = omp_get_num_threads()
    ompmyrank = 0
!$  ompmyrank = omp_get_thread_num()

    ompblock=0
    ompblock0=0
    do i=1,nzlma
      j=mod(i-1,ompnprocs)
      if ( j == ompmyrank ) then
        ompblock = ompblock + 1
      else if ( j < ompmyrank ) then
        ompblock0 = ompblock0 + 1
      end if
    end do

    ompnzlma1 = ompblock0 + 1
    ompnzlma2 = ompnzlma1 + ompblock - 1

! ---

!$omp master
    if ( .not.allocated(omplns1) ) then
      allocate( omplns1(0:np_grid-1,0:ompnprocs-1) )
      allocate( omplns2(0:np_grid-1,0:ompnprocs-1) )
    end if
    omplns1=0
    omplns2=0
    if ( allocated(ompMJJ1) ) deallocate(ompMJJ1)
    if ( allocated(ompMJJ2) ) deallocate(ompMJJ2)
    allocate( ompMJJ1(0:nzlma,0:ompnprocs-1)     ) ; ompMJJ1=0
    allocate( ompMJJ2(0:nzlma,0:ompnprocs-1)     ) ; ompMJJ2=0
!$omp end master
!$omp barrier

    omplns1(:,ompmyrank)=-1
    omplns2(:,ompmyrank)=-2
    do i=0,np_grid-1
       do j=1,lma_nsend(i)
          k = sendmap(j,i)
          if ( ompnzlma1 <= k .and. k <= ompnzlma2 ) then
            omplns2(i,ompmyrank) = j
          end if
       end do ! j
       do j=lma_nsend(i),1,-1
          k = sendmap(j,i)
          if ( ompnzlma1 <= k .and. k <= ompnzlma2 ) then
            omplns1(i,ompmyrank) = j
          end if
       end do ! j
    end do ! i

! ---

    ompblock=0
    ompblock0=0
    do i=1,mm
      j=mod(i-1,ompnprocs)
      if ( j == ompmyrank ) then
        ompblock = ompblock + 1
      else if ( j < ompmyrank ) then
        ompblock0 = ompblock0 + 1
      end if
    end do

    ompn1 = Igrid(1,0) + ompblock0
    ompn2 = ompn1 + ompblock - 1

! ---

    ompMJJ1(:,ompmyrank)=-1
    ompMJJ2(:,ompmyrank)=-2
    do k=1,nzlma
      do j=1,MJJ(k)
        i=JJP(j,k)
        if ( ompn1 <= i .and. i<= ompn2 ) then
          ompMJJ1(k,ompmyrank) = j
          exit
        end if
      end do
      do j=MJJ(k),1,-1
        i=JJP(j,k)
        if ( ompn1 <= i .and. i<= ompn2 ) then
          ompMJJ2(k,ompmyrank) = j
          exit
        end if
      end do
    end do

!$omp end parallel

    init_done = .true.
  
    call write_border( 0, " z_init_op_ps_nloc2_hp(end)" )

  end subroutine z_init_op_ps_nloc2_hp


  subroutine z_op_ps_nloc2_hp( tpsi, htpsi, n_in, k_in, s_in )
    use ps_nloc2_variables, only: MJJ, JJP, Mlma, z_uVk, z_uVunk, iuV, &
                                  z_uVunk0, flag_backup_uVunk_ps_nloc2, &
                                  nrlma_xyz, num_2_rank, lma_nsend, nzlma, &
                                  z_sbufnl3, z_rbufnl3, sendmap, recvmap
    use rgrid_module, only: Igrid, dV
    use rgrid_mol_module, only: iswitch_eqdiv
    use parallel_module, only: comm_grid
    implicit none
    complex(8),intent(in) :: tpsi(:,:)
    complex(8),intent(inout) :: htpsi(:,:)
    integer,optional,intent(in) :: n_in, k_in, s_in
    include 'mpif.h'
    integer :: i,ib,i1,i2,i3,j,lma,m,ML0,n,nb,k,s,i0
    integer :: ierr,nreq
    integer :: irank,jrank,istatus(MPI_STATUS_SIZE,512),ireq(512)
    complex(8) :: tmp_sum
    real(8) :: c, ttmp(2), ttmp1(2)

    if ( Mlma <= 0 ) return
    if ( .not.init_done ) call z_init_op_ps_nloc2_hp

    n=1 ; if ( present(n_in) ) n=n_in
    k=1 ; if ( present(k_in) ) k=k_in
    s=1 ; if ( present(s_in) ) s=s_in

    if ( flag_backup_uVunk_ps_nloc2 ) then
      call z_op_ps_nloc2_with_restore_uVunk( htpsi, n, k, s )
      return
    end if

    nb = size( tpsi, 2 )
    i0 = Igrid(1,0)

    !call watchb_omp( ttmp )

    do ib = 1, nb
      do lma = ompnzlma1, ompnzlma2
        tmp_sum=(0.0d0,0.0d0)
        do j = 1, MJJ(lma)
          i = JJP(j,lma) - i0 + 1
          tmp_sum = tmp_sum + conjg( z_uVk(j,lma,k) )*tpsi(i,ib)
        end do
        z_uVunk(lma,ib) = iuV(lma)*dV*tmp_sum
      end do ! lma
    end do ! ib

    !call watchb_omp( ttmp, time_nlpp(1,1) )

    select case( iswitch_eqdiv )
    case default

      do i = 1, 6

  !$omp barrier
          
        select case(i)
        case(1,3,5)
          j=i+1
          do ib = 1, nb
            z_uVunk0(ompnzlma1:ompnzlma2,ib) = z_uVunk(ompnzlma1:ompnzlma2,ib)
          end do
        case(2,4,6)
          j=i-1
        end select

  !$omp barrier

        do m = 1, nrlma_xyz(i)
          !call watchb_omp( ttmp1 )
          nreq=0
          irank=num_2_rank(m,i)
          jrank=num_2_rank(m,j)
          if ( irank >= 0 ) then
            do ib = 1, nb
              i2 = omplns1(irank,ompmyrank) - 1 + lma_nsend(irank)*(ib-1)
              do i1 = omplns1(irank,ompmyrank), omplns2(irank,ompmyrank)
                i2 = i2 + 1
                z_sbufnl3(i2,m,i) = z_uVunk0(sendmap(i1,irank),ib)
              end do
            end do
  !$omp barrier
  !$omp master
            nreq=nreq+1
            call MPI_Isend(z_sbufnl3(1,m,i),lma_nsend(irank)*nb &
            ,MPI_COMPLEX16,irank,1,comm_grid,ireq(nreq),ierr)
  !$omp end master
          end if
          !call watchb_omp( ttmp1, time_nlpp(1,4) )

  !$omp master
          if ( jrank >= 0 ) then
            nreq=nreq+1
            call MPI_Irecv(z_rbufnl3(1,m,i),lma_nsend(jrank)*nb &
            ,MPI_COMPLEX16,jrank,1,comm_grid,ireq(nreq),ierr)
          end if
          !call watchb_omp( ttmp1, time_nlpp(1,5) )
          call MPI_Waitall(nreq,ireq,istatus,ierr)
  !$omp end master

  !$omp barrier
          !call watchb_omp( ttmp1, time_nlpp(1,6) )
          if ( jrank >= 0 ) then
            do ib = 1, nb
              i2 = omplns1(jrank,ompmyrank) - 1 + lma_nsend(jrank)*(ib-1)
              do i1 = omplns1(jrank,ompmyrank), omplns2(jrank,ompmyrank)
                i2 = i2 + 1
                i3 = recvmap(i1,jrank)
                z_uVunk(i3,ib) = z_uVunk(i3,ib) + z_rbufnl3(i2,m,i)
              end do
            end do
          end if
          !call watchb_omp( ttmp1, time_nlpp(1,7) )

        end do ! m

      end do ! i

    case( 2 )

!$omp barrier
      call comm_eqdiv_ps_nloc2_mol( z_uVunk )

    end select

!$omp barrier
    !call watchb_omp( ttmp, time_nlpp(1,2) )

    do ib = 1, nb
      do lma = 1, nzlma
        do j = ompMJJ1(lma,ompmyrank), ompMJJ2(lma,ompmyrank)
          i = JJP(j,lma) - i0 + 1
          htpsi(i,ib) = htpsi(i,ib) + z_uVk(j,lma,k)*z_uVunk(lma,ib)
        end do
      end do
    end do

!$omp barrier
    !call watchb_omp( ttmp, time_nlpp(1,3) )

    return 
      
  end subroutine z_op_ps_nloc2_hp


  subroutine z_op_ps_nloc2_with_restore_uVunk( htpsi, n, k, s )
    use rgrid_module, only: Igrid
    use ps_nloc2_variables, only: z_uVk, z_uVunk, JJP, restore_uVunk_ps_nloc2
    implicit none
    complex(8),intent(inout) :: htpsi(:,:)
    integer,intent(in) :: n, k, s
    integer :: i,j,ib,nb,lma,i0,nzlma

!    call write_border( 1, " z_op_ps_nloc2_with_restore_uVUnk(start)" )

    nb = size( htpsi, 2 )
    i0 = Igrid(1,0)
    nzlma = size( JJP, 2 )

    call restore_uVunk_ps_nloc2( z_uVunk(:,1:nb), n, k, s )

    do ib = 1, nb
      do lma = 1, nzlma
        do j = ompMJJ1(lma,ompmyrank), ompMJJ2(lma,ompmyrank)
          i = JJP(j,lma) - i0 + 1
          htpsi(i,ib) = htpsi(i,ib) + z_uVk(j,lma,k)*z_uVunk(lma,ib)
        end do
      end do
    end do

!    call write_border( 1, " z_op_ps_nloc2_with_restore_uVUnk(end)" )

  end subroutine z_op_ps_nloc2_with_restore_uVunk


  subroutine z_init_op_ps_nloc2
    use rgrid_module, only: Igrid
    use ps_nloc2_variables, only: MJJ, JJP
    implicit none
    integer :: ML_0,ML_1,lma,i,j,n,nzlma

    call write_border( 0, " z_init_op_ps_nloc2(start)" )

    ML_0 = Igrid(1,0)
    ML_1 = Igrid(2,0)
    nzlma = size(MJJ)

! --- inverse map ( grid label i --> (j,lma) )

    if ( allocated(i2j) ) deallocate(i2j)
    if ( allocated(n_i2j) ) deallocate(n_i2j)

    allocate( n_i2j(ML_0:ML_1) ); n_i2j=0

    do lma=1,nzlma
      do j=1,MJJ(lma)
        i=JJP(j,lma)
        n_i2j(i)=n_i2j(i)+1
      end do
    end do

    n=maxval( n_i2j )
    allocate( i2j(2,n,ML_0:ML_1) ); i2j=0

    n_i2j(:)=0
    do lma=1,nzlma
      do j=1,MJJ(lma)
        i=JJP(j,lma)
        n_i2j(i)=n_i2j(i)+1
        i2j(1,n_i2j(i),i)=j
        i2j(2,n_i2j(i),i)=lma
      end do
    end do

    init_done = .true.

    call write_border( 0, " z_init_op_ps_nloc2_hp(end)" )

  end subroutine z_init_op_ps_nloc2


  subroutine z_op_ps_nloc2( tpsi, htpsi, n, k, s )
    use rgrid_module, only: Igrid, dV
    use rgrid_mol_module, only: iswitch_eqdiv
    use ps_nloc2_variables, only: Mlma, MJJ, JJP, z_uVk, z_uVunk, z_uVunk0, iuV, &
                                  num_2_rank, lma_nsend, sendmap, recvmap, &
                                  z_sbufnl1, z_rbufnl1, nrlma_xyz
    use parallel_module, only: comm_grid
    implicit none
    real(8),intent(in) :: tpsi(:,n:)
    real(8),intent(inout) :: htpsi(:,n:)
    integer,intent(in) :: n, k, s
    include 'mpif.h'
    integer :: i0,i,ib,j,jb,i1,i2,i3,m,lma,nb,ierr,nreq,lmani,lmanj
    integer :: irank,jrank,istatus(MPI_STATUS_SIZE,512),ireq(512),nzlma
    integer :: nb_0_omp,nb_1_omp,nzlma_0_omp,nzlma_1_omp,m_0,m_1,n_0,n_1
    real(8) :: et0,et1,ttmp(2)

    if ( Mlma <= 0 ) return
    if ( .not.init_done ) call z_init_op_ps_nloc2

!$omp barrier
    !call watchb_omp( ttmp )

    nb = size( tpsi, 2 )
    i0 = Igrid(1,0)
    nzlma = size(MJJ)

    call calc_range_omp(nb,nzlma,nb_0_omp,nb_1_omp,nzlma_0_omp,nzlma_1_omp)

    !call watchb_omp( ttmp, time_nlpp(1,4) )

    do ib = nb_0_omp, nb_1_omp
      do lma = nzlma_0_omp, nzlma_1_omp
        z_uVunk(lma,ib)=0.0d0
        do j = 1, MJJ(lma)
          i = JJP(j,lma) - i0 + 1
          z_uVunk(lma,ib) = z_uVunk(lma,ib) + conjg( z_uVk(j,lma,k) )*tpsi(i,ib)
        end do
        z_uVunk(lma,ib) = iuV(lma)*dV*z_uVunk(lma,ib)
      end do
    end do

    !call watchb_omp( ttmp, time_nlpp(1,1) )

    select case( iswitch_eqdiv )
    case default

      do i = 1, 6

!$omp barrier

        select case(i)
        case(1,3,5)
          j=i+1
          do ib = nb_0_omp, nb_1_omp
          do lma = nzlma_0_omp, nzlma_1_omp
            z_uVunk0(lma,ib) = z_uVunk(lma,ib)
          end do
          end do
        case(2,4,6)
          j=i-1
        end select

!$omp barrier

        do m = 1, nrlma_xyz(i)
          irank = num_2_rank(m,i)
          jrank = num_2_rank(m,j)
          nreq=0 
          if ( irank >= 0 ) then
            lmani = lma_nsend(irank)
            call calc_range_omp(nb,lmani,m_0,m_1,n_0,n_1)
            do ib = m_0, m_1
            do i1 = n_0, n_1
              i2 = i1 + (ib-1)*lmani
              z_sbufnl1(i2) = z_uVunk0(sendmap(i1,irank),ib)
            end do
            end do
!$omp barrier
!$omp master
            nreq=nreq+1
            call MPI_Isend(z_sbufnl1,lmani*nb,MPI_COMPLEX16,irank,1,comm_grid,ireq(nreq),ierr)
!$omp end master
          end if
!$omp master
          if ( jrank >= 0 ) then
            lmanj = lma_nsend(jrank)
            nreq=nreq+1
            call MPI_Irecv(z_rbufnl1,lmanj*nb,MPI_COMPLEX16,jrank,1,comm_grid,ireq(nreq),ierr)
          end if
          call MPI_Waitall(nreq,ireq,istatus,ierr)
!$omp end master
!$omp barrier
          if ( jrank >= 0 ) then
            lmanj = lma_nsend(jrank)
            call calc_range_omp(nb,lmanj,m_0,m_1,n_0,n_1)
            do ib = m_0, m_1
            do i1 = n_0, n_1
              i2 = i1 + (ib-1)*lmanj
              i3 = recvmap(i1,jrank)
              z_uVunk(i3,ib) = z_uVunk(i3,ib) + z_rbufnl1(i2)
            end do
            end do
          end if

        end do ! m

      end do ! i

    case( 2 )

!$omp barrier
      call comm_eqdiv_ps_nloc2_mol( z_uVunk )

    end select

!$omp barrier
    !call watchb_omp( ttmp, time_nlpp(1,2) )

    do ib = 1, nb
      do lma = 1, nzlma
!$omp do
        do j = 1, MJJ(lma)
          i = JJP(j,lma) - i0 + 1
          htpsi(i,ib) = htpsi(i,ib) + z_uVk(j,lma,k)*z_uVunk(lma,ib)
        end do
!$omp end do
       end do
    end do

!$omp barrier
    !call watchb_omp( ttmp, time_nlpp(1,3) )

  end subroutine z_op_ps_nloc2


  subroutine calc_range_omp( nb, nz, nb_0, nb_1, nz_0, nz_1 )
    !$ use omp_lib
    implicit none
    integer,intent(in)  :: nb, nz
    integer,intent(out) :: nb_0,nb_1,nz_0,nz_1
    integer :: mp,ip,k0,k1,id,a,i,j
    real(8) :: et0,et1
!$  et0=omp_get_wtime()
    mp=1
!$  mp=omp_get_num_threads()
    ip=0
!$  ip=omp_get_thread_num()
    k0=gcd(nb,mp)
    k1=mp/k0
    a=-1
    loop_j : do j=0,k0-1
    do i=0,k1-1
      a=a+1
      if ( a == ip ) then
        nb_0 = j*(nb/k0) + 1
        nb_1 = nb_0 + (nb/k0) - 1
        nz_0 = i*(nz/k1) + 1
        nz_1 = nz_0 + (nz/k1) - 1
        if ( i == k1-1 ) nz_1 = nz
        exit loop_j
      end if
    end do
    end do loop_j
!    et1=omp_get_wtime()
!    write(*,'(1x,8i4)') ip,mp,nb_0,nb_1,nz_0,nz_1
!    write(*,*) "time=",et1-et0
  end subroutine calc_range_omp

  function gcd(m0,n0)
    implicit none
    integer :: gcd,m0,n0
    integer :: m,n,mtmp,loop
    if ( m0 >= n0 ) then
      m=m0
      n=n0
    else
      m=n0
      n=m0
    end if
    do loop=1,10000
      if ( n == 0 ) exit
      mtmp = n
      n = mod(m,n)
      m = mtmp
    end do
    gcd = m
  end function gcd


  subroutine comm_eqdiv_ps_nloc2_mol( uVunk )
    use ps_nloc2_variables, only: lma_nsend, z_sbufnl, z_rbufnl, sendmap, recvmap, &
                                  z_uVunk
    use parallel_module, only: nprocs_g, myrank_g, comm_grid
    implicit none
    complex(8),intent(inout) :: uVunk(:,:)
    integer :: nreq,irank,m,i1,i2,i3,ib,ierr,nb,nzlma
    include 'mpif.h'
    integer :: istatus(MPI_STATUS_SIZE,512),ireq(512)
    nzlma=size(uVunk,1)
    nb=size(uVunk,2)
    nreq=0
    do irank = 0, nprocs_g-1
      m = lma_nsend(irank)
      if ( irank == myrank_g .or. m <= 0 ) cycle
      do ib = 1, nb
!$omp do
        do i1 = 1, m
          i2 = i1 + (ib-1)*m
          z_sbufnl(i2,irank) = z_uVunk(sendmap(i1,irank),ib)
        end do
!$omp end do
      end do !ib
!$omp master
      nreq = nreq + 1
      call MPI_Isend(z_sbufnl(1,irank),m*nb,MPI_COMPLEX16,irank,1,comm_grid,ireq(nreq),ierr)
      nreq=nreq+1
      call MPI_Irecv(z_rbufnl(1,irank),m*nb,MPI_COMPLEX16,irank,1,comm_grid,ireq(nreq),ierr)
!$omp end master
    end do !irank
!$omp master
    if ( nreq > 0 ) call MPI_Waitall(nreq,ireq,istatus,ierr)
!$omp end master
!$omp barrier
    do irank = 0, nprocs_g-1
      m = lma_nsend(irank)
      if ( irank == myrank_g .or. m <= 0 ) cycle
      do ib = 1, nb
!$OMP do
        do i1 = 1, m
          i2 = i1 + (ib-1)*m
          i3 = recvmap(i1,irank)
          z_uVunk(i3,ib) = z_uVunk(i3,ib) + z_rbufnl(i2,irank)
        end do
!$OMP end do
      end do !ib
    end do !irank
  end subroutine comm_eqdiv_ps_nloc2_mol


end module z_ps_nloc2_op_module
