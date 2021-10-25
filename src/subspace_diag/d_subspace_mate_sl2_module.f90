module d_subspace_mate_sl2_module

  implicit none

  private
  public :: d_subspace_mate_sl2

  real(8) :: dV=1.0d0
  logical :: init_done=.false.

  integer,parameter :: nblk_min = 4
  real(8),allocatable :: uhu(:,:)

contains

  subroutine init( u )
    use rsdft_allreduce_module, only: rsdft_allreduce
    implicit none
    real(8),intent(in) :: u(:)
    real(8) :: c
    c = sum( abs(u)**2 )
    call rsdft_allreduce( c, 'grid' )
    dV = 1.0d0/c
    init_done = .true.
  end subroutine init


  recursive subroutine mate_sub(  u, hu, m1,m2, n1,n2, nblk )
    implicit none
    real(8),intent(in) :: u(:,:)
    real(8),intent(in) :: hu(:,:)
    integer,intent(in) :: m1,m2,n1,n2,nblk
    integer :: mg,n,ns,ne,nn,m,ms,me,mm,nblkh
    integer :: ib1,ib2,jb1,jb2,ld
    real(8),parameter :: z0=0.0d0

    mg = size( u, 1 )
    ld = size( uhu, 1 )

    do ns = n1, n2, nblk
      ne = min( ns + nblk - 1, n2 )
      nn = ne - ns + 1
      jb1 = ns - n1 + 1
      jb2 = ne - n1 + 1

      do ms = m1, m2, nblk
        me = min( ms + nblk - 1, m2 )
        mm = me - ms + 1
        ib1 = ms - m1 + 1
        ib2 = me - m1 + 1

        if ( ms == ns ) then

          if ( mm <= nblk_min ) then
            do n = ns, ne
              call DGEMV('T',mg,ne-n+1,dV,u(1,n),mg,hu(1,n-n1+1),1,z0,uhu(n,n),1)
            end do
          else
            nblkh = max( nblk/2, nblk_min )
            call mate_sub( u, hu(:,jb1:jb2), ms,me, ns,ne, nblkh )
          end if

        else if ( ms > ne ) then

          call DGEMM('T','N',mm,nn,mg,dV,u(1,ms),mg,hu(1,jb1),mg,z0,uhu(ms,ns),ld)

        end if

      end do ! ms
    end do ! ns

    return
  end subroutine mate_sub


  subroutine d_subspace_mate_sl2( k, s, u )
    use parallel_module, only: myrank, MB_d, get_range_parallel,myrank_b,np_band,myrank_g
    use sl_variables, only: sl0, sl1, Dsub
    use hamiltonian_module, only: hamiltonian !, op_kinetic, op_nonlocal
    use rsdft_allreduce_module, only: rsdft_allreduce
    use rsdft_mpi_module, only: rsdft_ar=>rsdft_allreduce
    use sl_tools_module, only: slinfo, distribute_matrix
    use watch_module, only: watchb    
    ! use localpot_module, only: Vloc
    implicit none
    integer,intent(in) :: k, s
    real(8),intent(in) :: u(:,:)
    real(8),allocatable :: hu(:,:)
    real(8),parameter :: z0=0.0d0
    integer :: nb,mg,m,n,i,j,blk_size,ib,ms_next,nh,mh
    integer :: nb_0,nb_1,ns,ne,nn,ms,me,mm,ib1,ib2,jb1,jb2
    real(8),allocatable :: echk2(:,:)
    type(slinfo) :: sl
    real(8) :: ttmp(2), tt(2,8),ttmin(2,4),ttmax(2,4)
    logical :: disp_on
    integer,allocatable :: m_chk(:)
    integer :: num_blocks_1side, num_blocks_total, num_blocks_myrnk
    integer :: iblock
    logical :: np_band_is_odd
    character(3),parameter :: ON_OFF='off'

    call write_border( 1, ' d_subspace_mate_sl2(start)' )
    call check_disp_switch( disp_on, 0 )

    tt=0.0d0; call watchb( ttmp, barrier=ON_OFF )

    if ( .not.init_done ) call init( u(:,1) )

    mg = size( u, 1 )
    nb = sl0%nband

    call get_range_parallel( nb_0, nb_1, 'band' )

    blk_size = sl0%nbsize

    np_band_is_odd = ( mod( np_band, 2 ) == 1 )

    num_blocks_1side = nb/blk_size
    num_blocks_total = ((num_blocks_1side+1)*num_blocks_1side)/2
    num_blocks_myrnk = (num_blocks_total + np_band - 1)/ np_band - 1
    ! if ( disp_on ) then
    !   write(*,*) "nband=",nb
    !   write(*,*) "num_blocks_1side=",num_blocks_1side
    !   write(*,*) "num_blocks_total=",num_blocks_total
    !   write(*,*) "num_blocks_1side*blk_size=",num_blocks_1side*blk_size
    !   write(*,*) "num_blocks_myrnk(excluding diagonal part)=",num_blocks_myrnk
    !   write(*,*) "np_band_is_odd =",np_band_is_odd
    ! end if

    allocate( hu(size(u,1),blk_size) ); hu=z0

    allocate( echk2(nb,nb) ); echk2=z0

    call watchb( ttmp, tt(:,1), barrier=ON_OFF )

    do ib = nb_0, nb_1, blk_size

      ns = ib
      ne = min( ns + blk_size - 1, nb_1 )
      nn = ne - ns + 1

      call watchb( ttmp )

      hu=z0
      do ib1 = ns, ne, MB_d
        ib2 = min( ib1 + MB_d - 1, ne )
        jb1 = ib1 - ns + 1
        jb2 = ib2 - ns + 1
        call hamiltonian( u(:,ib1:ib2), hu(:,jb1:jb2), ib1,k,s  )
        !call op_kinetic( u(:,ib1:ib2), hu(:,jb1:jb2), k, Vloc(:,s)  )
      end do

      call watchb( ttmp, tt(:,2) )

      !do ib1 = ns, ne, MB_d
      !  ib2 = min( ib1 + MB_d - 1, ne )
      !  jb1 = ib1 - ns + 1
      !  jb2 = ib2 - ns + 1
      !  call op_nonlocal( u(:,ib1:ib2), hu(:,jb1:jb2), ib1,k,s  )
      !end do

      call watchb( ttmp, tt(:,3) )

      ! Diagonal part
      !
      ms = ns
      me = ne
      mm = nn

      allocate( uhu(ms:me,ns:ne) ); uhu=z0
      call mate_sub( u, hu, ms,me ,ns,ne ,blk_size )

      echk2(ms:me,ns:ne) = uhu
      deallocate( uhu )
      !
      ! -------------

      call watchb( ttmp, tt(:,4) )

      !do ms = ne+1, nb, blk_size
      !  me = min( ms + blk_size - 1, nb )
      !  mm = me - ms + 1
      !  allocate( uhu(ms:me,ns:ne) ); uhu=z0
      !  call DGEMM('T','N',mm,nn,mg,dV,u(1,ms),mg,hu,mg,z0,uhu,mm)
      !  echk2(ms:me,ns:ne) = uhu
      !  deallocate( uhu )
      !end do !ms
      !call watchb( ttmp, tt(:,4) )
      !cycle

      do iblock = 1, num_blocks_myrnk-1
        ms_next = ne+1+(iblock-1)*blk_size
        ms = mod( ms_next, nb )
        me = ms + blk_size - 1
        mm = me - ms + 1
        ! call flush_barrier(6)
        ! if ( myrank_g == 0 ) write(*,'(1x,"iblock,ms,me",3i8,3x,i8)') iblock,ms,me,myrank_b
        ! call flush_barrier(6)
        if ( ms_next <= nb ) then
          allocate( uhu(ms:me,ns:ne) ); uhu=z0
          call DGEMM('T','N',mm,nn,mg,dV,u(1,ms),mg,hu,mg,z0,uhu,mm)
          echk2(ms:me,ns:ne) = uhu
        else
          allocate( uhu(ns:ne,ms:me) ); uhu=z0
          call DGEMM('T','N',nn,mm,mg,dV,hu,mg,u(1,ms),mg,z0,uhu,nn)
          echk2(ns:ne,ms:me) = uhu
        end if
        deallocate( uhu )
      end do

      if ( num_blocks_myrnk <= 0 ) cycle

      if ( np_band_is_odd ) then

        ms_next = ne+1+(num_blocks_myrnk-1)*blk_size
        ms = mod( ms_next, nb )
        me = ms + blk_size - 1
        mm = me - ms + 1
        ! call flush_barrier(6)
        ! if ( myrank_g == 0 ) write(*,'(1x,"iblock,ms,me",3i8,3x,i8)') iblock,ms,me,myrank_b
        ! call flush_barrier(6)
        if ( ms_next <= nb ) then
          allocate( uhu(ms:me,ns:ne) ); uhu=z0
          call DGEMM('T','N',mm,nn,mg,dV,u(1,ms),mg,hu,mg,z0,uhu,mm)
          echk2(ms:me,ns:ne) = uhu
        else
          allocate( uhu(ns:ne,ms:me) ); uhu=z0
          call DGEMM('T','N',nn,mm,mg,dV,hu,mg,u(1,ms),mg,z0,uhu,nn)
          echk2(ns:ne,ms:me) = uhu
        end if
        deallocate( uhu )

      else

        if ( myrank_b <= np_band/2-1 ) then

          nh = ns + (nn+1)/2
          nn = ne - nh + 1
          ms = mod( ne+1+(num_blocks_myrnk-1)*blk_size, nb )
          me = ms + blk_size - 1
          mm = me - ms + 1
          ! call flush_barrier(6)
          ! if ( myrank_g == 0 ) write(*,'("(a)ms,me,mm,nh,ne,nn",6i8,3x,i4)') ms,me,mm,nh,ne,nn,myrank_b
          ! call flush_barrier(6)
          allocate( uhu(ms:me,nh:ne) ); uhu=z0
          call DGEMM('T','N',mm,nn,mg,dV,u(1,ms),mg,hu(1,nh-ns+1),mg,z0,uhu,mm)
          echk2(ms:me,nh:ne) = uhu
          deallocate( uhu )

        else

          ms = mod( ne+1+(num_blocks_myrnk-1)*blk_size, nb )
          mh = ms + nn/2 - 1
          mm = mh - ms + 1
          ! call flush_barrier(6)
          ! if ( myrank_g == 0 ) write(*,'("(b)ns,ne,nn,ms,mh,mm",6i8,3x,i4)') ns,ne,nn,ms,mh,mm,myrank_b
          ! call flush_barrier(6)
          allocate( uhu(ns:ne,ms:mh) ); uhu=z0
          call DGEMM('T','N',nn,mm,mg,dV,hu,mg,u(1,ms),mg,z0,uhu,nn)
          echk2(ns:ne,ms:mh) = uhu
          deallocate( uhu )

        end if

      end if

      call watchb( ttmp, tt(:,5) )

    end do !ns

    m = count(echk2/=z0)

    call watchb( ttmp, barrier=ON_OFF )

    call rsdft_allreduce( echk2, 'g' )
    call rsdft_allreduce( echk2, 'b' )

    call watchb( ttmp, tt(:,6), barrier=ON_OFF )

    ! do j = 1, nb
    !   do i = 1, nb
    !     if ( echk2(i,j) == 0.0d0 ) echk2(i,j)=echk2(j,i)
    !   end do
    ! end do

    ! allocate( m_chk(0:np_band-1) ); m_chk=0
    ! m_chk(myrank_b)=m
    ! call rsdft_allreduce( m_chk, 'b' )
    ! do j = 1, nb
    !   do i = 1, j-1
    !     echk2(i,j)=0.0d0
    !   end do
    ! end do
    ! n = count(echk2/=0.0d0)
    ! if ( disp_on ) then
    !   write(*,'(5i8,f30.15)') n,nb_0,nb_1,nb_1-nb_0+1,blk_size,sum(echk2**2)
    !   do i = 0, np_band-1
    !     write(*,*) "# of me in rank",i," =", m_chk(i)
    !   end do
    !   write(*,*) "sum of me =",sum(m_chk)
    ! end if
    ! deallocate( m_chk )

    call watchb( ttmp, tt(:,7), barrier=ON_OFF )

    if ( allocated(Dsub) ) deallocate(Dsub)

    if ( sl1%uplo /= '' ) then

      sl%nprow  = sl1%nprow
      sl%npcol  = sl1%npcol
      sl%myrow  = sl1%myrow
      sl%mycol  = sl1%mycol
      sl%mbsize = sl1%mbsize
      sl%nbsize = sl1%nbsize

      allocate( Dsub(sl1%ldr,sl1%ldc) ); Dsub=z0
      
      call distribute_matrix( sl, echk2, Dsub )

    end if

    call watchb( ttmp, tt(:,8), barrier=ON_OFF )

    ! ttmin(:,1:4)=tt(:,2:5); call rsdft_ar( ttmin, op_in='min' )
    ! ttmax(:,1:4)=tt(:,2:5); call rsdft_ar( ttmax, op_in='max' )

    ! if ( disp_on ) then
    !   write(*,'(1x,"time_mate(1)=",2f12.5)') tt(:,1)
    !   write(*,'(1x,"time_mate(2)=",2f12.5,2x,2f12.5)') ttmin(:,1), ttmax(:,1)
    !   write(*,'(1x,"time_mate(3)=",2f12.5,2x,2f12.5)') ttmin(:,2), ttmax(:,2)
    !   write(*,'(1x,"time_mate(4)=",2f12.5,2x,2f12.5)') ttmin(:,3), ttmax(:,3)
    !   write(*,'(1x,"time_mate(5)=",2f12.5,2x,2f12.5)') ttmin(:,4), ttmax(:,4)
    !   write(*,'(1x,"time_mate(6)=",2f12.5)') tt(:,6)
    !   write(*,'(1x,"time_mate(7)=",2f12.5)') tt(:,7)
    !   write(*,'(1x,"time_mate(8)=",2f12.5)') tt(:,8)
    ! end if

    call write_border( 1, ' d_subspace_mate_sl2(end)' )

  end subroutine d_subspace_mate_sl2

end module d_subspace_mate_sl2_module
