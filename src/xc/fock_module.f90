module fock_module

  use xc_hybrid_module, only: iflag_hybrid, alpha_hf &
                             ,FKMB_0,FKMB_1,FKBZ_0,FKBZ_1,FOCK_0,FOCK_1 &
                             ,occ_factor,gamma_hf
  use array_bound_module, only: ML_0,ML_1,MB,MB_0,MB_1,MBZ,MBZ_0,MBZ_1 &
                               ,MSP_0,MSP_1
  use wf_module, only: unk, occ, hunk
  use fock_fft_module, only: Fock_FFT_Double, Fock_FFT_DoubleComplex
  use parallel_module
  use watch_module
  use rsdft_mpi_module
  use hartree_mol_module, only: timer_reset_hartree_mol, timer_result_hartree_mol

  implicit none

  private
  public :: op_fock
  public :: UpdateWF_fock

  integer :: SYStype=0
  logical :: init_done = .false.

  interface op_fock
    module procedure d_op_fock, z_op_fock
  end interface

  real(8),allocatable :: occ_hf(:,:,:)
  real(8),allocatable :: d_unk_hf(:,:,:,:)
  complex(8),allocatable :: z_unk_hf(:,:,:,:)

contains

  subroutine init_fock
    use var_sys_parameter, only: use_real8_wf
    use io_ctrl_parameters, only: wf_available
    use io_read_wf_simple_module, only: read_wf_simple
    implicit none
    integer :: ng0,ng1,ng,ns0,ns1
    call get_range_parallel( ng0, ng1, 'g' ); ng=ng1-ng0+1
    call get_range_parallel( ns0, ns1, 's' )
    allocate( occ_hf(FKMB_0:FKMB_1,FKBZ_0:FKBZ_1,ns0:ns1) )
    occ_hf=0.0d0
    if ( use_real8_wf() ) then
      allocate( d_unk_hf(ng,FKMB_0:FKMB_1,FKBZ_0:FKBZ_1,ns0:ns1) )
      d_unk_hf=0.0d0
      if ( wf_available() ) call read_wf_simple( d_wf_out=d_unk_hf, occ_out=occ_hf )
    else
      allocate( z_unk_hf(ng,FKMB_0:FKMB_1,FKBZ_0:FKBZ_1,ns0:ns1) )
      z_unk_hf=(0.0d0,0.0d0)
      if ( wf_available() ) call read_wf_simple( z_wf_out=z_unk_hf, occ_out=occ_hf )
    end if
    init_done = .true.
  end subroutine init_fock


  subroutine d_op_fock( tpsi, hpsi, n,k,s )
    implicit none
    real(8),intent(in) :: tpsi(:,n:)
    real(8),intent(inout) :: hpsi(:,n:)
    integer,intent(in) :: n,k,s
    integer :: ib,ib1,ib2

    if ( iflag_hybrid == 0 ) return
    if ( .not.init_done ) call init_fock

    ib1 = n
    ib2 = ib1 + size(tpsi,2) - 1

    if ( iflag_hybrid == 2 ) then

      do ib = ib1, ib2
        hpsi(:,ib) = hpsi(:,ib) + hunk(:,ib,k,s)
      end do

    else if ( iflag_hybrid > 0 ) then

      do ib = ib1, ib2
        call Fock_Double( tpsi(:,ib), hpsi(:,ib), n,k,s )
      end do

    end if
  end subroutine d_op_fock

  subroutine z_op_fock( tpsi, hpsi, n,k,s )
    implicit none
    complex(8),intent(in) :: tpsi(:,n:)
    complex(8),intent(inout) :: hpsi(:,n:)
    integer,intent(in) :: n,k,s
    integer :: ib,ib1,ib2

    if ( iflag_hybrid == 0 ) return
    if ( .not.init_done ) call init_fock

    ib1 = n
    ib2 = ib1 + size(tpsi,2) - 1

    if ( iflag_hybrid == 2 ) then

      do ib = ib1, ib2
        hpsi(:,ib) = hpsi(:,ib) + hunk(:,ib,k,s)
      end do

    else if ( iflag_hybrid > 0 ) then

      do ib = ib1, ib2
        call Fock_DoubleComplex( tpsi(:,ib), hpsi(:,ib), n,k,s )
      end do

    end if
  end subroutine z_op_fock


  subroutine Fock_Double( psi, tpsi, n,k,s )
    implicit none
    real(8),intent(in) :: psi(:)
    real(8),intent(inout) :: tpsi(:)
    integer,intent(in) :: n,k,s
    complex(8),allocatable :: trho(:),tvht(:)
    real(8) :: c1,c2
    integer :: m1,m2,i,ng

    ng = size( psi, 1 )

    allocate( trho(ng) ); trho=(0.0d0,0.0d0)
    allocate( tvht(ng) ); tvht=(0.0d0,0.0d0)

    do m1 = FKMB_0, FKMB_1, 2

      m2 = min( m1+1, FKMB_1 )

      if ( abs(occ_hf(m1,k,s)) < 1.d-10 ) cycle

      if ( m1 < m2 ) then
        do i = 1, ng
          trho(i) = dcmplx( d_unk_hf(i,m1,k,s), d_unk_hf(i,m2,k,s) )*psi(i)
        end do
      else if ( m1 == m2 ) then
        do i = 1, ng
          trho(i) = d_unk_hf(i,m1,k,s)*psi(i)
        end do
      end if

      call Fock_FFT_Double( trho, tvht )

      c1 = alpha_hf*(2.0d0*occ_factor)*occ_hf(m1,k,s)
      c2 = alpha_hf*(2.0d0*occ_factor)*occ_hf(m2,k,s)

      if ( m1 < m2 ) then
        do i = 1, ng
          tpsi(i) = tpsi(i) - c1*real( tvht(i) )*d_unk_hf(i,m1,k,s) &
                            - c2*aimag( tvht(i) )*d_unk_hf(i,m2,k,s)
        end do
      else if ( m1 == m2 ) then
        do i = 1, ng
          tpsi(i) = tpsi(i) - c1*real( tvht(i) )*d_unk_hf(i,m1,k,s)
        end do
      end if

    end do ! m1

    deallocate( tvht, trho ) 

    return
  end subroutine Fock_Double

  subroutine Fock_DoubleComplex( psi, tpsi, n,k,s )
    use rsdft_mpi_module, only: rsdft_allreduce_sum
    implicit none
    complex(8),intent(in) :: psi(:)
    complex(8),intent(inout) :: tpsi(:)
    integer,intent(in) :: n,k,s
    complex(8),allocatable :: trho(:),tvht(:),tphi(:)
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    real(8) :: c
    integer :: m,q,t,ng

    ng = size( psi )

    allocate( trho(ng) ); trho=zero
    allocate( tvht(ng) ); tvht=zero
    allocate( tphi(ng) ); tphi=zero

! ---

    if ( gamma_hf == 0 ) then

      do t = FOCK_0, FOCK_1

        if ( t == 1 ) then

          do q = FKBZ_0, FKBZ_1
          do m = FKMB_0, FKMB_1

            if ( abs(occ_hf(m,q,s)) < 1.d-10 ) cycle

            trho(:) = conjg( z_unk_hf(:,m,q,s) )*psi(:)

            call Fock_FFT_DoubleComplex( trho, tvht, k,q,t )

            c = occ_factor*occ_hf(m,q,s)*alpha_hf

            tphi(:) = tphi(:) - c*tvht(:)*z_unk_hf(:,m,q,s)

          end do ! m
          end do ! q

        else ! [ t /= 1 ]

          do q = FKBZ_0, FKBZ_1
          do m = FKMB_0, FKMB_1

            if ( abs(occ_hf(m,q,s)) < 1.d-10 ) cycle

            trho(:) = z_unk_hf(:,m,q,s)*psi(:)

            call Fock_FFT_DoubleComplex( trho, tvht, k,q,t )

            c = occ_factor*occ_hf(m,q,s)*alpha_hf

            tphi(:) = tphi(:) - c*tvht(:)*conjg( z_unk_hf(:,m,q,s) )

          end do ! m
          end do ! q

        end if ! [ t ]

      end do ! t

    else ! [ gamma_hf /= 0 ]

      q = k

      do m = FKMB_0, FKMB_1

        if ( abs(occ_hf(m,q,s)) < 1.d-10 ) cycle

        trho(:) = conjg( z_unk_hf(:,m,q,s) )*psi(:)

        call Fock_FFT_DoubleComplex( trho, tvht, k,q,1 )

        c = alpha_hf*(2.0d0*occ_factor)*occ_hf(m,q,s)

        tphi(:) = tphi(:) - c*tvht(:)*z_unk_hf(:,m,q,s)

      end do ! m

    end if ! [ gamma_hf ]

    call rsdft_allreduce_sum( tphi, comm_fkmb )
    tpsi(:) = tpsi(:) + tphi(:)

    deallocate( tphi )
    deallocate( tvht ) 
    deallocate( trho )

    return
  end subroutine Fock_DoubleComplex


  subroutine UpdateWF_Fock( SYStype_in )
    implicit none
    integer,optional,intent(in) :: SYStype_in
    integer :: s,k,n,m,i_occ,i_orb,ierr
    type(time) :: tt
    logical :: disp_on

    if ( present(SYStype_in) ) SYStype = SYStype_in

    call write_border( 0, " UpdateWF_Fock(start)" )
    call start_timer( 'UpdateWF_Fock', tt )

    if ( .not.init_done ) call init_fock

    if ( allocated(d_unk_hf) ) then
      d_unk_hf(:,:,:,:) = 0.0d0
      occ_hf(:,:,:) = 0.0d0
      do s = MSP_0, MSP_1
      do k = MBZ_0, MBZ_1
      do n = MB_0 , MB_1
        d_unk_hf(:,n,k,s) = unk(:,n,k,s)
        occ_hf(n,k,s) = occ(n,k,s)
      end do !n
      end do !k
      end do !s
    else if ( allocated(z_unk_hf) ) then
      z_unk_hf(:,:,:,:) = (0.0d0,0.0d0)
      occ_hf(:,:,:) = 0.0d0
      do s = MSP_0, MSP_1
      do k = MBZ_0, MBZ_1
      do n = MB_0 , MB_1
        z_unk_hf(:,n,k,s) = unk(:,n,k,s)
        occ_hf(n,k,s) = occ(n,k,s)
      end do !n
      end do !k
      end do !s
    end if

    call write_border( 1, " UpdateWF_Fock(allreduce1)" )

    if ( allocated(d_unk_hf) ) then
      do s = MSP_0, MSP_1
      do k = MBZ_0, MBZ_1
        call rsdft_allreduce_sum( d_unk_hf(:,:,k,s), comm_band )
      end do ! k
      end do ! s
    else if ( allocated(z_unk_hf) ) then
      do s = MSP_0, MSP_1
      do k = MBZ_0, MBZ_1
        call rsdft_allreduce_sum( z_unk_hf(:,:,k,s), comm_band )
      end do ! k
      end do ! s
    end if

    call write_border( 1, " UpdateWF_Fock(allreduce2)" )

    if ( allocated(d_unk_hf) ) then
      do s = MSP_0, MSP_1
        call rsdft_allreduce_sum( d_unk_hf(:,:,:,s), comm_bzsm )
      end do ! s
    else if ( allocated(z_unk_hf) ) then
      do s = MSP_0, MSP_1
        call rsdft_allreduce_sum( z_unk_hf(:,:,:,s), comm_bzsm )
      end do ! s
    end if

    call write_border( 1, " UpdateWF_Fock(allreduce3)" )

    do s = MSP_0, MSP_1
    do k = MBZ_0, MBZ_1
      call rsdft_allreduce_sum( occ_hf(:,k,s), comm_band )
    end do ! k
    end do ! s

    call write_border( 1, " UpdateWF_Fock(allreduce4)" )

    do s = MSP_0, MSP_1
      call rsdft_allreduce_sum( occ_hf(:,:,s), comm_bzsm )
    end do ! s

    call end_timer( 'UpdateWF_Fock_allreduce' )

    hunk(:,:,:,:) = (0.0d0,0.0d0)
    ! ct_fock_fft(:)=0.0d0
    ! et_fock_fft(:)=0.0d0

    if ( allocated(d_unk_hf) ) then
      do s = MSP_0, MSP_1
        call Fock_4_Double( s )
      end do
    else if ( allocated(z_unk_hf) ) then
      do s = MSP_0, MSP_1
        if ( gamma_hf == 1 ) then
          call Fock_4( s )
        else
          call Fock_5( s )
        end if
      end do !s
    else
      call stop_program('Work-array for Hybdir-DFT (?_unk_hf) is not allocated')
    end if

    call end_timer( 'UpdateWF_fock_whole' )

    iflag_hybrid = 2

    ! call check_disp_switch( disp_on, 0 )
    ! if ( disp_on ) then
    !   write(*,'(1x,"iflag_hybrid=",i2)') iflag_hybrid
    !   write(*,*) "time(fock_fft1)=",ct_fock_fft(1),et_fock_fft(1)
    !   write(*,*) "time(fock_fft2)=",ct_fock_fft(2),et_fock_fft(2)
    !   write(*,*) "time(fock_fft3)=",ct_fock_fft(3),et_fock_fft(3)
    !   write(*,*) "time(fock_fft4)=",ct_fock_fft(4),et_fock_fft(4)
    !   write(*,*) "time(fock_fft5)=",ct_fock_fft(5),et_fock_fft(5)
    !   write(*,*) "time(fock_fft6)=",ct_fock_fft(6),et_fock_fft(6)
    !   write(*,*) "time(fock_fft7)=",ct_fock_fft(7),et_fock_fft(7)
    !   write(*,*) "time(fock_fft8)=",ct_fock_fft(8),et_fock_fft(8)
    ! end if

    call write_border( 0, " UpdateWF_fock(end)" )

  end subroutine UpdateWF_Fock


  subroutine Fock_4( s )
    implicit none
    integer,intent(in) :: s
    complex(8),allocatable :: trho(:),tvht(:)
    real(8) :: c
    integer :: m,n,i,j,k,nwork,iwork,ierr,nwork_0,nwork_1,ng
    integer,allocatable :: mapwork(:,:)

    call write_border( 1, " Fock_4(start)" )

    k = 1

! ---

    nwork=0
    i=0
    do n=1,MB
    do m=1,n
      if ( abs(occ(n,k,s))<1.d-10 .and. abs(occ(m,k,s))<1.d-10 ) cycle
      i=i+1
      j=mod(i-1,np_band)
      if ( j == myrank_b ) nwork=nwork+1
    end do
    end do

    allocate( mapwork(2,nwork) ); mapwork=0

    nwork=0
    i=0
    do n=1,MB
    do m=1,n
      if ( abs(occ(n,k,s))<1.d-10 .and. abs(occ(m,k,s))<1.d-10 ) cycle
      i=i+1
      j=mod(i-1,np_band)
      if ( j == myrank_b ) then
        nwork=nwork+1
        mapwork(1,nwork)=m
        mapwork(2,nwork)=n
      end if
    end do
    end do

    ir_fkmb(:)=0
    id_fkmb(:)=0
    do i=1,nwork
      j=mod(i-1,np_fkmb)
      ir_fkmb(j)=ir_fkmb(j)+1
    end do
    do j=0,np_fkmb-1
      id_fkmb(j)=sum(ir_fkmb(0:j))-ir_fkmb(j)
    end do

    nwork_0 = id_fkmb(myrank_f) + 1
    nwork_1 = id_fkmb(myrank_f) + ir_fkmb(myrank_f)

! ---

    ng = size( z_unk_hf, 1 )

    allocate( trho(ng) ); trho=(0.0d0,0.0d0)
    allocate( tvht(ng) ); tvht=(0.0d0,0.0d0)

! ---

    do iwork = nwork_0, nwork_1

      m = mapwork(1,iwork)
      n = mapwork(2,iwork)

      trho(:) = conjg( z_unk_hf(:,m,k,s) )*z_unk_hf(:,n,k,s)

      call Fock_FFT_DoubleComplex( trho, tvht, k,k,1 )

      if ( abs(occ(m,k,s)) >= 1.d-10 ) then
        c = alpha_hf*(2.0d0*occ_factor)*occ(m,k,s)
        hunk(:,n,k,s) = hunk(:,n,k,s) - c*tvht(:)*z_unk_hf(:,m,k,s)
      end if

      if ( m == n ) cycle

      if ( abs(occ(n,k,s)) >= 1.d-10 ) then
        c = alpha_hf*(2.0d0*occ_factor)*occ(n,k,s)
        hunk(:,m,k,s) = hunk(:,m,k,s) - c*conjg( tvht(:) )*z_unk_hf(:,n,k,s)
      end if

    end do !iwork

    deallocate( tvht ) 
    deallocate( trho )

    deallocate( mapwork )

! ---

    call rsdft_allreduce_sum( hunk(:,:,k,s), comm_band )
    call rsdft_allreduce_sum( hunk(:,:,k,s), comm_fkmb )

    call write_border( 1, " Fock_4(end)" )

    return
  end subroutine Fock_4

  subroutine Fock_4_Double( s )
    implicit none
    integer,intent(in) :: s
    complex(8),allocatable :: trho(:),tvht(:)
    real(8),parameter :: tol=1.d-10
    real(8) :: c
    integer :: m1,n1,m2,n2,m,n,i,j,k,nwork,iwork1,iwork2,ierr
    integer :: nwork_0,nwork_1,ng
    integer,allocatable :: mapwork(:,:)

    call write_border( 1, " Fock_4_Double(start)" )

    k = 1

! ---

    nwork=0
    i=0
    do n=1,MB
    do m=1,n
      if ( abs(occ(n,k,s)) < tol .and. abs(occ(m,k,s)) < tol ) cycle
      i=i+1
      j=mod(i-1,np_band)
      if ( j == myrank_b ) nwork=nwork+1
    end do
    end do

    allocate( mapwork(2,nwork) ); mapwork=0

    nwork=0
    i=0
    do n=1,MB
    do m=1,n
      if ( abs(occ(n,k,s)) < tol .and. abs(occ(m,k,s)) < tol ) cycle
      i=i+1
      j=mod(i-1,np_band)
      if ( j == myrank_b ) then
        nwork=nwork+1
        mapwork(1,nwork)=m
        mapwork(2,nwork)=n
      end if
    end do
    end do

    ir_fkmb(:)=0
    id_fkmb(:)=0
    do i=1,nwork
      j=mod(i-1,np_fkmb)
      ir_fkmb(j)=ir_fkmb(j)+1
    end do
    do j=0,np_fkmb-1
      id_fkmb(j)=sum(ir_fkmb(0:j))-ir_fkmb(j)
    end do

    nwork_0 = id_fkmb(myrank_f) + 1
    nwork_1 = id_fkmb(myrank_f) + ir_fkmb(myrank_f)

! ---

    ng = size( d_unk_hf, 1 )

    allocate( trho(ng) ); trho=(0.0d0,0.0d0)
    allocate( tvht(ng) ); tvht=(0.0d0,0.0d0)

! ---

    do iwork1 = nwork_0, nwork_1, 2

      iwork2 = min( iwork1+1, nwork_1 )

      m1 = mapwork(1,iwork1)
      n1 = mapwork(2,iwork1)
      m2 = mapwork(1,iwork2)
      n2 = mapwork(2,iwork2)

      if ( iwork1 < iwork2 ) then
        trho(:) = dcmplx( d_unk_hf(:,m1,k,s)*d_unk_hf(:,n1,k,s), &
                          d_unk_hf(:,m2,k,s)*d_unk_hf(:,n2,k,s) )
      else
        trho(:) = d_unk_hf(:,m1,k,s)*d_unk_hf(:,n1,k,s)
      end if

      call Fock_FFT_Double( trho, tvht )

      if ( abs(occ(m1,k,s)) >= tol ) then
        c = alpha_hf*(2.0d0*occ_factor)*occ(m1,k,s)
        hunk(:,n1,k,s) = hunk(:,n1,k,s) - c*real( tvht(:) )*d_unk_hf(:,m1,k,s)
      end if

      if ( m1 /= n1 .and. abs(occ(n1,k,s)) >= tol ) then
        c = alpha_hf*(2.0d0*occ_factor)*occ(n1,k,s)
        hunk(:,m1,k,s) = hunk(:,m1,k,s) - c*real( tvht(:) )*d_unk_hf(:,n1,k,s)
      end if

      if ( iwork1 < iwork2 ) then

        if ( abs(occ(m2,k,s)) >= tol ) then
          c = alpha_hf*(2.0d0*occ_factor)*occ(m2,k,s)
          hunk(:,n2,k,s) = hunk(:,n2,k,s) - c*aimag( tvht(:) )*d_unk_hf(:,m2,k,s)
        end if

        if ( m2 /= n2 .and. abs(occ(n2,k,s)) >= tol ) then
          c = alpha_hf*(2.0d0*occ_factor)*occ(n2,k,s)
          hunk(:,m2,k,s) = hunk(:,m2,k,s) - c*aimag( tvht(:) )*d_unk_hf(:,n2,k,s)
        end if

      end if

    end do ! iwork

    deallocate( tvht ) 
    deallocate( trho )

    deallocate( mapwork )

! ---

    call rsdft_allreduce_sum( hunk(:,:,k,s), comm_band )
    call rsdft_allreduce_sum( hunk(:,:,k,s), comm_fkmb )

    call write_border( 1, " Fock_4_Double(end)" )

    return
  end subroutine Fock_4_Double


  subroutine Fock_5( s )
    implicit none
    integer,intent(in) :: s
    complex(8),allocatable :: trho(:),tvht(:)
    real(8) :: c,ctt(0:3),ett(0:3)
    integer :: m,n,q,k,i,j,a,b,nwork,iwork,ierr,nwork_0,nwork_1
    integer :: ng
    integer,allocatable :: mapnk(:,:),mapwork(:,:)

    call write_border( 1, " Fock_5(start)" )

    call watch(ctt(0),ett(0))
  
! ---

    allocate( mapnk(2,MB*MBZ) ); mapnk=0

    i=0
    do k=1,MBZ
    do n=1,MB
      i=i+1
      mapnk(1,i)=n
      mapnk(2,i)=k
    end do
    end do

! ---

    nwork=0
    a=0
    do j=1,MB*MBZ
    do i=1,j
      m=mapnk(1,i)
      q=mapnk(2,i)
      n=mapnk(1,j)
      k=mapnk(2,j)
      if ( abs(occ_hf(n,k,s))<1.d-10 .and. abs(occ_hf(m,q,s))<1.d-10 ) cycle
      a=a+1
      b=mod(a-1,np_band*np_bzsm)
      if ( b == myrank_k+myrank_b*np_bzsm ) nwork=nwork+1
    end do
    end do

!    if ( disp_switch_parallel ) then
!       write(*,*) "total # of me    =",a
!       write(*,*) "# of me of rank0 =",nwork
!    end if

    allocate( mapwork(2,nwork) ); mapwork=0

    nwork=0
    a=0
    do j=1,MB*MBZ
    do i=1,j
      m=mapnk(1,i)
      q=mapnk(2,i)
      n=mapnk(1,j)
      k=mapnk(2,j)
      if ( abs(occ_hf(n,k,s))<1.d-10 .and. abs(occ_hf(m,q,s))<1.d-10 ) cycle
      a=a+1
      b=mod(a-1,np_band*np_bzsm)
      if ( b == myrank_k+myrank_b*np_bzsm ) then
        nwork=nwork+1
        mapwork(1,nwork)=i
        mapwork(2,nwork)=j
      end if
    end do ! i
    end do ! j

    ir_fkmb(:)=0
    id_fkmb(:)=0
    do i=1,nwork
      j=mod(i-1,np_fkmb)
      ir_fkmb(j)=ir_fkmb(j)+1
    end do
    do j=0,np_fkmb-1
      id_fkmb(j)=sum(ir_fkmb(0:j))-ir_fkmb(j)
    end do

    nwork_0 = id_fkmb(myrank_f) + 1
    nwork_1 = id_fkmb(myrank_f) + ir_fkmb(myrank_f)

! ---

    ng = size( z_unk_hf, 1 )

    allocate( trho(ng) ); trho=(0.0d0,0.0d0)
    allocate( tvht(ng) ); tvht=(0.0d0,0.0d0)

    call watch(ctt(1),ett(1))

! ---

    do iwork = nwork_0, nwork_1

      a = mapwork(1,iwork)
      b = mapwork(2,iwork)
      m = mapnk(1,a)
      q = mapnk(2,a)
      n = mapnk(1,b)
      k = mapnk(2,b)

! - normal part -

      trho(:) = conjg( z_unk_hf(:,m,q,s) )*z_unk_hf(:,n,k,s)

      call Fock_FFT_DoubleComplex( trho, tvht, k,q,1 )

      if ( abs(occ_hf(m,q,s)) >= 1.d-10 ) then
        c = alpha_hf*occ_factor*occ_hf(m,q,s)
        hunk(:,n,k,s) = hunk(:,n,k,s) - c*tvht(:)*z_unk_hf(:,m,q,s)
      end if

      if ( a /= b .and. abs(occ_hf(n,k,s)) >= 1.d-10 ) then
        c = alpha_hf*occ_factor*occ_hf(n,k,s)
        hunk(:,m,q,s) = hunk(:,m,q,s) - c*conjg( tvht(:) )*z_unk_hf(:,n,k,s)
      end if

! - time-reversal part -

      trho(:) = z_unk_hf(:,m,q,s)*z_unk_hf(:,n,k,s)

      call Fock_FFT_DoubleComplex( trho, tvht, k,q,2 )

      if ( abs(occ_hf(m,q,s)) >= 1.d-10 ) then
        c = alpha_hf*occ_factor*occ_hf(m,q,s)
        hunk(:,n,k,s) = hunk(:,n,k,s) - c*tvht(:)*conjg( z_unk_hf(:,m,q,s) )
      end if

      if ( a /= b .and. abs(occ_hf(n,k,s)) >= 1.d-10 ) then
        c = alpha_hf*occ_factor*occ_hf(n,k,s)
        hunk(:,m,q,s) = hunk(:,m,q,s) - c*tvht(:)*conjg( z_unk_hf(:,n,k,s) )
      end if

    end do ! iwork

    call watch(ctt(2),ett(2))

! ---

    call write_border( 1, " Fock_5(allreduce1)" )

    m=size( hunk,1 )*size( hunk,2 )
    do k=1,MBZ
      call rsdft_allreduce_sum( hunk(:,:,k,s), comm_band )
    end do

    call write_border( 1, " Fock_5(allreduce2)" )

    m=size( hunk,1 )*size( hunk,2 )*size( hunk,3 )
    call rsdft_allreduce_sum( hunk(:,:,:,s), comm_bzsm )

    call write_border( 1, " Fock_5(allreduce3)" )

    m=size( hunk,1 )*size( hunk,2 )*size( hunk,3 )
    call rsdft_allreduce_sum( hunk(:,:,:,s), comm_fkmb )

    call watch(ctt(3),ett(3))

    ! ct_focK_fft(6) = ctt(1) - ctt(0)
    ! et_focK_fft(6) = ett(1) - ett(0)
    ! ct_focK_fft(7) = ctt(2) - ctt(1)
    ! et_focK_fft(7) = ett(2) - ett(1)
    ! ct_focK_fft(8) = ctt(3) - ctt(2)
    ! et_focK_fft(8) = ett(3) - ett(2)

! ---

    deallocate( tvht ) 
    deallocate( trho )

    deallocate( mapwork )
    deallocate( mapnk )

! ---

    call write_border( 1, " Fock_5(end)" )

    return
  end subroutine Fock_5

end module fock_module
