module atomopt_ef_module

  use lattice_module
  use atom_module, only: atom, construct_atom, aa_atom &
                        ,write_coordinates_atom, shift_aa_coordinates_atom &
                        ,write_xyz_atom
  use force_module
  use pseudopot_module, only: pselect
  use eion_module
  use ps_local_module
  use ps_local_mol_module
  use ps_pcc_mol_module
  use ps_nloc2_mol_module
  use strfac_module
  use ps_pcc_module
  use ps_nloc3_module
  use ps_nloc2_module
  use ps_nloc_mr_module
  use scf_module
  use atomopt_diis_module, only: calc_coef_diis, calc_coef_diis_b
  use atomopt_io_module, only: flag_continue_atomopt &
                             , read_atomopt_io, write_atomopt_io
  use density_module, only: rho, normalize_density
  use ps_initrho_module, only: construct_ps_initrho
  use hartree_module, only: calc_hartree
  use xc_module, only: calc_xc

  implicit none

  private
  public :: atomopt_ef

  logical :: disp_sw, disp_scf
  integer :: SYStype
  integer :: NiterSCF

contains


  subroutine atomopt_ef( SYStype_in, fmax_tol, ncycle, okstep_in, max_linmin_in, max_alpha_in, NiterSCF_in )

    implicit none
    integer,intent(in) :: SYStype_in
    real(8),intent(in) :: fmax_tol
    integer,intent(in) :: ncycle
    real(8),intent(in) :: okstep_in
    integer,intent(in) :: max_linmin_in
    real(8),intent(in) :: max_alpha_in
    integer,optional,intent(in) :: NiterSCF_in
    type(atom) :: ion
    type(lattice) :: aa, bb
    integer :: max_loop, max_linmin, linmin, linmin_start
    integer :: ishape(1), ishape2(2), i1,i2, LWORK=0
    integer :: n,il,iu,m,a,ierr,loop,i,j,ip,loop_start
    integer,allocatable :: iwork(:),ifail(:)
    real(8) :: etot0, etot, fmax
    real(8) :: dxdg,dxHdx,c1,c2,c,dmax,dtmp
    real(8) :: vl,vu,tol,alpha,okstep,max_alpha
    real(8) :: aa_inv(3,3),da(3),da_tmp(3)
    real(8),parameter :: one=1.0d0, zero=0.0d0
    real(8),allocatable :: g(:),x(:),Hessian(:,:),Htmp(:,:)
    real(8),allocatable :: Vtmp(:,:),Wtmp(:,:)
    real(8),allocatable :: history(:,:),g0(:),x0(:),d(:)
    real(8),allocatable :: dx(:),dg(:),Hdx(:)
    real(8),allocatable :: w(:),z(:,:),work(:)
    real(8),allocatable :: info_linmin(:), xyzf_linmin(:,:,:)
    real(8),parameter :: eig_threshold=0.02d0

    integer :: ndiis, mdiis, nrhodiis, mrhodiis
    real(8),allocatable :: xdiis(:,:), ediis(:,:), coef_diis(:)
    real(8),allocatable :: xtmp(:),etmp(:),gtmp(:),gdiis(:,:)
    real(8),allocatable :: rhodiis(:,:,:)
    real(8),allocatable :: rho_ini(:,:), rho_ini_old(:,:)
    logical :: flag_diis, flag_sd , flag_forced_sd, flag_3pts
    logical :: flag_exit, flag_check_linmin_continue
    real(8) :: alpha0,alpha1,alpha2,alpha3,alpha4
    real(8) :: etot1,etot2,etot3,etot4
    character(40) :: msg
    real(8) :: factor, tau, fmax_min, etot_min
    integer :: icount, istep_golden_search
    integer :: myrank
    include 'mpif.h'

    call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr )

    call write_border( 0, "atomopt_ef(start)" )
    call check_disp_switch( disp_sw, 0 )
    disp_scf = disp_sw
    call check_disp_length( i, 0 )
    if ( i < 2 ) then
       call check_disp_switch( .false., 1 )
       disp_scf = .false.
    end if

! ---

    SYStype = SYStype_in

    NiterSCF = 50 ; if ( present(NiterSCF_in) ) NiterSCF=NiterSCF_in

    max_loop = ncycle

    loop_start = 1

    okstep = okstep_in

    ndiis = 0
    if ( disp_sw .and. ndiis > 0 ) write(*,*) "DIIS is used: ndiis=",ndiis
    mdiis = 0

    nrhodiis = nint(max_alpha_in)
    mrhodiis = 0
    if ( disp_sw ) write(*,*) "nrhodiis=",nrhodiis

    max_linmin = max_linmin_in

    max_alpha = 0.0d0 !max_alpha_in

    linmin_start = 1

    flag_check_linmin_continue = .false.

! ---

    call get_aa_lattice( aa )
    call get_inverse_lattice( aa%LatticeVector, aa_inv )

    call construct_atom( ion )

    do a=1,ion%natom
       ion%xyz(1:3,a) = matmul( aa%LatticeVector(1:3,1:3), ion%aaa(1:3,a) )
    end do

    if ( .not.flag_continue_atomopt() ) then

      if ( disp_sw ) then
        write(*,*) "Initial configuration"
        do a=1,ion%natom
          write(*,'(1x,3f20.15)') ion%aaa(1:3,a)
        end do
      end if

      call scf( etot, ierr ); if ( ierr  == -1 ) goto 999
      call calc_force( ion%natom, ion%force, fmax )

      if ( fmax <= fmax_tol ) goto 900

   end if

! ---

    n = 3*ion%natom
    allocate( g(n)  ) ; g=0.0d0
    allocate( x(n)  ) ; x=0.0d0
    allocate( d(n)  ) ; d=0.0d0
    allocate( g0(n) ) ; g0=0.0d0
    allocate( x0(n) ) ; x0=0.0d0
    allocate( Hessian(n,n) ) ; Hessian=0.0d0
    allocate( Htmp(n,n) ) ; Htmp=0.0d0
    allocate( Vtmp(n,n) ) ; Vtmp=0.0d0
    allocate( Hdx(n)    ) ; Hdx=0.0d0
    allocate( dx(n)     ) ; dx=0.0d0
    allocate( dg(n)     ) ; dg=0.0d0

    ishape(1) = n
    ishape2(1:2) = (/ 3, ion%natom /)

    allocate( history(6,0:max_loop) ) ; history=0.0d0

    allocate( info_linmin(size(history,1)) ); info_linmin=0.0d0
    allocate( xyzf_linmin(3,ion%natom,2) ); xyzf_linmin=0.0d0

    if ( flag_continue_atomopt() ) then

      call read_atomopt_io( &
           loop_start, &
           history, &
           x, x0, g, g0, &
           Hessian )

      flag_check_linmin_continue = .true.

    else

      history(1,0) = etot
      history(2,0) = fmax
      history(3,0) = ierr ; if (ierr == -2) history(3,0)=NiterSCF
      history(4,0) = sum( history(3,0:0) )
      history(5:,0) = 0.0d0

    end if

! --- DIIS

    if ( ndiis > 0 ) then
      allocate( xdiis(n,ndiis)   ); xdiis=0.0d0
      allocate( ediis(n,ndiis)   ); ediis=0.0d0
      allocate( gdiis(n,ndiis)   ); gdiis=0.0d0
      allocate( coef_diis(ndiis) ); coef_diis=0.0d0
      allocate( xtmp(n) ); xtmp=0.0d0
      allocate( etmp(n) ); etmp=0.0d0
      allocate( gtmp(n) ); gtmp=0.0d0
    end if

    if ( nrhodiis > 0 ) then
      allocate( xdiis(n,nrhodiis) ); xdiis=0.0d0
      allocate( coef_diis(nrhodiis) ); coef_diis=0.0d0
      allocate( rhodiis(size(rho,1),size(rho,2),nrhodiis) ); rhodiis=0.0d0
      allocate( rho_ini(size(rho,1),size(rho,2)) ); rho_ini=0.0d0
      allocate( rho_ini_old(size(rho,1),size(rho,2)) ); rho_ini_old=0.0d0
    end if

    !if ( nrhodiis > 0 ) then
    !  if ( .not.flag_continue_atomopt() ) then
    !    xdiis(:,1)=reshape( ion%xyz, ishape )
    !    rhodiis(:,:,1)=rho(:,:)
    !    mrhodiis=1
    !  end if
    !end if

! --- LAPACK

    if ( LWORK == 0 ) then
      allocate( w(n) ); w=0.0d0
      allocate( work(1) )
      call dsyev('V','U',n+1,Htmp,n+1,w,work,-1,ierr)
      LWORK=nint( work(1) )
      deallocate( work )
      allocate( work(LWORK) ); work=0.0d0
    end if

! ----------------------------------- loop


    do loop = loop_start, max_loop

      if ( disp_sw ) write(*,'(a60," loop=",i4)') repeat("-",60),loop

      alpha = 1.0d0
      flag_diis = .false.
      flag_sd = .false.
      flag_forced_sd = .false.
      linmin_start = 1

! ---

      call write_coordinates_atom( 197, 3 )
      call write_atomopt_io( loop, history(:,0:loop-1), &
                             x, x0, g, g0, Hessian )
! ---

      if ( loop == 1 ) then

        flag_sd = .true.

        Hessian=0.0d0
        do i=1,n
          Hessian(i,i) = one
        end do

        x(:) = reshape( ion%xyz(:,:)  , ishape )
        g(:) = reshape(-ion%force(:,:), ishape )
        d(:) = -g(:)

      else

        dx(:) = x(:) - x0(:)
        dg(:) = g(:) - g0(:)

        dxdg = sum( dx*dg )
        if ( disp_sw ) write(*,*) "dxdg=",dxdg

! ---

        if ( .false. ) then
        !if ( loop >= 3 ) then
          etot1 = history(1,loop-1)
          etot2 = history(1,loop-2)
          etot3 = history(1,loop-3)
          if ( etot1 > etot2 .and. etot2 > etot3 ) then
            if ( disp_sw ) write(*,*) "Climbing PES (SD is applied)"
            flag_forced_sd = .true.
          end if
          open(297,file='trajectory_ef.xyz',status='old')
          do j=1,loop-4
            read(297,*)
            read(297,*)
            do i=1,size(ion%xyz,2)
              read(297,*)
            end do
            read(297,*)
            do i=1,size(ion%xyz,2)
              read(297,*)
            end do
          end do
          read(297,*) n
          read(297,*) msg
          if ( disp_sw ) then
            write(*,*) "label, natom=",msg, n
            write(*,*) "etot(loop-3)=",history(1,loop-3)
            write(*,*) "etot(loop-2)=",etot2
            write(*,*) "etot(loop-1)=",etot1
          end if
          do i=1,size(ion%xyz,2)
            read(297,*) ion%xyz(:,i)
          end do
          read(297,*)
          do i=1,size(ion%force,2)
            read(297,*) ion%force(:,i)
          end do
          close(297)
          do i=1,size(ion%xyz,2)
            i2 = 3*i
            i1 = i2 - 3 + 1
            x(i1:i2) = ion%xyz(:,i)
            g(i1:i2) =-ion%force(:,i)
          end do
        end if

! ---

        if ( dxdg <= 0.0d0 .or. flag_forced_sd ) then

          if ( disp_sw ) write(*,*) "Steepest descent (Restart EF-opt process)"

          flag_sd = .true.
          flag_forced_sd = .false.

          d(:)  =-g(:)
          x0(:) = x(:)
          g0(:) = g(:)

          Hessian=0.0d0
          do i=1,size(Hessian,1)
            Hessian(i,i)=1.0d0
          end do

        else

          flag_sd = .false.

          call dgemv( 'N', n, n, one, Hessian, n, dx, 1, zero, Hdx, 1 )
          dxHdx = sum( dx*Hdx )
          c1 = 1.0d0/dxdg
          c2 = 1.0d0/dxHdx
          do j=1,n
          do i=1,n
            Hessian(i,j) = Hessian(i,j) + dg(i)*dg(j)*c1 - Hdx(i)*Hdx(j)*c2
          end do
          end do

          Htmp(:,:) = Hessian(:,:)
          call dsyev('V','U',n,Htmp,n,w,work,size(work),ierr)

          do i=1,n
            if ( w(i) < eig_threshold ) then
              if ( disp_sw ) write(*,*) "i,w(i)",i,w(i),"  ---> replaced to",eig_threshold
              w(i) = eig_threshold
            else
              exit
            end if
          end do

          Vtmp=0.0d0
          do i=1,n
            Vtmp(i,i)=1.0d0/w(i)
          end do
          Vtmp=matmul( Vtmp,transpose(Htmp) )
          Vtmp=matmul( Htmp,Vtmp )

! --- DIIS

          if ( ndiis > 0 ) then
            mdiis = min( mdiis+1, ndiis )
            do i=mdiis,2,-1
              xdiis(:,i) = xdiis(:,i-1)
              ediis(:,i) = ediis(:,i-1)
              gdiis(:,i) = gdiis(:,i-1)
            end do
            xdiis(:,1) = x(:)
            ediis(:,1) = x(:)-x0(:)
            gdiis(:,1) = g(:)
            if ( mdiis > 1 .or. mdiis == ndiis ) then
              call calc_coef_diis( coef_diis(1:mdiis), ediis(:,1:mdiis) )
              xtmp = matmul( xdiis(:,1:mdiis), coef_diis(1:mdiis) )
              etmp = matmul( ediis(:,1:mdiis), coef_diis(1:mdiis) )
              gtmp = matmul( gdiis(:,1:mdiis), coef_diis(1:mdiis) )
              flag_diis = .true.
            end if
            !dmax=0.0d0
            !do a=1,ion%natom
            !  i2 = 3*a
            !  i1 = i2 - 3 + 1
            !  dtmp = sqrt(sum((xtmp(i1:i2)-x0(i1:i2))**2))
            !  dmax = max(dmax,dtmp)
            !end do
            !write(*,*) "Maximum displacement(DIIS)",dmax
          end if

          if ( flag_diis ) then
            call DGEMV('N',n,n,-1.0d0,Vtmp,n,gtmp(:),1,0.0d0,d,1)
          else
            call DGEMV('N',n,n,-1.0d0,Vtmp,n,g(:),1,0.0d0,d,1)
          end if

        end if

      end if

! ---

      x0(:) = x(:)
      g0(:) = g(:)

! ---

      do a=1,ion%natom
        i2 = a*3
        i1 = i2 - 3 + 1 
        da = matmul( aa_inv, d(i1:i2) )
        do i=1,3
          da_tmp(1) = abs(da(i))
          da_tmp(2) = abs(da(i)+1.0d0)
          da_tmp(3) = abs(da(i)-1.0d0)
          j = minloc( da_tmp, 1 )
          select case(j)
          case(1)
          case(2); da(i)=da(i)+1.0d0
          case(3); da(i)=da(i)-1.0d0
          end select
        end do
        d(i1:i2) = matmul( aa%LatticeVector, da )
      end do

! ---

      tau = ( sqrt(5.0d0) - 1.0d0 )*0.5d0

      if ( flag_check_linmin_continue ) then
        call read_linmin_ef( linmin_start )
        flag_check_linmin_continue=.false.
        if ( disp_sw ) write(*,*) "linmin_start=",linmin_start
      end if

      if ( linmin_start == 1 ) then
        info_linmin=1.0d100
        xyzf_linmin=0.0d0
        alpha1=0.0d0; alpha2=0.0d0; alpha3=0.0d0; alpha4=0.0d0
        etot1 =0.0d0; etot2 =0.0d0; etot3 =0.0d0; etot4 =0.0d0
        flag_3pts = .false.
        flag_exit = .false.
        istep_golden_search=0
      else
        call read_linmin_ef
      end if

      do linmin = linmin_start, max_linmin+1 !---------- Line Minimization

        if ( myrank == 0 ) call write_linmin_ef( linmin )

        if ( .not.flag_3pts ) then

          if ( linmin == 1 ) then
            alpha1 = 0.0d0
            etot1  = etot
            alpha  = 1.0d0
          else if ( linmin == 2 ) then
            if ( etot1 < etot ) then
              factor = 0.3d0
              alpha3 = alpha
              etot3  = etot
              alpha  = alpha*factor
            else if ( etot1 >= etot ) then
              factor = 2.7d0
              alpha2 = alpha
              etot2  = etot
              alpha  = alpha*factor
            end if
          else if ( linmin >= 3 ) then
            if ( factor == 0.3d0 ) then
              if ( etot1 > etot ) then
                flag_3pts = .true.
                alpha1 = alpha1; etot1 = etot1
                alpha4 = alpha3; etot4 = etot3
                alpha2 = 0.0d0 ; etot2 = 0.0d0
                alpha3 = 0.0d0 ; etot3 = 0.0d0
              else
                alpha3 = alpha
                etot3  = etot
                alpha  = alpha*factor
              end if
            else if ( factor == 2.7d0 ) then
              if ( etot2 < etot ) then
                flag_3pts = .true.
                alpha1 = alpha1 ; etot1 = etot1
                alpha4 = alpha  ; etot4 = etot
                alpha2 = 0.0d0  ; etot2 = 0.0d0
                alpha3 = 0.0d0  ; etot3 = 0.0d0
              else
                alpha1 = alpha2
                etot1  = etot2
                alpha2 = alpha
                etot2  = etot
                alpha  = alpha*factor
              end if
            end if
          end if

        end if !.not.flag_3pts

        if ( flag_3pts ) then

          select case( istep_golden_search )
          case( 0 )
            alpha = alpha4 - tau*(alpha4-alpha1)
            istep_golden_search = istep_golden_search + 1
          case( 1 )
            alpha2 = alpha
            etot2  = etot
            alpha  = alpha1 + tau*(alpha4-alpha1)
            istep_golden_search = istep_golden_search + 1
          case( 2: )
            if ( alpha3 == 0.0d0 .and. etot3 == 0.0d0 ) then
              alpha3 = alpha
              etot3  = etot
            else
              alpha2 = alpha
              etot2  = etot
            end if
            if ( etot1 >= etot2 .and. etot2 <= etot3 ) then
              alpha4 = alpha3
              etot4  = etot3
              alpha3 = alpha2
              etot3  = etot2
              alpha2 = 0.0d0
              etot2  = 0.0d0
              alpha  = alpha4 - tau*(alpha4-alpha1)
            else if ( etot2 >= etot3 .and. etot3 <= etot4 ) then
              alpha1 = alpha2
              etot1  = etot2
              alpha2 = alpha3
              etot2  = etot3
              alpha3 = 0.0d0
              etot3  = 0.0d0
              alpha  = alpha1 + tau*(alpha4-alpha1)
            end if
            istep_golden_search = istep_golden_search + 1
          end select

        end if !flag_3pts

        if ( flag_sd .and. disp_sw ) then
          write(*,'(a10,"linmin=",i4)') repeat("-",10),linmin-1
          write(*,*) alpha1,etot1
          write(*,*) alpha2,etot2
          write(*,*) alpha3,etot3
          if ( flag_3pts ) write(*,*) alpha4,etot4
        end if

        if ( flag_exit .or. linmin==max_linmin ) exit

        dmax=0.0d0
        do a=1,ion%natom
          i2 = 3*a
          i1 = i2 - 3 + 1
          dtmp = sqrt(sum(d(i1:i2)**2))*alpha
          dmax = max(dmax,dtmp)
        end do
        !if ( disp_sw ) then
        !  write(*,*) "Maximum displacement(bohr):",dmax
        !end if
        if ( dmax > okstep ) then
          alpha=alpha*okstep/dmax
          !if ( disp_sw ) then
          !  write(*,*) "Maxmimum displacement is limited to",okstep
          !  write(*,*) "alpha is changed: alpha=",alpha
          !end if
        end if

! ---

        if ( flag_diis ) then
          do a=1,ion%natom
            i2 = 3*a
            i1 = i2 - 3 + 1
            ion%xyz(:,a) = xtmp(i1:i2) + alpha*d(i1:i2)
          end do
        else
          do a=1,ion%natom
            i2 = 3*a
            i1 = i2 - 3 + 1
            ion%xyz(:,a) = x0(i1:i2) + alpha*d(i1:i2)
          end do
        end if

! ---

        if ( nrhodiis > 0 ) then
          mrhodiis = mrhodiis + 1
          if ( mrhodiis > nrhodiis ) then
            mrhodiis = nrhodiis
            do i = 1, nrhodiis-1
              xdiis(:,i)=xdiis(:,i+1)
              rhodiis(:,:,i)=rhodiis(:,:,i+1)
            end do
          end if

          xdiis(:,mrhodiis)=reshape(ion%xyz,ishape)
do i=1,mrhodiis
if(disp_sw)write(*,*) i,sum(xdiis(:,i)**2),sum(rhodiis(:,:,i)**2)
end do
          if ( mrhodiis >= 3 ) then

            aa_atom(:,:) = matmul( aa_inv, ion%xyz )
            call construct_ps_initrho( rho_ini )

            call calc_coef_diis_b( coef_diis(1:mrhodiis-1), xdiis(:,1:mrhodiis) )

            rho=0.0d0
            do i=1,mrhodiis-1
              do a=1,ion%natom
                i2=3*a
                i1=i2-2
                aa_atom(:,a) = matmul( aa_inv, xdiis(i1:i2,i) )
              end do
              call construct_ps_initrho( rho_ini_old )
              rho(:,:) = rho(:,:) + coef_diis(i)*(rhodiis(:,:,i)-rho_ini_old(:,:))
            end do
            rho(:,:) = rho(:,:) + rho_ini(:,:)
            where( rho < 0.0d0 )
              rho=0.0d0
            end where
            call normalize_density( rho )
            call calc_hartree(lbound(rho,1),ubound(rho,1),size(rho,2),rho)
            call calc_xc
          end if
        end if

! ---

        aa_atom(:,:) = matmul( aa_inv, ion%xyz )
        call shift_aa_coordinates_atom( aa_atom )        

! ---

        !if ( disp_sw ) then
        !  write(*,*) "Next trial configuration"
        !end if
        dmax=0.0d0
        do a=1,ion%natom
          i2 = 3*a
          i1 = i2 - 3 + 1
          dtmp = sqrt(sum(d(i1:i2)**2))*alpha
          !if ( disp_sw ) then
          !  write(*,'(i4,2x,3f10.5,2x,3es14.5,2x,es14.5)') &
          !       a,aa_atom(:,a),ion%xyz(:,a),dtmp
          !end if
          dmax=max(dmax,dtmp)
        end do
        !if ( disp_sw ) then
        !  write(*,*) "Maximum displacement(bohr):",dmax
        !end if

        call write_coordinates_atom( 97, 3 )

        call scf( etot, ierr ); if ( ierr == -1 ) goto 999
        call calc_force( ion%natom, ion%force, fmax )

        write(msg,*) loop, etot
        call write_trajectory( ion%xyz, ion%force, msg )

        if ( nrhodiis > 0 ) rhodiis(:,:,mrhodiis)=rho(:,:)

        history(1,loop) = etot
        history(2,loop) = fmax
        icount=ierr
        if ( ierr == -1 ) then
          icount=0
        else if ( ierr == -2 ) then
          icount=NiterSCF
        end if
        history(3,loop) = history(3,loop) + icount
        history(4,loop) = sum( history(3,0:loop) )
        history(5,loop) = alpha
        history(6,loop) = sum(d*g)

        if ( disp_sw ) then
          do i=loop,loop
            write(*,'("linmin ",i4,f20.10,es14.5,i4,i6,2es14.5)') &
                 i, history(1:2,i), icount, nint(history(4,i)),history(5:6,i)
          end do
        end if

! ---

        flag_exit=.false.
        etot_min = info_linmin(1)
        fmax_min = info_linmin(2)
        if ( abs(etot-etot_min) < 1.d-7 ) then
          if ( disp_sw ) write(*,'("conv(etot)",3g20.10)') etot_min, etot, etot-etot_min
          flag_exit = .true.
        end if
        if ( abs(fmax-fmax_min) < fmax_tol ) then
          if ( disp_sw ) write(*,'("conv(fmax)",3g20.10)') fmax_min, fmax, fmax-fmax_min
          flag_exit=.true.
        end if

        if ( etot < info_linmin(1) ) then
          info_linmin(1:6) = history(1:6,loop)
          xyzf_linmin(:,:,1) = ion%xyz
          xyzf_linmin(:,:,2) = ion%force
        end if

        if ( .not.flag_sd .or. fmax <= fmax_tol ) flag_exit=.true.

      end do !linmin

      history(:,loop) = info_linmin(:)
      fmax = info_linmin(2)
      ion%xyz = xyzf_linmin(:,:,1)
      ion%force = xyzf_linmin(:,:,2)

      if ( disp_sw ) then
        do i=0,loop
          write(*,'("History",i4,f20.10,es14.5,i4,i6,2es14.5)') &
               i, history(1:2,i), nint(history(3:4,i)),history(5:6,i)
        end do
      end if

      if ( fmax <= fmax_tol ) exit

      x(:) =  reshape( ion%xyz(:,:)  , ishape )
      g(:) = -reshape( ion%force(:,:), ishape )

      linmin_start = 1

    end do ! loop

! ---

900 if ( disp_sw ) then
       write(*,'(1x,"etot    :",f25.15)') etot
       write(*,'(1x,"fmax/tol:",es12.5," /",es12.5)') fmax,fmax_tol
    end if

999 call check_disp_switch( disp_sw, 1 )
    call write_border( 0, "atomopt_ef(end)" )

    if ( allocated(xdiis) ) deallocate(xdiis)
    if ( allocated(coef_diis) ) deallocate(coef_diis)
    if ( allocated(rhodiis) ) deallocate(rhodiis)
    deallocate( work )
    deallocate( w )
    deallocate( history )
    deallocate( dg )
    deallocate( dx )
    deallocate( Hdx )
    deallocate( Htmp )
    deallocate( Hessian )
    deallocate( x0 )
    deallocate( g0 )
    deallocate( x )
    deallocate( g )

  contains

    subroutine write_linmin_ef(linmin)
      implicit none
      integer,intent(in) :: linmin
      integer,parameter :: u=298
      open(u,file='linmin_ef.dat',form='unformatted')
      write(u) linmin
      write(u) info_linmin
      write(u) xyzf_linmin
      write(u) alpha1,alpha2,alpha3,alpha4,alpha
      write(u) etot1,etot2,etot3,etot4,etot
      write(u) factor
      write(u) flag_3pts
      write(u) flag_exit
      write(u) istep_golden_search
      write(u) flag_sd
      close(u)
    end subroutine write_linmin_ef

    subroutine read_linmin_ef( linmin_start )
      implicit none
      integer,optional,intent(inout) :: linmin_start
      integer,parameter :: u=298
      integer :: linmin
      logical :: flag
      inquire( FILE='linmin_ef.dat', EXIST=flag )
      if ( .not.flag ) return
      if ( present(linmin_start) ) then
        open(u,file='linmin_ef.dat',form='unformatted')
        read(u) linmin_start
        close(u)
        return
      end if
      open(u,file='linmin_ef.dat',form='unformatted')
      read(u) linmin
      read(u) info_linmin
      read(u) xyzf_linmin
      read(u) alpha1,alpha2,alpha3,alpha4,alpha
      read(u) etot1,etot2,etot3,etot4,etot
      read(u) factor
      read(u) flag_3pts
      read(u) flag_exit
      read(u) istep_golden_search
      read(u) flag_sd
      close(u)
    end subroutine read_linmin_ef

  end subroutine atomopt_ef


  subroutine scf( etot, ierr_out )
    implicit none
    real(8),intent(out) :: etot
    integer,intent(out) :: ierr_out

    select case(SYStype)
    case default

       call calc_eion

       call construct_strfac
       call construct_ps_local
       call construct_ps_pcc
       call destruct_strfac

       select case(pselect)
       case(2)
          call prep_ps_nloc2
       case(3)
          call prep_ps_nloc3
       case(5)
          call prep_ps_nloc_mr
       end select

    case(1)

       call calc_eion

       call construct_ps_local_mol
       call construct_ps_pcc_mol
       call prep_ps_nloc2_mol

    end select

    call calc_scf( ierr_out, NiterSCF, Etot_out=etot )

  end subroutine scf


  subroutine write_trajectory( xyz, force, msg )
    implicit none
    real(8),intent(in) :: xyz(:,:), force(:,:)
    character(*),intent(in) :: msg
    integer,parameter :: u=297
    integer :: i, myrank
    include 'mpif.h'
    call MPI_Comm_rank(MPI_COMM_WORLD,myrank,i)
    if ( myrank == 0 ) then
      !
      open(u,file='trajectory_ef.xyz',position='append')
      call write_xyz_atom( u, xyz, "xyz_id_etot :"//msg )
      close(u)
      !
      open(u,file='trajectory_ef.force',position='append')
      write(u,*) size(force,2)
      write(u,*) "force_id_etot:",msg
      do i=1,size(force,2)
        write(u,*) force(:,i)
      end do
      close(u)
      !
    end if
  end subroutine write_trajectory


end module atomopt_ef_module
