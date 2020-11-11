module atomopt_bfgs_module

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
  use watch_module
  use atomopt_io_module, only: flag_continue_atomopt, read_bfgs_atomopt_io, write_bfgs_atomopt_io

  implicit none

  private
  public :: atomopt_bfgs

  integer :: SYStype
  logical :: disp_sw, disp_scf
  integer :: NiterSCF

contains


  subroutine atomopt_bfgs( SYStype_in, fmax_tol, ncycle, okstep_in, max_linmin_in, max_alpha_in, NiterSCF_in )

    implicit none
    integer,intent(in) :: SYStype_in
    real(8),intent(in) :: fmax_tol
    integer,intent(in) :: ncycle
    real(8),intent(in) :: okstep_in
    integer,intent(in) :: max_linmin_in
    real(8),intent(in) :: max_alpha_in
    integer,optional,intent(in) :: NiterSCF_in
    type(atom) :: ion
    type(lattice) :: aa
    integer,parameter :: max_loop_default=5
    integer :: np, max_loop,max_loop_linmin, linmin_start
    integer :: i1,i2,j1,j2,i,j,n,linmin,loop_start
    integer :: a, ierr, ip, jp, loop, icount,max_linmin
    integer :: istep_golden_search
    real(8) :: etot, fmax, tmp, okstep, max_alpha
    real(8) :: dxdg,dgHdg,dxg,dgHg,c,alpha,rdg
    real(8) :: aa_inv(3,3), da(3), da_tmp(3)
    real(8),allocatable :: g(:,:,:),x(:,:,:),Hg(:,:,:),Hdg(:,:)
    real(8),allocatable :: history(:,:),r(:,:)
    real(8),allocatable :: dx(:,:),dg(:,:),d(:,:)
    real(8),allocatable :: H(:,:), V(:,:), VT(:,:),W(:,:)
    real(8),allocatable :: HSR1(:,:)
    real(8),allocatable :: foraa(:,:)

    integer :: imax, myrank
    real(8) :: dd,dmax,dtmp,tau,factor
    real(8) :: tt_start(2),tt(2,0:10)
    logical :: flag_check(2)
    logical :: flag_3pts, flag_check_linmin_continue
    logical :: flag_exit, flag_sd
    real(8) :: Armijo, Wolfe2
    real(8) :: delta, sigma
    real(8),allocatable :: xyzf_linmin(:,:,:)
    real(8),allocatable :: info_linmin(:)
    real(8) :: alpha1,alpha2,alpha3,alpha4
    real(8) :: etot1,etot2,etot3,etot4
    real(8) :: etot_min, fmax_min
    include 'mpif.h'
    character(40) :: msg

    tt=0.0d0; call watchb( tt_start ); tt(:,0)=tt_start

    call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr )

    call write_border( 0, "atomopt_bfgs(start)" )
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

    max_linmin = max_linmin_in

    max_alpha = max_alpha_in

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

      !if ( disp_sw ) then
      !  write(*,*) "Initial configuration(in reduced coordinates)"
      !  do a=1,ion%natom
      !    write(*,'(1x,3f20.15)') ion%aaa(1:3,a)
      !  end do
      !end if

      call scf( etot, ierr ) ; if ( ierr == -1 ) goto 999

      call calc_force( ion%natom, ion%force, fmax )

      if ( fmax <= fmax_tol ) goto 900

      !allocate( foraa(3,ion%natom) ); foraa=0.0d0
      !foraa = matmul( transpose(aa_inv), ion%force )

    end if

! ---

    np = 1
    allocate( g(3,ion%natom,0:np)  ) ; g=0.0d0
    allocate( x(3,ion%natom,0:np)  ) ; x=0.0d0
    allocate( Hg(3,ion%natom,0:np) ) ; Hg=0.0d0
    allocate( Hdg(3,ion%natom) ) ; Hdg=0.0d0
    allocate( dx(3,ion%natom) ) ; dx=0.0d0
    allocate( dg(3,ion%natom) ) ; dg=0.0d0
    allocate( d(3,ion%natom) ) ; d=0.0d0
    allocate( r(3,ion%natom) ) ; r=0.0d0

    allocate( history(7,0:max_loop+1) ) ; history=0.0d0

    allocate( xyzf_linmin(3,ion%natom,2) ); xyzf_linmin=0.0d0
    allocate( info_linmin(size(history,1)) ); info_linmin=0.0d0

    if ( flag_continue_atomopt() ) then

      call read_bfgs_atomopt_io( &
           loop_start, &
           history, &
           x(:,:,1), x(:,:,0), g(:,:,1), g(:,:,0), &
           H )

      flag_check_linmin_continue = .true.

    else

      history(1,0) = etot
      history(2,0) = fmax
      history(3,0) = ierr ; if ( ierr == -2 ) history(3,0)=NiterSCF
      history(4,0) = sum(history(3,0:0))
      history(5:,0) = 0.0d0

    end if

    n = 3*ion%natom
    allocate( H(n,n)  ); H=0.0d0
    !allocate( V(n,n)  ); V=0.0d0
    !allocate( VT(n,n) ); VT=0.0d0
    !allocate( W(n,n)  ); W=0.0d0
    !allocate( HSR1(n,n) ); HSR1=0.0d0

    do i=1,size(H,1)
      H(i,i) = 1.0d0
    end do
    !do i=1,size(HSR1,1)
    !  HSR1(i,i) = 1.0d0
    !end do

! ---

    x(:,:,0) = ion%xyz(:,:)
    g(:,:,0) =-ion%force(:,:)

    do loop = loop_start, max_loop

      if ( disp_sw ) write(*,'(a50," loop",i2)') repeat("-",50),loop

      call watchb( tt(:,0) )

! ---

      flag_sd = .false.
      linmin_start = 1

! ---

      call write_coordinates_atom( 197, 3 )
      call write_bfgs_atomopt_io( loop, history(:,0:loop-1), &
           x(:,:,1), x(:,:,0), g(:,:,1), g(:,:,0), H )

      if ( loop == 1 ) then

        flag_sd = .true.

        x(:,:,1) = x(:,:,0)
        g(:,:,1) = g(:,:,0)
        d(:,:)   =-g(:,:,0)

      else

        flag_sd = .false.

        dx(:,:) = x(:,:,1) - x(:,:,0)
        dg(:,:) = g(:,:,1) - g(:,:,0)

        dxdg = sum( dx(:,:)*dg(:,:) )
        c = 1.0d0/dxdg

        n=3*ion%natom
        call DGEMV('N',n,n,1.0d0,H,n,dg(:,:),1,0.0d0,Hdg,1)
        dgHdg = sum(dg*Hdg)

        j=0
        do j2=1,ion%natom
        do j1=1,3
          j=j+1
          i=0
          do i2=1,ion%natom
          do i1=1,3
            i=i+1
            H(i,j) = H(i,j) - ( Hdg(i1,i2)*dx(j1,j2) + dx(i1,i2)*Hdg(j1,j2) )*c &
                 + ( 1.0d0 + dgHdg*c )*dx(i1,i2)*dx(j1,j2)*c
          end do
          end do
        end do
        end do

!        j=0
!        do j2=1,ion%natom
!        do j1=1,3
!          j=j+1
!          i=0
!          do i2=1,ion%natom
!          do i1=1,3
!            i=i+1
!            V(i,j) = -dg(i1,i2)*dx(j1,j2)*c
!          end do
!          end do
!        end do
!        end do
!        do i=1,size(V,1)
!          V(i,i) = V(i,i) + 1.0d0
!        end do
!        VT = transpose( V )
!        W = matmul( VT, H )
!        H = matmul( W , V )
!        j=0
!        do j2=1,ion%natom
!        do j1=1,3
!          j=j+1
!          i=0
!          do i2=1,ion%natom
!          do i1=1,3
!            i=i+1
!            H(i,j) = H(i,j) + dx(i1,i2)*dx(j1,j2)*c
!          end do
!          end do
!        end do
!        end do

! ---
!        n=3*ion%natom
!        call DGEMV('N',n,n,1.0d0,HSR1,n,dg(:,:),1,0.0d0,Hdg,1)
!        r(:,:) = dx(:,:) - Hdg(:,:)
!        rdg = sum(r*dg)
!        c = 1.0d0/rdg
!        j=0
!        do j2=1,ion%natom
!        do j1=1,3
!          j=j+1
!          i=0
!          do i2=1,ion%natom
!          do i1=1,3
!            i=i+1
!            HSR1(i,j) = HSR1(i,j) + r(i1,i2)*r(j1,j2)*c
!          end do
!          end do
!        end do
!        end do
! ---

        n=3*ion%natom
        call DGEMV('N',n,n,-1.0d0,H,n,g(:,:,1),1,0.0d0,d,1)
!        call DGEMV('N',n,n,-1.0d0,HSR1,n,g(:,:,1),1,0.0d0,d,1)

        !dmax=0.0d0
        !do i=1,ion%natom
        !  dd = sum( d(:,i)**2 )*alpha**2
        !  if ( dmax < dd ) then
        !    dmax = dd
        !    imax = i
        !  end if
        !end do

        c = sum(g(:,:,1)*d)
        if ( c > 0.0d0 ) then
          if ( disp_sw ) write(*,*) "------------------ Steepest Decent"
          flag_sd = .true.
          d(:,:)=-g(:,:,1)
          H=0.0d0
          do i=1,size(H,1)
            H(i,i) = 1.0d0
          end do
!          HSR1=0.0d0
!          do i=1,size(HSR1,1)
!            HSR1(i,i) = 1.0d0
!          end do
        end if

      end if !(loop/=1)

      do a=1,ion%natom
        da = matmul( aa_inv, d(:,a) )
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
        d(:,a) = matmul( aa%LatticeVector, da )
      end do

      call watchb( tt(:,0), tt(:,1) )

! ---

      tau = ( sqrt(5.0d0) - 1.0d0 )*0.5d0

      if ( flag_check_linmin_continue ) then
        call read_linmin_bfgs( linmin_start )
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
        call read_linmin_bfgs
      end if

      do linmin = linmin_start, max_linmin+1 !------------ Line-Mimization loop (start)

        if ( disp_sw ) write(*,'(a30," linmin",i2)') repeat("-",30),linmin

        call watchb( tt(:,0) )

        if ( myrank == 0 ) call write_linmin_bfgs( linmin )

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

        if ( disp_sw ) then
          write(*,'(a10,"linmin=",i4," (3pts:",l1,")")') repeat("-",10),linmin-1,flag_3pts
          write(*,*) alpha1,etot1
          write(*,*) alpha2,etot2
          write(*,*) alpha3,etot3
          if ( flag_3pts ) write(*,*) alpha4,etot4
        end if

        !if ( flag_3pts .and. istep_golden_search == 0 ) then
        !  alpha2=0.0d0; etot2=0.0d0
        !  alpha3=0.0d0; etot3=0.0d0
        !end if

        if ( flag_exit .or. linmin==max_linmin+1 ) exit

        dmax=0.0d0
        do a=1,ion%natom
          dtmp = sqrt(sum(d(:,a)**2))*alpha
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
      
        ion%xyz = x(:,:,1) + alpha*d(:,:)
        aa_atom(:,:) = matmul( aa_inv, ion%xyz )

        call shift_aa_coordinates_atom( aa_atom )
        !if ( disp_sw ) then
        !  write(*,*) "Next trial configuration"
        !  do a=1,ion%natom
        !    write(*,'(i2,2x,3f10.5,2x,3f10.5,2x,f10.5)') &
        !         a,aa_atom(:,a),ion%xyz(:,a),sqrt(sum(d(:,a)**2))*alpha
        !  end do
        !end if
        call write_coordinates_atom( 97, 3 )

        call watchb( tt(:,0), tt(:,2) )

        call scf( etot, ierr ) ; if ( ierr == -1 ) goto 999
        call calc_force( ion%natom, ion%force, fmax )
    
        write(msg,*) loop, etot
        call write_trajectory( ion%xyz, ion%force, msg )

        history(1,loop) = etot
        history(2,loop) = fmax
        icount = ierr
        if ( ierr == -1 ) then
          icount = 0
        else if ( ierr == -2 ) then
          icount=NiterSCF
        end if
        history(3,loop) = history(3,loop) + icount
        history(4,loop) = sum(history(3,0:loop))
        c=sum(d*g(:,:,0))
        history(5,loop) = c
        history(6,loop) = alpha
        history(7,loop) = linmin
        if ( disp_sw ) then
          do i=loop,loop
            write(*,'("linmin",i4,5x,es15.8,es14.5,i4,i6,2es14.5,1x,i6)') &
                 i, history(1:2,i), icount, nint(history(4,i)),history(5:6,i),nint(history(7,i))
          end do
          write(*,'("etime:",4f10.5)') tt(2,1),tt(2,2),tt(2,3),tt(2,3)
        end if

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
          info_linmin = history(:,loop)
          xyzf_linmin(:,:,1) = ion%xyz(:,:)
          xyzf_linmin(:,:,2) = ion%force(:,:)
        end if

        !delta = 1.0d-4
        !sigma = 0.99d0
        !c = sum(g(:,:,1)*d(:,:))
        !write(*,*) etot,fmax,alpha
        !write(*,*) "linmin,alpha=",loop_linmin,alpha
        !write(*,*) etot,history(icount,1)
        !write(*,*) fmax,history(icount,2)
        !write(*,*) "sum(g*d)",c
        !write(*,*) "Check Armijo condition",etot,"<=?",history(icount,1)+delta*alpha*c
        !write(*,*) "Check Wolfe2 condition",-sum(ion%force*d),">=?",sigma*c
        !write(*,*) "Check WolfeS condition",abs(sum(ion%force*d)),"<=?",sigma*abs(c)
        !write(*,*) "delta, sigma",delta,sigma
        !flag_check(:)=.false.
        !tmp = history(1,loop)+delta*alpha*c
        !if ( etot <= tmp ) flag_check(1)=.true.
        !tmp = abs(sum(ion%force*d))
        !if ( tmp <= sigma*abs(c) ) flag_check(2)=.true.

        if ( fmax <= fmax_tol ) flag_exit=.true.

        call watchb( tt(:,0), tt(:,3) )

      end do !linmin

      history(:,loop) = info_linmin(:)
      fmax = info_linmin(2)
      ion%xyz = xyzf_linmin(:,:,1)
      ion%force = xyzf_linmin(:,:,2)

      call watchb( tt(:,0) )

      x(:,:,0) = x(:,:,1)
      g(:,:,0) = g(:,:,1)
      x(:,:,1) = ion%xyz(:,:)
      g(:,:,1) =-ion%force(:,:)

      call watchb( tt(:,0), tt(:,4) )

      if ( disp_sw ) then
        write(*,'(a7,2x,a15,a14,2a6,2a14,1x,a6)') "History","Etot","Fmax","iter","iter" &
                                          ,"(d,g)","alpha","linmin"
        do i=0,loop
          write(*,'(i4,5x,es15.8,es14.5,2i6,2es14.5,1x,i6)') &
               i, history(1:2,i), nint(history(3:4,i)),history(5:6,i),nint(history(7,i))
        end do
        write(*,'("etime:",4f10.5)') tt(2,1),tt(2,2),tt(2,3),tt(2,3)
      end if

      linmin_start = 1

      if ( fmax <= fmax_tol ) exit

    end do ! loop    

! ---

900 if ( disp_sw ) then
       write(*,'(1x,"etot    :",f25.15)') etot
       write(*,'(1x,"fmax/tol:",es12.5," /",es12.5)') fmax,fmax_tol
    end if

999 call check_disp_switch( disp_sw, 1 )
    call write_border( 0, "atomopt_bfgs(end)" )

    if ( allocated(history) ) then
       deallocate( history )
       deallocate( dg )
       deallocate( dx )
       deallocate( Hdg )
       deallocate( Hg )
       deallocate( x )
       deallocate( g )
    end if

  contains

    subroutine write_linmin_bfgs( linmin )
      implicit none
      integer,intent(in) :: linmin
    end subroutine write_linmin_bfgs

    subroutine read_linmin_bfgs( linmin )
      implicit none
      integer,optional,intent(out) :: linmin
      if ( present(linmin) ) linmin=0
    end subroutine read_linmin_bfgs

  end subroutine atomopt_bfgs


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
      open(u,file='trajectory_bfgs.xyz',position='append')
      call write_xyz_atom( u, xyz, "xyz_id_etot :"//msg )
      close(u)
      !
      open(u,file='trajectory_bfgs.force',position='append')
      write(u,*) size(force,2)
      write(u,*) "force_id_etot:",msg
      do i=1,size(force,2)
        write(u,*) force(:,i)
      end do
      close(u)
      !
    end if
  end subroutine write_trajectory


end module atomopt_bfgs_module
