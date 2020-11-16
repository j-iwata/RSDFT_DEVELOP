module atomopt_rf_module

  use lattice_module
  use atom_module, only: atom, construct_atom, aa_atom &
                        ,write_coordinates_atom, shift_aa_coordinates_atom
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
  use atomopt_io_module, only: flag_continue_atomopt, write_atomopt_io, read_atomopt_io
  use efield_module, only: sawtooth_efield

  implicit none

  private
  public :: atomopt_rf

  logical :: disp_sw, disp_scf
  integer :: SYStype
  integer :: NiterSCF

contains


  subroutine atomopt_rf( SYStype_in, fmax_tol, ncycle, okstep_in, max_linmin_in, max_alpha_in, NiterSCF_in )

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
    integer,parameter :: np=20
    integer :: ishape(1), ishape2(2), i1,i2, LWORK=0
    integer :: n,il,iu,m,a,ierr,loop,i,j,icount,ip
    integer :: max_loop,linmin,max_linmin
    integer,allocatable :: iwork(:),ifail(:)
    real(8) :: etot0, etot, fmax, okstep,max_alpha
    real(8) :: dxdg,dxHdx,c1,c2,c
    real(8) :: vl,vu,tol,alpha,dmax,dtmp
    real(8) :: aa_inv(3,3),da(3),da_tmp(3)
    real(8),parameter :: one=1.0d0, zero=0.0d0
    real(8),allocatable :: g(:),x(:),Hessian(:,:),Htmp(:,:)
    real(8),allocatable :: history(:,:),g0(:),x0(:),d(:)
    real(8),allocatable :: dx(:),dg(:),Hdx(:)
    real(8),allocatable :: w(:),z(:,:),work(:)
    integer :: loop_start
    real(8),allocatable :: info_linmin(:), xyzf_linmin(:,:,:)
    real(8) :: alpha1,alpha2,alpha3,alpha4
    real(8) :: etot1,etot2,etot3,etot4
    logical :: flag_sd

    call write_border( 0, "atomopt_rf(start)" )
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

      call scf( etot, ierr ) ; if ( ierr == -1 ) goto 999
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
    allocate( Htmp(n+1,n+1) ) ; Htmp=0.0d0
    allocate( Hdx(n)    ) ; Hdx=0.0d0
    allocate( dx(n)     ) ; dx=0.0d0
    allocate( dg(n)     ) ; dg=0.0d0

    ishape(1) = n
    ishape2(1:2) = (/ 3, ion%natom /)

    allocate( history(6,0:max_loop*np) ) ; history=0.0d0

    allocate( info_linmin(size(history,1)) ); info_linmin=0.0d0
    allocate( xyzf_linmin(3,ion%natom,2) ); xyzf_linmin=0.0d0

    if ( flag_continue_atomopt() ) then

      call read_atomopt_io( &
           loop_start, &
           history, &
           x, x0, g, g0, &
           Hessian )

    else

      history(1,0) = etot
      history(2,0) = fmax
      history(3,0) = ierr ; if ( ierr == -2 ) history(3,0)=NiterSCF
      history(4,0) = sum( history(3,0:0) )
      history(5,0) = 0.0d0

    end if

! --- for LAPACK

    tol=1.d-10
    vl=0.0d0
    vu=0.0d0
    il=1
    iu=1
    m=iu-il+1
    allocate( w(n+1)         ) ; w=0.0d0
    allocate( z(n+1,m)       ) ; z=0.0d0
    allocate( work(8*(n+1))  ) ; work=0.0d0
    allocate( iwork(5*(n+1)) ) ; iwork=0
    allocate( ifail(n+1)     ) ; ifail=0

! ---

    do loop = loop_start, max_loop

      if ( disp_sw ) write(*,'(a60," loop=",i4)') repeat("-",60),loop

      alpha = 1.0d0

      call write_coordinates_atom( 197, 3 )
      call write_atomopt_io( loop,history(:,0:loop-1),x,x0,g,g0,Hessian )

      if ( loop == 1 ) then

        flag_sd = .true.

        Hessian=0.0d0
        do i=1,n
          Hessian(i,i) = one
        end do

        x(:) = reshape( ion%xyz(:,:)  , ishape )
        g(:) = reshape(-ion%force(:,:), ishape )
        d(:) = -g(:)

        x0(:) = x(:)
        g0(:) = g(:)

      else

        dx(:) = x(:) - x0(:)
        dg(:) = g(:) - g0(:)

        dxdg = sum( dx*dg )
        if ( disp_sw ) write(*,*) "dxdg=",dxdg

        if ( dxdg < 0.0d0 ) then

          if ( disp_sw ) write(*,*) "Steepest Descent"

          flag_sd = .true.

          d(:) = -g(:)

          Hessian = 0.0d0
          do i = 1, size(Hessian,1)
            Hessian(i,i) = 1.0d0
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

          Htmp=0.0d0
          Htmp(1:n,1:n) = Hessian(:,:)
          Htmp(1:n,n+1) = g(:)
          Htmp(n+1,1:n) = g(:)

!---
          if ( LWORK == 0 ) then
            call dsyevx('V','I','U',n+1,Htmp,n+1,vl,vu,il,iu,tol,m,w,z,n+1 &
                 ,work,-1,iwork,ifail,ierr)
            LWORK=nint(work(1))
            deallocate( work )
            allocate( work(LWORK) ); work=0.0d0
          end if
          call dsyevx('V','I','U',n+1,Htmp,n+1,vl,vu,il,iu,tol,m,w,z,n+1 &
               ,work,size(work),iwork,ifail,ierr)
          d(:) = z(1:n,1)/z(n+1,1)
!---
!          if ( LWORK == 0 ) then
!            call dsyev('V','U',n+1,Htmp,n+1,w,work,-1,ierr)
!            LWORK=nint( work(1) )
!            deallocate( work )
!            allocate( work(LWORK) ); work=0.0d0
!          end if
!          call dsyev('V','U',n+1,Htmp,n+1,w,work,size(work),ierr)
!          d(:) = Htmp(1:n,1)/Htmp(n+1,1)
!---

        end if

        x0(:) = x(:)
        g0(:) = g(:)

      end if

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

      info_linmin=0.0d0
      xyzf_linmin=0.0d0
      alpha1=0.0d0; alpha2=0.0d0; alpha3=0.0d0; alpha4=0.0d0
      etot1 =0.0d0; etot2 =0.0d0; etot3 =0.0d0; etot4 =0.0d0

      do linmin = 1, max_linmin !---------- line minimization

        select case(linmin)
        case( 1 )
          alpha  = 1.0d0
          alpha1 = alpha
        case(2)
          alpha  = (max_alpha-1.0d0)*0.3d0
          alpha2 = alpha
        case(3)
          alpha  = max_alpha
          alpha3 = alpha
        case(4)
          if ( etot1 >= etot2 .and. etot2 <= etot3 ) then
            alpha  = alpha1 + (alpha3-alpha2)
            alpha4 = alpha
         else if ( etot1 <= etot2 .and. etot1 <= etot3 ) then
            alpha = alpha1*0.3d0
            alpha4=alpha1 ; etot4=etot1
            alpha1=alpha  ; etot1=0.0d0 
          else if ( etot3 <= etot2 .and. etot3 <= etot1 ) then
            alpha = alpha3*1.7d0
            alpha4=alpha3; etot4=etot3
            alpha3=alpha ; etot3=0.0d0
          else
            write(*,'("alpha1,apha2,alpha3,alpha4",4f10.5)') alpha1,alpha2,alpha3,alpha4
            write(*,'("etot1,apha2,etot3,etot4",4f10.5)') etot1,etot2,etot3,etot4
            call stop_program('zzz')
          end if
        case(5:)
          
        end select

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
          if ( disp_sw ) then
            write(*,*) "Maxmimum displacement is limited to",okstep
            write(*,*) "alpha is changed: alpha=",alpha
          end if
        end if

        x(:) = x0(:) + alpha*d(:)
        ion%xyz(:,:) = reshape( x, ishape2 )

! ---
        aa_atom(:,:) = matmul( aa_inv, ion%xyz )
        call shift_aa_coordinates_atom( aa_atom )
        if ( .false. ) then
        !if ( disp_sw ) then
          write(*,*) "Next trial configuration"
          do a=1,ion%natom
            i2 = 3*a
            i1 = i2 - 3 + 1
            write(*,'(i2,2x,3f10.5,2x,3es14.5,2x,es14.5)') &
                 a,aa_atom(:,a),ion%xyz(:,a),sqrt(sum(d(i1:i2)**2))*alpha
          end do
        end if
          
        call write_coordinates_atom( 97, 3 )

        call scf( etot, ierr ) ; if ( ierr == -1 ) goto 999
        call calc_force( ion%natom, ion%force, fmax )

        history(1,loop) = etot
        history(2,loop) = fmax
        history(3,loop) = ierr ; if ( ierr == -2 ) ierr=NiterSCF
        history(4,loop) = sum( history(3,0:loop) )
        history(5,loop) = alpha
        c=0.0d0
        do a = 1, ion%natom
          i2 = 3*a
          i1 = i2 - 3 + 1
          c = c - sum( d(i1:i2)*ion%force(:,a) )
        end do
        history(6,loop) = c

        if ( disp_sw ) then
          do i=loop,loop
            write(*,'("linmin ",i4,f20.10,es14.5,i4,i6,f10.5,es14.5)') &
                 i, history(1:2,i), nint(history(3:4,i)),history(5:6,i)
          end do
        end if

        if ( etot < info_linmin(1) ) then
          info_linmin(1:6) = history(1:6,loop)
          xyzf_linmin(:,:,1) = ion%xyz
          xyzf_linmin(:,:,2) = ion%force
        end if

        select case(linmin)
        case(1); etot1=etot
        case(2); etot2=etot
        case(3); etot3=etot
        end select

        if ( linmin >= 4 ) then
          if ( etot1 == 0.0d0 ) then
            etot1=etot
          else if ( etot3 == 0.0d0 ) then
            etot3=etot
          else if ( etot4 == 0.0d0 ) then
            etot4=etot
          end if
        end if

        !write(*,'("alpha1,alpha2,alpha3,alpha4",4f12.7)') alpha1,alpha2,alpha3,alpha4
        !write(*,'("etot1 ,etot2 ,etot3 ,etot4 ",4f12.7)') etot1,etot2,etot3,etot4

        if ( .not.flag_sd .or. fmax <= fmax_tol ) exit

      end do ! linmin

      history(:,loop) = info_linmin(:)
      fmax = info_linmin(2)
      ion%xyz   = xyzf_linmin(:,:,1)
      ion%force = xyzf_linmin(:,:,2)

      if ( disp_sw ) then
        do i=0,loop
          write(*,'("History",i4,f20.10,es14.5,i4,i6,f10.5,es14.5)') &
               i, history(1:2,i), nint(history(3:4,i)),history(5:6,i)
        end do
      end if

      if ( fmax <= fmax_tol ) exit

      x(:) =  reshape( ion%xyz(:,:), ishape )
      g(:) = -reshape( ion%force(:,:), ishape )

    end do ! loop

! ---

900 if ( disp_sw ) then
       write(*,'(1x,"etot    :",f25.15)') etot
       write(*,'(1x,"fmax/tol:",es12.5," /",es12.5)') fmax,fmax_tol
    end if

999 call check_disp_switch( disp_sw, 1 )
    call write_border( 0, "atomopt_rf(end)" )

    deallocate( ifail )
    deallocate( iwork )
    deallocate( work )
    deallocate( z )
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

  end subroutine atomopt_rf


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

    call sawtooth_efield( Vion )

    call calc_scf( ierr_out, NiterSCF, Etot_out=etot )

  end subroutine scf


end module atomopt_rf_module
