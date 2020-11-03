module atomopt_diis_module

  use lattice_module
  use atom_module, only: atom, construct_atom, aa_atom, write_coordinates_atom
  use force_module
  use atomopt_io_module, only: flag_continue_atomopt, read_diis_atomopt_io, write_diis_atomopt_io
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
  use density_module, only: rho, normalize_density
  use hartree_module, only: calc_hartree
  use xc_module, only: calc_xc

  implicit none

  private
  public :: atomopt_diis
  public :: calc_coef_diis

  logical :: disp_sw
  integer :: SYStype
  real(8) :: alpha=1.0d0
  integer :: NiterSCF

contains


  subroutine atomopt_diis( SYStype_in, fmax_tol, ncycle, okstep_in, NiterSCF_in )

    implicit none
    integer,intent(in) :: SYStype_in
    real(8),intent(in) :: fmax_tol
    integer,intent(in) :: ncycle
    real(8),intent(in) :: okstep_in
    integer,optional,intent(in) :: NiterSCF_in
    type(atom) :: ion
    type(lattice) :: aa
    integer :: max_loop
    integer :: np,loop_start
    integer :: a, ierr, ip, jp, loop,i,j
    real(8) :: etot0, etot, fmax, tmp, okstep
    real(8) :: aa_inv(3,3)
    real(8),allocatable :: g(:,:,:),x(:,:,:),d(:,:)
    real(8),allocatable :: rho_diis(:,:,:), coef_diis(:)
    real(8),allocatable :: history(:,:)
    logical :: do_rho_diis

    call write_border( 0, "atomopt_diis(start)" )
    call check_disp_switch( disp_sw, 0 )
    call check_disp_length( i, 0 )
    if ( i < 2 ) call check_disp_switch( .false., 1 )

! ---

    SYStype = SYStype_in

    NiterSCF = 50 ; if ( present(NiterSCF_in) ) NiterSCF=NiterSCF_in

    max_loop = ncycle

    loop_start = 1

    okstep = okstep_in

    do_rho_diis = .false.

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

    np = 20
    allocate( g(3,ion%natom,0:np) ); g=0.0d0
    allocate( x(3,ion%natom,0:np) ); x=0.0d0
    allocate( d(3,ion%natom)      ); d=0.0d0
    allocate( coef_diis(0:np)     ); coef_diis=0.0d0
    if ( do_rho_diis ) then
      allocate( rho_diis(size(rho,1),size(rho,2),0:np) ); rho_diis=0.0d0
    end if
    allocate( history(5,0:max_loop) ) ; history=0.0d0

    if ( flag_continue_atomopt() ) then

      if ( do_rho_diis ) then
        call read_diis_atomopt_io( &
             loop_start, &
             history(:,0:loop_start-1), &
             ip, x, g, rho_diis )
      else
        call read_diis_atomopt_io( &
             loop_start, &
             history(:,0:loop_start-1), &
             ip, x, g )
      end if

      fmax = history(2,ip)

    else

      ip = 0

      history(1,0) = etot
      history(2,0) = fmax
      history(3,0) = ierr ; if ( ierr == -2 ) history(3,0)=NiterSCF 
      history(4,0) = sum( history(3,0:ip) )
      history(5,0) = 0.0d0

      x(:,:,ip) = ion%xyz(:,:)
      g(:,:,ip) = ion%force(:,:)

      if ( do_rho_diis ) rho_diis(:,:,ip) = rho(:,:)

    end if

! ---

    do loop = loop_start, max_loop

      if ( disp_sw ) write(*,'(a60," loop=",i4)') repeat("-",60),loop

! ---

      call write_coordinates_atom( 197, 3 )
      if ( do_rho_diis ) then
        call write_diis_atomopt_io( &
             loop, &
             history(:,0:loop-1), &
             ip, x(:,:,0:ip), g(:,:,0:ip), rho_diis(:,:,0:ip) )
      else
        call write_diis_atomopt_io( &
             loop, &
             history(:,0:loop-1), &
             ip, x(:,:,0:ip), g(:,:,0:ip) )
      end if

! ---

      alpha = min( 1.0d0, 1.d-2/fmax )
      !alpha = 1.d-2/fmax
      !alpha = 1.0d0

      if ( ip > 0 ) then

        call diis_private( size(g(:,:,0)), ip+1, g, x, coef_diis )

        d(:,:)=0.0d0
        do i=0,ip
          d(:,:) = d(:,:) + g(:,:,i)*coef_diis(i)
        end do
        call for_safety_move( d, alpha, okstep, aa%LatticeVector, aa_inv ) 
        ion%xyz(:,:) = alpha*d(:,:)
        do i=0,ip
          ion%xyz(:,:) = ion%xyz(:,:) + x(:,:,i)*coef_diis(i)
        end do

        if ( do_rho_diis ) then
          rho=0.0d0
          do i=0,ip
            rho(:,:) = rho(:,:) + coef_diis(i)*rho_diis(:,:,i)
          end do
          call normalize_density( rho )
          call calc_hartree(lbound(rho,1),ubound(rho,1),size(rho,2),rho)
          call calc_xc
        end if

      else

        call for_safety_move( g(:,:,ip), alpha, okstep, aa%LatticeVector, aa_inv ) 
        ion%xyz(:,:) = x(:,:,ip) + alpha*g(:,:,ip)

      end if

! ---

       aa_atom(:,:) = matmul( aa_inv, ion%xyz )

! ---

       call write_coordinates_atom( 97, 3 )

       call scf( etot, ierr ) ; if ( ierr == -1 ) goto 999
       call calc_force( ion%natom, ion%force, fmax )

       history(1,loop) = etot
       history(2,loop) = fmax
       history(3,loop) = ierr ; if ( ierr == -2 ) history(3,loop)=NiterSCF
       history(4,loop) = sum( history(3,0:loop) )
       history(5,loop) = alpha

       if ( disp_sw ) then
          do i=0,loop
             write(*,'(1x,i4,f20.10,es14.5,i4,i6,f10.5)') &
             i,history(1:2,i),nint(history(3:4,i)),history(5,i)
          end do
       end if

       if ( fmax <= fmax_tol ) goto 900

       ip = ip + 1
       if ( ip > np ) then
          do jp=1,np-1
             x(:,:,jp) = x(:,:,jp+1)
             g(:,:,jp) = g(:,:,jp+1)
             if ( do_rho_diis ) rho_diis(:,:,jp) = rho_diis(:,:,jp+1)
          end do
          ip = np
       end if

       x(:,:,ip) = ion%xyz(:,:)
       g(:,:,ip) = ion%force(:,:)

       if ( do_rho_diis ) rho_diis(:,:,ip) = rho(:,:)

    end do ! loop    

! ---

900 if ( disp_sw ) then
       write(*,'(1x,"etot    :",f25.15)') etot
       write(*,'(1x,"fmax/tol:",es12.5," /",es12.5)') fmax,fmax_tol
    end if

999 call check_disp_switch( disp_sw, 1 )
    call write_border( 0, "atomopt_diis(end)" )

    if ( allocated(history) ) deallocate( history )
    if ( allocated(coef_diis) ) deallocate( coef_diis )
    if ( allocated(rho_diis) ) deallocate( rho_diis )
    if ( allocated(d) ) deallocate( d )
    if ( allocated(x) ) deallocate( x )
    if ( allocated(g) ) deallocate( g )

  end subroutine atomopt_diis


  subroutine diis_private( m,n,g,x,coef_out )
    implicit none
    integer,intent(in)  :: m,n
    real(8),intent(in)  :: g(m,n),x(m,n)
    real(8),intent(out) :: coef_out(n)
    real(8),allocatable :: gg(:,:),xx(:)
    integer,allocatable :: ipiv(:)
    integer :: info,i,j

    allocate( gg(n+1,n+1) ) ; gg=0.0d0
    allocate( xx(n+1)     ) ; xx=0.0d0
    allocate( ipiv(n+1)   ) ; ipiv=0

    do j=1,n
       do i=1,n
          gg(i,j) = sum( g(:,i)*g(:,j) )
       end do
    end do

    gg(1:n,n+1) = 1.0d0
    gg(n+1,1:n) = 1.0d0
    xx(n+1) = 1.0d0

    call DGESV( n+1, 1, gg, n+1, ipiv, xx, n+1, info )

    coef_out(:) = xx(1:n)

    deallocate( ipiv )
    deallocate( xx   )
    deallocate( gg   )

  end subroutine diis_private


  subroutine calc_coef_diis( coef, f )
    implicit none
    real(8),intent(out) :: coef(:)
    real(8),intent(in)  :: f(:,:)
    real(8),allocatable :: A(:,:),x(:)
    integer,allocatable :: ipiv(:)
    integer :: info,i,j,n

    n = size( f, 2 )
    allocate( A(n+1,n+1) ); A=0.0d0
    allocate( x(n+1)     ); x=0.0d0
    allocate( ipiv(n+1)  ); ipiv=0

    do j=1,n
       do i=1,n
          A(i,j) = sum( f(:,i)*f(:,j) )
       end do
    end do

    A(1:n,n+1) = 1.0d0
    A(n+1,1:n) = 1.0d0
    x(n+1) = 1.0d0

    call DGESV( n+1, 1, A, n+1, ipiv, x, n+1, info )

    coef(:) = x(1:n)

    deallocate( ipiv )
    deallocate( x )
    deallocate( A )

  end subroutine calc_coef_diis


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


  subroutine for_safety_move( d, alpha, okstep, aa, aa_inv )
    implicit none
    real(8),intent(inout) :: d(:,:)
    real(8),intent(inout) :: alpha
    real(8),intent(in) :: okstep
    real(8),intent(in) :: aa(:,:), aa_inv(:,:)
    integer :: a,i,j
    real(8) :: da(3), da_tmp(3), dmax, dtmp
    dmax=0.0d0
    do a = 1, size(d,2)
      da = matmul( aa_inv, d(:,a) )
      do i = 1, 3
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
      d(:,a) = matmul( aa, da(:) )
      dtmp = sqrt(sum(da**2))*alpha
      dmax = max(dtmp,dmax)
    end do
    if ( dmax > okstep ) then
      alpha=alpha*okstep/dmax
      if ( disp_sw ) then
        write(*,*) "Maxmimum displacement is limited to",okstep
        write(*,*) "alpha is changed: alpha=",alpha
      end if
    end if
  end subroutine for_safety_move


end module atomopt_diis_module
