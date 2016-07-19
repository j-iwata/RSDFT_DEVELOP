MODULE atomopt_diis_module

  use lattice_module
  use atom_module, only: atom, construct_atom, aa_atom, write_coordinates_atom
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

  implicit none

  PRIVATE
  PUBLIC :: atomopt_diis

  logical :: disp_sw
  integer :: SYStype
  real(8) :: beta=0.5d0

CONTAINS


  SUBROUTINE atomopt_diis( SYStype_in, fmax_tol )

    implicit none
    integer,intent(IN) :: SYStype_in
    real(8),intent(IN) :: fmax_tol
    type(atom) :: ion
    type(lattice) :: aa
    integer,parameter :: max_loop=50
    integer :: np
    integer :: a, ierr, ip, jp, loop,i,j
    real(8) :: etot0, etot, fmax, tmp
    real(8) :: aa_inv(3,3)
    real(8),allocatable :: g(:,:,:),x(:,:,:)
    real(8),allocatable :: history(:,:)

    call write_border( 0, "atomopt_diis(start)" )
    call check_disp_switch( disp_sw, 0 )
    call check_disp_switch( .false., 1 )

    SYStype = SYStype_in

    call get_aa_lattice( aa )
    call get_inverse_lattice( aa%LatticeVector, aa_inv )

    call construct_atom( ion )

    do a=1,ion%natom
       ion%xyz(1:3,a) = matmul( aa%LatticeVector(1:3,1:3), ion%aaa(1:3,a) )
    end do

    if ( disp_sw ) then
       write(*,*) "Initial configuration"
       do a=1,ion%natom
          write(*,'(1x,3f20.15)') ion%aaa(1:3,a)
       end do
    end if

    call scf( etot )
    call calc_force( ion%natom, ion%force, fmax )

    if ( fmax <= fmax_tol ) goto 900

! ---

    np = 3
    allocate( g(3,ion%natom,0:np) ) ; g=0.0d0
    allocate( x(3,ion%natom,0:np) ) ; x=0.0d0

    allocate( history(0:max_loop,2) ) ; history=0.0d0

    history(0,1) = etot
    history(0,2) = fmax

! ---

    ip = 0

    x(:,:,ip) = ion%xyz(:,:)
    g(:,:,ip) = ion%force(:,:)

    do loop=1,max_loop

       if ( disp_sw ) write(*,'(a60," ICY=",i4)') repeat("-",60),loop

! ---

       beta = min( 1.0d0, 1.d-2/fmax )
       !beta = 1.d-2/fmax

       if ( ip > 0 ) then
          call diis( size(g(:,:,0)), ip+1, g, x, ion%xyz )
       else
          ion%xyz(:,:) = x(:,:,ip) + beta*g(:,:,ip)
       end if

! ---

!       do a=1,ion%natom
!          ion%xyz(1:3,a) = matmul( aa%LatticeVector(1:3,1:3), ion%aaa(1:3,a) )
!       end do

       aa_atom(:,:) = matmul( aa_inv, ion%xyz )

! ---

       call write_coordinates_atom( 97, 3 )

       call scf( etot )
       call calc_force( ion%natom, ion%force, fmax )

       if ( fmax <= fmax_tol ) goto 900

       history(loop,1) = etot
       history(loop,2) = fmax

       if ( disp_sw ) then
          do i=0,loop
             write(*,*) i, (history(i,j),j=1,2)
          end do
       end if

       ip = ip + 1
       if ( ip > np ) then
          do jp=1,np-1
             x(:,:,jp) = x(:,:,jp+1)
             g(:,:,jp) = g(:,:,jp+1)
          end do
          ip = np
       end if

       x(:,:,ip) = ion%xyz(:,:)
       g(:,:,ip) = ion%force(:,:)

    end do ! loop    

! ---

900 if ( disp_sw ) then
       write(*,'(1x,"etot    :",f25.15)') etot
       write(*,'(1x,"fmax/tol:",es12.5," /",es12.5)') fmax,fmax_tol
    end if

    call check_disp_switch( disp_sw, 1 )
    call write_border( 0, "atomopt_diis(end)" )

    if ( allocated(x) ) deallocate( x )
    if ( allocated(g) ) deallocate( g )

  END SUBROUTINE atomopt_diis


  SUBROUTINE diis( m,n,g,x,xout )
    implicit none
    integer,intent(IN)  :: m,n
    real(8),intent(IN)  :: g(m,n),x(m,n)
    real(8),intent(OUT) :: xout(m)
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

    xout(:) = matmul( x(:,1:n),xx(1:n) ) + beta*matmul( g(:,1:n),xx(1:n) )

    deallocate( ipiv )
    deallocate( xx   )
    deallocate( gg   )

  END SUBROUTINE diis


  SUBROUTINE scf( etot )
    implicit none
    real(8),intent(OUT) :: etot
    integer :: ierr

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

    call calc_scf( .false., ierr, 50, Etot_out=etot )

  END SUBROUTINE scf


END MODULE atomopt_diis_module
