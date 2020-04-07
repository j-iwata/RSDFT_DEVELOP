MODULE atomopt_rf_module

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
  use vector_tools_module, only: vinfo

  implicit none

  PRIVATE
  PUBLIC :: atomopt_rf

  logical :: disp_sw, disp_scf
  integer :: SYStype
  integer :: NiterSCF

CONTAINS


  SUBROUTINE atomopt_rf( v, SYStype_in, fmax_tol, NiterSCF_in )

    implicit none
    type(vinfo),intent(IN) :: v(2)
    integer,intent(IN) :: SYStype_in
    real(8),intent(IN) :: fmax_tol
    integer,optional,intent(IN) :: NiterSCF_in

    type(atom) :: ion
    type(lattice) :: aa, bb

    integer,parameter :: max_loop=50, np=20
    integer :: ishape(1), ishape2(2)
    integer :: n,il,iu,m,a,ierr,loop,i,j,icount,ip
    integer,allocatable :: iwork(:),ifail(:)
    real(8) :: etot0, etot, fmax
    real(8) :: dxdg,dxHdx
    real(8) :: vl,vu,tol
    real(8) :: aa_inv(3,3)
    real(8),parameter :: one=1.0d0, zero=0.0d0
    real(8),allocatable :: g(:),x(:),Hessian(:,:),Htmp(:,:)
    real(8),allocatable :: history(:,:),g0(:),x0(:)
    real(8),allocatable :: dx(:),dg(:),Hdx(:)
    real(8),allocatable :: w(:),z(:,:),work(:)

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

! ---

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

    call scf( v, etot, ierr ) ; if ( ierr == -1 ) goto 999
    call calc_force( ion%natom, ion%force, fmax )

    if ( fmax <= fmax_tol ) goto 900

! ---

    n = 3*ion%natom
    allocate( g(n)  ) ; g=0.0d0
    allocate( x(n)  ) ; x=0.0d0
    allocate( g0(n) ) ; g0=0.0d0
    allocate( x0(n) ) ; x0=0.0d0
    allocate( Hessian(n+1,n+1) ) ; Hessian=0.0d0
    allocate( Htmp(n+1,n+1) ) ; Htmp=0.0d0
    allocate( Hdx(n)    ) ; Hdx=0.0d0
    allocate( dx(n)     ) ; dx=0.0d0
    allocate( dg(n)     ) ; dg=0.0d0

    ishape(1) = n
    ishape2(1:2) = (/ 3, ion%natom /)

    allocate( history(0:max_loop*np,3) ) ; history=0.0d0

    history(0,1) = etot
    history(0,2) = fmax
    history(0,3) = ierr

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

    icount = 0

    do loop=1,max_loop

       do i=1,n
          Hessian(i,i) = one
       end do

       x(:) = reshape( ion%xyz(:,:)  , ishape )
       g(:) = reshape(-ion%force(:,:), ishape )
       ion%xyz(:,:) = ion%xyz(:,:) + ion%force(:,:)

       do ip=1,np

          if ( disp_sw ) write(*,'(a60," ICY=",2i4)') repeat("-",60),loop,ip

          aa_atom(:,:) = matmul( aa_inv, ion%xyz )
          call shift_aa_coordinates_atom( aa_atom )

          if ( disp_sw ) then
             write(*,*) "Next trial configuration"
             do a=1,ion%natom
!                write(*,'(1x,3f20.15)') aa_atom(1:3,a)
                write(*,*) aa_atom(1:3,a)
             end do
          end if

          call write_coordinates_atom( 97, 3 )

          call scf( v, etot, ierr ) ; if ( ierr == -1 ) goto 999
          call calc_force( ion%natom, ion%force, fmax )

          if ( fmax <= fmax_tol ) goto 900

          icount = icount + 1
          history(icount,1) = etot
          history(icount,2) = fmax
          history(icount,3) = ierr

          if ( disp_sw ) then
             do i=0,icount
                write(*,'(1x,i4,f20.10,es14.5,i4)') &
                     i, (history(i,j),j=1,2), nint(history(i,3))
             end do
          end if

          x(:) = reshape( ion%xyz(:,:)  , ishape )
          g(:) = reshape(-ion%force(:,:), ishape )

          dx(:) = x(:) - x0(:)
          dg(:) = g(:) - g0(:)

          dxdg = sum( dx*dg )
          if ( disp_sw ) write(*,*) "dxdg=",dxdg
          !if ( dxdg < 0.0d0 ) then
          !   ion%xyz(:,:)   = reshape( x0, ishape2 )
          !   ion%force(:,:) = reshape( -g0, ishape2 )
          !   exit
          !end if

          x0(:) = x(:)
          g0(:) = g(:)

          call dgemv( 'N', n, n, one, Hessian, n+1, dx, 1, zero, Hdx, 1 )

          dxHdx = sum( dx*Hdx )

          do j=1,n
          do i=1,n
             Hessian(i,j) = Hessian(i,j) &
                  + dg(i)*dg(j)/dxdg - Hdx(i)*Hdx(j)/dxHdx
          end do
          end do

          Hessian(1:n,n+1) = g(:)
          Hessian(n+1,1:n) = g(:)

          Htmp(:,:) = Hessian(:,:)

          call dsyevx('V','I','U',n+1,Htmp,n+1,vl,vu,il,iu,tol,m,w,z,n+1 &
               ,work,size(work),iwork,ifail,ierr)

          dx(:) = z(1:n,1)/z(n+1,1)

          i=0
          do a=1,ion%natom
             do j=1,3
                i=i+1
                ion%xyz(j,a) = x(i) + dx(i)
             end do
             if ( disp_sw ) write(*,*) a, sqrt(sum(dx(i-2:i)**2))
          end do

       end do ! ip

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

  END SUBROUTINE atomopt_rf


  SUBROUTINE scf( v, etot, ierr_out )
    implicit none
    type(vinfo),intent(IN) :: v(2)
    real(8),intent(OUT) :: etot
    integer,intent(OUT) :: ierr_out

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

    call calc_scf( v, disp_scf, ierr_out, NiterSCF, Etot_out=etot )

  END SUBROUTINE scf


END MODULE atomopt_rf_module
