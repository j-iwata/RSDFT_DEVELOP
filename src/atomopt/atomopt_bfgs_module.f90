MODULE atomopt_bfgs_module

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

  implicit none

  PRIVATE
  PUBLIC :: atomopt_bfgs

  integer :: SYStype
  logical :: disp_sw, disp_scf
  integer :: NiterSCF

CONTAINS


  SUBROUTINE atomopt_bfgs( SYStype_in, fmax_tol, NiterSCF_in )

    implicit none
    integer,intent(IN) :: SYStype_in
    real(8),intent(IN) :: fmax_tol
    integer,optional,intent(IN) :: NiterSCF_in
    type(atom) :: ion
    type(lattice) :: aa
    integer,parameter :: max_loop=5
    integer :: np
    integer :: a, ierr, ip, jp, loop,i,j, icount
    real(8) :: etot, fmax
    real(8) :: dxdg,dgHdg,dxg,dgHg
    real(8) :: aa_inv(3,3)
    real(8),allocatable :: g(:,:,:),x(:,:,:),Hg(:,:,:),Hdg(:,:,:)
    real(8),allocatable :: history(:,:)
    real(8),allocatable :: dx(:,:),dg(:,:)

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

    call scf( etot, ierr ) ; if ( ierr == -1 ) goto 999
    call calc_force( ion%natom, ion%force, fmax )

    if ( fmax <= fmax_tol ) goto 900

! ---

    np = 20
    allocate( g(3,ion%natom,0:np)  ) ; g=0.0d0
    allocate( x(3,ion%natom,0:np)  ) ; x=0.0d0
    allocate( Hg(3,ion%natom,0:np) ) ; Hg=0.0d0
    allocate( Hdg(3,ion%natom,np)  ) ; Hdg=0.0d0
    allocate( dx(3,ion%natom) ) ; dx=0.0d0
    allocate( dg(3,ion%natom) ) ; dg=0.0d0

    allocate( history(0:max_loop*np,3) ) ; history=0.0d0

    history(0,1) = etot
    history(0,2) = fmax
    history(0,3) = ierr

! ---

    icount=0

    do loop=1,max_loop

       x(:,:,0)  = ion%xyz(:,:)
       g(:,:,0)  =-ion%force(:,:)
       Hg(:,:,0) = g(:,:,0)

       do ip=1,np

          if ( disp_sw ) write(*,'(a60," ICY=",2i4)') repeat("-",60),loop,ip

          ion%xyz(:,:) = x(:,:,ip-1) - Hg(:,:,ip-1)

          if ( disp_sw ) then
             do a=1,ion%natom
                write(*,*) a, sqrt(sum(Hg(:,a,ip-1)**2))
             end do
          end if

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

          call scf( etot, ierr ) ; if ( ierr == -1 ) goto 999
          call calc_force( ion%natom, ion%force, fmax )

          if ( fmax <= fmax_tol ) goto 900

          icount = icount + 1
          history(icount,1) = etot
          history(icount,2) = fmax
          history(icount,3) = ierr

          if ( ierr == -2 ) then
             ion%xyz(:,:)   =  x(:,:,ip-1)
             ion%force(:,:) = -g(:,:,ip-1)
             exit
          end if

          if ( disp_sw ) then
             do i=0,icount
                write(*,'(1x,i4,f20.10,es14.5,i4)') &
                     i, (history(i,j),j=1,2), nint(history(i,3))
             end do
          end if

          x(:,:,ip) = ion%xyz(:,:)
          g(:,:,ip) =-ion%force(:,:)

          Hg(:,:,ip) = g(:,:,ip)
          if ( ip == 1 ) Hdg(:,:,1) = Hg(:,:,1) - Hg(:,:,0)
          do jp=1,ip
             dx(:,:) = x(:,:,jp) - x(:,:,jp-1)
             dg(:,:) = g(:,:,jp) - g(:,:,jp-1)
             dxdg = sum( dx*dg )
             if ( disp_sw ) write(*,*) jp,dxdg
             dgHdg = sum( dg(:,:)*Hdg(:,:,jp) )
             dxg = sum( dx(:,:)*g(:,:,ip) )
             dgHg = sum( dg(:,:)*Hg(:,:,ip) )
             Hg(:,:,ip) = Hg(:,:,ip) + ( dxdg+dgHdg )/dxdg**2*dxg*dx(:,:) &
                  - ( Hdg(:,:,ip)*dxg + dx(:,:)*dgHg )/dxdg
             if ( jp == ip-1 ) Hdg(:,:,ip) = Hg(:,:,ip) - Hg(:,:,ip-1)
          end do

       end do ! ip

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

  END SUBROUTINE atomopt_bfgs


  SUBROUTINE scf( etot, ierr_out )
    implicit none
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

    call calc_scf( ierr_out, NiterSCF, Etot_out=etot )

  END SUBROUTINE scf


END MODULE atomopt_bfgs_module
