module atomopt_ef_module

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

  private
  public :: atomopt_ef

  logical :: disp_sw, disp_scf
  integer :: SYStype
  integer :: NiterSCF

contains


  subroutine atomopt_ef( SYStype_in, fmax_tol, ncycle, NiterSCF_in )

    implicit none
    integer,intent(in) :: SYStype_in
    real(8),intent(in) :: fmax_tol
    integer,intent(in) :: ncycle
    integer,optional,intent(in) :: NiterSCF_in
    type(atom) :: ion
    type(lattice) :: aa, bb
    integer :: max_loop
    integer :: ishape(1), ishape2(2), i1,i2, LWORK=0
    integer :: n,il,iu,m,a,ierr,loop,i,j,icount,ip
    integer,allocatable :: iwork(:),ifail(:)
    real(8) :: etot0, etot, fmax
    real(8) :: dxdg,dxHdx,c1,c2,dmax,dtmp
    real(8) :: vl,vu,tol,alpha
    real(8) :: aa_inv(3,3),da(3),da_tmp(3)
    real(8),parameter :: one=1.0d0, zero=0.0d0
    real(8),allocatable :: g(:),x(:),Hessian(:,:),Htmp(:,:)
    real(8),allocatable :: Vtmp(:,:),Wtmp(:,:)
    real(8),allocatable :: history(:,:),g0(:),x0(:),d(:)
    real(8),allocatable :: dx(:),dg(:),Hdx(:)
    real(8),allocatable :: w(:),z(:,:),work(:)
    real(8),parameter :: eig_threshold=0.02d0

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

    n = 3*ion%natom
    allocate( g(n)  ) ; g=0.0d0
    allocate( x(n)  ) ; x=0.0d0
    allocate( d(n)  ) ; d=0.0d0
    allocate( g0(n) ) ; g0=0.0d0
    allocate( x0(n) ) ; x0=0.0d0
    allocate( Hessian(n,n) ) ; Hessian=0.0d0
    allocate( Htmp(n,n) ) ; Htmp=0.0d0
    allocate( Vtmp(n,n) ) ; Vtmp=0.0d0
    allocate( Wtmp(n,n) ) ; Wtmp=0.0d0
    allocate( Hdx(n)    ) ; Hdx=0.0d0
    allocate( dx(n)     ) ; dx=0.0d0
    allocate( dg(n)     ) ; dg=0.0d0

    ishape(1) = n
    ishape2(1:2) = (/ 3, ion%natom /)

    allocate( history(0:max_loop,4) ) ; history=0.0d0

    history(0,1) = etot
    history(0,2) = fmax
    history(0,3) = ierr

! --- LAPACK

    if ( LWORK == 0 ) then
      allocate( w(n) ); w=0.0d0
      allocate( work(1) )
      call dsyev('V','U',n+1,Htmp,n+1,w,work,-1,ierr)
      LWORK=nint( work(1) )
      deallocate( work )
      allocate( work(LWORK) ); work=0.0d0
    end if

    icount = 0

    do loop=1,max_loop

      if ( disp_sw ) write(*,'(a60," ICY=",i4)') repeat("-",60),loop

      alpha = 1.0d0

      if ( loop == 1 ) then

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

        if ( dxdg <= 0.0d0 ) then

          d(:)  = -g(:)
          x0(:) = x(:)
          g0(:) = g(:)

          Hessian=0.0d0
          do i=1,size(Hessian,1)
            Hessian(i,i)=1.0d0
          end do

        else

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
              write(*,*) "i,w(i)",i,w(i),"  ---> replaced to",eig_threshold
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

          Wtmp=matmul( Vtmp, Hessian )
          write(*,*) "1sum(Wtmp**2)",sum(Wtmp**2)
          Wtmp=matmul( Hessian, Vtmp )
          write(*,*) "2sum(Wtmp**2)",sum(Wtmp**2)

          x0(:) = x(:)
          g0(:) = g(:)
        
          call DGEMV('N',n,n,-1.0d0,Vtmp,n,g(:),1,0.0d0,d,1)

        end if

      end if

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

      if ( disp_sw ) then
        write(*,*) "Next trial configuration"
      end if
      dmax=0.0d0
      do a=1,ion%natom
        i2 = 3*a
        i1 = i2 - 3 + 1
        dtmp = sqrt(sum(d(i1:i2)**2))*alpha
        if ( disp_sw ) then
          write(*,'(i2,2x,3f10.5,2x,3es14.5,2x,es14.5)') &
               a,aa_atom(:,a),ion%xyz(:,a),dtmp
        end if
        dmax=max(dmax,dtmp)
      end do
      if ( disp_sw ) then
        write(*,*) "Maximum displacement(bohr):",dmax
      end if

      if ( dmax > 2.0d0 ) then
        alpha=alpha/dmax
      end if

      x(:) = x0(:) + alpha*d(:)

      ion%xyz(:,:) = reshape( x, ishape2 )

      aa_atom(:,:) = matmul( aa_inv, ion%xyz )
      call shift_aa_coordinates_atom( aa_atom )        
          
      call write_coordinates_atom( 97, 3 )

      call scf( etot, ierr ) ; if ( ierr == -1 ) goto 999
      call calc_force( ion%natom, ion%force, fmax )

      icount = icount + 1
      history(icount,1) = etot
      history(icount,2) = fmax
      history(icount,3) = ierr
      history(icount,4) = alpha

      if ( disp_sw ) then
        do i=0,icount
          write(*,'(1x,i4,f20.10,es14.5,i4,es14.5)') &
               i, (history(i,j),j=1,2), nint(history(i,3)),history(i,4)
        end do
      end if

      if ( fmax <= fmax_tol ) goto 900

      g(:) = -reshape( ion%force(:,:), ishape )

    end do ! loop

! ---

900 if ( disp_sw ) then
       write(*,'(1x,"etot    :",f25.15)') etot
       write(*,'(1x,"fmax/tol:",es12.5," /",es12.5)') fmax,fmax_tol
    end if

999 call check_disp_switch( disp_sw, 1 )
    call write_border( 0, "atomopt_ef(end)" )

!    deallocate( ifail )
!    deallocate( iwork )
    deallocate( work )
!    deallocate( z )
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


end module atomopt_ef_module
