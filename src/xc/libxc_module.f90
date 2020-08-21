module libxc_module

#ifdef _LIBXC_
  use xc_f90_lib_m
#endif
  use io_tools_module
  use basic_type_factory
  use physical_type_methods
  use xc_variables, only: xcene, xcpot
  use grid_module,only: grid
  use gradient_module

  implicit none

  private
  public :: init_libxc
  public :: finalize_libxc
  public :: calc_libxc
  public :: calc_fxc_libxc

  character(len=256) :: func_string(2)
  integer :: func_id(2)
  integer :: num_func
#ifdef _LIBXC_
  type(xc_f90_func_t) :: xc_func(2)
  type(xc_f90_func_info_t) :: xc_info(2)
#endif
  logical :: flag_init = .false.

contains

#ifdef _LIBXC_
  subroutine read_libxc
    implicit none
    call IOTools_readStringKeyword( "LIBXC", func_string )
    func_id(1) = xc_f90_functional_get_number( func_string(1) )
    func_id(2) = xc_f90_functional_get_number( func_string(2) )
  end subroutine read_libxc
#endif

  subroutine init_libxc( spin )

    implicit none
    integer,intent(in) :: spin
    character(len=120) :: s1, s2
    logical :: disp_sw
    integer :: i, j
    real(8) :: coef

    if ( flag_init ) return
    flag_init = .true.

    call write_border( 0, " init_libxc(start)" )
    call check_disp_switch( disp_sw, 0 )

#ifndef _LIBXC_
    call stop_program( "LIBXC is not linked to this code!" )
#else

    call read_libxc
    num_func = count( func_id > 0 )
    if ( disp_sw ) then
       write(*,*) "func_id =",func_id
       write(*,*) "num_func=",num_func
    end if
    if ( all(func_id<=0) ) call stop_program( "stop@libxc(1)" )

    do i=1,num_func

       if ( spin == 2 ) then
          call xc_f90_func_init(xc_func(i),func_id(i),XC_POLARIZED)
       else
          call xc_f90_func_init(xc_func(i),func_id(i),XC_UNPOLARIZED)
       end if

       xc_info(i) = xc_f90_func_get_info(xc_func(i))

       if ( disp_sw ) then

          select case( xc_f90_func_info_get_kind(xc_info(i)) )
          case( XC_EXCHANGE )
             write(*,'(a)') "Exchange"
          case( XC_CORRELATION )
             write(*,'(a)') "Correlation"
          case( XC_EXCHANGE_CORRELATION )
             write(*,'(a)') "Exchange-correlation"
          case( XC_KINETIC )
             write(*,'(a)') "Kinetic"
          end select

          s1 = xc_f90_func_info_get_name( xc_info(i) )

          select case( xc_f90_func_info_get_family(xc_info(i)) )
          case( XC_FAMILY_LDA )
             write(s2,'(a)') "LDA"
          case( XC_FAMILY_GGA )
             write(s2,'(a)') "GGA"
             !if ( spin == 2 ) call stop_program( "spin==2 is not implemented")
          case( XC_FAMILY_HYB_GGA )
             write(s2,'(a)') "Hybrid GGA"
             !call xc_f90_hyb_gga_xc_hse_set_par( xc_func(i), 0.25d0, 0.11d0 )
             !call xc_f90_hyb_exx_coef( xc_func(i), coef )
             !write(*,*) "Ratio of EXX =",coef
             !if ( spin == 2 ) call stop_program( "spin==2 is not implemented")
          case( XC_FAMILY_MGGA )
             write(s2,'(a)') "MGGA"
          case( XC_FAMILY_HYB_MGGA )
             write(s2,'(a)') "Hybrid MGGA"
          case( XC_FAMILY_LCA )
             write(s2,'(a)') "LCA"
          end select

          write(*,'(4a)') trim(s1), ' (', trim(s2), ')'

          j=0
          s1 = xc_f90_func_reference_get_ref( xc_f90_func_info_get_references(xc_info(i),j) )
          do while( j >= 0 )
             write(*,'(a,i1,2a)') "[", j, "] ", trim(s1)
             s1 = xc_f90_func_reference_get_ref( xc_f90_func_info_get_references(xc_info(i),j) )
          end do

       end if

    end do ! i

    call write_border( 0, " init_libxc(end)" )
#endif
  end subroutine init_libxc


  subroutine finalize_libxc
    implicit none
    integer :: i
#ifdef _LIBXC_
    do i=1,num_func
       call xc_f90_func_end( xc_func(i) )
    end do
#endif
  end subroutine finalize_libxc


  subroutine calc_libxc( rgrid, density, ene, pot )
    implicit none
    type( grid ),intent(in) :: rgrid
    type( GSArray ),intent(in) :: density
    type( xcene ),optional :: ene
    type( xcpot ),optional :: pot
    integer :: m,n,nn,s
    real(8),allocatable :: exc(:), vxc(:,:)
    real(8),allocatable :: vrho(:,:), sigma(:,:), vsigma(:,:)
    real(8) :: s1(2)
#ifdef _LIBXC_
    call write_border( 1, " calc_libxc(start)" )

    m  = density%g_range%size
    n  = density%s_range%size
    nn = pot%xc%s_range%size

    allocate( exc(m)      ) ; exc=0.0d0
    allocate( vxc(m,n)    ) ; vxc=0.0d0
    allocate( vrho(1,1)   ) ; vrho=0.0d0
    allocate( sigma(1,1)  ) ; sigma=0.0d0
    allocate( vsigma(1,1) ) ; vsigma=0.0d0

    call calc_libxc_a( 1, density%val, exc, vxc, rgrid, density )

    if ( present(pot) ) then

       if ( allocated(pot%x%val) ) pot%x%val(:,:) = vxc(:,1:nn)

       pot%xc%val(:,:) = vxc(:,1:nn)

    end if

    if ( present(ene) ) then

       s1(1)=0.0d0
       do s=1,n
          s1(1) = s1(1) + sum( density%val(:,s)*exc(:) )
       end do

       call dSpatialIntegral( s1(1:1) )

       ene%Ex = s1(1)

       ene%Exc = ene%Ex

    end if

! ---

    if ( num_func >= 2 ) then

       call calc_libxc_a( 2, density%val, exc, vxc, rgrid, density )

       if ( present(pot) ) then

          if ( allocated(pot%c%val) ) pot%c%val(:,:) = vxc(:,1:nn)

          pot%xc%val(:,:) = pot%xc%val(:,:) + vxc(:,1:nn)

       end if

       if ( present(ene) ) then

          s1(2)=0.0d0
          do s=1,n
             s1(2) = s1(2) + sum( density%val(:,s)*exc(:) )
          end do

          call dSpatialIntegral( s1(2:2) )

          ene%Ec = s1(2)

          ene%Exc = ene%Exc + ene%Ec

       end if

    end if

    deallocate( vsigma )
    deallocate( sigma  )
    deallocate( vrho   )
    deallocate( vxc    )
    deallocate( exc    )

    call write_border( 1, " calc_libxc(end)" )
#endif
  end subroutine calc_libxc

#ifdef _LIBXC_
  subroutine calc_libxc_a( i, rho, exc, vxc, rgrid, density )
    implicit none
    integer,intent(in)  :: i
    real(8),intent(in)  :: rho(:,:)
    real(8),intent(out) :: exc(:), vxc(:,:)
    type(grid),intent(in) :: rgrid
    type(GSArray),intent(in) :: density
    real(8),allocatable :: vrho(:,:)
    real(8),allocatable :: sigma(:,:)
    real(8),allocatable :: vsigma(:,:)
    !real(8),allocatable :: a(:,:), b(:,:)
    real(8),allocatable :: a(:), b(:)
    integer(8) :: np
    integer :: ns,j,s
    type(gradient) :: grad

    np = size( rho, 1 )
    ns = size( rho, 2 )

    !allocate( a(ns,np) ) ; a=0.0d0
    !do j=1,np
    !   a( 1,j) = rho(j, 1)
    !   a(ns,j) = rho(j,ns)
    !end do
    allocate( a(ns*np) ); a=0.0d0
    do s=1,ns
       do j=1,np
          a(j+np*(s-1)) = rho(j, s)
       end do
    end do

    select case( xc_f90_func_info_get_family(xc_info(i)) )
    case( XC_FAMILY_LDA )

       !allocate( b(ns,np) ) ; b=0.0d0
       allocate( b(ns*np) ) ; b=0.0d0

       !call xc_f90_lda_exc_vxc( xc_func(i), np, a(1,1), exc(1), b(1,1) )
       call xc_f90_lda_exc_vxc( xc_func(i), np, a, exc, b )

       !do j=1,np
       !   vxc(j, 1) = b( 1,j)
       !   vxc(j,ns) = b(ns,j)
       !end do
       do s=1,ns
          do j=1,np
             vxc(j,s) = b(j+np*(s-1))
          end do
       end do

       deallocate( b )

    case( XC_FAMILY_GGA )

       allocate( sigma(2*ns-1,np)  ) ; sigma=0.0d0
       allocate( vsigma(2*ns-1,np) ) ; vsigma=0.0d0
       allocate( vrho(ns,np)       ) ; vrho=0.0d0

       call construct_gradient( rgrid, density, grad )

       sigma(1,:) = grad%gg(:)

       !call xc_f90_gga_exc_vxc( xc_func(i), np, a(1,1), sigma(1,1), &
       !                         exc(1), vrho(1,1), vsigma(1,1) )
       call xc_f90_gga_exc_vxc( xc_func(i), np, a, sigma, &
                                exc, vrho, vsigma )

       call calc_gga_pot( rgrid, vrho, vsigma, grad%gx, grad%gy, grad%gz, vxc )

       deallocate( vrho   )
       deallocate( vsigma )
       deallocate( sigma  )

    case( XC_FAMILY_HYB_GGA )

       allocate( sigma(2*ns-1,np)  ) ; sigma=0.0d0
       allocate( vsigma(2*ns-1,np) ) ; vsigma=0.0d0
       allocate( vrho(ns,np)       ) ; vrho=0.0d0

       call construct_gradient( rgrid, density, grad )

       !call xc_gga_exc_vxc( xc_func(i), np, a(1,1), sigma(1,1), &
       !                         exc(1), vrho(1,1), vsigma(1,1) )
       call xc_f90_gga_exc_vxc( xc_func(i), np, a, sigma, &
                                exc, vrho, vsigma )

       call calc_gga_pot( rgrid, vrho, vsigma, grad%gx, grad%gy, grad%gz, vxc )

       deallocate( vrho   )
       deallocate( vsigma )
       deallocate( sigma  )

    case( XC_FAMILY_MGGA ); goto 900
    case( XC_FAMILY_HYB_MGGA ); goto 900
    case( XC_FAMILY_LCA ); goto 900
    end select

    deallocate( a )

    return

900 call stop_program("stop@calc_libxc_a")

  end subroutine calc_libxc_a


  subroutine calc_gga_pot( rgrid, vrho, vsigma, gx, gy, gz, vxc )
    implicit none
    type(grid),intent(in) :: rgrid
    real(8),intent(in)  :: vrho(:,:), vsigma(:,:)
    real(8),intent(in)  :: gx(:), gy(:), gz(:)
    real(8),intent(out) :: vxc(:,:)
    real(8),allocatable :: f(:),g(:),h(:)
    integer :: m,n,i

    m = size( vxc, 1 )
    n = size( vxc, 2 )

    do i=1,n
       vxc(:,i) = vrho(i,:)
    end do

    allocate( f(m) ) ; f=0.0d0
    allocate( g(m) ) ; g=0.0d0
    allocate( h(m) ) ; h=0.0d0

    f(:) = vsigma(1,:)

    g(:) = f(:)*gx(:)
    call calc_xyz_gradient( 1, rgrid, g, h )
    vxc(:,1) = vxc(:,1) - 2.0d0*h(:)

    g(:) = f(:)*gy(:)
    call calc_xyz_gradient( 2, rgrid, g, h )
    vxc(:,1) = vxc(:,1) - 2.0d0*h(:)

    g(:) = f(:)*gz(:)
    call calc_xyz_gradient( 3, rgrid, g, h )
    vxc(:,1) = vxc(:,1) - 2.0d0*h(:)

    deallocate( h, g, f )

  end subroutine calc_gga_pot
#endif


  subroutine calc_fxc_libxc( rho, fxc )
    implicit none
    real(8),intent(in)  :: rho(:,:)
    real(8),intent(out) :: fxc(:,:)
    real(8),allocatable :: a(:,:), b(:,:)
    integer :: np4,ns,i,j,s
    integer(8) :: np

    fxc=0.0d0

#ifdef _LIBXC_

    np = size( rho, 1 )
    ns = size( rho, 2 )

    allocate( a(ns,np) ) ; a=0.0d0

    do j=1,np
       a( 1,j) = rho(j, 1)
       a(ns,j) = rho(j,ns)
    end do

    allocate( b(2*ns-1,np) ) ; b=0.0d0

    do i=1,num_func

       select case( xc_f90_func_info_get_family(xc_info(i)) )
       case( XC_FAMILY_LDA )

          call xc_f90_lda_fxc( xc_func(i), np, a(1,1), b(1,1) )

          do s=1,2*ns-1
             do j=1,np
                fxc(j,s) = fxc(j,s) + b(s,j)
             end do
          end do

       case default

          call stop_program( "fxc of LDA_X & LDA_C_PZ is only available" )

       end select

    end do

    deallocate( b )
    deallocate( a )

#endif
  end subroutine calc_fxc_libxc


end module libxc_module
