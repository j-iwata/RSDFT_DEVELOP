MODULE linear_response_module

  use momentum_module
  use wf_module
  use inner_product_module, only: calc_inner_product
  use rgrid_module
  use hamiltonian_module
  use density_module, only: rho
  use hartree_module
  use libxc_module, only: calc_fxc_libxc
  use parallel_module, only: comm_grid, comm_band, comm_bzsm, comm_spin
  use rsdft_mpi_module, only: rsdft_allreduce_sum
  use symmetry_module, only: sym_rho

  implicit none

  PRIVATE
  PUBLIC :: calc_dielectric_constant

  logical :: disp_sw

CONTAINS


  SUBROUTINE calc_dielectric_constant
    implicit none
    complex(8),parameter :: zero=(0.0d0,0.0d0)
#ifdef _DRSDFT_
    real(8),allocatable :: psi_bar(:,:,:,:)
    real(8),allocatable :: dlt_psi(:,:,:,:)
    real(8),allocatable :: work(:,:)
#else
    complex(8),allocatable :: psi_bar(:,:,:,:)
    complex(8),allocatable :: dlt_psi(:,:,:,:)
    complex(8),allocatable :: work(:,:)
#endif
    integer,parameter :: max_loop_scf=100
    integer :: mg,mb,mk,ms,n,k,s,mg0,mb0,mk0,ms0
    integer :: mvb,n0,k0,s0,loop_scf
    real(8) :: rtmp(8)
    real(8) :: Eexternal, epsilon, epsilon0, polarization
    real(8),allocatable :: dlt_rho(:,:),dlt_vhxc(:),fxc(:,:),rwork(:,:)

    call write_border( 0, "calc_dielectric_constant(start)" )
    call check_disp_switch( disp_sw, 0 )

    mg  = size( unk, 1 )
    mb  = size( unk, 2 )
    mk  = size( unk, 3 )
    ms  = size( unk, 4 )
    mg0 = ML_0_WF
    mb0 = MB_0_WF
    mk0 = MK_0_WF
    ms0 = MS_0_WF

    mvb = nint( sum(occ)*0.5d0 )

    allocate( psi_bar(mg,mb,mk,ms) ) ; psi_bar=zero
    allocate( dlt_psi(mg,mb,mk,ms) ) ; psi_bar=zero
    allocate( work(mg,1)           ) ; work=zero
    allocate( dlt_rho(mg,MS_WF)    ) ; dlt_rho=0.0d0
    allocate( dlt_vhxc(mg)         ) ; dlt_vhxc=0.0d0
    allocate( fxc(mg,MS_WF)        ) ; fxc=0.0d0
    allocate( rwork(mg,MS_WF)      ) ; rwork=0.0d0

    call set_initial_wf( "random", psi_bar )

    do s=1,ms
    do k=1,mk
    do n=1,mb
       n0=mb0+n-1
       k0=mk0+k-1
       s0=ms0+s-1
       if ( n0 > mvb ) cycle
       call op_momentum( "z", k0, unk(:,n0:n0,k0,s0), work )
       call ortho_valence( unk(:,1:mvb,k0,s0), work(:,1) )
       call ortho_valence( unk(:,1:mvb,k0,s0), psi_bar(:,n,k,s) )
       call solve_sternheimer( n0,k0,s0, work(:,1), psi_bar(:,n:n,k,s) )
       call ortho_valence( unk(:,1:mvb,k0,s0), psi_bar(:,n,k,s) )
    end do
    end do
    end do

! ---

    Eexternal = 0.1d0
    call set_initial_wf( "random", dlt_psi )

    do s=1,ms
    do k=1,mk
    do n=1,mb
       n0=mb0+n-1
       k0=mk0+k-1
       s0=ms0+s-1
       if ( n0 > mvb ) cycle
       work(:,1) = -Eexternal*psi_bar(:,n,k,s)
       call ortho_valence( unk(:,1:mvb,k0,s0), dlt_psi(:,n,k,s) )
       call solve_sternheimer( n0,k0,s0, work(:,1), dlt_psi(:,n:n,k,s) )
       call ortho_valence( unk(:,1:mvb,k0,s0), dlt_psi(:,n,k,s) )
    end do
    end do
    end do

! ---

    dlt_rho(:,:)=0.0d0
    do s=1,ms
    do k=1,mk
    do n=1,mb
       n0=mb0+n-1
       k0=mk0+k-1
       s0=ms0+s-1
       if ( n0 > mvb ) cycle
#ifdef _DRSDFT_
       dlt_rho(:,s0) = dlt_rho(:,s0) + occ(n0,k0,s0) &
            *unk(:,n0,k0,s0)*dlt_psi(:,n,k,s)
#else
       dlt_rho(:,s0) = dlt_rho(:,s0) + occ(n0,k0,s0) &
            *2.0d0*real( conjg(unk(:,n0,k0,s0))*dlt_psi(:,n,k,s) )
#endif
    end do
    end do
    end do
    call rsdft_allreduce_sum( dlt_rho, comm_band )
    call rsdft_allreduce_sum( dlt_rho, comm_bzsm )
    call rsdft_allreduce_sum( dlt_rho, comm_spin )
    call sym_rho( ML_0_WF,ML_1_WF,MS_WF,MS_0_WF,MS_1_WF,dlt_rho )

! ---

    call calc_fxc_libxc( rho, fxc )

    epsilon=0.0d0

    do loop_scf=1,max_loop_scf

       epsilon0=epsilon

       polarization=0.0d0
       do s=1,ms
       do k=1,mk
       do n=1,mb
          n0=mb0+n-1
          k0=mk0+k-1
          s0=ms0+s-1
          if ( n0 > mvb ) cycle
#ifdef _DRSDFT_
          polarization = polarization - occ(n0,k0,s0) &
               *sum( psi_bar(:,n,k,s)*dlt_psi(:,n,k,s) )
#else
          polarization = polarization - occ(n0,k0,s0) &
               *2.0d0*real( sum(conjg(psi_bar(:,n,k,s))*dlt_psi(:,n,k,s)) )
#endif
       end do
       end do
       end do
       polarization=polarization*dV/( dV*Ngrid(0) )
       call rsdft_allreduce_sum( polarization, comm_grid )
       call rsdft_allreduce_sum( polarization, comm_band )
       call rsdft_allreduce_sum( polarization, comm_bzsm )
       call rsdft_allreduce_sum( polarization, comm_spin )

       epsilon = 1.0d0 + 4.0d0*acos(-1.0d0)*polarization/Eexternal
       if ( disp_sw ) then
          write(*,'(1x,i5,2x,2f20.10)') loop_scf, epsilon, epsilon-epsilon0
       end if

       call calc_hartree( dlt_rho, dlt_vhxc )
       dlt_vhxc(:) = dlt_vhxc(:) + fxc(:,1)*dlt_rho(:,1)

       do s=1,ms
       do k=1,mk
       do n=1,mb
          n0=mb0+n-1
          k0=mk0+k-1
          s0=ms0+s-1
          if ( n0 > mvb ) cycle
          work(:,1) = -Eexternal*psi_bar(:,n,k,s) - dlt_vhxc(:)*unk(:,n0,k0,s0)
          call ortho_valence( unk(:,1:mvb,k0,s0), work(:,1) )
          call solve_sternheimer( n0,k0,s0, work(:,1), dlt_psi(:,n:n,k,s) )
          call ortho_valence( unk(:,1:mvb,k0,s0), dlt_psi(:,n,k,s) )
       end do
       end do
       end do

       rwork(:,:)=0.0d0
       do s=1,ms
       do k=1,mk
       do n=1,mb
          n0=mb0+n-1
          k0=mk0+k-1
          s0=ms0+s-1
          if ( n0 > mvb ) cycle
#ifdef _DRSDFT_
          rwork(:,s0) = rwork(:,s0) + occ(n0,k0,s0) &
               *unk(:,n0,k0,s0)*dlt_psi(:,n,k,s)
#else
          rwork(:,s0) = rwork(:,s0) + occ(n0,k0,s0) &
               *2.0d0*real( conjg(unk(:,n0,k0,s0))*dlt_psi(:,n,k,s) )
#endif
       end do
       end do
       end do
       call rsdft_allreduce_sum( rwork, comm_band )
       call rsdft_allreduce_sum( rwork, comm_bzsm )
       call rsdft_allreduce_sum( rwork, comm_spin )
       call sym_rho( ML_0_WF,ML_1_WF,MS_WF,MS_0_WF,MS_1_WF,rwork )
       do s=1,MS_WF
          dlt_rho(:,s) = dlt_rho(:,s) + 0.3d0*( rwork(:,s)-dlt_rho(:,s) )
       end do

    end do ! max_loop_scf

! ---

    deallocate( rwork    )
    deallocate( dlt_vhxc )
    deallocate( dlt_rho  )
    deallocate( work     )
    deallocate( dlt_psi  )
    deallocate( psi_bar  )

    call write_border( 0, "calc_dielectric_constant(end)" )

  END SUBROUTINE calc_dielectric_constant


  SUBROUTINE ortho_valence( u, v )
    implicit none
#ifdef _DRSDFT_
    real(8),intent(IN)    :: u(:,:)
    real(8),intent(INOUT) :: v(:)
    real(8) :: uv
#else
    complex(8),intent(IN)    :: u(:,:)
    complex(8),intent(INOUT) :: v(:)
    complex(8) :: uv
#endif
    integer :: n
    do n=1,size(u,2)
       call calc_inner_product( u(:,n), v(:), uv )
       v(:) = v(:) - u(:,n)*uv
    end do
  END SUBROUTINE ortho_valence


  SUBROUTINE solve_sternheimer( n, k, s, rhs, psi )
    implicit none
    integer,intent(IN)       :: n, k, s
#ifdef _DRSDFT_
    real(8),intent(IN)    :: rhs(:)
    real(8),intent(INOUT) :: psi(:,:)
    real(8),allocatable :: p(:,:), r(:,:), Ap(:,:), xmin(:)
    real(8) :: r0r0,rr,pAp,bb,ztmp,alpha,beta
#else
    complex(8),intent(IN)    :: rhs(:)
    complex(8),intent(INOUT) :: psi(:,:)
    complex(8),allocatable :: p(:,:), r(:,:), Ap(:,:), xmin(:)
    complex(8) :: r0r0,rr,pAp,bb,ztmp,alpha,beta
#endif
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    integer,parameter :: maxcg=1000
    integer :: m, icg
    real(8),parameter :: tol=1.d-8
    real(8) :: err,errmin !,alpha,beta

    call calc_inner_product( rhs, rhs, bb )

    m =size( psi, 1 )
    allocate( r(m,1)  ) ; r=zero
    allocate( p(m,1)  ) ; p=zero
    allocate( Ap(m,1) ) ; Ap=zero
    allocate( xmin(m) ) ; xmin=zero

    call op_matrix( n,k,s, psi, r )
    r(:,1) = rhs - r(:,1)
    p = r

    call calc_inner_product( r(:,1), r(:,1), r0r0 )

    err = sqrt( abs(r0r0)/abs(bb) )
    errmin = err
    if ( err <= tol ) goto 900

    xmin = psi(:,1)

    do icg=1,maxcg

       call op_matrix( n,k,s, p, Ap )
       call calc_inner_product( p(:,1), Ap(:,1), pAp )

!------
!(1)
       alpha = r0r0/pAp
!(2)
!       call calc_inner_product( r, p, ztmp )
!       alpha = ztmp/pAp
!------

       psi = psi + alpha*p

       r = r - alpha*Ap

       call calc_inner_product( r(:,1), r(:,1), rr )

       err = sqrt( abs(rr)/abs(bb) )
       if ( err < errmin ) then
          errmin = err
          xmin = psi(:,1)
       end if
       if ( err <= tol ) exit

!------
!(1)
       beta = rr/r0r0
!(2)
!       call calc_inner_product( r, Ap, ztmp )
!       beta = -ztmp/pAp
!------

!       if ( mod(icg,100)==0 ) then
!       if ( abs(rr) > abs(r0r0) ) then
!          write(*,*) icg, err
!          write(*,*) alpha,rr,beta
!       end if

       p = r + beta*p

       r0r0 = rr

    end do ! icg

900 continue

    if ( err > tol ) then
!       write(*,*) "solve sternheimer",icg, err, errmin
       if ( err > errmin ) psi(:,1) = xmin
    end if

    deallocate( Ap )
    deallocate( p  )
    deallocate( r  )

  END SUBROUTINE solve_sternheimer

  SUBROUTINE op_matrix( n, k, s, x, Ax )
    implicit none
    integer,intent(IN)     :: n,k,s
#ifdef _DRSDFT_
    real(8),intent(IN)  :: x(:,:)
    real(8),intent(OUT) :: Ax(:,:)
#else
    complex(8),intent(IN)  :: x(:,:)
    complex(8),intent(OUT) :: Ax(:,:)
#endif
    integer :: n1,n2
    n1=Igrid(1,0)
    n2=Igrid(2,0)
    call hamiltonian( k, s, x, Ax, n1, n2, 1, 1 )
    Ax = Ax - esp(n,k,s)*x
  END SUBROUTINE op_matrix


END MODULE linear_response_module
