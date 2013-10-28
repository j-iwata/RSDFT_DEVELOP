MODULE sseig

  use global_variables, only: disp_switch
  use parallel_module
  use iter_lin_solvers, only: z_dot,comm_rhs,solver_timers,solver_info &
                             ,myrank_r,nprocs_r,seed_val,z_mgs_qr,z_matmat_cn,z_matmat_nx &
                             ,current_solver_info,current_solver_timers
  use hamiltonian_module
  use timer_module

  implicit none

  PRIVATE
  PUBLIC :: z_sseig_one_circle,my_ccurve_id,ccurve_id_rank,numeig,num_basis &
           ,n_ccurve_local,n_ccurve,rho_scale,ssrho,n_rhs_set_local &
           ,unit_ss_band_cor,unit_ss_band_val,ssgamma,SS_IO,L_fixed,opt &
           ,solver_timers_array,sseig_timers_array,solver_info_array,my_rhs_set_id &
           ,first_time,eigvec,eigval,ssres,inputV,eigvec_merged,eigval_merged &
           ,cmyrank,unit_ss_band_wf,file_ss_wf,file_ss_wf_split &
           ,old_idx,current_sseig_timers,n_rhs_set,current_ccurve_local_id &
           ,disp_sseig_result,z_rayleigh_ritz,z_normalize_eigenvectors,z_eval_residual &
           ,opt_N,out_iter_max,out_iter_init,resval_merged
  
  type optss
     integer :: N = 32
     integer :: M = 16
     integer :: L = 8
     real(8) :: delta = 1d-10
  end type optss

  type sseig_timers
     type(timer) :: total_time
     type(timer) :: ls_time
     type(timer) :: post_time
     type(timer) :: svd_time
     type(timer) :: orth_time
     type(timer) :: rr_time
     type(timer) :: mate_time
     type(timer) :: dsyev_time
     type(timer) :: rot_time
  end type sseig_timers

  type(sseig_timers),save :: current_sseig_timers

  integer :: n_ccurve, L_fixed = 32
  
  integer,allocatable :: n_ccurve_local(:)
  integer,allocatable :: my_ccurve_id(:)
  integer,allocatable :: ccurve_id_rank(:)

  integer :: n_rhs_set = 1

  integer,allocatable :: n_rhs_set_local(:)
  integer,allocatable :: my_rhs_set_id(:)
  
  integer,allocatable :: numeig(:), num_basis(:)

  real(8),allocatable :: ssgamma(:),ssrho(:)
  real(8),allocatable :: ssres(:,:)

  type(optss),save :: opt
  

  real(8),allocatable :: eigval(:,:), eigval_merged(:)
  real(8),allocatable :: resval_merged(:)
  complex(8),allocatable :: eigvec(:,:,:), inputV(:,:,:)
  complex(8),allocatable :: eigvec_merged(:,:)
  integer,allocatable :: old_idx(:)
  logical :: first_time = .true.
  
  real(8) :: rho_scale = 1.0d0

  integer :: current_ccurve_local_id, current_rhs_local_id

  type(solver_timers),allocatable :: solver_timers_array(:,:)
  type(solver_info),allocatable :: solver_info_array(:,:,:)
  type(sseig_timers),allocatable :: sseig_timers_array(:)

  integer :: unit_ss_band_val = 321
  integer :: unit_ss_band_cor = 1212
  integer :: unit_ss_band_wf  = 503

  character(len=5) :: cmyrank
  character(len=9) :: file_ss_wf = 'ss_wf.dat'
  character(len=32) :: file_ss_wf_split

  integer :: SS_IO = 0

  integer :: opt_N,out_iter_max,out_iter_init

CONTAINS


  SUBROUTINE z_sseig_one_circle(gamma,rho,opt,lambda,X,residual,LSfun,V,evnum)

    implicit none

    INTERFACE
       SUBROUTINE LSfun(Z,V,X)
         complex(8),intent(in) :: Z(:)
         complex(8),intent(in) :: V(:,:)
         complex(8),intent(inout) :: X(:,:,:)
       END SUBROUTINE LSfun
    END INTERFACE

    real(8),intent(in)     :: gamma
    real(8),intent(in)     :: rho
    type(optss),intent(in) :: opt

    real(8),intent(out)    :: lambda(:)
    complex(8),intent(out) :: X(:,:)      ! eigvec
    real(8),intent(out)    :: residual(:)
    complex(8),intent(in)  :: V(:,:)      ! inputV
    integer,intent(inout)  :: evnum

    complex(8),allocatable :: theta(:),omega(:),Y(:,:,:),Y_tmp(:,:,:)
    real(8) :: c,s,PI2
    integer :: i,count,ierr
    integer :: n_rhs_blk, n_rhs_blk_here, start_L, end_L, offset
    
    allocate( theta(opt%N), omega(opt%N) )

    PI2 = 2.d0*acos(-1.d0)
    do i = 1,opt%N
       c = cos( PI2*( (i-1)+0.5d0 )/opt%N )
       s = sin( PI2*( (i-1)+0.5d0 )/opt%N )
       theta(i) = dcmplx( c,s )
       omega(i) = gamma + rho * theta(i)
    end do

    seed_val = dble( omega(opt%N/2) )

    n_rhs_blk = (opt%L-1)/n_rhs_set + 1

    start_L = n_rhs_blk*((n_rhs_set-1)/nprocs_r+1)*myrank_r+1
    end_L   = min(n_rhs_blk*((n_rhs_set-1)/nprocs_r+1)*(myrank_r+1),opt%L)

    allocate(  Y( size(X,1), end_L-start_L+1, opt%N )  )

    do i = 1, n_rhs_set_local(myrank_r)

       offset = start_L + n_rhs_blk*(i-1)
       current_rhs_local_id = i

       n_rhs_blk_here = min(offset+n_rhs_blk-1,opt%L) - offset + 1

       allocate(  Y_tmp(size(X,1), n_rhs_blk_here, opt%N)  )

       if ( disp_switch ) write(*,'(a40," LSfun")') repeat("-",60)
       call tic(current_sseig_timers%ls_time,comm_grid)
       call LSfun(omega,V(:,offset:offset+n_rhs_blk_here-1),Y_tmp)
       call toc(current_sseig_timers%ls_time,comm_grid)

       Y(:,offset-start_L+1:offset-start_L+n_rhs_blk_here,:) = Y_tmp

       deallocate( Y_tmp )

       solver_timers_array(current_rhs_local_id,current_ccurve_local_id) &
            = current_solver_timers
       solver_info_array(:,current_rhs_local_id,current_ccurve_local_id) &
            = current_solver_info(:)

    end do

    if ( disp_switch ) write(*,'(a40," z_sseig_post_linear_system")') repeat("-",60)
    call tic(current_sseig_timers%post_time,comm_grid)
    call z_sseig_post_linear_system &
         (gamma,rho,opt,theta,omega,Y,lambda,X,residual,evnum)
    call toc(current_sseig_timers%post_time,comm_grid)

    deallocate( Y )
    deallocate( omega, theta )

  END SUBROUTINE z_sseig_one_circle


  SUBROUTINE z_sseig_post_linear_system &
       (gamma,rho,opt,theta,omega,Y,lambda,X,residual,evnum)
    implicit none
    real(8),intent(in)     :: gamma
    real(8),intent(in)     :: rho
    type(optss),intent(in) :: opt
    complex(8),intent(in)  :: theta(:),omega(:),Y(:,:,:)
    real(8),intent(out)    :: lambda(:)
    complex(8),intent(out) :: X(:,:)
    real(8),intent(out)    :: residual(:)
    integer,intent(inout)  :: evnum
    
    integer :: n1,n2,k
    complex(8),allocatable :: dummy(:,:)

    n1 = id_grid(myrank_g) + 1
    n2 = id_grid(myrank_g) + ir_grid(myrank_g)

    call z_calc_S( Y, theta, opt%N, opt%M, opt%L, X )

    call tic(current_sseig_timers%svd_time,comm_grid)
    call z_SVD_hankel( X, opt%delta, k )
    call toc(current_sseig_timers%svd_time,comm_grid)

    evnum = max(min(k,size(lambda,1)),evnum)

    allocate( dummy(evnum,evnum) )

    call tic(current_sseig_timers%orth_time,comm_grid)
    call z_MGS_QR( X(:,1:evnum),dummy )
    call z_MGS_QR( X(:,1:evnum),dummy )
    call toc(current_sseig_timers%orth_time,comm_grid)

    deallocate(dummy)

    call tic(current_sseig_timers%rr_time,comm_grid)
    call z_rayleigh_ritz( X(:,1:evnum),gamma,lambda(1:evnum), n1,n2 )
    call toc(current_sseig_timers%rr_time,comm_grid)

    call z_normalize_eigenvectors( X(:,1:evnum) )

    call z_eval_residual( lambda(1:evnum),X(:,1:evnum),residual(1:evnum), n1,n2 )

  END SUBROUTINE z_sseig_post_linear_system


!  SUBROUTINE calc_omega( gamma, rho, N, theta, omega )
!    complex(8),intent(in)  :: gamma 
!    real(8),intent(in)     :: rho
!    integer,intent(in)     :: N
!    complex(8),intent(out) :: theta(:),omega(:)
!    integer :: i
!    real(8) :: PI2,c,s
!    PI2 = 2.d0*acos(-1.d0)
!    do i = 1,N
!       c = cos( PI2*( (i-1)+0.5D0 )/N )
!       s = sin( PI2*( (i-1)+0.5D0 )/N )
!       theta(i) = dcmplx( c,s )
!       omega(i) = gamma + rho * theta(i)
!    end do
!  END SUBROUTINE calc_omega


  SUBROUTINE z_calc_S( Y, theta, N, M, L, S )
    implicit none
    complex(8),intent(in)  :: Y(:,:,:),theta(:)
    integer,intent(in)     :: N,M,L
    complex(8),intent(out) :: S(:,:)
    integer :: i,k,j,start_L,end_L,ierr,n_rhs_blk
    integer,allocatable  :: ircnt(:),idisp(:)
    complex(8),allocatable :: S_local(:,:)

    allocate( ircnt(0:nprocs_r-1) )
    allocate( idisp(0:nprocs_r-1) )

    n_rhs_blk = (opt%L-1)/n_rhs_set + 1

    do i=0,nprocs_r-1
       start_L = n_rhs_blk*((n_rhs_set-1)/nprocs_r+1)*i+1
       end_L   = min(n_rhs_blk*((n_rhs_set-1)/nprocs_r+1)*(i+1),opt%L)
       ircnt(i) = (end_L - start_L + 1)*(size(Y,1))
       idisp(i) = (start_L-1)*(size(Y,1))
    end do

    start_L = n_rhs_blk*((n_rhs_set-1)/nprocs_r+1)*myrank_r+1
    end_L   = min(n_rhs_blk*((n_rhs_set-1)/nprocs_r+1)*(myrank_r+1),opt%L)

    k=size(Y,1)
    allocate( S_local(k,end_L-start_L+1) )

    do k = 1,M
       S_local(:,:) = (0d0,0d0)
       do j = 1,N
          S_local(:,:) = S_local(:,:) + theta(j)**k * Y(:,:,j)
       end do
       S_local(:,:) = S_local(:,:) / dble(N)     
       call mpi_allgatherv( S_local, ircnt(myrank_r), mpi_complex16 &
            , S(1,(k-1)*L+1), ircnt, idisp, mpi_complex16, comm_rhs, ierr )
    end do

    deallocate( S_local )
    deallocate( idisp )
    deallocate( ircnt )

  END SUBROUTINE z_calc_S


  SUBROUTINE z_SVD_hankel( S, delta, k )
    implicit none
    complex(8),intent(in) :: S(:,:)
    real(8),intent(in) :: delta
    integer,intent(out) :: k

    integer :: ms,bs,ldwork,info,i,j
    real(8),allocatable :: SIGMA(:)
    real(8) :: max_sigma
    real(8),allocatable :: H(:,:),work(:)
    complex(8),allocatable :: Hc(:,:)
    real(8),allocatable :: U(:,:),VT(:,:)

    ms = size(S,1)
    bs = size(S,2)
    ldwork = 5*bs

    allocate(H(bs,bs), SIGMA(bs), work(ldwork))
    allocate(U(bs,bs), VT(bs,bs))
    allocate(Hc(bs,bs))

    call z_matmat_CN(S,S,Hc)

    H = dble(Hc)

    call DGESVD('A','A',bs,bs,H,bs,SIGMA,U,bs,VT,bs,work,ldwork,info)

    SIGMA = dsqrt(SIGMA)
    max_sigma = SIGMA(1)
    do i=1,bs
       if ( SIGMA(i) < delta*max_sigma ) then
          exit
       end if
    end do
    k = i - 1

    deallocate(H, SIGMA, work)
    deallocate(U, VT)
    deallocate(Hc)
  END SUBROUTINE z_SVD_hankel


  SUBROUTINE z_rayleigh_ritz( Q,gamma,lambda, n1,n2 )
    implicit none
    real(8),intent(in)       :: gamma
    complex(8),intent(inout) :: Q(:,:)
    integer,intent(in)       :: n1,n2
    real(8),intent(out)      :: lambda(:)

    integer :: ms,bs,i,lwork,info
    complex(8),parameter :: zero=(0.d0,0.d0)
    complex(8),allocatable :: AQ(:,:),A1(:,:),work(:),QQ(:,:)
    real(8),allocatable :: rwork(:)
    intrinsic size

    ms    = size(Q,1)
    bs    = size(Q,2)
    lwork = 2*bs-1

    allocate( AQ(ms,bs),A1(bs,bs),work(lwork),QQ(ms,bs),rwork(3*bs-2) )

    call tic(current_sseig_timers%mate_time,comm_grid)
    do i=1,bs
#ifndef _DRSDFT_
       call hamiltonian(1,1,Q(:,i),AQ(:,i),n1,n2,1,1)
#endif
    end do
    call z_matmat_CN(Q,AQ,A1)
    call toc(current_sseig_timers%mate_time,comm_grid)

    do i=1,bs
       A1(i,i) = A1(i,i) - gamma
    end do 

    call tic(current_sseig_timers%dsyev_time,comm_grid)
    call ZHEEV('V','U',bs,A1,bs,lambda,work,lwork,rwork,info)
    call toc(current_sseig_timers%dsyev_time,comm_grid)

    lambda(:) = lambda(:) + gamma

    call tic(current_sseig_timers%rot_time,comm_grid)
    call z_matmat_NX( 'N',zero,Q,A1,QQ )
    Q(:,:) = QQ(:,:)
    call toc(current_sseig_timers%rot_time,comm_grid)

    deallocate( AQ,A1,work,QQ,rwork )

  END SUBROUTINE z_rayleigh_ritz


  SUBROUTINE z_normalize_eigenvectors(X)
    implicit none
    complex(8),intent(inout) :: X(:,:)
    integer :: blk_size,i,j
    blk_size = size(X,2)
    do i=1,blk_size
       X(:,i) = X(:,i) / sqrt( z_dot(X(:,i),X(:,i)) )
    end do
  END SUBROUTINE z_normalize_eigenvectors

  
  SUBROUTINE z_eval_residual(lambda,X,residual,n1,n2)
    implicit none
    real(8),intent(in) :: lambda(:)
    complex(8),intent(in) :: X(:,:)
    real(8),intent(out) :: residual(:)
    integer,intent(in) :: n1,n2

    integer :: ms,bs,i,j,ierr
    real(8),allocatable :: tmp(:)
    complex(8),allocatable :: r(:),tmpAX(:,:)

    ms = size(X,1)
    bs = size(residual)

    allocate( r(ms),tmpAX(ms,1),tmp(bs) )

    do i=1,bs
#ifndef _DRSDFT_
       call hamiltonian(1,1,X(:,i),tmpAX(:,1),n1,n2,1,1)
#endif
       r = tmpAX(:,1) - lambda(i)*X(:,i)
       residual(i) = 0d0
       do j=1,ms
          tmp(i) = abs(r(j))
          residual(i) = residual(i) + tmp(i) * tmp(i)
       end do
       tmp(i) = residual(i)
    end do

    call mpi_allreduce(tmp,residual,bs,mpi_real8,mpi_sum,comm_grid,ierr)

    residual = dsqrt(residual)

    deallocate( r,tmpAX,tmp )

  END SUBROUTINE z_eval_residual


  SUBROUTINE disp_sseig_result( timers, gamma, rho, num_basis, eigval, res )    
    implicit none
    type(sseig_timers),intent(in) :: timers
    real(8),intent(in) :: gamma, rho, eigval(:), res(:)
    integer,intent(in) :: num_basis
    integer :: m
        
    write(*,'(3x,a)') ' --- SSEIG TIME --- '
    write(*,'(1x,a,1x,2f10.3)') 'total               ', timers%total_time%ct, timers%total_time%et
    write(*,'(1x,a,1x,2f10.3)') '  |--solver         ', timers%ls_time%ct, timers%ls_time%et
    write(*,'(1x,a,1x,2f10.3)') '  |--post           ', timers%post_time%ct, timers%post_time%et
    write(*,'(1x,a,1x,2f10.3)') '       |--SVD       ', timers%svd_time%ct, timers%svd_time%et
    write(*,'(1x,a,1x,2f10.3)') '       |--orth      ', timers%orth_time%ct, timers%orth_time%et
    write(*,'(1x,a,1x,2f10.3)') '       |--RR        ', timers%rr_time%ct, timers%rr_time%et
    write(*,'(1x,a,1x,2f10.3)') '            |--mate ', timers%mate_time%ct, timers%mate_time%et
    write(*,'(1x,a,1x,2f10.3)') '            |--zheev', timers%dsyev_time%ct, timers%dsyev_time%et
    write(*,'(1x,a,1x,2f10.3)') '            |--rot  ', timers%rot_time%ct, timers%rot_time%et
    write(*,*) 'num basis = ', num_basis
    write(*,'(1x,"center = ",e14.7,3x,"radius = ",e14.7)') gamma,rho
    write(*,'(1x,"interval = [",e14.7,",",e14.7,"]")') gamma-rho,gamma+rho
    write(*,'(1x,a4,3x,a25,3x,a14,3x,a21)') 'band','lambda','res'
    do m=1, num_basis
       write(*,'(1x,I4,3x,f25.15,3x)',advance='no') m,eigval(m)
       write(*,*) res(m)
    end do
  END SUBROUTINE disp_sseig_result


END MODULE sseig
