MODULE iter_lin_solvers

  use parallel_module
  use hamiltonian_module
  use timer_module

  implicit none

  PRIVATE
  PUBLIC :: solver_info,solver_timers,seed_val,nprocs_r,myrank_r &
           ,comm_rhs,z_dot, myrank_c,nprocs_c,unit_kry_conv,comm_ccurve &
           ,unroll_step,orth_impl,iterative_method_args &
           ,solver,z_mgs_qr,z_matmat_cn,z_matmat_nx,disp_lin_result &
           ,current_solver_info,current_solver_timers,tol_iter

  type itargs
     real(8) :: epsmax = 1D-10 ! input
     real(8) :: epsmin = 1D-10 ! input
     integer :: imax = 3000 ! input
     integer :: iter ! output
     real(8) :: residual ! output
  end type itargs

  type solver_info
     real(8) :: shift_r
     real(8) :: shift_i
     integer :: iter
     real(8) :: rec_res
     real(8) :: true_res
  end type solver_info

  type solver_timers
     type(timer) :: total_time
     type(timer) :: shift_time
     type(timer) :: hpsi_time
     type(timer) :: lin_orth_time
     type(timer) :: scalar_time
     type(timer) :: input_time
     type(timer) :: gemm_time
     type(timer) :: out_time
  end type solver_timers

  type(itargs), save :: iterative_method_args

  type(solver_info),save, allocatable :: sol_info

  real(8) :: seed_val
  integer :: shift_count
  type(timer),save :: shift_time, hpsi_time, lin_orth_time &
                    , scalar_time, input_time, gemm_time, out_time
  type(solver_timers),save :: current_solver_timers
  type(solver_info),allocatable :: current_solver_info(:)

  integer :: comm_ccurve, comm_rhs, nprocs_c, nprocs_r, myrank_c, myrank_r  
  integer :: orth_impl = 0, unroll_step = 16

  integer :: unit_kry_conv = 123;

  complex(8),allocatable :: X_seed(:,:)

  real(8) :: tol_iter(4)

CONTAINS


  SUBROUTINE solver(Z,B,X)
    implicit none
    complex(8),intent(in)    :: Z(:)      ! complex energy
    complex(8),intent(in)    :: B(:,:)    ! rhs vectors
    complex(8),intent(inout) :: X(:,:,:)  ! solution vectors
    real(8) :: epsmax,epsmin
    integer :: imax,i,j
    
    epsmax = iterative_method_args%epsmax
    epsmin = iterative_method_args%epsmin
    imax   = iterative_method_args%imax

    if ( disp_switch_parallel ) write(*,*) "epsmax,epsmin,imax=",epsmax,epsmin,imax

    call clear_timer(current_solver_timers%total_time)
    call tic(current_solver_timers%total_time,comm_grid)
    if ( unroll_step > 0 ) then
       call shifted_block_CG_rQ_unroll_iter(B,Z,epsmax,epsmin,imax,X)
    else
       call shifted_block_CG_rQ_naive(B,Z,epsmax,epsmin,imax,X)
    end if
    call toc(current_solver_timers%total_time,comm_grid)

    if ( nprocs_r == 1 .and. nprocs_c == 1 .and. myrank_g == 0 ) then
       write(*,*) 'Result of shifted block CG rQ method'
       call disp_lin_result(current_solver_info, size(Z), current_solver_timers)
    end if

  END SUBROUTINE solver


  SUBROUTINE shifted_block_CG_rQ_naive(B,sigma,epsmax,epsmin,imax,X)
    implicit none
    integer, intent(in)     :: imax
    real(8), intent(in)     :: epsmax,epsmin
    complex(8), intent(in)  :: B(:,:)    ! rhs vectors
    complex(8), intent(in)  :: sigma(:)  ! complex energy
    complex(8), intent(out) :: X(:,:,:)  ! solution vectors

    integer :: i, j, k, n, n1, n2, L, m, ierr, maxidx, conv_count, iter
    integer, allocatable :: ipiv(:),iter_array(:)    
    real(8) :: res_true, res_temp
    real(8), allocatable :: res_seed(:), Bnorm_inv(:), residual(:,:)
    complex(8) :: zero = (0d0,0d0), one = (1d0,0d0)
    complex(8), allocatable :: &
           Q(:,:), P(:,:,:), P_seed(:,:), ts_tmp(:,:) &
         , AP(:,:), PAP(:,:), PAP_old(:,:), alpha(:,:) &
         , delta(:,:), rho(:,:), rho_old(:,:), eye(:,:), blk_tmp1(:,:) &
         , blk_tmp2(:,:), work(:) &
         , xi1(:,:,:), xi1_hat(:,:,:), xi2(:,:,:), xitld_hat(:,:,:)

    call clear_timer(current_solver_timers%shift_time)
    call clear_timer(current_solver_timers%hpsi_time)
    call clear_timer(current_solver_timers%lin_orth_time)
    call clear_timer(current_solver_timers%scalar_time)
    call clear_timer(current_solver_timers%input_time)
    call clear_timer(current_solver_timers%gemm_time)
    call clear_timer(current_solver_timers%out_time)

    n = size(X,1) ! # of grid points
    L = size(X,2) ! # of rhs vectors (blocks?)
    m = size(X,3) ! # of complex energies

    n1 = id_grid(myrank_g)+1
    n2 = id_grid(myrank_g)+ir_grid(myrank_g)

    if ( .not.allocated(X_seed) ) then
       allocate( X_seed(n,L) )
       X_seed(1:n,1:L) = (0d0,0d0)
    end if

    allocate( Q(n,L),P_seed(n,L),P(n,L,m),AP(n,L),ts_tmp(n,L) &
             ,PAP(L,L),PAP_old(L,L),alpha(L,L),delta(L,L),rho(L,L) &
             ,rho_old(L,L),eye(L,L),blk_tmp1(L,L),blk_tmp2(L,L) &
             ,xi1(L,L,m),xi1_hat(L,L,m),xi2(L,L,m),xitld_hat(L,L,m) &
             ,ipiv(L),work(L),Bnorm_inv(L),res_seed(L) &
             ,iter_array(m),residual(L,m) )

    X(1:n,1:L,1:m) = (0d0,0d0)
    eye(1:L,1:L) = (0d0,0d0)
    do j=1,L
       eye(j,j) = (1d0,0d0)
    end do

    Q(1:n,1:L) = B(1:n,1:L)
    call z_MGS_QR(Q,delta)
    rho(1:L,1:L) = delta(1:L,1:L)
    P_seed(1:n,1:L) = Q(1:n,1:L)
    do k=1,m
       P(1:n,1:L,k) = Q(1:n,1:L)
       xi1(1:L,1:L,k) = rho(1:L,1:L)
       xi2(1:L,1:L,k) = eye(1:L,1:L)
    end do

    xitld_hat(1:L,1:L,1:m) = (0d0,0d0)
    call norm2_as_block(B,Bnorm_inv,n,L)
    Bnorm_inv(1:L) = 1.d0 / Bnorm_inv(1:L)
    PAP(1:L,1:L) = (0d0,0d0)
    residual(1:L,1:m) = 1.d0

    do iter=1,imax

       call tic(current_solver_timers%hpsi_time,comm_grid)
!       call hamiltonian(1,1,P_seed,AP,n1,n2,1,L)
       do i=1,L
#ifndef _DRSDFT_
          call hamiltonian(1,1,P_seed(1,i),AP(1,i),n1,n2,i,i)
#endif
       end do
       call toc(current_solver_timers%hpsi_time,comm_grid)

       AP(1:n,1:L) = seed_val*P_seed(1:n,1:L) - AP(1:n,1:L)

       PAP_old(1:L,1:L) = PAP(1:L,1:L)
       call z_matmat_CN(P_seed,AP,PAP)
       call inv(PAP,alpha,ipiv,work)
       call z_matmat_NX('N',zero,alpha,delta,blk_tmp1)
       call z_matmat_NX('N',one,P_seed,blk_tmp1,X_seed)

       rho_old(1:L,1:L) = rho(1:L,1:L)
       
       call z_matmat_NX('N',one,AP,-alpha,Q)

       call tic(current_solver_timers%lin_orth_time,comm_grid)
       if ( orth_impl == 1) then
          call z_MGS_QR(Q,rho)
       else if ( orth_impl == 2 ) then
          call z_TS_QR(Q,rho)
       else
          call z_CGS_QR_Takahashi(Q,rho)
       end if
       call toc(current_solver_timers%lin_orth_time,comm_grid)

       call z_matmat_NX('N',zero,rho,delta,blk_tmp1)
       delta(1:L,1:L) = blk_tmp1(1:L,1:L)

       call z_matmat_NX('C',zero,P_seed,rho,ts_tmp)
       P_seed(1:n,1:L) = Q(1:n,1:L) + ts_tmp(1:n,1:L)

       call norm2_as_block_serial(delta,res_seed,L,L)
       res_seed(1:L) = res_seed(1:L)*Bnorm_inv(1:L)

       call tic(current_solver_timers%shift_time,comm_grid)
       do k=1,m
          if ( maxval(residual(1:L,k)) > epsmin ) then             
             call tic(current_solver_timers%scalar_time,comm_grid)
             xi2(1:L,1:L,k) = xi1(1:L,1:L,k)
             blk_tmp1(1:L,1:L) = eye(1:L,1:L) - xitld_hat(1:L,1:L,k)
             call z_matmat_NX('N',zero,rho_old,blk_tmp1,blk_tmp2)
             call z_matmat_NX('N',zero,blk_tmp2,PAP_old,blk_tmp1)
             call z_matmat_NX('C',zero,blk_tmp1,rho_old,blk_tmp2)
             call z_matmat_NX('N',zero,blk_tmp2,alpha,blk_tmp1)
             call inv(eye-(seed_val-sigma(k))*alpha+blk_tmp1 &
                     ,xitld_hat(:,:,k),ipiv,work)
             call z_matmat_NX('N',zero,xitld_hat(:,:,k),xi2(:,:,k),xi1_hat(:,:,k))
             call z_matmat_NX('N',zero,rho,xi1_hat(:,:,k),xi1(:,:,k))
             call z_matmat_NX('N',zero,alpha,xi1_hat(:,:,k),blk_tmp1)
             call toc(current_solver_timers%scalar_time,comm_grid)

             call tic(current_solver_timers%gemm_time,comm_grid)
             call z_matmat_NX('N',one,P(:,:,k),blk_tmp1,X(:,:,k))
             call toc(current_solver_timers%gemm_time,comm_grid)

             call tic(current_solver_timers%scalar_time,comm_grid)
             call z_matmat_NX('N',zero,alpha,xitld_hat(:,:,k),blk_tmp1)
             call z_matmat_NX('N',zero,blk_tmp1,PAP,blk_tmp2)
             call z_matmat_NX('C',zero,blk_tmp2,rho,blk_tmp1)
             call toc(current_solver_timers%scalar_time,comm_grid)

             call tic(current_solver_timers%input_time,comm_grid)
             call zcopy(L*n,Q,1,ts_tmp,1)
             call toc(current_solver_timers%input_time,comm_grid)

             call tic(current_solver_timers%gemm_time,comm_grid)
             call z_matmat_NX('N',one,P(:,:,k),blk_tmp1,ts_tmp)
             call toc(current_solver_timers%gemm_time,comm_grid)

             call tic(current_solver_timers%out_time,comm_grid)
             call zcopy(L*n,ts_tmp,1,P(:,:,k),1)
             call toc(current_solver_timers%out_time,comm_grid)

             call norm2_as_block_serial(xi1(:,:,k),residual(:,k),L,L)
             residual(:,k) = residual(:,k)*Bnorm_inv(:)
             shift_count = shift_count + 1
          end if

          if ( maxval(residual(1:L,k)) > epsmin ) then
             iter_array(k) = iter
          end if

       end do ! k
       call toc(current_solver_timers%shift_time,comm_grid)
       
       if ( myrank_g == 0 ) then
          write(unit_kry_conv,'(I8,1x)',advance='no') iter
          do k=1,m
             write(unit_kry_conv,'(1pe8.1,1x)',advance='no') maxval(residual(:,k))
          end do
          write(unit_kry_conv,*)
       end if

       conv_count = 0
       do k=1,m
          if ( iter_array(k) < iter ) then
             conv_count = conv_count + 1
          end if
       end do
       if ( maxval(residual) < epsmax .or. conv_count > m/2 ) then
     !  if ( maxval(residual) < epsmax ) then
          exit
       end if

    end do ! iter

    do i=1,m
       res_temp = 0d0
       maxidx = maxloc(residual(1:L,i),1)
#ifndef _DRSDFT_
       call hamiltonian(1,1,X(:,maxidx,i),AP(:,maxidx),n1,n2,1,1)
#endif
       AP(:,maxidx) = B(:,maxidx) - (sigma(i)*X(:,maxidx,i) - AP(:,maxidx))
          
       res_true = dsqrt(dble(z_dot(AP(:,maxidx),AP(:,maxidx))))
       res_true = res_true*Bnorm_inv(maxidx)

       current_solver_info(i)%shift_r  = dble(sigma(i))
       current_solver_info(i)%shift_i  = dimag(sigma(i))
       current_solver_info(i)%iter     = iter_array(i)
       current_solver_info(i)%rec_res  = residual(maxidx,i)
       current_solver_info(i)%true_res = res_true
    end do

  END SUBROUTINE shifted_block_CG_rQ_naive


  SUBROUTINE shifted_block_CG_rQ( B,sigma,epsmax,epsmin,imax,X )
    implicit none
    integer, intent(in)     :: imax
    real(8), intent(in)     :: epsmax,epsmin
    complex(8), intent(in)  :: B(:,:), sigma(:)
    complex(8), intent(out) :: X(:,:,:)

    integer :: i, j, k, n, n1, n2, L, m, iter, ierr, maxidx, conv_count, is, ie
    integer, allocatable :: ipiv(:),iter_array(:)    
    real(8) :: res_true, res_temp
    real(8), allocatable :: res_seed(:), Bnorm_inv(:), residual(:,:)
    complex(8) :: zero = (0d0,0d0), one = (1d0,0d0)
    complex(8), allocatable :: &
           Q(:,:,:), P(:,:,:), P_seed(:,:), ts_tmp(:,:) &
         , AP(:,:), PAP(:,:), PAP_old(:,:), alpha(:,:), delta(:,:) &
         , rho(:,:), rho_old(:,:), eye(:,:), blk_tmp1(:,:) &
         , blk_tmp2(:,:), work(:) &
         , xi1(:,:,:), xi1_hat(:,:,:), xi2(:,:,:), xitld_hat(:,:,:) &
         , XP_coef(:,:,:,:), ts_tmp_dble(:,:), blk_tmp_dble1(:,:) &
         , blk_tmp_dble2(:,:), Q_cache(:,:,:)
    logical, allocatable :: conv_flag(:)
    integer :: unroll_grid = 128

    call clear_timer(current_solver_timers%shift_time)
    call clear_timer(current_solver_timers%hpsi_time)
    call clear_timer(current_solver_timers%lin_orth_time)
    call clear_timer(current_solver_timers%scalar_time)
    call clear_timer(current_solver_timers%input_time)
    call clear_timer(current_solver_timers%gemm_time)
    call clear_timer(current_solver_timers%out_time)

    n = size(X,1)
    L = size(X,2)
    m = size(X,3)

    n1 = id_grid(myrank_g) + 1
    n2 = id_grid(myrank_g) + ir_grid(myrank_g)

    if ( .not.allocated(X_seed) ) then
       allocate( X_seed(n,L) )
       X_seed(1:n,1:L) = (0d0,0d0)
    end if

    allocate( Q(n,L,unroll_step+1),P_seed(n,L),P(n,L,m),AP(n,L) &
             ,ts_tmp(n,L),PAP(L,L),PAP_old(L,L),alpha(L,L) &
             ,delta(L,L),rho(L,L),rho_old(L,L),eye(L,L),blk_tmp1(L,L) &
             ,blk_tmp2(L,L) &
             ,xi1(L,L,m),xi1_hat(L,L,m),xi2(L,L,m),xitld_hat(L,L,m) &
             ,ipiv(L),work(L),Bnorm_inv(L),res_seed(L), iter_array(m) &
             ,XP_coef(L,2*L,unroll_step+1,m),ts_tmp_dble(unroll_grid,2*L) &
             ,Q_cache(unroll_grid,L,unroll_step+1) &
             ,blk_tmp_dble1(L,2*L),blk_tmp_dble2(L,2*L),conv_flag(m) &
             ,residual(L,m) )

    X(:,:,:)    = (0d0,0d0)
    eye(:,:)    = (0d0,0d0)
    do j=1,L
       eye(j,j) = (1d0,0d0)
    end do

    Q(:,:,unroll_step+1) = B
    call z_MGS_QR( Q(:,:,unroll_step+1), delta )

    rho(:,:) = delta(:,:)
    P_seed(:,:) = Q(:,:,unroll_step+1)
    
    do k=1,m
       P(:,:,k) = Q(:,:,unroll_step+1)
       xi1(:,:,k) = rho(:,:)
       xi2(:,:,k) = eye(:,:)
    end do
    xitld_hat(:,:,:) = (0d0,0d0)
    call norm2_as_block(B,Bnorm_inv,n,L)
    Bnorm_inv(:) = 1d0 / Bnorm_inv(:)
    PAP(:,:) = (0d0,0d0)
    residual(:,:) = 1d0
    conv_flag(:)  = .false.
    iter_array(:) = unroll_step+1

    do iter=1,imax

       call tic(current_solver_timers%hpsi_time,comm_grid)
#ifndef _DRSDFT_
       call hamiltonian(1,1,P_seed,AP,n1,n2,1,L)
#endif
       call toc(current_solver_timers%hpsi_time,comm_grid)

       AP(:,:) = seed_val*P_seed(:,:) - AP(:,:)
       PAP_old(:,:) = PAP(:,:)
       call z_matmat_CN(P_seed,AP,PAP)
       call inv(PAP,alpha,ipiv,work)
       call z_matmat_NX('N',zero,alpha,delta,blk_tmp1)
       call z_matmat_NX('N',one,P_seed,blk_tmp1,X_seed)

       rho_old(:,:) = rho(:,:)
       
       Q(:,:,modulo(iter-1,unroll_step+1)+1) &
            = Q(:,:,modulo(iter-2,unroll_step+1)+1)

       call z_matmat_NX('N',one,AP,-alpha,Q(:,:,modulo(iter-1,unroll_step+1)+1))

       call tic(current_solver_timers%lin_orth_time,comm_grid)
       if ( orth_impl == 1) then
          call z_MGS_QR(Q(:,:,mod(iter-1,unroll_step+1)+1),rho)
       else if ( orth_impl == 2 ) then
          call z_TS_QR(Q(:,:,mod(iter-1,unroll_step+1)+1),rho)
       else
          call z_CGS_QR_Takahashi(Q(:,:,mod(iter-1,unroll_step+1)+1),rho)
       end if
       call toc(current_solver_timers%lin_orth_time,comm_grid)

       call z_matmat_NX('N',zero,rho,delta,blk_tmp1)
       delta(:,:) = blk_tmp1(:,:)

       call z_matmat_NX('C',zero,P_seed,rho,ts_tmp)
       P_seed(:,:) = Q(:,:,mod(iter-1,unroll_step+1)+1) + ts_tmp(:,:)

       call norm2_as_block_serial(delta,res_seed,L,L)
       res_seed(:) = res_seed(:)*Bnorm_inv(:)

       call tic(current_solver_timers%shift_time,comm_grid)
       do k=1,m
          if ( .not. conv_flag(k) ) then
             call tic(current_solver_timers%scalar_time,comm_grid)
             xi2(:,:,k) = xi1(:,:,k)
             blk_tmp1(:,:) = eye(:,:) - xitld_hat(:,:,k)
             call z_matmat_NX('N',zero,rho_old,blk_tmp1,blk_tmp2)
             call z_matmat_NX('N',zero,blk_tmp2,PAP_old,blk_tmp1)
             call z_matmat_NX('C',zero,blk_tmp1,rho_old,blk_tmp2)
             call z_matmat_NX('N',zero,blk_tmp2,alpha,blk_tmp1)

             blk_tmp1(:,:) = eye(:,:) + blk_tmp1(:,:)
             blk_tmp1(:,:) = blk_tmp1(:,:) + (sigma(k) - seed_val)*alpha(:,:)
             call inv(blk_tmp1,xitld_hat(:,:,k),ipiv,work)
             call z_matmat_NX('N',zero,xitld_hat(:,:,k),xi2(:,:,k),xi1_hat(:,:,k))
             call z_matmat_NX('N',zero,rho,xi1_hat(:,:,k),xi1(:,:,k))
             call z_matmat_NX('N',zero,alpha,xi1_hat(:,:,k),blk_tmp1)
             XP_coef(:,1:L,mod(iter-1,unroll_step+1)+1,k) = blk_tmp1(:,:)
             
             call z_matmat_NX('N',zero,alpha,xitld_hat(:,:,k),blk_tmp1)
             call z_matmat_NX('N',zero,blk_tmp1,PAP,blk_tmp2)
             call z_matmat_NX('C',zero,blk_tmp2,rho,blk_tmp1)
             XP_coef(:,L+1:2*L,mod(iter-1,unroll_step+1)+1,k) = blk_tmp1(:,:)
             
             call norm2_as_block_serial(xi1(:,:,k),residual(:,k),L,L)
             residual(:,k) = residual(:,k)*Bnorm_inv(:)
             call toc(current_solver_timers%scalar_time,comm_grid)
          end if
       end do
       call toc(current_solver_timers%shift_time,comm_grid)
       
       if ( mod(iter,unroll_step+1) == 0 ) then
          call tic(current_solver_timers%shift_time,comm_grid)
          call tic(current_solver_timers%scalar_time,comm_grid)
          do k=1,m
             if ( .not. conv_flag(k) ) then
                do j=unroll_step,1,-1
                   call z_matmat_NX('N',zero,XP_coef(:,L+1:2*L,j,k) &
                                    ,XP_coef(:,:,j+1,k),blk_tmp_dble1)
                   XP_coef(:,1:L,j,k) = XP_coef(:,1:L,j,k) + blk_tmp_dble1(:,1:L)
                   XP_coef(:,L+1:2*L,j,k) = blk_tmp_dble1(:,L+1:2*L)
                end do
             end if
          end do
          call toc(current_solver_timers%scalar_time,comm_grid)

          do i=1,n,unroll_grid
             is = i
             ie = min(i+unroll_grid-1,n)
             Q_cache(1:ie-is+1,:,:) = Q(is:ie,:,:)
             do k=1,m
                if ( .not. conv_flag(k) ) then
                   call tic(current_solver_timers%input_time,comm_grid)
                   ts_tmp_dble(1:ie-is+1,1:L) = X(is:ie,:,k)
                   ts_tmp_dble(1:ie-is+1,L+1:2*L) = Q_cache(1:ie-is+1,:,unroll_step+1)
                   call toc(current_solver_timers%input_time,comm_grid)
                   call tic(current_solver_timers%gemm_time,comm_grid)
                   do j=unroll_step,1,-1
                      call z_matmat_NX('N',one,Q_cache(1:ie-is+1,:,j) &
                           ,XP_coef(:,:,j+1,k),ts_tmp_dble(1:ie-is+1,:))
                   end do
                   call z_matmat_NX('N',one,P(is:ie,:,k),XP_coef(:,:,1,k) &
                        ,ts_tmp_dble(1:ie-is+1,:))
                   call toc(current_solver_timers%gemm_time,comm_grid)
                   call tic(current_solver_timers%out_time,comm_grid)
                   X(is:ie,:,k) = ts_tmp_dble(1:ie-is+1,1:L)
                   P(is:ie,:,k) = ts_tmp_dble(1:ie-is+1,L+1:2*L)
                   call toc(current_solver_timers%out_time,comm_grid)
                end if
             end do
          end do
          call toc(current_solver_timers%shift_time,comm_grid)
             
          conv_count = 0
          do k=1,m
             conv_flag(k) = maxval(residual(:,k)) < epsmin
             if ( maxval(residual(:,k)) > epsmin ) then
                iter_array(k) = iter + unroll_step + 1
             end if
             
             if ( iter_array(k) < iter ) then
                conv_count = conv_count + 1
             end if
          end do

          if ( maxval(residual) < epsmax .or. conv_count > m/2 ) then
             exit
          end if          
       end if
    end do
    
    do i=1,m
       res_temp = 0d0
       maxidx = maxloc(residual(:,i),1)
#ifndef _DRSDFT_
       call hamiltonian(1,1,X(:,maxidx,i),AP(:,maxidx),n1,n2,1,1)
#endif
       AP(:,maxidx) = B(:,maxidx) - (sigma(i)*X(:,maxidx,i) - AP(:,maxidx))
       
       res_true = dsqrt(dble(z_dot(AP(:,maxidx),AP(:,maxidx))))
       res_true = res_true*Bnorm_inv(maxidx)

       current_solver_info(i)%shift_r = dble(sigma(i))
       current_solver_info(i)%shift_i = dimag(sigma(i))
       current_solver_info(i)%iter = iter_array(i)
       current_solver_info(i)%rec_res = residual(maxidx,i)
       current_solver_info(i)%true_res = res_true
    end do
  END SUBROUTINE shifted_block_CG_rQ


  SUBROUTINE shifted_block_CG_rQ_unroll_iter(B,sigma,epsmax,epsmin,imax,X)
    implicit none
    integer, intent(in) :: imax
    real(8), intent(in) :: epsmax,epsmin
    complex(8), intent(in) :: B(:,:), sigma(:)
    complex(8), intent(out) :: X(:,:,:)

    integer :: i, j, k, n, n1, n2, L, m, iter, ierr, maxidx, conv_count
    integer, allocatable :: ipiv(:),iter_array(:)    
    real(8) :: res_true, res_temp
    real(8), allocatable :: res_seed(:), Bnorm_inv(:), residual(:,:)
    complex(8) :: zero = (0d0,0d0), one = (1d0,0d0)
    complex(8), allocatable :: &
           Q(:,:,:), P(:,:,:), P_seed(:,:), ts_tmp(:,:) &
         , AP(:,:), PAP(:,:), PAP_old(:,:), alpha(:,:) &
         , delta(:,:), rho(:,:), rho_old(:,:), eye(:,:) &
         , blk_tmp1(:,:), blk_tmp2(:,:), work(:) &
         , xi1(:,:,:), xi1_hat(:,:,:), xi2(:,:,:), xitld_hat(:,:,:) &
         , XP_coef(:,:,:,:), ts_tmp_dble(:,:), blk_tmp_dble1(:,:), blk_tmp_dble2(:,:)
    logical, allocatable :: conv_flag(:)

    call clear_timer(current_solver_timers%shift_time)
    call clear_timer(current_solver_timers%hpsi_time)
    call clear_timer(current_solver_timers%lin_orth_time)
    call clear_timer(current_solver_timers%scalar_time)
    call clear_timer(current_solver_timers%input_time)
    call clear_timer(current_solver_timers%gemm_time)
    call clear_timer(current_solver_timers%out_time)

    n = size(X,1)
    L = size(X,2)
    m = size(X,3)

    n1 = id_grid(myrank_g)+1
    n2 = id_grid(myrank_g)+ir_grid(myrank_g)

    if ( .not.allocated(X_seed) ) then
       allocate( X_seed(n,L) )
       X_seed(1:n,1:L) = (0d0,0d0)
    end if

    allocate(Q(n,L,unroll_step+1),P_seed(n,L),P(n,L,m),AP(n,L) &
         ,ts_tmp(n,L),PAP(L,L),PAP_old(L,L),alpha(L,L) &
         ,delta(L,L),rho(L,L),rho_old(L,L),eye(L,L),blk_tmp1(L,L),blk_tmp2(L,L) &
         ,xi1(L,L,m),xi1_hat(L,L,m),xi2(L,L,m),xitld_hat(L,L,m) &
         ,ipiv(L),work(L),Bnorm_inv(L),res_seed(L), iter_array(m) &
         ,XP_coef(L,2*L,unroll_step+1,m),ts_tmp_dble(n,2*L) &
         ,blk_tmp_dble1(L,2*L),blk_tmp_dble2(L,2*L),conv_flag(m),residual(L,m))

    X = (0d0,0d0)
    eye = (0d0,0d0)
    do j=1,L
       eye(j,j) = (1d0,0d0)
    end do

    Q(:,:,unroll_step+1) = B
    call z_MGS_QR(Q(:,:,unroll_step+1),delta)

    rho = delta
    P_seed = Q(:,:,unroll_step+1)
    
    do k=1,m
       P(:,:,k) = Q(:,:,unroll_step+1)
       xi1(:,:,k) = rho
       xi2(:,:,k) = eye
    end do
    xitld_hat = (0d0,0d0)
    call norm2_as_block(B,Bnorm_inv,n,L)
    Bnorm_inv = 1d0 / Bnorm_inv
    PAP = (0d0,0d0)
    residual = 1d0
    conv_flag = .false.
    iter_array = unroll_step + 1

    do iter=1,imax
       call tic(current_solver_timers%hpsi_time,comm_grid)
!       call hamiltonian(1,1,P_seed,AP,n1,n2,1,L)
       do i=1,L
#ifndef _DRSDFT_
          call hamiltonian(1,1,P_seed(1,i),AP(1,i),n1,n2,i,i)
#endif
       end do
       call toc(current_solver_timers%hpsi_time,comm_grid)
       AP = seed_val*P_seed - AP
       PAP_old = PAP
       call z_matmat_CN(P_seed,AP,PAP)
       call inv(PAP,alpha,ipiv,work)
       call z_matmat_NX('N',zero,alpha,delta,blk_tmp1)
       call z_matmat_NX('N',one,P_seed,blk_tmp1,X_seed)

       rho_old = rho
       
       Q(:,:,modulo(iter-1,unroll_step+1)+1) = Q(:,:,modulo(iter-2,unroll_step+1)+1)

       call z_matmat_NX('N',one,AP,-alpha,Q(:,:,modulo(iter-1,unroll_step+1)+1))

       call tic(current_solver_timers%lin_orth_time,comm_grid)
       if ( orth_impl == 1) then
          call z_MGS_QR(Q(:,:,mod(iter-1,unroll_step+1)+1),rho)
       else if ( orth_impl == 2 ) then
          call z_TS_QR(Q(:,:,mod(iter-1,unroll_step+1)+1),rho)
       else
          call z_CGS_QR_Takahashi(Q(:,:,mod(iter-1,unroll_step+1)+1),rho)
       end if
       call toc(current_solver_timers%lin_orth_time,comm_grid)

       call z_matmat_NX('N',zero,rho,delta,blk_tmp1)
       delta = blk_tmp1

       call z_matmat_NX('C',zero,P_seed,rho,ts_tmp)
       P_seed = Q(:,:,mod(iter-1,unroll_step+1)+1) + ts_tmp

       call norm2_as_block_serial(delta,res_seed,L,L)
       res_seed = res_seed*Bnorm_inv

       call tic(current_solver_timers%shift_time,comm_grid)
       do k=1,m
          if ( .not. conv_flag(k) ) then
             call tic(current_solver_timers%scalar_time,comm_grid)
             xi2(:,:,k) = xi1(:,:,k)
             blk_tmp1 = eye - xitld_hat(:,:,k)
             call z_matmat_NX('N',zero,rho_old,blk_tmp1,blk_tmp2)
             call z_matmat_NX('N',zero,blk_tmp2,PAP_old,blk_tmp1)
             call z_matmat_NX('C',zero,blk_tmp1,rho_old,blk_tmp2)
             call z_matmat_NX('N',zero,blk_tmp2,alpha,blk_tmp1)
             call inv(eye - (seed_val - sigma(k))*alpha + blk_tmp1,xitld_hat(:,:,k),ipiv,work)
             call z_matmat_NX('N',zero,xitld_hat(:,:,k),xi2(:,:,k),xi1_hat(:,:,k))
             call z_matmat_NX('N',zero,rho,xi1_hat(:,:,k),xi1(:,:,k))
             call z_matmat_NX('N',zero,alpha,xi1_hat(:,:,k),blk_tmp1)
             XP_coef(:,1:L,mod(iter-1,unroll_step+1)+1,k) = blk_tmp1
             
             call z_matmat_NX('N',zero,alpha,xitld_hat(:,:,k),blk_tmp1)
             call z_matmat_NX('N',zero,blk_tmp1,PAP,blk_tmp2)
             call z_matmat_NX('C',zero,blk_tmp2,rho,blk_tmp1)
             XP_coef(:,L+1:2*L,mod(iter-1,unroll_step+1)+1,k) = blk_tmp1
             
             call norm2_as_block_serial(xi1(:,:,k),residual(:,k),L,L)
             residual(:,k) = residual(:,k)*Bnorm_inv
             call toc(current_solver_timers%scalar_time,comm_grid)
          end if
       end do
       call toc(current_solver_timers%shift_time,comm_grid)

       if ( myrank_g == 0 ) then
          write(unit_kry_conv,'(I8,1x)',advance='no') iter
          do k=1,m
             write(unit_kry_conv,'(1pe8.1,1x)',advance='no') maxval(residual(:,k))
          end do
          write(unit_kry_conv,*)
       end if
       
       if ( mod(iter,unroll_step+1) == 0 ) then
          call tic(current_solver_timers%shift_time,comm_grid)
          call tic(current_solver_timers%scalar_time,comm_grid)
          do k=1,m
             if ( .not. conv_flag(k) ) then
                do j=unroll_step,1,-1
                   blk_tmp_dble2 = XP_coef(:,:,j+1,k)
                   call z_matmat_NX('N',zero,XP_coef(:,L+1:2*L,j,k),blk_tmp_dble2 &
                        ,blk_tmp_dble1)
                   XP_coef(:,1:L,j,k) = XP_coef(:,1:L,j,k) + blk_tmp_dble1(:,1:L)
                   XP_coef(:,L+1:2*L,j,k) = blk_tmp_dble1(:,L+1:2*L)
                end do
             end if
          end do
          call toc(current_solver_timers%scalar_time,comm_grid)

          do k=1,m
             if ( .not. conv_flag(k) ) then
                call tic(current_solver_timers%input_time,comm_grid)
                call zcopy(L*n,X(:,:,k),1,ts_tmp_dble(:,1:L),1)
                call zcopy(L*n,Q(:,:,unroll_step+1),1,ts_tmp_dble(:,L+1:2*L),1)
                call toc(current_solver_timers%input_time,comm_grid)
                call tic(current_solver_timers%gemm_time,comm_grid)
                do j=unroll_step,1,-1
                   call z_matmat_NX('N',one,Q(:,:,j),XP_coef(:,:,j+1,k),ts_tmp_dble)
                end do
                call z_matmat_NX('N',one,P(:,:,k),XP_coef(:,:,1,k),ts_tmp_dble)
                call toc(current_solver_timers%gemm_time,comm_grid)
                call tic(current_solver_timers%out_time,comm_grid)
                call zcopy(L*n,ts_tmp_dble(:,1:L),1,X(:,:,k),1)
                call zcopy(L*n,ts_tmp_dble(:,L+1:2*L),1,P(:,:,k),1)
                call toc(current_solver_timers%out_time,comm_grid)
             end if
          end do
          call toc(current_solver_timers%shift_time,comm_grid)
             
          conv_count = 0
          do k=1,m
             conv_flag(k) = maxval(residual(:,k)) < epsmin
             if ( maxval(residual(:,k)) > epsmin ) then
                iter_array(k) = iter + unroll_step + 1
             end if
             
             if ( iter_array(k) < iter ) then
                conv_count = conv_count + 1
             end if
          end do

          if ( maxval(residual) < epsmax .or. conv_count > m/2 ) then
             exit
          end if          
       end if
    end do
    
    do i=1,m
       res_temp = 0d0
       maxidx = maxloc(residual(:,i),1)
#ifndef _DRSDFT_
       call hamiltonian(1,1,X(:,maxidx,i),AP(:,maxidx),n1,n2,1,1)
#endif
       AP(:,maxidx) = B(:,maxidx) - (sigma(i)*X(:,maxidx,i) - AP(:,maxidx))
       
       res_true = dsqrt(dble(z_dot(AP(:,maxidx),AP(:,maxidx))))
       res_true = res_true*Bnorm_inv(maxidx)

       current_solver_info(i)%shift_r = dble(sigma(i))
       current_solver_info(i)%shift_i = dimag(sigma(i))
       current_solver_info(i)%iter = iter_array(i)
       current_solver_info(i)%rec_res = residual(maxidx,i)
       current_solver_info(i)%true_res = res_true
    end do
  END SUBROUTINE shifted_block_CG_rQ_unroll_iter


  SUBROUTINE shifted_CG_as_block(B,sigma,epsmax,epsmin,imax,X,iter,residual)    
    implicit none
    integer, intent(in) :: imax
    real(8), intent(in) :: epsmax,epsmin
    complex(8), intent(in) :: B(:,:), sigma(:)
    integer, intent(out) :: iter
    real(8), intent(out) :: residual(:,:)
    complex(8), intent(out) :: X(:,:,:)
    
    integer :: i, j, k, n, n1, n2, L, m, ierr, maxidx, conv_count

    integer, allocatable :: iter_array(:)
    real(8) :: res_true, res_temp
    real(8), allocatable :: Bnorm(:), res_seed(:)
    complex(8), allocatable :: R(:,:), RR(:), RR0(:), P(:,:,:), P_seed(:,:) &
         , AP(:,:), alpha(:,:), beta(:,:), pi0(:,:), pi1(:,:), pi2(:,:)
    complex(8) :: alpha_tmp, beta_tmp

    shift_count = 0
    call clear_timer(shift_time)
    call clear_timer(hpsi_time)
    call clear_timer(lin_orth_time)

    n = size(X,1)
    L = size(X,2)
    m = size(X,3)

    n1 = id_grid(myrank_g) + 1
    n2 = id_grid(myrank_g) + ir_grid(myrank_g)

    if ( .not.allocated(X_seed) ) then
       allocate( X_seed(n,L) )
       X_seed(1:n,1:L) = (0d0,0d0)
    end if

    allocate(Bnorm(L),R(n,L),RR(L),RR0(L),P_seed(n,L),P(n,L,m),AP(n,L) &
         ,alpha(L,imax),beta(L,imax),pi0(L,m),pi1(L,m),pi2(L,m),res_seed(L))
    allocate(iter_array(m))
    iter_array = 0

!---Initialization---

    X = (0D0,0D0)
    R = B
    P_seed = (0D0,0D0)
    P = (0D0,0D0)
    call z_dot_as_block(R,R,RR,n,L)
    alpha(:,1) = (1D0,0D0)
    pi0 = (1D0,0D0)
    pi1 = (1D0,0D0)
    beta(:,1) = (0D0,0D0)
    call norm2_as_block(R,Bnorm,n,L)
    residual = (1D0,0D0)
    Bnorm = 1D0/Bnorm    
    
!---Iteration---

    do iter=2,imax
       if ( mod(iter,100) == 0 .and. myrank_g == 0 ) then
          write(*,*) 'shifted  iter = ',iter
       end if
       
       !---Seed system---
       do i=1,n
          P_seed(i,:) = R(i,:) + beta(:,iter-1)*P_seed(i,:)
       end do

       call tic(current_solver_timers%hpsi_time,comm_grid)
#ifndef _DRSDFT_
       call hamiltonian(1,1,P_seed,AP,n1,n2,1,L)
#endif
       call toc(current_solver_timers%hpsi_time,comm_grid)

       AP = seed_val*P_seed - AP

       call z_dot_as_block(P_seed,AP,alpha(:,iter),n,L)
       alpha(:,iter) = RR / alpha(:,iter)
       do i=1,n
          X_seed(i,:) = X_seed(i,:) + alpha(:,iter)*P_seed(i,:)
       end do
       !----------------
       
       !---Shifted systems---
       call tic(current_solver_timers%shift_time,comm_grid)
       do i=1,m
          if ( maxval(residual(:,i)) > epsmin ) then
             do j=1,L
                pi2(j,i) = (1d0 - alpha(j,iter)*(seed_val - sigma(i)))*pi1(j,i) + &
                     (beta(j,iter-1)/alpha(j,iter-1))*alpha(j,iter)*(pi1(j,i) - pi0(j,i))
                beta_tmp = (pi0(j,i)/pi1(j,i))*(pi0(j,i)/pi1(j,i))*beta(j,iter-1)
                alpha_tmp = (pi1(j,i)/pi2(j,i))*alpha(j,iter)
                P(:,j,i) = (1d0/pi1(j,i))*R(:,j) + beta_tmp*P(:,j,i)
                X(:,j,i) = X(:,j,i) + alpha_tmp*P(:,j,i)
                shift_count = shift_count + 1
             end do
          end if
          if ( maxval(residual(:,i)) > epsmax ) then
             iter_array(i) = iter
          end if
       end do
       pi0 = pi1
       pi1 = pi2
       call toc(current_solver_timers%shift_time,comm_grid)
       
       do i=1,n
          R(i,:) = R(i,:) - alpha(:,iter)*AP(i,:)
       end do
       RR0 = RR
       call z_dot_as_block(R,R,RR,n,L)
       beta(:,iter) = RR / RR0
       call norm2_as_block(R,res_seed,n,L)
       res_seed = res_seed*Bnorm
       
       !---------------------
       
       do i=1,m
          do j=1,L
             if ( residual(j,i) > epsmin ) then
                residual(j,i) = abs((1d0/pi2(j,i)))*res_seed(j)
             end if
          end do
       end do

       !---------------------
       !---Convergence test---
       conv_count = 0
       do i=1,m
          if ( iter_array(i) < iter ) then
             conv_count = conv_count + 1
          end if
       end do

       if ( maxval(residual) <= epsmax .or. conv_count > m/2 ) THEN
          exit
       end if
       !----------------------
    end do

    if ( myrank_g == 0 ) then
       write(*,*) 'Result of shifted CG method'
       write(*,'(a5,1x,a10,1x,a10,1x,a6,1x,a15,1x,a4,1x,a15)') &
            'id','real','imag','iter','rec res','rhs','true res'
    end if

    do i=1,m
       res_temp = 0d0
       maxidx = maxloc(residual(:,i),1)
#ifndef _DRSDFT_
       call hamiltonian(1,1,X(:,maxidx,i),AP(:,maxidx),n1,n2,1,1)
#endif
       AP(:,maxidx) = B(:,maxidx) - (sigma(i)*X(:,maxidx,i) - AP(:,maxidx))
             
       res_true = dsqrt(dble(z_dot(AP(:,maxidx),AP(:,maxidx))))
       res_true = res_true*Bnorm(maxidx)
       if ( myrank_g == 0 ) then
          write(*,'(I5,1x,1pe10.2,1x,1pe10.2,1x,I6,1x,1pe15.7,1x,I4,1x,1pe15.7)') &
          i,real(sigma(i)),aimag(sigma(i)),iter_array(i),residual(maxidx,i),maxidx,res_true
       end if
    end do
  END SUBROUTINE shifted_CG_as_block  


  FUNCTION d_dot(x,y)
    implicit none
    real(8), intent(in) :: x(:),y(:)
    real(8) :: d_dot
    integer :: i,ierr
    real(8) :: tmp
    tmp = 0d0
    do i=1,size(x)
       tmp = tmp + x(i)*y(i)
    end do
    call mpi_allreduce(tmp,d_dot,1,mpi_real8,mpi_sum,comm_grid,ierr)
  END FUNCTION d_dot


  FUNCTION z_dot(x,y)
    complex(8), intent(in) :: x(:),y(:)
    complex(8) :: z_dot
    integer :: i,ierr
    complex(8) :: tmp
    complex(8) :: zdotc
    tmp = (0d0,0d0)
    do i=1,size(x)
       tmp = tmp + dconjg(x(i))*y(i)
    end do
    call mpi_allreduce(tmp,z_dot,1,mpi_complex16,mpi_sum,comm_grid,ierr)
  END FUNCTION z_dot


  SUBROUTINE z_dot_as_block(X,Y,dot,n,L)
    implicit none
    complex(8), intent(in) :: X(:,:),Y(:,:)
    integer, intent(in) :: n,L
    complex(8), intent(out) :: dot(:)
    integer :: i, ierr
    complex(8) :: temp(L)
    dot(:) = (0.D0,0.D0)
    do i=1,n
       dot(:) = dot(:) + dconjg(X(i,:))*Y(i,:)
    end do
    temp = dot
    call mpi_allreduce(temp,dot,L,mpi_complex16,mpi_sum,comm_grid,ierr)
  END SUBROUTINE z_dot_as_block


  SUBROUTINE norm2_as_block(V,norm,n,L)
    implicit none
    complex(8), intent(in) :: V(:,:)
    integer, intent(in) :: n,L
    real(8), intent(out) :: norm(:)
    integer :: i,ierr
    real(8) :: temp(L)
    norm(:) = 0D0
    do i=1,n
       temp(:) = abs(V(i,:))
       norm(:) = norm(:) + temp(:)*temp(:)
    end do
    temp(:) = norm(:)
    call mpi_allreduce(temp,norm,L,mpi_real8,mpi_sum,comm_grid,ierr)
    norm(:) = dsqrt( dble( norm(:) ) )
  END SUBROUTINE norm2_as_block


  SUBROUTINE norm2_as_block_serial(V,norm,n,L)
    implicit none
    complex(8), intent(in) :: V(:,:)
    integer, intent(in) :: n,L
    real(8), intent(out) :: norm(:)
    integer :: i,ierr
    real(8) :: temp(L)
    norm = 0D0
    do i=1,n
       temp = abs(V(i,:))
       norm = norm + temp*temp
    end do
    norm = dsqrt(dble(norm))
  END SUBROUTINE norm2_as_block_serial


  SUBROUTINE z_matmat_NX(trans,beta,X,Y,mat)
    character, intent(in) :: trans
    complex(8), intent(in) :: beta,X(:,:),Y(:,:)
    complex(8), intent(out) :: mat(:,:)
    integer :: ierr,m,n,k
    complex(8) :: one
    one = (1d0,0d0)
    m = size(X,1)
    n = size(mat,2)
    k = size(X,2)
    call ZGEMM('N',trans,m,n,k,one,X,m,Y,k,beta,mat,m)
  END SUBROUTINE z_matmat_NX


  SUBROUTINE d_matmat_TN(X,Y,mat)
    real(8), intent(in) :: X(:,:),Y(:,:)
    real(8), intent(out) :: mat(:,:)
    integer :: ierr,ms,bsX,bsY
    real(8) :: zero,one
    real(8),allocatable :: tmp(:,:)
    zero = 0d0; one = 1d0
    ms = size(X,1)
    bsX = size(X,2)
    bsY = size(Y,2)
    allocate(tmp(bsX,bsY))
    call DGEMM('T','N',bsX,bsY,ms,one,X,ms,Y,ms,zero,tmp,bsX)
    call mpi_allreduce(tmp,mat,bsX*bsY,mpi_real8,mpi_sum,comm_grid,ierr)
    deallocate(tmp)
  END SUBROUTINE d_matmat_TN


  SUBROUTINE z_matmat_CN(X,Y,mat)
    complex(8), intent(in) :: X(:,:),Y(:,:)
    complex(8), intent(out) :: mat(:,:)
    integer :: ierr,ms,bsX,bsY
    complex(8) :: zero,one
    complex(8),allocatable :: tmp(:,:)
    zero = (0d0,0d0); one = (1d0,0d0)
    ms = size(X,1)
    bsX = size(X,2)
    bsY = size(Y,2)
    allocate(tmp(bsX,bsY))
    call ZGEMM('C','N',bsX,bsY,ms,one,X,ms,Y,ms,zero,tmp,bsX)
    call mpi_allreduce(tmp,mat,bsX*bsY,mpi_complex16,mpi_sum,comm_grid,ierr)
    deallocate(tmp)
  END SUBROUTINE z_matmat_CN


  SUBROUTINE inv(X,mat,ipiv,work)
    complex(8), intent(in) :: X(:,:)
    complex(8), intent(out) :: mat(:,:)
    integer, intent(out) :: ipiv(:)
    complex(8), intent(out) :: work(:)
    integer :: bs,info
    bs = size(X,1)
    mat = X
    call ZGETRF(bs,bs,mat,bs,ipiv,info)
    call ZGETRI(bs,mat,bs,ipiv,work,bs,info)
  END SUBROUTINE inv


  SUBROUTINE d_LAPACK_QR(Q,R)
    real(8), intent(inout) :: Q(:,:)
    real(8), intent(out) :: R(:,:)
    integer :: i,j,ms,bs,lwork,info
    real(8),allocatable :: TAU(:),work(:)
    ms = size(Q,1)
    bs = size(Q,2)
    allocate(TAU(min(ms,bs)), work(bs))
    call DGEQRF(ms,bs,Q,ms,TAU,work,bs,info)
    do i=1,bs
       do j=1,bs
          if ( i <= j ) then
             R(i,j) = Q(i,j)
          else
             R(i,j) = 0d0
          end if
       end do
    end do
    call DORGQR(ms,bs,bs,Q,ms,TAU,work,bs,info)
    deallocate(TAU , work)
  END SUBROUTINE d_LAPACK_QR


  SUBROUTINE z_LAPACK_QR(Q,R)
    complex(8), intent(inout) :: Q(:,:)
    complex(8), intent(out) :: R(:,:)
    integer :: i,j,ms,bs,lwork,info
    complex(8),allocatable :: TAU(:),work(:)
    ms = size(Q,1)
    bs = size(Q,2)
    allocate(TAU(min(ms,bs)), work(bs))
    call ZGEQRF(ms,bs,Q,ms,TAU,work,bs,info)
    do i=1,bs
       do j=1,bs
          if ( i <= j ) then
             R(i,j) = Q(i,j)
          else
             R(i,j) = (0d0,0d0)
          end if
       end do
    end do
    call ZUNGQR(ms,bs,bs,Q,ms,TAU,work,bs,info)
    deallocate(TAU , work)
  END SUBROUTINE z_LAPACK_QR


  SUBROUTINE z_CGS_QR_naive(Q,R)
    ! QR factorization via classic Gram-Schmidt method
    complex(8),intent(out) :: Q(:,:),R(:,:)
    integer i,j
    R = (0d0,0d0)
    do i=1,size(Q,2)
       do j=1,i-1
          R(j,i) = z_dot(Q(:,j),Q(:,i))
          Q(:,i) = Q(:,i) - R(j,i)*Q(:,j)
       end do
       R(i,i) = dsqrt(dble(z_dot(Q(:,i),Q(:,i))))
       Q(:,i) = Q(:,i) / R(i,i)
    end do
  END SUBROUTINE z_CGS_QR_naive


  SUBROUTINE d_MGS_QR(Q,R)
    ! QR factorization via modified Gram-Schmidt method
    real(8),intent(out) :: Q(:,:),R(:,:)
    integer i,j,k,s1,s2
    s1 = size(Q,1)
    s2 = size(Q,2)
    R = 0d0
    do i=1,s2
       R(i,i) = dsqrt(dble(d_dot(Q(:,i),Q(:,i))))
       do k=1,s1
          Q(k,i) = Q(k,i) / R(i,i)
       end do
       do j=i+1,s2
          R(i,j) = d_dot(Q(:,i),Q(:,j))
          do k=1,s1
             Q(k,j) = Q(k,j) - R(i,j)*Q(k,i)
          end do
       end do
    end do
  END SUBROUTINE d_MGS_QR


  SUBROUTINE z_MGS_QR(Q,R)
    ! QR factorization via modified Gram-Schmidt method
    complex(8),intent(out) :: Q(:,:),R(:,:)
    integer i,j,k,s1,s2
    s1 = size(Q,1)
    s2 = size(Q,2)
    R(:,:) = (0d0,0d0)
    do i=1,s2
       R(i,i) = dsqrt(dble(z_dot(Q(:,i),Q(:,i))))
       do k=1,s1
          Q(k,i) = Q(k,i) / R(i,i)
       end do
       do j=i+1,s2
          R(i,j) = z_dot(Q(:,i),Q(:,j))
          do k=1,s1
             Q(k,j) = Q(k,j) - R(i,j)*Q(k,i)
          end do
       end do
    end do
  END SUBROUTINE z_MGS_QR

#ifdef TEST
  subroutine tic(arg,comm)
    type(timer), intent(out) :: arg
    integer, intent(in) :: comm
    integer :: ierr
    call mpi_barrier(comm,ierr)
    call cpu_time(arg%tmpct)
    arg%tmpet = mpi_wtime()
  end subroutine tic

  subroutine toc(arg,comm)
    type(timer), intent(inout) :: arg
    integer, intent(in) :: comm
    real(8) :: ct, et
    integer :: ierr
    call mpi_barrier(comm,ierr)
    call cpu_time(ct)
    et = mpi_wtime()
    arg%ct = arg%ct + ct - arg%tmpct
    arg%et = arg%et + et - arg%tmpet
  end subroutine toc

  subroutine clear_timer(arg)
    type(timer), intent(out) :: arg
    arg%ct = 0d0
    arg%et = 0d0
    arg%tmpct = 0d0
    arg%tmpet = 0d0
  end subroutine clear_timer
#endif

  SUBROUTINE z_CGS_QR_Takahashi(Q,R)
    implicit none
    complex(8),intent(inout) :: Q(:,:)
    complex(8),intent(out) :: R(:,:)
    integer :: blk_size,min_blk,i,j
    blk_size = size(Q,2)
    if ( blk_size == 1 ) then
       call z_MGS_QR(Q,R)
       return
    end if
    min_blk = 2
    do j=1,blk_size
       do i=j+1,blk_size
          R(i,j) = (0d0,0d0)
       end do
    end do
    call z_CGS_QR_Takahashi_sub(Q,R,1,blk_size,1,blk_size,blk_size,min_blk)
  END SUBROUTINE z_CGS_QR_Takahashi

  RECURSIVE SUBROUTINE z_CGS_QR_Takahashi_sub(Q,R,mm1,mm2,nn1,nn2,MBLK,min_blk)
    implicit none
    complex(8),intent(inout) :: Q(:,:)
    complex(8),intent(out) :: R(:,:)
    integer,intent(in) :: mm1,mm2,nn1,nn2,MBLK,min_blk
    integer :: n,ns,ne,m,ms,me,m1,m2,mm,nn,MBLKH,ierr
    integer :: n1,n2,ML0,i,s1
    real(8) :: c,d
    complex(8) :: one,zero
    complex(8),allocatable :: utmp(:),vtmp(:),utmp2(:,:),vtmp2(:,:)
    one = (1d0,0d0)
    zero = (0d0,0d0)
    n1  = id_grid(myrank_g)+1
    n2  = id_grid(myrank_g)+ir_grid(myrank_g)
    ML0 = n2-n1+1
    s1 = size(Q,1)
    do ms=mm1,mm2,MBLK
       me=min(ms+MBLK-1,mm2)
       mm=me-ms+1
       do ns=nn1,nn2,MBLK
          ne=min(ns+MBLK-1,nn2)
          ne=min(ne,me-1)
          nn=ne-ns+1
          if ( nn<=0 ) cycle
          if ( ms>=ne+1 ) then
             allocate( utmp2(ns:ne,ms:me),vtmp2(ns:ne,ms:me) )
             call zgemm('C','N',nn,mm,ML0,one,Q(:,ns),ML0,Q(:,ms),ML0,zero,utmp2,nn)
             call mpi_allreduce(utmp2,vtmp2,nn*mm,mpi_complex16,mpi_sum,comm_grid,ierr)
             R(ns:ne,ms:me) = vtmp2
             call zgemm('N','N',ML0,mm,nn,-one,Q(:,ns),ML0,vtmp2,nn,one,Q(:,ms),ML0)
             deallocate( vtmp2,utmp2 )
             if ( ms==ne+1 ) then
                d=0.d0
                do i=1,s1
                   d=d+abs(Q(i,ms))**2
                end do
                call mpi_allreduce(d,c,1,mpi_real8,mpi_sum,comm_grid,ierr)
                R(ms,ms) = sqrt(c)
                c=1.d0/sqrt(c)
                do i=1,s1
                   Q(i,ms)=c*Q(i,ms)
                end do
             end if
          else if ( mm<=min_blk ) then
             allocate( utmp(min_blk),vtmp(min_blk) )
             do m=ms,me
                n=min(m-1,ne)
                if ( n-ns+1>0 ) then
                   call zgemv('C',ML0,n-ns+1,one,Q(:,ns),ML0,Q(:,m),1,zero,utmp,1)
                   call mpi_allreduce(utmp,vtmp,n-ns+1,mpi_complex16,mpi_sum,comm_grid,ierr)
                   R(ns:n,m) = vtmp(1:n-ns+1) 
                   call zgemv('N',ML0,n-ns+1,-one,Q(:,ns),ML0,vtmp,1,one,Q(:,m),1)
                end if
                if ( m==1 .or. (n==m-1 .and. m/=ns) ) then
                   d=0.d0
                   do i=1,s1
                      d=d+abs(Q(i,m))**2
                   end do
                   call mpi_allreduce(d,c,1,mpi_real8,mpi_sum,comm_grid,ierr)
                   R(m,m) = sqrt(c)
                   c=1.d0/sqrt(c)
                   do i=1,s1
                      Q(i,m)=c*Q(i,m)
                   end do
                end if
             end do
             deallocate( vtmp,utmp )
          else
             MBLKH=max(MBLK/2,min_blk)
             call z_CGS_QR_Takahashi_sub(Q,R,ms,me,ns,ne,MBLKH,min_blk)
          end if
       end do
    end do
  END SUBROUTINE z_CGS_QR_Takahashi_sub


  SUBROUTINE z_TS_QR(Q,R)
    complex(8),intent(inout) :: Q(:,:)
    complex(8),intent(inout) :: R(:,:)
    integer :: i,m
    integer,allocatable :: id(:)
    complex(8),allocatable :: Q_mm(:,:)
    if ( nprocs_g > 1 ) then
       allocate(id(nprocs_g))
       do i = 1, nprocs_g
          id(i) = i-1
       end do
       m = size(Q,2)
       allocate(Q_mm(m,m))
       call z_TS_QR_sub(Q,Q_mm,R,id,0)
    else
       call z_LAPACK_QR(Q,R)
    end if
  END SUBROUTINE z_TS_QR


  RECURSIVE SUBROUTINE z_TS_QR_sub(Q,Q_half,R,id,depth)
    complex(8),intent(inout) :: Q(:,:),Q_half(:,:)
    complex(8),intent(inout) :: R(:,:)
    integer,intent(inout) :: id(:)
    integer,intent(in) :: depth
    integer :: i, m, irank, my_group_id, partner_rank, ierr
    logical :: is_factor
    integer :: istatus(mpi_status_size)
    complex(8),allocatable :: Q_tmp(:,:), Q_2mm(:,:)

    m = size(Q,2)

    if ( depth /= ceiling(log(dble(nprocs_g))/log(2d0)) ) then
       is_factor = mod(id(myrank_g+1),2) == 0
       my_group_id = int(id(myrank_g+1)/2)
       if ( is_factor ) then
          partner_rank = -1
          do irank = 1, nprocs_g, 2**depth
             if ( int(id(irank)/2) == my_group_id .and. mod(id(irank),2) == 1 ) then
                partner_rank = irank - 1
             end if
          end do
          if ( partner_rank == -1 ) then
             do i = 1, nprocs_g, 2**(depth+1)
                id(i) = int((i-1)/(2**(depth+1)))
             end do
             call z_TS_QR_sub(Q,Q_half,R,id,depth+1)
             return
          end if

          call z_LAPACK_QR(Q,R) 

          call mpi_recv(Q_half, m*m, mpi_complex16, partner_rank &
               , mpi_any_tag, comm_grid, istatus, ierr)
          allocate(Q_2mm(2*m,m))
          Q_2mm(1:m,:) = R
          Q_2mm(m+1:2*m,:) = Q_half
          do i = 1, nprocs_g, 2**(depth+1)
             id(i) = int((i-1)/(2**(depth+1)))
          end do
          call z_TS_QR_sub(Q_2mm,Q_half,R,id,depth+1)
          Q_half = Q_2mm(1:m,:)
          Q_2mm(1:m,:) = R
          call mpi_send(Q_2mm, 2*m*m, mpi_complex16, partner_rank, 1, comm_grid, ierr)
       else
          call z_LAPACK_QR(Q,R) 
          allocate(Q_2mm(2*m,m))
          do irank = 1, nprocs_g, 2**depth
             if ( int(id(irank)/2) == my_group_id .and. mod(id(irank),2) == 0 ) then
                partner_rank = irank - 1
             end if
          end do
          call mpi_send(R, m*m, mpi_complex16, partner_rank, 1, comm_grid, ierr)
          call mpi_recv(Q_2mm, 2*m*m, mpi_complex16, partner_rank &
               , mpi_any_tag, comm_grid, istatus, ierr)
          R = Q_2mm(1:m,:)
          Q_half = Q_2mm(m+1:2*m,:)
       end if
       allocate(Q_tmp(size(Q,1),size(Q,2)))
       Q_tmp = Q
       call z_matmat_NX('N',(0d0,0d0),Q_tmp,Q_half,Q)
    else
       call z_LAPACK_QR(Q,R) 
    end if

  END SUBROUTINE z_TS_QR_sub


  SUBROUTINE disp_lin_result(sinfo,info_size,stimer)
    type(solver_info),intent(in) :: sinfo(:)
    integer,intent(in) :: info_size
    type(solver_timers),intent(in) :: stimer
    integer :: m
      
    write(*,'(a5,1x,a10,1x,a10,1x,a6,1x,a15,1x,a15)') &
         'id','real','imag','iter','rec res','true res'
    do m=1,info_size
       write(*,'(I5,1x,1pe10.2,1x,1pe10.2,1x,I6,1x,1pe15.7,1x,1pe15.7)') &
            m, sinfo(m)%shift_r, sinfo(m)%shift_i, sinfo(m)%iter &
            , sinfo(m)%rec_res, sinfo(m)%true_res
    end do
    write(*,'(3x,a)') ' --- LINEAR SOLVER TIME --- '
    write(*,'(1x,a,1x,2f10.3)') 'total           ',stimer%total_time%ct,stimer%total_time%et
    write(*,'(1x,a,1x,2f10.3)') '  |--hpsi       ',stimer%hpsi_time%ct,stimer%hpsi_time%et
    write(*,'(1x,a,1x,2f10.3)') '  |--orth       ',stimer%lin_orth_time%ct,stimer%lin_orth_time%et
    write(*,'(1x,a,1x,2f10.3)') '  |--shift      ',stimer%shift_time%ct,stimer%shift_time%et
    write(*,'(1x,a,1x,2f10.3)') '       |--scalar',stimer%scalar_time%ct,stimer%scalar_time%et
    write(*,'(1x,a,1x,2f10.3)') '       |--input ',stimer%input_time%ct,stimer%input_time%et
    write(*,'(1x,a,1x,2f10.3)') '       |--gemm  ',stimer%gemm_time%ct,stimer%gemm_time%et
    write(*,'(1x,a,1x,2f10.3)') '       |--out   ',stimer%out_time%ct,stimer%out_time%et
    call flush(6)
  END SUBROUTINE disp_lin_result
    
END MODULE iter_lin_solvers
