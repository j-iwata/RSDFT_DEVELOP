MODULE apply_sseig_module

  use global_variables, only: disp_switch
  use parallel_module
  use bc_module, only: allocate_bcset
  use ps_nloc2_module, only: allocate_ps_nloc2
  use sseig, only: n_ccurve,ssgamma,ssrho,rho_scale,opt,inputV,L_fixed,n_ccurve_local &
                  ,ssres,eigvec,eigval,n_rhs_set,n_rhs_set_local,my_rhs_set_id &
                  ,solver_timers_array,solver_info_array,sseig_timers_array,first_time &
                  ,my_ccurve_id,numeig,num_basis,current_sseig_timers,ccurve_id_rank &
                  ,current_ccurve_local_id,eigval_merged,eigvec_merged,old_idx &
                  ,eigvec_merged_old,unit_ss_band_cor,unit_ss_band_wf,cmyrank,file_ss_wf &
                  ,unit_ss_band_val,file_ss_wf_split,SS_IO,z_sseig_one_circle,opt &
                  ,disp_sseig_result, z_rayleigh_ritz,z_normalize_eigenvectors,z_eval_residual &
                  ,opt_N,out_iter_max,out_iter_init
  use iter_lin_solvers, only: comm_rhs,comm_ccurve,nprocs_c,nprocs_r,myrank_c,myrank_r &
                             ,iterative_method_args,current_solver_info &
                             ,solver_timers,solver_info,solver &
                             ,disp_lin_result,z_mgs_qr,z_matmat_cn,tol_iter
  use timer_module

  implicit none

  PRIVATE
  PUBLIC :: apply_sseig,total_numeig

  integer :: total_numeig

CONTAINS


  SUBROUTINE apply_sseig(dk_bz)
    real(8),intent(IN) :: dk_bz(3)
    integer :: ev_size,sscount
    integer :: i,j,k,m,n1,n2,ierr,count,i_cc,j_cc
    integer :: i_rhs,j_rhs,tmp_disp
    real(8) :: rannum,randd(2),eigrestol,left,right
    integer :: MB_d_tmp
    real(8),allocatable :: tmp_mat_r(:,:), res_merged(:) &
         , eigval_merged_tmp(:), res_merged_tmp(:)
    complex(8),allocatable :: tmp_mat_c(:,:),dummy(:,:),eigvec_merged_tmp(:,:)
    integer :: istatus(mpi_status_size), out_iter
    type(solver_info),allocatable :: sinfo(:)
    type(solver_timers) :: stimer
    type(timer) :: apply_sseig_time, merge_time, write_wf_time
    type(timer) :: tot_time
    integer :: a1,a2,a3,b1,b2,b3,l,n
    integer :: k1,num_candidate,ntmp
    real(8),parameter :: restol_ji=1.d-2
    complex(8) :: ovlp1,ovlp2
    real(8) :: norm1,norm2,c 

    if ( disp_switch ) then
       write(*,'(a60," apply_sseig")') repeat("-",60)
       call flush(6)
    end if

    eigrestol = 1d-2

    n1 = id_grid(myrank_g)+1
    n2 = id_grid(myrank_g)+ir_grid(myrank_g)

    if ( myrank_c == 0 .and. myrank_r == 0 .and. myrank_g == 0 ) then
       write(*,*) 'closed curve id, interval'
       do i=1,n_ccurve
          write(*,*) i,'[',ssgamma(i)-ssrho(i)/rho_scale &
               ,',',ssgamma(i)+ssrho(i)/rho_scale,']'
       end do
       call flush(6)
    end if

    opt%L   = L_fixed
    opt%N   = opt_N
!    opt%N   = 32
!    opt%N   = 64
!    opt%N   = 128
    opt%M   = opt%N / 2

    ev_size = opt%M * opt%L

    if ( disp_switch ) then
       write(*,*) "ev_size,myrank=",ev_size,myrank
       call flush(6)
    end if

    MB_d_tmp = MB_d
    MB_d     = opt%L
    call allocate_bcset
    call allocate_ps_nloc2(MB_d)

    if ( .not.allocated(inputV) ) then
       allocate( inputV(n1:n2,opt%L,n_ccurve_local(myrank_c)) )
       inputV=0.d0
    end if
    if ( .not.allocated(eigvec)  ) then
       allocate( eigvec(n1:n2,ev_size,n_ccurve_local(myrank_c)) )
       eigvec=0.d0
    end if
    if ( .not.allocated(eigval)  ) then
       allocate( eigval(ev_size,n_ccurve_local(myrank_c)) )
       eigval=0.d0
    end if
    if ( .not.allocated(ssres)  ) then
       allocate( ssres(ev_size,n_ccurve_local(myrank_c)) )
       ssres=0.d0
    end if

    do i=0,nprocs_r-1
       n=(n_rhs_set-1)/nprocs_r+1
       n_rhs_set_local(i) = min( n*(i+1), opt%L ) - n*i
    end do

    if ( allocated(my_rhs_set_id) ) deallocate(my_rhs_set_id)
    allocate( my_rhs_set_id(n_rhs_set_local(myrank_r)) )

    count = 0
    do i=1,n_rhs_set
       if ( (i-1)/((n_rhs_set-1)/nprocs_r+1) == myrank_r ) then
          count = count + 1
          my_rhs_set_id(count) = i
       end if
    end do

    if ( allocated(current_solver_info) ) deallocate(current_solver_info)
    allocate( current_solver_info(opt%N) )

    if ( allocated(sinfo) ) deallocate(sinfo)
    allocate( sinfo(opt%N) )

    m = n_rhs_set_local(myrank_r)
    n = n_ccurve_local(myrank_c)

    if ( allocated(solver_timers_array) ) deallocate(solver_timers_array)
    allocate( solver_timers_array(m,n) )

    if ( allocated(solver_info_array) ) deallocate(solver_info_array)
    allocate( solver_info_array(opt%N,m,n) )

    if ( allocated(sseig_timers_array) ) deallocate(sseig_timers_array)
    allocate( sseig_timers_array(n) )


    num_candidate = 0

    call clear_timer(tot_time)

    do out_iter=1,out_iter_max

       if ( DISP_SWITCH ) then
          write(*,'(a40," out_iter=",i2)') repeat("-",40),out_iter
          call flush(6)
       end if

       call tic(tot_time,mpi_comm_world)

       do i=1,n_ccurve

          call mpi_barrier(comm_ccurve,ierr)

          do j=1,n_ccurve_local(myrank_c)

             if ( my_ccurve_id(j) == i ) then

                if ( out_iter <= out_iter_init ) then

                   if ( myrank_r == 0 .and. myrank_g == 0 ) then
                      write(*,*) i,'input vectors are random vectors'
                      call flush(6)
                   end if
                   iterative_method_args%epsmax = tol_iter(1)
                   iterative_method_args%epsmin = tol_iter(2)
                   do k = 1,opt%L
                      do m = n1,n2
                         call random_number(randd)
                         randd = 2d0*randd - 1d0
                         inputV(m,k,j) = dcmplx( randd(1),randd(2) )
                      end do
!                      do k1=1,k-1
!                         ovlp1=sum( conjg(inputV(:,k1,j))*inputV(:,k,j) )
!                         call mpi_allreduce(ovlp1,ovlp2,1,mpi_complex16,mpi_sum,comm_grid,ierr)
!                         inputV(:,k,j)=inputV(:,k,j)-inputV(:,k1,j)*ovlp2
!                      end do
!                      norm1=sum(abs(inputV(:,k,j))**2)
!                      call mpi_allreduce(norm1,norm2,1,mpi_real8,mpi_sum,comm_grid,ierr)
!                      c=1.d0/sqrt(norm2)
!                      inputV(:,k,j)=c*inputV(:,k,j)
                   end do ! k

                else ! out_iter > out_iter_init

                   if ( myrank_r == 0 .and. myrank_g == 0 ) then
                      write(*,*) i,'input vectors are approximate eigenvectors'
                      call flush(6)
                   end if
                   iterative_method_args%epsmax = tol_iter(3)
                   iterative_method_args%epsmin = tol_iter(4)
                   inputV(:,:,j)=0.d0
!tmp
                   do k=1,opt%L
                      do m=1,num_basis(i)
                         call random_number(randd)
                         inputV(:,k,j)=inputV(:,k,j)+dcmplx(randd(1),randd(2))*eigvec(:,m,j)
                      end do
                   end do
                   cycle
                   if ( num_candidate == 0 ) then
                      do k = 1,opt%L
                         do m = n1,n2
                            call random_number(randd)
                            randd = 2d0*randd - 1d0
                            inputV(m,k,j) = dcmplx( randd(1),randd(2) )
                         end do
!                         do m=1,num_basis(i)
!                            if ( numeig(i)/=0 .and. ssres(m,i) >= eigrestol ) cycle
!                            call random_number(rannum)
!                            inputV(:,k,j)=inputV(:,k,j)+rannum*eigvec(:,m,j)
!                         end do
                      end do ! k
                   else
                      do k = 1,opt%L
!                         do m = n1,n2
!                            call random_number(randd)
!                            randd = 2d0*randd - 1d0
!                            inputV(m,k,j) = dcmplx( randd(1),randd(2) )
!                         end do
                         do m=1,num_candidate
                            call random_number(randd)
                            randd=2.d0*randd-1.d0
                            inputV(:,k,j)=inputV(:,k,j)+dcmplx(randd(1),randd(2))*eigvec_merged(:,m)
                         end do
                      end do ! k
!                      do k = 1,opt%L
!                         do m = n1,n2
!                            call random_number(rannum)
!                            rannum = 2d0*rannum - 1d0
!                            inputV(:,k,j) = rannum
!                         end do
!                      end do
!                      call ortho_sub(n2-n1+1,opt%L,inputV(n1,1,j),num_candidate,eigvec_merged)
                   end if

                end if

             end if

          end do ! j [1:n_ccurve_local(myrank_c)]

       end do ! i [1:n_ccurve]

       eigval(:,:)=0.d0
       ssres(:,:)=0.d0

       call mpi_barrier(comm_ccurve,ierr)
     

       if ( myrank_c == 0 .and. myrank_r == 0 .and. myrank_g == 0 ) then
          write(*,*) 'START SSEIG'
          write(*,*) 'opt%N, opt%M, opt%L',opt%N,opt%M,opt%L
          call flush(6)
       end if
     

       do i=1,n_ccurve_local(myrank_c)

          current_ccurve_local_id = i
          num_basis(i) = ev_size
          n = my_ccurve_id(i)

          call clear_timer(current_sseig_timers%total_time)
          call clear_timer(current_sseig_timers%ls_time)
          call clear_timer(current_sseig_timers%post_time)
          call clear_timer(current_sseig_timers%svd_time)
          call clear_timer(current_sseig_timers%orth_time)
          call clear_timer(current_sseig_timers%rr_time)
          call clear_timer(current_sseig_timers%mate_time)
          call clear_timer(current_sseig_timers%dsyev_time)
          call clear_timer(current_sseig_timers%rot_time)

          if ( disp_switch ) then
             write(*,'(a40," z_sseig_one_circle",i4)') repeat("-",40),i
             call flush(6)
          end if

          call tic(current_sseig_timers%total_time,comm_grid)
          call z_sseig_one_circle( ssgamma(n),ssrho(n),opt,eigval(:,i) &
               ,eigvec(:,:,i),ssres(:,i),solver,inputV(:,:,i),num_basis(i) )
          call toc(current_sseig_timers%total_time,comm_grid)
          if ( nprocs_r == 1 .and. nprocs_c == 1 .and. myrank_g == 0 ) then
             if ( disp_switch ) then
                write(*,'(a40," disp_sseig_result",i4)') repeat("-",40),i
                call flush(6)
             end if
             call disp_sseig_result( current_sseig_timers,ssgamma(n) &
             ,ssrho(n),num_basis(i), eigval(:,i), ssres(:,i) )
          end if

          sseig_timers_array(i) = current_sseig_timers

          numeig(n) = 0
          do m=1,num_basis(i)
             if ( ssres(m,i) < eigrestol ) numeig(n) = numeig(n) + 1
          end do

       end do ! i


       do i=1,n_ccurve
          if ( myrank_c == 0 ) then
             if ( ccurve_id_rank(i) /= myrank_c ) then
                call mpi_recv(numeig(i),1,mpi_integer,ccurve_id_rank(i) &
                     ,mpi_any_tag,comm_ccurve,istatus,ierr)
             end if
          else
             if ( ccurve_id_rank(i) == myrank_c ) then
                call mpi_send(numeig(i),1,mpi_integer,0,0,comm_ccurve,ierr)
             end if
          end if
       end do
       call mpi_bcast(numeig,n_ccurve,mpi_integer,0,comm_ccurve,ierr)

     
       if ( nprocs_r > 1 .or. nprocs_c > 1 ) then
          call mpi_barrier(comm_ccurve,ierr)
          do i_cc = 1, n_ccurve
          do j_cc = 1, n_ccurve_local(myrank_c)
             if ( my_ccurve_id(j_cc) == i_cc ) then
                call mpi_barrier(comm_rhs,ierr)
                do i_rhs = 1, n_rhs_set
                do j_rhs = 1, n_rhs_set_local(myrank_r)
                   if ( my_rhs_set_id(j_rhs) == i_rhs .and. myrank_g == 0 ) then
                      write(*,*) 'Result of shifted block CG rQ method'
                      write(*,'(1x,"closed curve id = ",I3,3x,"rhs set id = ",I3)') i_cc,i_rhs
                      write(*,'(1x,"myrank_c = ",I3,3x,"myrank_r = ",I3)') myrank_c,myrank_r
                      call disp_lin_result(solver_info_array(:,j_rhs,j_cc) &
                           , opt%N, solver_timers_array(j_rhs,j_cc))
                      call flush(6)
                   end if
                end do
                call mpi_barrier(comm_rhs,ierr)
                end do
                if ( myrank_r == 0 .and. myrank_g == 0 ) then
                   write(*,'(1x,a,I3)') 'closed curve id = ',i_cc
                   write(*,'(1x,a,I3)') 'myrank_c = ',myrank_c
                   call disp_sseig_result(sseig_timers_array(j_cc) &
                        ,ssgamma(j_cc), ssrho(j_cc), num_basis(j_cc) &
                        ,eigval(:,j_cc), ssres(:,j_cc))
                   call flush(6)
                end if
             end if
          end do
          call mpi_barrier(comm_ccurve,ierr)
          end do
       end if ! ( nprocs_r > 1 .or. nprocs_c > 1 )

       if ( myrank_c == 0 .and. myrank_r == 0 .and. myrank_g == 0 ) then
          do i=1,n_ccurve
             write(*,*) 'numeig', i, numeig(i)
          end do
          call flush(6)
       end if

       ntmp = num_candidate + sum( numeig(1:n_ccurve) )

       if ( all(numeig(1:n_ccurve)==0) ) then
          cycle
       end if

       allocate( eigvec_merged_tmp(n1:n2,ntmp) )

       if ( allocated(eigvec_merged) ) then
          eigvec_merged_tmp(:,1:num_candidate) = eigvec_merged(:,:)
          deallocate( eigvec_merged )
          deallocate( eigval_merged )
       end if
       n=num_candidate
       do i=1,n_ccurve
          do m=1,num_basis(i)
             if ( ssres(m,i) < eigrestol ) then
                n=n+1
                eigvec_merged_tmp(:,n) = eigvec(:,m,i)
             end if
          end do
       end do
       if ( n /= ntmp ) stop "stop@apply_sseig(1)"

       allocate( dummy(ntmp,ntmp) )

       call z_MGS_QR(eigvec_merged_tmp,dummy)
       call z_MGS_QR(eigvec_merged_tmp,dummy)

       deallocate( dummy )
       allocate( eigval_merged_tmp(ntmp) )

       call z_rayleigh_ritz(eigvec_merged_tmp,0d0,eigval_merged_tmp,n1,n2)
       call z_normalize_eigenvectors(eigvec_merged_tmp)

       allocate( res_merged_tmp(ntmp) )

       call z_eval_residual(eigval_merged_tmp,eigvec_merged_tmp,res_merged_tmp,n1,n2)

       left=ssgamma(1)-ssrho(1)
       right=ssgamma(1)+ssrho(1)
       m=0
       n=0
       do i=1,ntmp
          if ( res_merged_tmp(i) < restol_ji ) then
             n=n+1
             if ( left <= eigval_merged_tmp(i) .and. eigval_merged_tmp(i) <= right ) m=m+1
          end if
       end do

       call toc(tot_time,mpi_comm_world)

       if ( disp_switch ) then
          write(*,'(1x,"m,n,ntmp=",3i6,2f15.10,4f10.3)') m,n,ntmp,left,right &
          ,current_sseig_timers%total_time%ct,current_sseig_timers%total_time%et &
          ,tot_time%ct,tot_time%et
          call flush(6)
       end if

       allocate( eigval_merged(n) )
       allocate( eigvec_merged(n1:n2,n) )
       num_candidate = n
       n=0
       do i=1,ntmp
          if ( res_merged_tmp(i) < restol_ji ) then
             n=n+1
             eigval_merged(n) = eigval_merged_tmp(i)
             eigvec_merged(:,n) = eigvec_merged_tmp(:,i)
             if ( disp_switch ) then
                write(*,*) i,eigval_merged_tmp(i),res_merged_tmp(i), "   ji"
                call flush(6)
             end if
          end if
       end do
       if ( n /= num_candidate ) stop "stop@apply_sseig(2)"

       deallocate( res_merged_tmp )
       deallocate( eigval_merged_tmp )
       deallocate( eigvec_merged_tmp )
10     continue

    end do ! out_iter


    stop

 !! merge eigenvectors from all closed curves

    if ( disp_switch ) write(*,'(a40," merge eigenvectors from all closed curves")')

    total_numeig = sum(numeig)
    if ( disp_switch ) write(*,*) "total_numeig=",total_numeig
    allocate( eigval_merged_tmp(total_numeig) )
    allocate( eigvec_merged_tmp(n1:n2,total_numeig) )
 
    tmp_disp = 1
    do i=1,n_ccurve
       do j=1,n_ccurve_local(myrank_c)
          if ( my_ccurve_id(j) == i ) then
             count = 0
             do m=1,num_basis(j)
                if ( ssres(m,j) < eigrestol ) then
                   eigval_merged_tmp(tmp_disp+count) = eigval(m,j)
                   eigvec_merged_tmp(:,tmp_disp+count) = eigvec(:,m,j)
                   count = count + 1
                end if
             end do
             if ( myrank_c /= 0 ) then
                call mpi_send(eigval_merged_tmp(tmp_disp),numeig(i) &
                     ,mpi_real8,0,0,comm_ccurve,ierr)
                call mpi_send(eigvec_merged_tmp(:,tmp_disp) &
                     ,numeig(i)*(n2-n1+1),mpi_complex16,0,0,comm_ccurve,ierr)
             end if
          end if
       end do
       if ( myrank_c == 0 .and. ccurve_id_rank(i) /= 0 ) then
          call mpi_recv(eigval_merged_tmp(tmp_disp),numeig(i) &
               ,mpi_real8,ccurve_id_rank(i),mpi_any_tag,comm_ccurve,istatus,ierr)
          call mpi_recv(eigvec_merged_tmp(:,tmp_disp),numeig(i)*(n2-n1+1) &
               ,mpi_complex16,ccurve_id_rank(i),mpi_any_tag,comm_ccurve,istatus,ierr)
       end if
       tmp_disp = tmp_disp + numeig(i)
    end do

    call mpi_bcast(eigval_merged_tmp,total_numeig,mpi_real8,0,comm_ccurve,ierr)
    call mpi_bcast(eigvec_merged_tmp,total_numeig*(n2-n1+1) &
                  ,mpi_complex16,0,comm_ccurve,ierr)

    if ( allocated(eigval_merged) ) deallocate(eigval_merged)
    if ( allocated(eigvec_merged) ) deallocate(eigvec_merged)

    if ( n_ccurve > 1 ) then
       allocate( res_merged_tmp(total_numeig) )
       allocate( dummy(total_numeig,total_numeig) )
       call z_MGS_QR(eigvec_merged_tmp,dummy)
       call z_MGS_QR(eigvec_merged_tmp,dummy)
       call z_rayleigh_ritz(eigvec_merged_tmp,0d0,eigval_merged_tmp,n1,n2)
       call z_normalize_eigenvectors(eigvec_merged_tmp)
       call z_eval_residual(eigval_merged_tmp,eigvec_merged_tmp,res_merged_tmp,n1,n2)
       count = 0
       do m=1,total_numeig
          if ( res_merged_tmp(m) < eigrestol ) then
             count = count + 1
          end if
       end do
       allocate( eigval_merged(count) )
       allocate( eigvec_merged(n1:n2,count) )
       allocate( res_merged(count) )
       count = 0
       do m=1,total_numeig
          if ( res_merged_tmp(m) < eigrestol ) then
             count = count + 1
             eigval_merged(count) = eigval_merged_tmp(m)
             eigvec_merged(:,count) = eigvec_merged_tmp(:,m)
             res_merged(count) = res_merged_tmp(m)
          end if
       end do
       total_numeig = count
       if ( disp_switch ) write(*,*) "total_numeig(=count)=",count
       deallocate( dummy )
       deallocate( res_merged_tmp )
    else ! n_ccurve==1
       allocate( eigval_merged(total_numeig) )
       allocate( eigvec_merged(n1:n2,total_numeig) )
       allocate( res_merged(total_numeig) )
       eigval_merged(:)   = eigval_merged_tmp(:)
       eigvec_merged(:,:) = eigvec_merged_tmp(:,:)
       res_merged(:) = 0.d0
       count = 1
       do m=1,num_basis(1)
          if ( ssres(m,1) < eigrestol ) then
             res_merged(count) = ssres(m,1)
             count = count + 1
          end if
       end do
       if ( disp_switch ) write(*,*) "numeig,count=",numeig(1),count
    end if

    deallocate( eigvec_merged_tmp )
    deallocate( eigval_merged_tmp )

!-- band connectivity check

    if ( disp_switch ) write(*,'(1x,a40," band connetivity")')

    if ( allocated(old_idx) ) deallocate(old_idx)
    allocate( old_idx(total_numeig) )
    old_idx(:) = 0

    if ( first_time ) then
       old_idx(:) = 0
       first_time = .false.
    else
       m=size(eigvec_merged_old,2)
       n=size(eigvec_merged,2)
       allocate( tmp_mat_r(m,n), tmp_mat_c(m,n) )
       call z_matmat_CN(eigvec_merged_old,eigvec_merged,tmp_mat_c)
       tmp_mat_r(:,:) = abs( tmp_mat_c(:,:) )
       do i=1,m
          old_idx( maxloc(tmp_mat_r(i,:)) ) = i
       end do
       if ( myrank_g == 0 .and. myrank_c == 0 .and. myrank_r == 0 ) then
          write(unit_ss_band_cor,'(I7,1x)',advance='no') (j,j=0,size(tmp_mat_r,2))
          write(unit_ss_band_cor,*)
          do i=1,size(tmp_mat_r,1)
             write(unit_ss_band_cor,'(I8,1x)',advance='no') i
             write(unit_ss_band_cor,'(1pe8.1,1x)',advance='no') &
                  (tmp_mat_r(i,j),j=1,size(tmp_mat_r,2))
             write(unit_ss_band_cor,*)
          end do
          call flush(unit_ss_band_cor)
       end if
       deallocate( tmp_mat_c, tmp_mat_r )
    end if

    if ( myrank_c == 0 .and. myrank_r == 0 .and. myrank_g == 0 ) then
       write(*,*) 'merged_value'
       write(*,'(1x,a5,2x,a25,3x,a14,3x,a9)') 'index','eigval','res','old_index'
       write(unit_ss_band_val,'(1x,i6,3f20.12)') total_numeig,dk_bz(1:3)
       count = 0
       do m=1,total_numeig
          write(*,'(1x,I4,3x,f25.15,3x,1pe14.7,3x,I4)') &
               m,eigval_merged(m),res_merged(m),old_idx(m)
          write(unit_ss_band_val,'(1x,I4,3x,f25.15,3x,1pe14.7,3x,I4)') &
               m,eigval_merged(m),res_merged(m),old_idx(m)
       end do
       call flush(unit_ss_band_val)
    end if


  ! SS_IO 0:no read no write, 1:write only, 2:read and write
    if ( SS_IO == 1 .or. SS_IO == 2 ) then
       if ( myrank == 0 ) then
          write(*,*) 'write wave functions'
       end if
       write(cmyrank,'(i5.5)') myrank
       file_ss_wf_split = trim(file_ss_wf)//"."//trim(adjustl(cmyrank))
       open(unit_ss_band_wf,file=file_ss_wf_split,form="unformatted",status="replace")
       write(unit_ss_band_wf) n1,n2,size(eigvec,2),size(eigvec,3),size(eigvec_merged,2)
       write(unit_ss_band_wf) numeig(:)
       write(unit_ss_band_wf) num_basis(:)
       write(unit_ss_band_wf) eigvec(:,:,:)
       write(unit_ss_band_wf) ssres(:,:)
       write(unit_ss_band_wf) eigvec_merged(:,:)
       close(unit_ss_band_wf)
    end if


    if ( allocated(eigvec_merged_old) ) deallocate(eigvec_merged_old)

    allocate( eigvec_merged_old(n1:n2,size(eigvec_merged,2)) )
    eigvec_merged_old = eigvec_merged

    deallocate(solver_timers_array)
    deallocate(solver_info_array)
    deallocate(sseig_timers_array)
    deallocate(my_rhs_set_id)
    deallocate(current_solver_info)

    MB_d = MB_d_tmp
    call allocate_bcset
    call allocate_ps_nloc2(MB_d)
  
  END SUBROUTINE apply_sseig
  
  SUBROUTINE ortho_sub(mm,na,a,nb,b)
    integer,intent(IN) :: mm,na,nb
    complex(8),intent(IN) :: b(mm,nb)
    complex(8),intent(INOUT) :: a(mm,na)
    complex(8),allocatable :: ba0(:),ba(:)
    integer :: ia,ib,ierr
    allocate( ba0(nb),ba(nb) )
    do ia=1,na
       do ib=1,nb
          ba0(ib) = sum( conjg(b(:,ib))*a(:,ia) )
       end do
       call mpi_allreduce(ba0,ba,nb,MPI_COMPLEX16,MPI_SUM,comm_grid,ierr)
       do ib=1,nb
          a(:,ia) = a(:,ia) - b(:,ib)*ba(ib)
       end do
    end do
    deallocate( ba,ba0 )
  END SUBROUTINE ortho_sub


END MODULE apply_sseig_module
