MODULE apply_sseig_module

  use global_variables, only: disp_switch
  use parallel_module
  use bc_module, only: allocate_bcset
  use ps_nloc2_module, only: allocate_ps_nloc2
  use sseig, only: n_ccurve,ssgamma,ssrho,rho_scale,opt,inputV,L_fixed,n_ccurve_local &
                  ,ssres,eigvec,eigval,n_rhs_set,n_rhs_set_local,my_rhs_set_id &
                  ,solver_timers_array,solver_info_array,sseig_timers_array,first_time &
                  ,my_ccurve_id,numeig,num_basis,current_sseig_timers,ccurve_id_rank &
                  ,current_ccurve_local_id,eigvec_merged,eigval_merged &
                  ,unit_ss_band_cor,unit_ss_band_wf,cmyrank,file_ss_wf &
                  ,unit_ss_band_val,file_ss_wf_split,SS_IO,z_sseig_one_circle,opt &
                  ,disp_sseig_result, z_rayleigh_ritz,z_normalize_eigenvectors,z_eval_residual &
                  ,opt_N,out_iter_max,out_iter_init,resval_merged
  use iter_lin_solvers, only: comm_rhs,comm_ccurve,nprocs_c,nprocs_r,myrank_c,myrank_r &
                             ,iterative_method_args,current_solver_info &
                             ,solver_timers,solver_info,solver &
                             ,disp_lin_result,z_mgs_qr,z_matmat_cn,tol_iter
  use timer_module
  use watch_module

  implicit none

  PRIVATE
  PUBLIC :: apply_sseig,total_numeig

  integer :: total_numeig
  integer,allocatable :: ir_c(:),id_c(:)

CONTAINS


  SUBROUTINE apply_sseig(dk_bz)
    implicit none
    real(8),intent(IN) :: dk_bz(3)
    integer :: ev_size,sscount
    integer :: i,j,k,m,n1,n2,ierr,count,i_cc,j_cc
    integer :: i_rhs,j_rhs,tmp_disp
    real(8) :: rannum,randd(2),eigrestol,left,right
    integer :: MB_d_tmp
    real(8),allocatable :: tmp_mat_r(:,:),eigval_merged_tmp(:),resval_merged_tmp(:)
    complex(8),allocatable :: tmp_mat_c(:,:),dummy(:,:),eigvec_merged_tmp(:,:)
    integer :: istatus(mpi_status_size), out_iter
    type(solver_info),allocatable :: sinfo(:)
    type(solver_timers) :: stimer
    type(timer) :: apply_sseig_time, merge_time, write_wf_time
    type(timer) :: tot_time
    integer :: a1,a2,a3,b1,b2,b3,l,n,irank_c
    integer :: k1,ntmp,tot_num_states,tot_num_states_old
    integer,allocatable :: num_states_ccurve(:)
    real(8),parameter :: restol_ji=1.d-2
    complex(8) :: ovlp1,ovlp2
    real(8) :: norm1,norm2,c 
    logical :: flag_end

    if ( disp_switch ) then
       write(*,'(a60," apply_sseig")') repeat("-",60)
       call flush(6)
    end if

    flag_end = .false.

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
       inputV=(0.d0,0.d0)
    end if
    if ( .not.allocated(eigvec)  ) then
       allocate( eigvec(n1:n2,ev_size,n_ccurve_local(myrank_c)) )
       eigvec=(0.d0,0.d0)
    end if
    if ( .not.allocated(eigval)  ) then
       allocate( eigval(ev_size,n_ccurve_local(myrank_c)) )
       eigval=0.d0
    end if
    if ( .not.allocated(ssres)  ) then
       allocate( ssres(ev_size,n_ccurve_local(myrank_c)) )
       ssres=0.d0
    end if

!--
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
!--

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


    allocate( num_states_ccurve(n_ccurve) ) ; num_states_ccurve=0
    tot_num_states=0
    tot_num_states_old=-100000


    total_numeig = 0

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
                   end do ! k

                else ! out_iter > out_iter_init

                   if ( myrank_r == 0 .and. myrank_g == 0 ) then
                      write(*,*) i,'input vectors are approximate eigenvectors'
                      call flush(6)
                   end if
                   iterative_method_args%epsmax = tol_iter(3)
                   iterative_method_args%epsmin = tol_iter(4)
                   inputV(:,:,j)=0.d0

                   do k=1,opt%L
                      do m=1,num_basis(j)
                         call random_number(randd)
                         inputV(:,k,j)=inputV(:,k,j)+dcmplx(randd(1),randd(2))*eigvec(:,m,j)
                      end do
                   end do
                   cycle

                   if ( total_numeig == 0 ) then
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
                         do m = n1,n2
                            call random_number(randd)
                            randd = 2d0*randd - 1d0
                            inputV(m,k,j) = dcmplx( randd(1),randd(2) )
                         end do
                      end do
!                      do k = 1,opt%L
!                         do m=1,total_numeig
!                            call random_number(randd)
!                            randd=2.d0*randd-1.d0
!                            inputV(:,k,j)=inputV(:,k,j)+dcmplx(randd(1),randd(2))*eigvec_merged(:,m)
!                         end do
!                      end do ! k
                      call ortho_sub(n2-n1+1,opt%L,inputV(n1,1,j),total_numeig,eigvec_merged)
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

       call global_watch(flag_end)
       if ( flag_end ) then
          if ( myrank == 0 ) write(*,*) "etime limit exceeded!"
          call mpi_finalize(ierr)
          stop "stop@apply_sseig"
       end if

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
             write(*,*) 'i,numeig(i)', i, numeig(i)
          end do
          call flush(6)
       end if


       ntmp = total_numeig + sum( numeig(1:n_ccurve) )


       if ( all(numeig(1:n_ccurve)==0) ) then
          if ( disp_switch ) write(*,*) "No new vectors( numeig(:)==0 )"
          exit
       end if


       allocate( eigvec_merged_tmp(n1:n2,ntmp) )
       eigvec_merged_tmp(:,:)=(0.d0,0.d0)

       n=0
       do irank_c=0,myrank_c-1
          do i=1,n_ccurve
             if ( ccurve_id_rank(i) == irank_c ) then
                n=n+numeig(i)
             end if
          end do
       end do

       do i=1,n_ccurve_local(myrank_c)
          do m=1,num_basis(i)
             if ( ssres(m,i) < eigrestol ) then
                n=n+1
                eigvec_merged_tmp(:,n) = eigvec(:,m,i)
             end if
          end do
       end do

       n=sum(numeig)
       call gather_sub(n2-n1+1,n,eigvec_merged_tmp)


       if ( allocated(eigvec_merged) ) then
          eigvec_merged_tmp(:,n+1:n+total_numeig) = eigvec_merged(:,:)
          deallocate( eigvec_merged )
          deallocate( eigval_merged )
          deallocate( resval_merged )
       end if


       allocate( dummy(ntmp,ntmp) )

       call z_MGS_QR(eigvec_merged_tmp,dummy)
       call z_MGS_QR(eigvec_merged_tmp,dummy)

       deallocate( dummy )

       allocate( eigval_merged_tmp(ntmp) )

       call z_rayleigh_ritz(eigvec_merged_tmp,0d0,eigval_merged_tmp,n1,n2)
       call z_normalize_eigenvectors(eigvec_merged_tmp)

       allocate( resval_merged_tmp(ntmp) )

       call z_eval_residual(eigval_merged_tmp,eigvec_merged_tmp,resval_merged_tmp,n1,n2)

       num_states_ccurve(:)=0
       m=0
       n=0
       do i=1,ntmp
          if ( resval_merged_tmp(i) < restol_ji ) then
             n=n+1
             do j=1,n_ccurve
                left =ssgamma(j)-ssrho(j)
                right=ssgamma(j)+ssrho(j)
                if ( left <= eigval_merged_tmp(i) .and. eigval_merged_tmp(i) <= right ) then
                   m=m+1
                   num_states_ccurve(j)=num_states_ccurve(j)+1
                end if
             end do
          end if
       end do

       call toc(tot_time,mpi_comm_world)

       if ( disp_switch ) then
          write(*,'(1x,"m,n,ntmp=",3i6,2f15.10,4f10.3)') m,n,ntmp,left,right &
          ,current_sseig_timers%total_time%ct,current_sseig_timers%total_time%et &
          ,tot_time%ct,tot_time%et
          write(*,'(1x,"num_states_ccurve:",50(i3,i5,2x))') ( j,num_states_ccurve(j),j=1,n_ccurve )
          call flush(6)
       end if

       allocate( resval_merged(n)       ) ; resval_merged=0.d0
       allocate( eigval_merged(n)       ) ; eigval_merged=0.d0
       allocate( eigvec_merged(n1:n2,n) ) ; eigvec_merged=0.d0
       total_numeig = n
       n=0
       do i=1,ntmp
          if ( resval_merged_tmp(i) < restol_ji ) then
             n=n+1
             resval_merged(n) = resval_merged_tmp(i)
             eigval_merged(n) = eigval_merged_tmp(i)
             eigvec_merged(:,n) = eigvec_merged_tmp(:,i)
             if ( disp_switch ) then
                write(*,*) n,i,eigval_merged_tmp(i),resval_merged_tmp(i), "   ji"
                call flush(6)
             end if
          end if
       end do
       if ( n /= total_numeig ) stop "stop@apply_sseig(2)"

       deallocate( resval_merged_tmp )
       deallocate( eigval_merged_tmp )
       deallocate( eigvec_merged_tmp )

       tot_num_states = sum( num_states_ccurve(1:n_ccurve) )
       if ( disp_switch ) write(*,*) "tot_num_states,old=",tot_num_states,tot_num_states_old
       if ( tot_num_states <= tot_num_states_old ) then
          exit
       else
          tot_num_states_old = tot_num_states
       end if

    end do ! out_iter


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

    deallocate(num_states_ccurve)
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

  SUBROUTINE gather_sub(mm,nn,f)
    implicit none
    integer,intent(IN) :: mm,nn
    complex(8),intent(INOUT) :: f(mm,nn)
    integer :: irank_c,n,i,ierr
    if ( .not.allocated(ir_c) ) then
       allocate( ir_c(0:nprocs_c-1) ) ; ir_c=0
       allocate( id_c(0:nprocs_c-1) ) ; id_c=0
    end if
    id_c(:)=0
    ir_c(:)=0
    do irank_c=0,nprocs_c-1
       do i=1,n_ccurve
          if ( ccurve_id_rank(i) == irank_c ) then 
             ir_c(irank_c)=ir_c(irank_c)+numeig(i)
          end if
       end do
    end do
    ir_c(:)=ir_c(:)*mm
    do irank_c=0,nprocs_c-1
       id_c(irank_c)=sum(ir_c(0:irank_c))-ir_c(irank_c)
    end do
    n=min(nn,id_c(myrank_c)/mm+1)
    call mpi_allgatherv(f(1,n),ir_c(myrank_c),MPI_COMPLEX16 &
         ,f,ir_c,id_c,MPI_COMPLEX16,comm_ccurve,ierr)
  END SUBROUTINE gather_sub

END MODULE apply_sseig_module
