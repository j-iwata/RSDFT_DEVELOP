MODULE prepare_sseig_module

  use parallel_module
  use iter_lin_solvers, only: orth_impl,unroll_step,nprocs_c,iterative_method_args &
                             ,comm_rhs,myrank_r,comm_ccurve,myrank_c,nprocs_r &
                             ,unit_kry_conv,tol_iter
  use atom_module, only: Natom
  use electron_module, only: Nband
  use sseig, only: rho_scale,ssrho,numeig,num_basis,my_ccurve_id,ccurve_id_rank &
                  ,n_ccurve_local,n_ccurve,opt,unit_ss_band_cor,SS_IO,n_rhs_set_local &
                  ,unit_ss_band_val,ssgamma,L_fixed,opt_N,out_iter_max,out_iter_init
  use timer_module, only: timer

  implicit none

CONTAINS


  SUBROUTINE prepare_sseig(unit,disp_switch_in)
    implicit none
    integer,intent(IN) :: unit
    logical,intent(IN) :: disp_switch_in
    integer :: i,j,m,n,s,k,n1,n2,ierr,numeig_per_circle, count
    real(8) :: left,right
    integer,allocatable :: cluster_id(:)
    real(8),allocatable :: cluster_center(:)
    integer :: iter,itermax,mindisidx
    real(8) :: step, mindis, tmpval0, tmpval
    logical :: change_flag
    integer,allocatable :: ir(:),id(:)
    integer :: max_read
    character(7) :: label
    integer :: eigvec_size2, eigvec_size3, eigvec_merged_size2
    type(timer) :: read_wf_time
    integer :: MI,MB

    MI = Natom
    MB = Nband

    if ( disp_switch_in ) write(*,'(a60," prepare_sseig")') repeat("-",60)
  
    left  = 0.d0
    right = 0.d0
    tol_iter(:) = 0.d0

    if ( myrank == 0 ) then
       max_read = max(2*MI,1000)
       do i=1,max_read
          read(unit,'(a7)') label
          if ( label=='# SSEIG' ) then
             read(unit,*) L_fixed, orth_impl, unroll_step
             read(unit,*) n_ccurve
             allocate( ssgamma(n_ccurve), ssrho(n_ccurve) )
             do j=1,n_ccurve
                read(unit,*) left,right
                ssgamma(j) = 0.5d0*(left + right)
                ssrho(j) = 0.5d0*(right - left)
                write(*,'(1x,i3,2f10.5,2f10.5)') j,left,right,ssgamma(j),ssrho(j)
             end do
             read(unit,*) SS_IO
             read(unit,*) opt_N,out_iter_max,out_iter_init
             read(unit,*) tol_iter(1:4)
             exit
          end if
          if ( i>=max_read ) stop "Label '# SSEIG' was not found."
       end do
       write(*,*) "L_fixed     =",L_fixed
       write(*,*) "orth_impl   =",orth_impl
       write(*,*) "unroll_step =",unroll_step
       write(*,*) "n_ccurve    =",n_ccurve
       do i=1,n_ccurve
          write(*,*) "circle id =",i
          write(*,*) "center of circle =",ssgamma(i)
          write(*,*) "radius of circle =",ssrho(i)
       end do
       write(*,*) "SS_IO =",SS_IO
       write(*,*) "opt_N=",opt_N
       write(*,*) "out_iter_max,out_iter_init=",out_iter_max,out_iter_init
       write(*,'(1x,"tol_iter=",4g20.10)') tol_iter(1:4)
    end if

    call mpi_bcast(L_fixed    ,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(orth_impl  ,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(unroll_step,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(n_ccurve   ,1,mpi_real8  ,0,mpi_comm_world,ierr)
    if ( myrank /= 0 ) then
       allocate( ssgamma(n_ccurve),ssrho(n_ccurve) )
    end if
    call mpi_bcast(ssgamma,n_ccurve,mpi_real8  ,0,mpi_comm_world,ierr)
    call mpi_bcast(ssrho  ,n_ccurve,mpi_real8  ,0,mpi_comm_world,ierr)
    call mpi_bcast(SS_IO  ,1       ,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(left ,1,mpi_real8  ,0,mpi_comm_world,ierr)
    call mpi_bcast(right,1,mpi_real8  ,0,mpi_comm_world,ierr)
    call mpi_bcast(opt_N,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(out_iter_max,1,mpi_integer,0,mpi_comm_world,ierr)
    if ( out_iter_init <= 0 ) then
       out_iter_init = out_iter_max
       if ( myrank == 0 ) write(*,*) "out_iter_init=",out_iter_init
    end if
    call mpi_bcast(out_iter_init,1,mpi_integer,0,mpi_comm_world,ierr)
    if ( tol_iter(1) == 0.d0 ) tol_iter(1)=1.d-10
    if ( tol_iter(2) == 0.d0 ) tol_iter(2)=tol_iter(1)
    if ( tol_iter(3) == 0.d0 ) tol_iter(3)=tol_iter(1)
    if ( tol_iter(4) == 0.d0 ) tol_iter(4)=tol_iter(2)
    if ( myrank == 0 ) write(*,'(1x,"tol_iter=",4g20.10)') tol_iter(1:4)
    call mpi_bcast(tol_iter,4,mpi_real8,0,mpi_comm_world,ierr)

! --- open files ---
    open(unit_kry_conv   ,file='ss_kry_conv.txt',status='replace')
    open(unit_ss_band_val,file='ss_band_val.txt',status='replace')
    open(unit_ss_band_cor,file='ss_band_cor.txt',status='replace')

    n1 = id_grid(myrank_g) + 1
    n2 = id_grid(myrank_g) + ir_grid(myrank_g)

!    allocate( ir(0:nprocs_b-1),id(0:nprocs_b-1) )
!    ir(0:nprocs_b-1)=ir_band(0:nprocs_b-1)*(n2-n1+1)
!    id(0:nprocs_b-1)=id_band(0:nprocs_b-1)*(n2-n1+1)
!    do s=MSP_0,MSP_1
!       do k=MBZ_0,MBZ_1
!          call mpi_allgatherv(unk(n1,MB_0,k,s),ir(id_class(myrank,4)) &
!         ,mpi_complex16,unk(n1,1,k,s),ir,id,mpi_complex16,comm_band,ierr)
!       end do
!    end do
!    deallocate( id,ir )

  ! START parallel setting
    nprocs_r = 1

    n=mod(myrank_b,nprocs_r)
    call mpi_comm_split(comm_band,n,myrank_b,comm_ccurve,ierr)
    call mpi_comm_size(comm_ccurve,nprocs_c,ierr)
    call mpi_comm_rank(comm_ccurve,myrank_c,ierr)

    if ( n_ccurve < nprocs_c ) stop "stop@prepare_sseig"

    call mpi_comm_split(comm_band,myrank_b/nprocs_r,myrank_b,comm_rhs,ierr)
    call mpi_comm_size(comm_rhs,nprocs_r,ierr)
    call mpi_comm_rank(comm_rhs,myrank_r,ierr)

    allocate( n_rhs_set_local(0:nprocs_r-1) )
  ! END parallel setting


  ! START parameter setting
    opt%delta = 1D-14
    iterative_method_args%epsmax = 1d-5
    iterative_method_args%epsmin = 1d-12
    iterative_method_args%imax   = 100*MB
  ! END parameter setting


  ! START set number of circles
    allocate( n_ccurve_local(0:nprocs_c-1) )
    allocate( ccurve_id_rank(n_ccurve)     )

    do j=0,nprocs_c-1
       count = 0
       do i=1,n_ccurve
          if ( mod(i-1,nprocs_c) == j ) then
             ccurve_id_rank(i) = j
             count = count + 1
          end if
       end do
       n_ccurve_local(j) = count
    end do

    if ( myrank == 0 ) then
       write(*,'(1x,a8,a16)') "rank_c","n_ccurve_local"
       do j=0,nprocs_c-1
          write(*,'(1x,i8,i16)') j,n_ccurve_local(j)
       end do
       write(*,*) "n_ccurve=",n_ccurve
       write(*,'(1x,a10,a10)') "ccurve_id","rank_c"
       do i=1,n_ccurve
          write(*,'(1x,i10,i10)') i,ccurve_id_rank(i)
       end do
    end if

    n=n_ccurve_local(myrank_c)
    allocate( my_ccurve_id(n) )

    count = 0
    do i=1,n_ccurve
       if ( mod(i-1,nprocs_c) == myrank_c ) then
          count = count + 1
          my_ccurve_id(count) = i
       end if
    end do

    if ( myrank_c == 0 .and. myrank_r == 0 .and. myrank_g == 0 ) then
       write(*,'(A,I3)') 'number of closed curves = ', n_ccurve
    end if

    call mpi_barrier(comm_ccurve,ierr)
  
  ! END set number of circles


    allocate( numeig(n_ccurve) )
    allocate( num_basis(n_ccurve_local(myrank_c)) )


  ! START set sseig module variables 'eigval' and 'eigvec'
    num_basis(:)=0
    numeig(:) = 0

    ssrho(:) = ssrho(:)*rho_scale

    if ( myrank_c == 0 .and. myrank_r == 0 .and. myrank_g == 0 ) then
       write(*,*) 'Information of closed curves'
       write(*,'(a4,2x,a8,2x,a14,1x,a14,2x,a7)') &
            'id', 'myrank_c', 'center', 'radius', 'num eig'
    end if
    do i=1,n_ccurve
       call mpi_barrier(comm_ccurve,ierr)
       if ( myrank == 0 ) then
          write(*,'(I4,2x,I8,2x,e14.7,1x,e14.7,2x,I7)') &
               i,ccurve_id_rank(i),ssgamma(i),ssrho(i),numeig(i)
       end if
    end do
    call mpi_barrier(comm_ccurve,ierr)
  ! END set sseig module variables eigval and eigvec

  END SUBROUTINE prepare_sseig


END MODULE prepare_sseig_module
