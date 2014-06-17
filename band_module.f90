MODULE band_module

  use watch_module

  implicit none

  PRIVATE
  PUBLIC :: read_band,band,band_sseig

  integer,parameter :: mnbk=100
  integer :: nbk,mb_band,mb2_band,maxiter_band
  integer :: nfki(mnbk),nskip_band
  real(8) :: ak(3,mnbk+1)
  real(8) :: esp_conv_tol=-1.d0

  integer,parameter :: unit_band_eigv=100
  integer,parameter :: unit_band_ovlp=110
  integer,parameter :: unit_band_dedk=120

CONTAINS


  SUBROUTINE read_band(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i,ierr
    integer,parameter :: max_read=10000
    character(8) :: cbuf,ckey
    include 'mpif.h'
    esp_conv_tol = 1.d-5
    maxiter_band = 50
    nskip_band=0
    if ( rank == 0 ) then
       rewind unit
       do i=1,max_read
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:4) == 'BAND' ) then
             backspace(unit)
             read(unit,*) cbuf,nbk,mb_band,mb2_band,esp_conv_tol,maxiter_band
             read(unit,*) ak(1,1:nbk+1)
             read(unit,*) ak(2,1:nbk+1)
             read(unit,*) ak(3,1:nbk+1)
             read(unit,*) nfki(1:nbk)
          else if ( ckey(1:8) == 'SKIPBAND' ) then
             backspace(unit)
             read(unit,*) cbuf,nskip_band
          end if
       end do
999    continue
       write(*,*) "nbk     =",nbk
       write(*,*) "mb_band =",mb_band
       write(*,*) "mb2_band=",mb2_band
       if ( esp_conv_tol < 0.d0 ) esp_conv_tol=1.d-5
       write(*,*) "esp_von_tol=",esp_conv_tol
       write(*,*) "maxiter_band=",maxiter_band
       write(*,*) "nskip_band=",nskip_band
    end if
    call mpi_bcast(nbk,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(ak,3*mnbk,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(nfki,mnbk,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(mb_band,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(mb2_band,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(esp_conv_tol,1,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(maxiter_band,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(nskip_band,1,mpi_integer,0,mpi_comm_world,ierr)
  END SUBROUTINE read_band


  SUBROUTINE band_sseig(disp_switch)

    use prepare_sseig_module
    use apply_sseig_module, only: apply_sseig,total_numeig
    use sseig

    use parallel_module
    use rgrid_module

    use array_bound_module
    use bz_module
    use aa_module, only: aa
    use bb_module, only: bb
    use kinetic_module, only: init_kinetic
    use ps_nloc2_module, only: prep_uvk_ps_nloc2,prep_rvk_ps_nloc2
    use momentum_module

    implicit none

    logical,intent(IN) :: disp_switch
    integer :: nktrj,i,j,k,n,ierr
    real(8) :: dak(3),sum0,sum1
    real(8),allocatable :: ktrj(:,:),pxyz(:,:)
    complex(8),allocatable :: eigvec_merged_old(:,:)
    integer,save :: total_numeig_old=0

#ifdef _DRSDFT_

    write(*,*) "Bnad calculation is not available for REAL8 version"

#else

    nktrj = sum( nfki(1:nbk) )
    if ( nktrj > 1 ) nktrj=nktrj+1

    allocate( ktrj(6,nktrj) )

    k=0
    do i=1,nbk
       dak(1:3) = ( ak(1:3,i+1) - ak(1:3,i) )/dble( nfki(i) )
       do j=1,nfki(i)
          k=k+1
          ktrj(1:3,k) = ak(1:3,i) + (j-1)*dak(1:3)
          ktrj(4:6,k) = matmul( bb(1:3,1:3),dak(1:3) )
       end do
    end do
    if ( nktrj > 1 ) then
       ktrj(1:3,k+1) = ak(1:3,nbk+1)
       ktrj(4:6,k+1) = 0.d0
    end if

    if ( disp_switch ) then
       write(*,*) "nktrj=",nktrj
       do k=1,nktrj
          write(*,'(1x,i4,2x,3f15.10,2x,3f15.10)') k,ktrj(1:6,k)
       end do
    end if

    if ( myrank == 0 ) then
       open(unit_band_eigv,file="band_eigv")
       open(unit_band_dedk,file="band_dedk")
       open(unit_band_ovlp,file="band_ovlp_2",form="unformatted")
       write(unit_band_eigv,'(1x,3f20.15)') bb(1:3,1)
       write(unit_band_eigv,'(1x,3f20.15)') bb(1:3,2)
       write(unit_band_eigv,'(1x,3f20.15)') bb(1:3,3)
    end if

    call prepare_sseig(1,disp_switch)

    do k=1,nktrj

       if ( k <= nskip_band ) then
          if ( DISP_SWITCH_PARALLEL ) write(*,*) "Band ",k," is skipped"
          cycle
       end if

       kbb(1:3,1) = ktrj(1:3,k)
       call init_kinetic(aa,bb,Nbzsm,kbb,disp_switch)
       call prep_uvk_ps_nloc2(1,Nbzsm,kbb)

       call apply_sseig( k, ktrj(4:6,k) )

       call prep_rvk_ps_nloc2(1,Nbzsm,kbb)

       eigvec_merged(:,:)=eigvec_merged(:,:)/sqrt(dV)

       if ( disp_switch ) then
          write(*,*) "k,total_numeig=",k,total_numeig
       end if

       allocate( pxyz(3,total_numeig) ) ; pxyz=0.d0

       do n=1,total_numeig
          sum0=sum( abs(eigvec_merged(:,n))**2 )*dV
          call mpi_allreduce(sum0,sum1,1,mpi_real8,mpi_sum,comm_grid,ierr)
          call calc_expectval_momentum(1,ML_0,ML_1,1,1,eigvec_merged(ML_0,n),pxyz(1,n))
          if ( disp_switch ) then
             write(*  ,'(1x,i5,2x,3g22.12)') n,eigval_merged(n),pxyz(3,n),sum1
          end if
       end do
       if ( myrank == 0 ) then
          write(unit_band_eigv,'(1x,2i6,3f20.12,i8)') k,total_numeig,kbb(1:3,1),0
          write(unit_band_dedk,'(1x,2i6,3f20.12)') k,total_numeig,kbb(1:3,1)
          do n=1,total_numeig
             write(unit_band_eigv,'(1x,i5,2x,2(1x,g22.12,1x,g15.5))') &
                  n,eigval_merged(n),resval_merged(n)
             write(unit_band_dedk,'(1x,i5,2x,2(g22.12,1x,3g15.5))') &
                  n,eigval_merged(n),pxyz(1:3,n)
          end do
       end if

       deallocate( pxyz )

       if ( .not.allocated(eigvec_merged_old) ) then
          allocate( eigvec_merged_old(ML_0:ML_1,total_numeig*2) )
          eigvec_merged_old(:,:)=(0.d0,0.d0)
       else
          call calc_overlap_2(ML_0,ML_1,total_numeig,total_numeig_old,eigvec_merged,eigvec_merged_old)
          n=size(eigvec_merged_old,2)
          if ( n < total_numeig ) then
             deallocate( eigvec_merged_old )
             allocate( eigvec_merged_old(ML_0:ML_1,total_numeig*2) )
             eigvec_merged_old=(0.d0,0.d0)
          end if
       end if
       eigvec_merged_old(:,1:total_numeig) = eigvec_merged(:,1:total_numeig)
       total_numeig_old = total_numeig

    end do ! k

    if ( myrank == 0 ) then
       close(unit_band_ovlp)
       close(unit_band_dedk)
       close(unit_band_eigv)
    end if

    deallocate( eigvec_merged_old )
    deallocate( ktrj )
#endif
  END SUBROUTINE band_sseig


  SUBROUTINE band(MBV_in,disp_switch)

    use parallel_module
    use rgrid_module

    use array_bound_module
    use bz_module
    use aa_module, only: aa
    use bb_module, only: bb
    use kinetic_module, only: init_kinetic
    use ps_nloc2_module, only: prep_uvk_ps_nloc2,prep_rvk_ps_nloc2
    use momentum_module
    use wf_module, only: esp,res,unk,gather_wf
    use cg_module
    use gram_schmidt_t_module
    use subspace_diag_sl_module
    use subspace_diag_la_module
    use esp_gather_module

    implicit none

    integer,intent(IN) :: MBV_in
    logical,intent(IN) :: disp_switch
    integer :: MBV,nktrj,i,j,k,s,n,ibz,ierr,iktrj,iter,Diter_band
    integer :: iktrj_0,iktrj_1,iktrj_2,iktrj_00,iktrj_tmp,ireq
    integer,allocatable :: ir_k(:),id_k(:)
    real(8) :: dak(3),sum0,sum1,max_err,max_err0
    real(8),allocatable :: ktrj(:,:),esp0(:,:,:),pxyz(:,:,:,:)
    real(8),allocatable :: kbb_tmp(:,:),esp_tmp(:,:,:),esp0_tmp(:,:,:)
    complex(8),allocatable :: unk0(:,:,:),unk00(:,:,:)
    logical :: disp_switch_parallel_bak,flag_end
    character(32) :: file_ovlp
    character(5) :: cc

    disp_switch_parallel_bak = disp_switch_parallel
!    disp_switch_parallel = .false.

    Diter_band = maxiter_band

    flag_end = .false.

    MBV=MBV_in
    if ( MBV < 1 .or. mb_band < MBV ) MBV=1
    if ( disp_switch_parallel ) write(*,*) "MBV=",MBV

    nktrj = sum( nfki(1:nbk) )
    if ( nktrj > 1 ) nktrj=nktrj+1

    allocate( ktrj(6,nktrj) )

    k=0
    do i=1,nbk
       dak(1:3) = ( ak(1:3,i+1) - ak(1:3,i) )/dble( nfki(i) )
       do j=1,nfki(i)
          k=k+1
          ktrj(1:3,k) = ak(1:3,i) + (j-1)*dak(1:3)
          ktrj(4:6,k) = matmul( bb(1:3,1:3),dak(1:3) )
       end do
    end do
    if ( nktrj > 1 ) then
       ktrj(1:3,k+1) = ak(1:3,nbk+1)
       ktrj(4:6,k+1) = 0.d0
    end if

    if ( disp_switch ) then
       write(*,*) "nktrj=",nktrj
       do k=1,nktrj
          write(*,'(1x,i4,2x,3f15.10,2x,3f15.10)') k,ktrj(1:6,k)
       end do
    end if

    allocate( ir_k(0:np_bzsm-1),id_k(0:np_bzsm-1) )
    id_k(:)=0
    ir_k(:)=0
    do k=1,nktrj
       i=mod(k-1,np_bzsm)
       ir_k(i)=ir_k(i)+1
    end do
    i=maxval(ir_k)
    ir_k(:)=i
    do i=0,np_bzsm-1
       id_k(i) = sum(ir_k(0:i)) - ir_k(i)
    end do
    if ( myrank == 0 ) then
       write(*,'(1x,3a6)') "rank","id_k","ir_k"
       do i=0,np_bzsm-1
          write(*,'(1x,3i6)') i,id_k(i),ir_k(i)
       end do
       write(*,*) "sum(ir_k),nktrj=",sum(ir_k),nktrj
    end if

    if ( MB < mb_band ) then
       call modify_mb
       call modify_arraysize
    end if

    allocate( esp0_tmp(MB,0:np_bzsm-1,MSP)   ) ; esp0_tmp=0.d0
    allocate( esp_tmp(MB,0:np_bzsm-1,MSP)    ) ; esp_tmp=0.d0
    allocate( kbb_tmp(3,0:np_bzsm-1)         ) ; kbb_tmp=0.d0
    allocate( pxyz(3,MB,0:np_bzsm-1,MSP)     ) ; pxyz=0.d0
    allocate( esp0(MB,MBZ,MSP)               ) ; esp0=0.d0
    allocate( unk0(ML_0:ML_1,MB,MSP_0:MSP_1) ) ; unk0=0.d0

!    k=MBZ_0

    if ( myrank == 0 ) then
       open(unit_band_eigv,file="band_eigv")
       open(unit_band_dedk,file="band_dedk")
       write(unit_band_eigv,'(1x,3f20.15)') bb(1:3,1)
       write(unit_band_eigv,'(1x,3f20.15)') bb(1:3,2)
       write(unit_band_eigv,'(1x,3f20.15)') bb(1:3,3)
    end if

    if ( np_bzsm == 1 ) then
       if ( myrank == 0 ) then
          open(unit_band_ovlp,file="band_ovlp")
       end if
    else
       write(cc,'(i5.5)') myrank
       file_ovlp = "band_ovlp_"//trim(adjustl(cc))
       if ( myrank_g == 0 .and. myrank_b == 0 ) then
          open(unit_band_ovlp,file=file_ovlp)
       end if
    end if

    iktrj_0 = id_k(myrank_k)+1
    iktrj_1 = id_k(myrank_k)+ir_k(myrank_k)
    iktrj_2 = id_k(myrank_k)+maxval(ir_k)

    loop_iktrj : do iktrj = iktrj_0, iktrj_2

       if ( iktrj <= nskip_band ) then
          if ( DISP_SWITCH_PARALLEL ) write(*,*) "Band ",iktrj," is skipped"
          cycle
       end if

       iktrj_00 = id_k(0) + iktrj - iktrj_0 + 1

       kbb(1:3,MBZ_0) = ktrj(1:3,min(iktrj,iktrj_1,nktrj))
       call mpi_allgather(kbb(1,MBZ_0),3,mpi_real8,kbb_tmp,3,mpi_real8,comm_bzsm,ierr)
       if ( myrank == 0 ) then
          write(*,*) "kbb_tmp"
          do ibz=0,np_bzsm-1
             write(*,'(1x,i8,3f10.5)') id_k(ibz)+iktrj-iktrj_0+1,kbb_tmp(1:3,ibz)
          end do
       endif

       call init_kinetic(aa,bb,Nbzsm,kbb,disp_switch)
       call prep_uvk_ps_nloc2(MBZ_0,MBZ_0,kbb(1,MBZ_0))

       esp0(:,:,:) = 1.d10

       do iter=1,Diter_band
          if ( disp_switch ) write(*,'(a60,"band_iter=",i3)') repeat("-",60),iter
          do s=MSP_0,MSP_1
             call conjugate_gradient(ML_0,ML_1,MB,MBZ_0,s,Ncg,iswitch_gs &
                  ,unk(ML_0,1,MBZ_0,s),esp(1,MBZ_0,s),res(1,MBZ_0,s))
             call gram_schmidt_t(1,MB,MBZ_0,s)
#ifdef _LAPACK_
             call subspace_diag_la(MBZ_0,s)
#else
             call subspace_diag_sl(MBZ_0,s,.false.)
#endif
          end do
          max_err = maxval( abs(  esp(1:mb2_band,MBZ_0,MSP_0:MSP_1) &
                                -esp0(1:mb2_band,MBZ_0,MSP_0:MSP_1) ) )
          if ( iktrj_1 < iktrj ) max_err=0.d0
          call mpi_allreduce(max_err,max_err0,1,MPI_REAL8,MPI_MAX,comm_spin,ierr)
          call mpi_allreduce(max_err0,max_err,1,MPI_REAL8,MPI_MAX,comm_bzsm,ierr)
          if ( disp_switch ) write(*,*) "iktrj,max_err,iter=",iktrj,max_err,iter
          if ( max_err < esp_conv_tol ) then
             exit
          else if ( iter == Diter_band ) then
             stop "band is not converged"
          end if
          esp0(:,:,:)=esp(:,:,:)
          call global_watch(flag_end)
          if ( flag_end ) then
             if ( myrank == 0 ) write(*,*) "etime limit !!!"
             exit loop_iktrj
          end if
       end do ! iter

       call esp_gather(MB,Nbzsm,MSP,esp)

       call gather_wf

#ifndef _DRSDFT_
       if ( iktrj_0 < iktrj ) then
          if ( np_bzsm == 1 ) then
             if ( myrank == 0 ) then
                write(unit_band_ovlp,'(1x,2i6,3f20.12)') iktrj,MB,kbb(1:3,MBZ_0)
             end if
          else
             if ( myrank_g == 0 .and. myrank_b == 0 ) then
                write(unit_band_ovlp,'(1x,2i6,3f20.12)') iktrj,MB,kbb(1:3,MBZ_0)
             end if
          end if
       end if
       do s=MSP_0,MSP_1
          if ( np_bzsm == 1 ) then
             if ( 1 < iktrj .and. iktrj <= iktrj_1 ) then
                call calc_overlap(ML_0,ML_1,MB_0,MB_1,MB,unk(:,1,MBZ_0,s),unk0(:,1,s))
             end if
          else
             if ( iktrj_0 < iktrj .and. iktrj <= iktrj_1 ) then
                call calc_overlap_kpara(ML_0,ML_1,MB_0,MB_1,MB,unk(:,1,MBZ_0,s),unk0(:,1,s))
             end if
             if ( iktrj == iktrj_0 ) then
                if ( .not.allocated(unk00) ) then
                   allocate( unk00(ML_0:ML_1,MB,MSP_0:MSP_1) )
                   unk00=0.d0
                end if
                unk00(:,:,s)=unk(:,:,MBZ_0,s)
             end if
             if ( iktrj == iktrj_1 ) then
                call send_wf_band(ML_1-ML_0+1,MB,unk00(ML_0,1,s))
                unk0(:,:,s) = unk(:,:,MBZ_0,s)
                call recv_wf_band(ML_1-ML_0+1,MB,unk(ML_0,1,MBZ_0,s))
                if ( myrank_k < np_bzsm-1 ) then
                   if ( myrank_g == 0 .and. myrank_b == 0 ) then
                      write(unit_band_ovlp,'(1x,2i6,3f20.12)') iktrj+1,MB,ktrj(1:3,iktrj+1)
                   end if
                   call calc_overlap_kpara(ML_0,ML_1,MB_0,MB_1,MB,unk(:,1,MBZ_0,s),unk0(:,1,s))
                end if
                if ( allocated(unk00) ) deallocate(unk00)
             end if
          end if
          unk0(:,:,s) = unk(:,:,MBZ_0,s)
       end do
#endif

       call prep_rvk_ps_nloc2(MBZ_0,MBZ_0,kbb(1,MBZ_0))

       pxyz(:,:,:,:)=0.d0
       do n=1,mb2_band
          do s=MSP_0,MSP_1
#ifndef _DRSDFT_
             call calc_expectval_momentum(MBZ_0,ML_0,ML_1,1,1,unk(ML_0,n,MBZ_0,1),pxyz(1,n,myrank_k,s))
#endif
          end do ! s
       end do ! n

       do s=MSP_0,MSP_1
          call mpi_allgather(pxyz(1,1,myrank_k,s),3*MB,MPI_REAL8 &
               ,pxyz(1,1,0,s),3*MB,MPI_REAL8,comm_bzsm,ierr)
       end do ! s
       call mpi_allgather(pxyz(1,1,0,MSP_0),3*MB*np_bzsm*(MSP_1-MSP_0+1),MPI_REAL8 &
            ,pxyz,3*MB*np_bzsm*(MSP_1-MSP_0+1),MPI_REAL8,comm_spin,ierr)

       do s=1,MSP
          call mpi_allgather(esp(1,MBZ_0,s),MB,mpi_real8,esp_tmp(1,0,s),MB,mpi_real8,comm_bzsm,ierr)
       end do
       do s=1,MSP
          esp0_tmp(:,myrank_k,s)=esp0(:,MBZ_0,s)
          call mpi_allgather(esp0_tmp(1,myrank_k,s),MB,mpi_real8,esp0_tmp(1,0,s),MB,mpi_real8,comm_bzsm,ierr)
       end do

       do ibz=0,np_bzsm-1
          iktrj_tmp = iktrj_00 + ibz
          if ( iktrj_tmp > nktrj ) exit
          if ( myrank == 0 ) then
             write(unit_band_eigv,'(1x,2i6,3f20.12,i8)') iktrj_tmp,mb2_band,kbb_tmp(1:3,ibz),MBV
             write(unit_band_dedk,'(1x,2i6,3f20.12)') iktrj_tmp,mb2_band,kbb_tmp(1:3,ibz)
          end if
          do n=1,mb2_band
             if ( disp_switch ) then
                write(*,'(1x,i5,2x,2(1x,g22.12,1x,g15.5))') n &
                     ,( esp_tmp(n,ibz,s),abs(esp_tmp(n,ibz,s)-esp0_tmp(n,ibz,s)),s=1,MSP )
             end if
             if ( myrank == 0 ) then
                write(unit_band_eigv,'(1x,i5,2x,2(1x,g22.12,1x,g15.5))') n &
                     ,( esp_tmp(n,ibz,s),abs(esp_tmp(n,ibz,s)-esp0_tmp(n,ibz,s)),s=1,MSP )
                write(unit_band_dedk,'(1x,i5,2x,2(g22.12,1x,3g15.5))') n &
                     ,( esp_tmp(n,ibz,s),pxyz(1:3,n,ibz,s),s=1,MSP )
             end if
          end do ! n
       end do ! ibz

    end do loop_iktrj

    if ( myrank == 0 ) then
       close(unit_band_dedk)
       close(unit_band_eigv)
    end if
    if ( np_bzsm == 1 ) then
       if ( myrank == 0 ) then
          close(unit_band_ovlp)
       end if
    else
       if ( myrank_g == 0 .and. myrank_b == 0 ) then
          close(unit_band_ovlp)
       end if
    end if

    deallocate( unk0 )
    deallocate( esp0 )
    deallocate( pxyz )
    deallocate( kbb_tmp )
    deallocate( esp_tmp )
    deallocate( esp0_tmp )
    deallocate( id_k )
    deallocate( ir_k )
    deallocate( ktrj )

    disp_switch_parallel = disp_switch_parallel_bak

  END SUBROUTINE band

  SUBROUTINE send_wf_band(mm,nn,f)
    use parallel_module
    implicit none
    integer,intent(IN) :: mm,nn
    complex(8),intent(IN) :: f(mm,nn)
    integer :: ierr,ireq,istatus(MPI_STATUS_SIZE)
    if ( myrank_k > 0 ) then
       call MPI_ISEND(f,mm*nn,MPI_COMPLEX16,myrank_k-1,1,comm_bzsm,ireq,ierr)
       call MPI_WAIT(ireq,istatus,ierr)
    end if
  END SUBROUTINE send_wf_band
  SUBROUTINE recv_wf_band(mm,nn,f)
    use parallel_module
    implicit none
    integer,intent(IN) :: mm,nn
    complex(8),intent(OUT) :: f(mm,nn)
    integer :: ierr,ireq,istatus(MPI_STATUS_SIZE)
    if ( myrank_k < np_bzsm-1 ) then
       call MPI_IRECV(f,mm*nn,MPI_COMPLEX16,myrank_k+1,1,comm_bzsm,ireq,ierr)
       call MPI_WAIT(ireq,istatus,ierr)
    end if
  END SUBROUTINE recv_wf_band

  SUBROUTINE write_data_band(iktrj,rank,mm,nn,unk)
    implicit none
    integer,intent(IN) :: iktrj,rank,mm,nn
    complex(8),intent(IN) :: unk(mm,nn)
    character(5) :: c1,c2
    character(32) :: filename
    write(c1,'(i5.5)') iktrj
    write(c2,'(i5.5)') rank
    filename = "wf_band."//trim(adjustl(c1))//trim(adjustl(c2))
    open(10+rank,file=filename,form='unformatted')
    write(10+rank) unk
    call flush(10+rank)
    close(10+rank)
  END SUBROUTINE write_data_band

  SUBROUTINE read_data_band(iktrj,rank,mm,nn,unk)
    implicit none
    integer,intent(IN) :: iktrj,rank,mm,nn
    complex(8),intent(OUT) :: unk(mm,nn)
    character(5) :: c1,c2
    character(32) :: filename
    write(c1,'(i5.5)') iktrj-1
    write(c2,'(i5.5)') rank
    filename = "wf_band."//trim(adjustl(c1))//trim(adjustl(c2))
    open(10+rank,file=filename,form='unformatted')
    read(10+rank) unk
    close(10+rank)
  END SUBROUTINE read_data_band

  SUBROUTINE calc_overlap(n1,n2,m1,m2,mm,f1,f0)
    use parallel_module
    implicit none
    integer,intent(IN) :: n1,n2,m1,m2,mm
    complex(8),intent(IN) :: f1(n1:n2,mm),f0(n1:n2,mm)
    real(8),allocatable :: sqovlp01(:,:),sqovlp10(:,:),sqovlp10_tmp(:,:)
    complex(8),allocatable :: zwork(:,:),ovlp01(:,:)
    complex(8),parameter :: zero=(0.d0,0.d0),one=(1.d0,0.d0)
    integer,parameter :: hyaku=100, jyuu=10
    integer :: nn,ierr,ib1,ib2,nib,i,j
    integer,allocatable :: indx(:),maxloc10(:,:),maxloc10_tmp(:,:)

    nn = n2 - n1 + 1

    allocate( maxloc10(mm,jyuu) ) ; maxloc10=0
    allocate( sqovlp10(mm,jyuu) ) ; sqovlp10=0.d0

    allocate( zwork(mm,hyaku) )
    allocate( ovlp01(mm,hyaku) )
    allocate( sqovlp01(mm,hyaku) )
    allocate( indx(mm) )

    do ib1=m1,m2,hyaku
       ib2=min(ib1+hyaku-1,m2,mm)
       nib=ib2-ib1+1
       if ( nib < 0 ) cycle
       call zgemm('C','N',mm,nib,nn,one,f0(n1,1),nn,f1(n1,ib1),nn,zero,zwork(1,1),mm)
       call mpi_allreduce(zwork,ovlp01,hyaku*mm,mpi_complex16,mpi_sum,comm_grid,ierr)
       sqovlp01(1:mm,1:nib) = abs( ovlp01(1:mm,1:nib) )**2
       do j=1,nib
          call indexx(mm,sqovlp01(1,j),indx)
          do i=mm,max(mm-jyuu+1,1),-1
             maxloc10(j-1+ib1,mm-i+1) = indx(i)
             sqovlp10(j-1+ib1,mm-i+1) = sqovlp01( indx(i),j )
          end do
       end do
    end do ! ib1

    deallocate( indx )
    deallocate( sqovlp01 )
    deallocate( ovlp01 )
    deallocate( zwork )

    allocate( maxloc10_tmp(mm,jyuu) ) ; maxloc10_tmp(:,:)=maxloc10(:,:)
    allocate( sqovlp10_tmp(mm,jyuu) ) ; sqovlp10_tmp(:,:)=sqovlp10(:,:)
    call mpi_reduce(maxloc10_tmp,maxloc10,mm*jyuu,mpi_integer,mpi_sum,0,comm_band,ierr)
    call mpi_reduce(sqovlp10_tmp,sqovlp10,mm*jyuu,mpi_real8,mpi_sum,0,comm_band,ierr)
    deallocate( sqovlp10_tmp,maxloc10_tmp )

    if ( myrank == 0 ) then
       do i=1,mm
          write(unit_band_ovlp,'(1x,i6,2x,10i6)') i,( maxloc10(i,j),j=1,jyuu )
       end do
       do i=1,mm
          write(unit_band_ovlp,'(1x,i6,2x,10g15.6)') i,( sqovlp10(i,j),j=1,jyuu )
       end do
    end if

    deallocate( sqovlp10 )
    deallocate( maxloc10 )

  END SUBROUTINE calc_overlap

  SUBROUTINE calc_overlap_kpara(n1,n2,m1,m2,mm,f1,f0)
    use parallel_module
    implicit none
    integer,intent(IN) :: n1,n2,m1,m2,mm
    complex(8),intent(IN) :: f1(n1:n2,mm),f0(n1:n2,mm)
    real(8),allocatable :: sqovlp01(:,:),sqovlp10(:,:),sqovlp10_tmp(:,:)
    complex(8),allocatable :: zwork(:,:),ovlp01(:,:)
    complex(8),parameter :: zero=(0.d0,0.d0),one=(1.d0,0.d0)
    integer,parameter :: hyaku=100, jyuu=10
    integer :: nn,ierr,ib1,ib2,nib,i,j
    integer,allocatable :: indx(:),maxloc10(:,:),maxloc10_tmp(:,:)

    nn = n2 - n1 + 1

    allocate( maxloc10(mm,jyuu) ) ; maxloc10=0
    allocate( sqovlp10(mm,jyuu) ) ; sqovlp10=0.d0

    allocate( zwork(mm,hyaku) )
    allocate( ovlp01(mm,hyaku) )
    allocate( sqovlp01(mm,hyaku) )
    allocate( indx(mm) )

    do ib1=m1,m2,hyaku
       ib2=min(ib1+hyaku-1,m2,mm)
       nib=ib2-ib1+1
       if ( nib < 0 ) cycle
       call zgemm('C','N',mm,nib,nn,one,f0(n1,1),nn,f1(n1,ib1),nn,zero,zwork(1,1),mm)
       call mpi_allreduce(zwork,ovlp01,hyaku*mm,mpi_complex16,mpi_sum,comm_grid,ierr)
       sqovlp01(1:mm,1:nib) = abs( ovlp01(1:mm,1:nib) )**2
       do j=1,nib
          call indexx(mm,sqovlp01(1,j),indx)
          do i=mm,max(mm-jyuu+1,1),-1
             maxloc10(j-1+ib1,mm-i+1) = indx(i)
             sqovlp10(j-1+ib1,mm-i+1) = sqovlp01( indx(i),j )
          end do
       end do
    end do ! ib1

    deallocate( indx )
    deallocate( sqovlp01 )
    deallocate( ovlp01 )
    deallocate( zwork )

    allocate( maxloc10_tmp(mm,jyuu) ) ; maxloc10_tmp(:,:)=maxloc10(:,:)
    allocate( sqovlp10_tmp(mm,jyuu) ) ; sqovlp10_tmp(:,:)=sqovlp10(:,:)
    call mpi_reduce(maxloc10_tmp,maxloc10,mm*jyuu,mpi_integer,mpi_sum,0,comm_band,ierr)
    call mpi_reduce(sqovlp10_tmp,sqovlp10,mm*jyuu,mpi_real8,mpi_sum,0,comm_band,ierr)
    deallocate( sqovlp10_tmp,maxloc10_tmp )

    if ( myrank_g == 0 .and. myrank_b == 0 ) then
       do i=1,mm
          write(unit_band_ovlp,'(1x,i6,2x,10i6)') i,( maxloc10(i,j),j=1,jyuu )
       end do
       do i=1,mm
          write(unit_band_ovlp,'(1x,i6,2x,10g15.6)') i,( sqovlp10(i,j),j=1,jyuu )
       end do
    end if

    deallocate( sqovlp10 )
    deallocate( maxloc10 )

  END SUBROUTINE calc_overlap_kpara

  SUBROUTINE calc_overlap_2(n1,n2,mm1,mm0,f1,f0)
    use parallel_module
    implicit none
    integer,intent(IN) :: n1,n2,mm1,mm0
    complex(8),intent(IN) :: f1(n1:n2,mm1),f0(n1:n2,mm0)
    real(8),allocatable :: sqovlp01(:,:),sqovlp10(:,:)
    complex(8),allocatable :: zwork(:,:),ovlp01(:,:)
    complex(8),parameter :: zero=(0.d0,0.d0),one=(1.d0,0.d0)
    integer,parameter :: hyaku=100
    integer :: nn,ierr,ib1,ib2,nib,i,j
    integer,allocatable :: indx(:),maxloc10(:,:)

    nn = n2 - n1 + 1

    allocate( maxloc10(mm1,mm0) ) ; maxloc10=0
    allocate( sqovlp10(mm1,mm0) ) ; sqovlp10=0.d0
    allocate( zwork(mm0,hyaku)  )
    allocate( ovlp01(mm0,hyaku) )
    allocate( sqovlp01(mm0,mm1) )
    allocate( indx(mm0)         )

    do ib1=1,mm1,hyaku
       ib2=min(ib1+hyaku-1,mm1)
       nib=ib2-ib1+1
       if ( nib < 0 ) cycle
       call zgemm('C','N',mm0,nib,nn,one,f0(n1,1),nn,f1(n1,ib1),nn,zero,zwork(1,1),mm0)
       call mpi_allreduce(zwork,ovlp01,hyaku*mm0,mpi_complex16,mpi_sum,comm_grid,ierr)
       sqovlp01(1:mm0,ib1:ib2) = abs( ovlp01(1:mm0,1:nib) )**2
       do j=ib1,ib2
          call indexx(mm0,sqovlp01(1,j),indx)
          do i=mm0,1,-1
             maxloc10(j,mm0-i+1) = indx(i)
             sqovlp10(j,mm0-i+1) = sqovlp01( indx(i),j )
          end do
       end do
    end do ! ib1

    deallocate( indx )
    deallocate( sqovlp01 )
    deallocate( ovlp01 )
    deallocate( zwork )

!    allocate( maxloc10_tmp(mm,jyuu) ) ; maxloc10_tmp(:,:)=maxloc10(:,:)
!    allocate( sqovlp10_tmp(mm,jyuu) ) ; sqovlp10_tmp(:,:)=sqovlp10(:,:)
!    call mpi_reduce(maxloc10_tmp,maxloc10,mm*jyuu,mpi_integer,mpi_sum,0,comm_band,ierr)
!    call mpi_reduce(sqovlp10_tmp,sqovlp10,mm*jyuu,mpi_real8,mpi_sum,0,comm_band,ierr)
!    deallocate( sqovlp10_tmp,maxloc10_tmp )

    if ( myrank == 0 ) then
       write(unit_band_ovlp) mm1,mm0
       write(unit_band_ovlp) maxloc10(:,:)
       write(unit_band_ovlp) sqovlp10(:,:)
    end if

    deallocate( sqovlp10 )
    deallocate( maxloc10 )

  END SUBROUTINE calc_overlap_2

#ifdef _OLD_VERSION_TEST_
  SUBROUTINE calc_overlap(n1,n2,m1,m2,mm,f1,f0,maxloc_ovlp)
    use parallel_module
    implicit none
    integer,intent(IN) :: n1,n2,m1,m2,mm
    complex(8),intent(IN) :: f1(n1:n2,mm),f0(n1:n2,mm)
    integer,intent(OUT) :: maxloc_ovlp(mm)
    real(8) :: xmax
    real(8),allocatable :: ovlp(:),ovlp2(:,:),ovlp_tmp(:,:)
    complex(8),allocatable :: ovlp0(:,:),ovlp1(:,:)
    complex(8),parameter :: zero=(0.d0,0.d0),one=(1.d0,0.d0)
    integer,parameter :: hyaku=100, jyuu=10
    integer :: nn,ierr,ib1,ib2,nib,i,j,k,l,imax,loop
    integer,allocatable :: indx(:),count_maxovlp(:),maxloc_tmp(:,:)
    integer,allocatable :: count_tmp(:,:)
    logical,allocatable :: msk(:),msk2(:)

    nn = n2 - n1 + 1

    maxloc_ovlp(:)=0

    allocate( ovlp_tmp(mm,jyuu) )
    allocate( ovlp(mm)        )
    allocate( ovlp0(mm,hyaku) )
    allocate( ovlp1(mm,hyaku) )
    allocate( ovlp2(mm,hyaku) )
    allocate( msk(mm) ) ; msk=.true.
    allocate( msk2(hyaku) ) ; msk2=.true.
    allocate( count_maxovlp(mm) )

    allocate( count_tmp(mm,jyuu) )
    allocate( maxloc_tmp(mm,jyuu) )
    allocate( indx(mm) )

    maxloc_tmp(:,:)=0
    ovlp_tmp(:,:)=0.d0

    do ib1=m1,m2,hyaku
       ib2=min(ib1+hyaku-1,m2,mm)
       nib=ib2-ib1+1
       if ( nib < 0 ) cycle
       call zgemm('C','N',mm,nib,nn,one,f0(n1,1),nn,f1(n1,ib1),nn,zero,ovlp0(1,1),mm)
       call mpi_allreduce(ovlp0,ovlp1,hyaku*mm,mpi_complex16,mpi_sum,comm_grid,ierr)
       ovlp2(1:mm,1:nib) = abs( ovlp1(1:mm,1:nib) )**2
       do k=1,nib
          call indexx(mm,ovlp2(1,k),indx)
          do i=mm,max(mm-jyuu+1,1),-1
             maxloc_tmp(k-1+ib1,mm-i+1) = indx(i)
             ovlp_tmp(k-1+ib1,mm-i+1) = ovlp2( indx(i),k )
          end do
       end do
    end do ! ib1

    count_tmp = 0
    do j=1,jyuu
       do i=1,mm
          if ( maxloc_tmp(i,j) == 0 ) cycle
          count_tmp( maxloc_tmp(i,j),j ) = count_tmp( maxloc_tmp(i,j),j ) + 1
       end do
    end do

    do loop=1,2

       do i=1,mm
          if ( count_tmp(i,1) == 1 ) then
             do j=1,mm
                if ( maxloc_tmp(j,1) == i ) maxloc_tmp(j,2:jyuu)=0
                do k=2,jyuu
                   if ( maxloc_tmp(j,k) == i ) maxloc_tmp(j,k)=0
                end do
             end do
          end if
       end do
       count_tmp(:,:)=0
       do j=1,jyuu
          do i=1,mm
             if ( maxloc_tmp(i,j) == 0 ) cycle
             count_tmp( maxloc_tmp(i,j),j ) = count_tmp( maxloc_tmp(i,j),j ) + 1
          end do
       end do

       do i=1,mm
          if ( count_tmp(i,1) == 0 .and. count_tmp(i,2) == 1 ) then
             do j=1,mm
                if ( maxloc_tmp(j,2) == i ) then
                   maxloc_tmp(j,1)=0
                   maxloc_tmp(j,3:jyuu)=0
                end if
                do k=3,jyuu
                   if ( maxloc_tmp(j,k) == i ) maxloc_tmp(j,k)=0
                end do
             end do
          end if
       end do
       count_tmp(:,:)=0
       do j=1,jyuu
          do i=1,mm
             if ( maxloc_tmp(i,j) == 0 ) cycle
             count_tmp( maxloc_tmp(i,j),j ) = count_tmp( maxloc_tmp(i,j),j ) + 1
          end do
       end do

    end do ! loop

    do i=1,mm
       if ( count(maxloc_tmp(i,:)/=0) > 1 ) then
          do j=1,jyuu
             l=maxloc_tmp(i,j)
             if ( l /= 0 ) then
                where( maxloc_tmp == l )
                   maxloc_tmp = 0
                end where
                maxloc_tmp(i,:)=0
                maxloc_tmp(i,j)=l
             end if
          end do
       end if
    end do
    count_tmp(:,:)=0
    do j=1,jyuu
       do i=1,mm
          if ( maxloc_tmp(i,j) == 0 ) cycle
          count_tmp( maxloc_tmp(i,j),j ) = count_tmp( maxloc_tmp(i,j),j ) + 1
       end do
    end do

    do i=1,mm
       if ( count(maxloc_tmp(i,:)/=0) /= 1 ) then
          if ( myrank == 0 ) then
             write(*,'(1x,"band connection failed",i5,10i4)') i,maxloc_tmp(i,1:10)
          end if
       end if
       do j=1,jyuu
          if ( maxloc_tmp(i,j) /= 0 ) then
             maxloc_ovlp(i) = maxloc_tmp(i,j)
             exit
          end if
       end do
    end do

    deallocate( indx )
    deallocate( maxloc_tmp )
    deallocate( count_tmp )
    deallocate( count_maxovlp )
    deallocate( msk2 )
    deallocate( msk )
    deallocate( ovlp2 )
    deallocate( ovlp1 )
    deallocate( ovlp0 )
    deallocate( ovlp  )
    deallocate( ovlp_tmp )

  END SUBROUTINE calc_overlap
#endif

  SUBROUTINE modify_mb

    use array_bound_module
    use electron_module
    use parallel_module
    use scalapack_module
    use subspace_diag_module, only: prep_subspace_diag
    use subspace_mate_sl_0_module
    use subspace_rotv_sl_0_module

    implicit none

    integer :: i,j

    if ( mb_band <= Nband ) return

    Nband = mb_band

    call prep_0_scalapack(Nband,myrank==0)

    ir_band(:)=0
    do i=0,Nband-1
       j=mod(i,np_band)
       ir_band(j) = ir_band(j) + 1
    end do
    do i=0,np_band-1
       id_band(i) = sum( ir_band(0:i) ) - ir_band(i)
    end do
    MB   = Nband
    MB_0 = id_band(myrank_b) + 1
    MB_1 = id_band(myrank_b) + ir_band(myrank_b)

    call prep_subspace_diag(Nband,myrank==0)
    call reset_subspace_mate_sl_0
    call reset_subspace_rotv_sl_0

  END SUBROUTINE modify_mb


  SUBROUTINE modify_arraysize

    use wf_module
    use array_bound_module

    implicit none

    complex(8),allocatable :: utmp(:,:,:,:)
    integer :: ml_old,mb_old,mbz_old,msp_old

    ml_old  = size( unk,1 )
    mb_old  = size( unk,2 )
    mbz_old = size( unk,3 )
    msp_old = size( unk,4 )

    allocate( utmp(ml_old,mb_old,mbz_old,msp_old) )
    utmp=unk

    call init_wf

    unk(ML_0:ML_0+ml_old-1,1:mb_old,MBZ_0:MBZ_0+mbz_old-1,MSP_0:MSP_0+msp_old-1) &
         = utmp(1:ml_old,1:mb_old,1:mbz_old,1:msp_old)

    deallocate( utmp )

  END SUBROUTINE modify_arraysize


END MODULE band_module
