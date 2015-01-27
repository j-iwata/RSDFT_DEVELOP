MODULE band_module

  use watch_module
  use parallel_module
  use rgrid_module
  use array_bound_module
  use bz_module
  use aa_module, only: aa
  use bb_module, only: bb
  use kinetic_module, only: init_kinetic
  use momentum_module
  use wf_module
  use cg_module
  use gram_schmidt_t_module
  use subspace_diag_module
  use esp_gather_module
  use electron_module
  use scalapack_module
  use subspace_mate_sl_0_module, only: reset_subspace_mate_sl_0
  use subspace_rotv_sl_0_module, only: reset_subspace_rotv_sl_0
  use band_variables, only: nfki,nbk,ak,nskip_band &
       ,esp_conv_tol, mb_band, mb2_band, maxiter_band &
       ,unit_band_eigv,unit_band_dedk,unit_band_ovlp,read_band
  use sweep_module, only: calc_sweep, init_sweep
  use pseudopot_module, only: pselect
  use ps_nloc2_module, only: prep_uvk_ps_nloc2, prep_rvk_ps_nloc2
  use ps_nloc3_module, only: prep_ps_nloc3, init_ps_nloc3
  use io_module, only: Init_IO
  use xc_hybrid_module, only: iflag_hybrid, prep_kq_xc_hybrid
  use fock_ffte_module, only: init_fock_ffte
  
  implicit none

  PRIVATE
  PUBLIC :: band

  integer :: unit = 1

CONTAINS


  SUBROUTINE band( MBV_in, disp_switch )

    implicit none

    integer,intent(IN) :: MBV_in
    logical,intent(IN) :: disp_switch
    integer :: MBV,nktrj,i,j,k,s,n,ibz,ierr,iktrj,iter,Diter_band
    integer :: iktrj_0,iktrj_1,iktrj_2,iktrj_00,iktrj_tmp,ireq
    integer,allocatable :: ir_k(:),id_k(:)
    real(8) :: dak(3),sum0,sum1,max_err,max_err0
    real(8),allocatable :: ktrj(:,:),pxyz(:,:,:,:)
    real(8),allocatable :: kbb_tmp(:,:),esp_tmp(:,:,:),esp0_tmp(:,:,:)
    complex(8),allocatable :: unk0(:,:,:),unk00(:,:,:)
    logical :: disp_switch_parallel_bak,flag_end
    character(32) :: file_ovlp
    character(5) :: cc

    call read_band( myrank, unit )

    disp_switch_parallel_bak = disp_switch_parallel
!    disp_switch_parallel = .false.

    Diter_band = maxiter_band

    flag_end = .false.

    MBV=MBV_in
    if ( MBV < 1 .or. mb_band < MBV ) MBV=1
    if ( disp_switch_parallel ) write(*,*) "MBV=",MBV

! ---

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

! ---

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

! ---

    if ( iflag_hunk == 0 ) call deallocate_work_wf

    if ( MB < mb_band ) then
       call modify_mb
       call modify_arraysize
    end if

! ---

    allocate( esp0_tmp(MB,0:np_bzsm-1,MSP)   ) ; esp0_tmp=0.d0
    allocate( esp_tmp(MB,0:np_bzsm-1,MSP)    ) ; esp_tmp=0.d0
    allocate( kbb_tmp(3,0:np_bzsm-1)         ) ; kbb_tmp=0.d0
    allocate( pxyz(3,MB,0:np_bzsm-1,MSP)     ) ; pxyz=0.d0
    allocate( unk0(ML_0:ML_1,MB,MSP_0:MSP_1) ) ; unk0=0.d0

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

    call init_sweep( 2, mb2_band, esp_conv_tol )

    call Init_IO( "band" )

    MBZ_1 = MBZ_0

    loop_iktrj : do iktrj = iktrj_0, iktrj_2

       if ( iktrj <= nskip_band ) then
          if ( DISP_SWITCH_PARALLEL ) write(*,*) "Band ",iktrj," is skipped"
          cycle
       end if

       iktrj_00 = id_k(0) + iktrj - iktrj_0 + 1

       if ( iktrj <= nktrj ) then
          kbb(1:3,MBZ_0) = ktrj(1:3,iktrj)
       else
          kbb(1:3,MBZ_0) = 0.0d0
       end if
       call mpi_allgather &
            (kbb(1,MBZ_0),3,mpi_real8,kbb_tmp,3,mpi_real8,comm_bzsm,ierr)
       if ( myrank == 0 ) then
          write(*,*) "kbb_tmp"
          do ibz=0,np_bzsm-1
             iktrj_tmp = id_k(ibz)+iktrj-iktrj_0+1
             write(*,'(1x,i8,3f10.5)') iktrj_tmp,kbb_tmp(1:3,ibz)
          end do
       endif

       call init_kinetic(aa,bb,Nbzsm,kbb,disp_switch)

       select case( pselect )
       case( 2 )
          call prep_uvk_ps_nloc2(MBZ_0,MBZ_0,kbb(1,MBZ_0))
       case( 3 )
          call init_ps_nloc3
          call prep_ps_nloc3
       end select

       if ( iflag_hybrid > 0 ) then
          if ( disp_switch ) write(*,*) "iflag_hybrid=",iflag_hybrid
          call prep_kq_xc_hybrid(Nbzsm,MBZ_0,MBZ_1,kbb,bb,disp_switch)
          call init_fock_ffte
       end if

! --- sweep ---

       call calc_sweep( Diter_band, iswitch_gs, ierr, disp_switch )

! ---

       if ( ierr == -1 ) then
          if ( myrank == 0 ) write(*,*) "etime limit !!!"
          exit loop_iktrj
       end if
       if ( ierr == -2 ) then
          write(*,*) "band is not converged"
          return
       end if

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
                   unk00=0.0d0
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
!                if ( allocated(unk00) ) deallocate(unk00)
             end if
          end if
          unk0(:,:,s) = unk(:,:,MBZ_0,s)
       end do ! s
#endif

       select case( pselect )
       case( 2 )
          call prep_rvk_ps_nloc2(MBZ_0,MBZ_0,kbb(1,MBZ_0))
       case( 3 )
!          call prep_rvk_ps_nloc3(MBZ_0,MBZ_0,kbb(1,MBZ_0))
       end select

       pxyz(:,:,:,:)=0.d0
       do n=1,mb2_band
          do s=MSP_0,MSP_1
#ifndef _DRSDFT_
             call calc_expectval_momentum(MBZ_0,ML_0,ML_1,1,1,unk(ML_0,n,MBZ_0,s),pxyz(1,n,myrank_k,s))
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
          esp0_tmp(:,myrank_k,s)=esp(:,MBZ_0,s)
          call mpi_allgather(esp0_tmp(1,myrank_k,s),MB,mpi_real8,esp0_tmp(1,0,s),MB,mpi_real8,comm_bzsm,ierr)
       end do

       do ibz=0,np_bzsm-1
          iktrj_tmp = iktrj_00 + id_k(ibz)
          if ( iktrj_tmp > nktrj ) exit
          if ( myrank == 0 ) then
             write(unit_band_eigv,'(1x,2i6,3f20.12,i8)') &
                  iktrj_tmp,mb2_band,kbb_tmp(1:3,ibz),MBV
             write(unit_band_dedk,'(1x,2i6,3f20.12)') &
                  iktrj_tmp,mb2_band,kbb_tmp(1:3,ibz)
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

    if ( allocated(unk00) ) deallocate( unk00 )
    deallocate( unk0 )
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
    implicit none
    integer,intent(IN) :: mm,nn
    complex(8),intent(OUT) :: f(mm,nn)
    integer :: ierr,ireq,istatus(MPI_STATUS_SIZE)
    if ( myrank_k < np_bzsm-1 ) then
       call MPI_IRECV(f,mm*nn,MPI_COMPLEX16,myrank_k+1,1,comm_bzsm,ireq,ierr)
       call MPI_WAIT(ireq,istatus,ierr)
    end if
  END SUBROUTINE recv_wf_band


  SUBROUTINE calc_overlap(n1,n2,m1,m2,mm,f1,f0)
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


  SUBROUTINE modify_mb

    implicit none

    integer :: i,j

    if ( mb_band <= Nband ) return

    Nband = mb_band

    call init_scalapack( Nband )

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

    call init_subspace_diag( Nband )
    call reset_subspace_mate_sl_0
    call reset_subspace_rotv_sl_0

  END SUBROUTINE modify_mb


  SUBROUTINE modify_arraysize

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
