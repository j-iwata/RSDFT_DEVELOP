MODULE band_sseig_module

  use watch_module
  use parallel_module
  use prepare_sseig_module
  use apply_sseig_module, only: apply_sseig,total_numeig
  use sseig
  use rgrid_module
  use array_bound_module
  use bz_module
  use aa_module, only: aa
  use bb_module, only: bb
  use kinetic_module, only: init_kinetic
  use ps_nloc2_module, only: prep_uvk_ps_nloc2,prep_rvk_ps_nloc2
  use momentum_module
  use band_variables, only: nfki,nbk,ak,nskip_band &
       ,unit_band_eigv,unit_band_dedk,unit_band_ovlp,read_band

  implicit none

  PRIVATE
  PUBLIC :: band_sseig

  integer,parameter :: unit=1

CONTAINS


  SUBROUTINE band_sseig(disp_switch)
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

    call read_band( myrank, unit )

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


  SUBROUTINE calc_overlap_2(n1,n2,mm1,mm0,f1,f0)
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

    if ( myrank == 0 ) then
       write(unit_band_ovlp) mm1,mm0
       write(unit_band_ovlp) maxloc10(:,:)
       write(unit_band_ovlp) sqovlp10(:,:)
    end if

    deallocate( sqovlp10 )
    deallocate( maxloc10 )

  END SUBROUTINE calc_overlap_2


END MODULE band_sseig_module
