MODULE subspace_mate_sl_0_module

  use rgrid_module, only: dV,zdV
  use parallel_module
  use hamiltonian_module
  use wf_module, only: unk,esp
  use scalapack_module
  use subspace_diag_variables
  use array_bound_module, only: ML_0,ML_1,MB_0,MB_1
  use watch_module

  implicit none

  PRIVATE
  PUBLIC :: subspace_mate_sl_0, reset_subspace_mate_sl_0

#ifdef _DRSDFT_
  real(8),allocatable :: hunk(:,:),vtmp2(:,:),wtmp2(:,:)
  character(1),parameter :: TRANSA='T', TRANSB='N'
#else
  complex(8),allocatable :: hunk(:,:),vtmp2(:,:),wtmp2(:,:)
  character(1),parameter :: TRANSA='C', TRANSB='N'
#endif

  integer,allocatable :: mat_indx(:,:),nblock(:)

CONTAINS

  SUBROUTINE reset_subspace_mate_sl_0
    if ( allocated(mat_indx) ) deallocate( mat_indx )
  END SUBROUTINE reset_subspace_mate_sl_0

  SUBROUTINE subspace_mate_sl_0(k,s)
    implicit none
    integer,intent(IN) :: k,s
    integer :: i,i0,i1,i2,ib1,ib2,j,j0,j1,j2,m,me,mm,ms,MB
    integer :: MBLK,MBLKH,ML0,n1,n2,mms,mme,nnn,nns,n,ne,nn,ns
    integer :: istatus(MPI_STATUS_SIZE,123),nreq,ireq(123),itag,ierr
    integer :: IPROW,IPCOL,iroot1,iroot2,mrnk,nrecv_me,nsend_me
    complex(8) :: ztmp
    real(8) :: ct0,ct1,et0,et1,ctt(2),ett(2)

    MB = MB_diag

    if ( .not.allocated(mat_indx) ) then
       if ( .not.allocated(nblock) ) allocate( nblock(0:nprocs-1) )
       nblock=0
       m=-1
       do i=1,MB,MBSIZE
          m=m+1
          mm=mod(m+NPROW,NPROW)
          n=-1
          do j=1,MB,MBSIZE
             n=n+1
             nn=mod(n+NPCOL,NPCOL)
             nblock(usermap(mm,nn,1))=nblock(usermap(mm,nn,1))+1 
          end do
       end do
       allocate( mat_indx(8,maxval(nblock)) ) ; mat_indx(:,:)=0
       nblock(:)=0
       m=-1
       do i0=1,MB,MBSIZE
          i1=min(i0+MBSIZE-1,MB)
          m=m+1
          mm=mod(m+NPROW,NPROW)
          n=-1
          do j0=1,MB,MBSIZE
             j1=min(j0+MBSIZE-1,MB)
             n=n+1
             nn=mod(n+NPCOL,NPCOL)
             nblock(usermap(mm,nn,1))=nblock(usermap(mm,nn,1))+1
             if ( usermap(mm,nn,1) == myrank ) then
                mat_indx(1,nblock(myrank))=i0
                mat_indx(2,nblock(myrank))=i1
                mat_indx(3,nblock(myrank))=j0
                mat_indx(4,nblock(myrank))=j1
             end if
          end do
       end do
       mat_indx(5,1)=1
       mat_indx(6,1)=mat_indx(5,1)+mat_indx(2,1)-mat_indx(1,1)
       mat_indx(7,1)=1
       mat_indx(8,1)=mat_indx(7,1)+mat_indx(4,1)-mat_indx(3,1)
       do i=2,nblock(myrank)
          if ( mat_indx(1,i) == mat_indx(1,i-1) ) then
             mat_indx(5,i) = mat_indx(5,i-1)
             mat_indx(6,i) = mat_indx(6,i-1)
             mat_indx(7,i) = mat_indx(8,i-1)+1
             mat_indx(8,i) = mat_indx(7,i)+mat_indx(4,i)-mat_indx(3,i)
          else
             mat_indx(5,i) = mat_indx(6,i-1)+1
             mat_indx(6,i) = mat_indx(5,i)+mat_indx(2,i)-mat_indx(1,i)
             mat_indx(7,i) = 1
             mat_indx(8,i)=mat_indx(7,i)+mat_indx(4,i)-mat_indx(3,i)
          end if
       end do
    end if

    UPLO = 'L'

    n1    = ML_0
    n2    = ML_1
    ML0   = ML_1-ML_0+1
    ctt(:)= 0.d0
    ett(:)= 0.d0

    NBLK1 = 4

!    MBLK  = min( MBSIZE,NBSIZE )
    MBLK = (MB+1)/2

    allocate( hunk(n1:n2,MBLK) )
    hunk=zero

    do ns=1,MB,MBLK

       ne=min(ns+MBLK-1,MB)
       nn=ne-ns+1

       call watch(ct0,et0)
       do ib1=ns,ne,MB_d
          ib2=min(ib1+MB_d-1,ne)
          call hamiltonian(k,s,unk(n1,ib1,k,s),hunk(n1,ib1-ns+1),n1,n2,ib1,ib2)
       end do
       call watch(ct1,et1) ; ctt(1)=ctt(1)+ct1-ct0 ; ett(1)=ett(1)+et1-et0

       do ms=1,MB,MBLK

          me=min(ms+MBLK-1,MB)
          mm=me-ms+1

          if ( ns > me ) cycle

          allocate( vtmp2(ms:me,ns:ne) ) ; vtmp2=zero
          allocate( wtmp2(ms:me,ns:ne) ) ; wtmp2=zero

          if ( ne <= ms ) then
#ifdef _DRSDFT_
             call dgemm('T','N',mm,nn,ML0, dV,unk(n1,ms,k,s),ML0,hunk(n1,1),ML0,zero,vtmp2,mm)
#else
             call zgemm('C','N',mm,nn,ML0,zdV,unk(n1,ms,k,s),ML0,hunk(n1,1),ML0,zero,vtmp2,mm)
#endif
          else
             call mate_sub(ML0,mm,MBLK,unk(n1,ms,k,s),hunk(n1,1),vtmp2) !"mm==nn" is assumed
          end if

          call watch(ct0,et0)
!!!!!!!!! If a single conetxt calls pd???, mpi_allreduce may be able to be replaced by mpi_reduce
          call mpi_allreduce(vtmp2,wtmp2,mm*nn,TYPE_MAIN,mpi_sum,comm_grid,ierr)
          call watch(ct1,et1) ; ctt(2)=ctt(2)+ct1-ct0 ; ett(2)=ett(2)+et1-et0

          do n=1,nblock(myrank)
             if ( ( ( ms <= mat_indx(1,n) .and. mat_indx(1,n) <= me ) .or. &
                    ( ms <= mat_indx(2,n) .and. mat_indx(2,n) <= me ) ) .and. &
                  ( ( ns <= mat_indx(3,n) .and. mat_indx(3,n) <= ne ) .or. &
                    ( ns <= mat_indx(4,n) .and. mat_indx(4,n) <= ne ) ) ) then
                i0=mat_indx(5,n)
                j0=mat_indx(7,n)
                do j=ns,ne
                   if ( j < mat_indx(3,n) .or. mat_indx(4,n) < j ) cycle
                do i=ms,me
                   if ( i < mat_indx(1,n) .or. mat_indx(2,n) < i ) cycle
                   Hsub(i0+i-mat_indx(1,n),j0+j-mat_indx(3,n))=wtmp2(i,j)
                end do
                end do
             end if
          end do ! n

          deallocate( wtmp2 )
          deallocate( vtmp2 )

       end do ! ms

    end do ! ns

    deallocate( hunk )

    if ( disp_switch_parallel ) then
       write(*,*) "time_hamil(mate)=",ctt(1),ett(1)
       write(*,*) "time_allrd(mate)=",ctt(2),ett(2)
    end if

  END SUBROUTINE subspace_mate_sl_0


  RECURSIVE SUBROUTINE mate_sub(m,n,nblk,a,b,c)
    implicit none
    integer,intent(IN)  :: m,n,nblk
#ifdef _DRSDFT_
    real(8),intent(IN)  :: a(m,n),b(m,n)
    real(8),intent(OUT) :: c(n,n)
    real(8),allocatable :: ctmp(:,:)
#else
    complex(8),intent(IN)  :: a(m,n),b(m,n)
    complex(8),intent(OUT) :: c(n,n)
    complex(8),allocatable :: ctmp(:,:)
#endif
    integer :: i,j,i0,i1,j0,j1,ni,nj,nblkh

    do i0=1,n,nblk
       i1=min(i0+nblk-1,n)
       ni=i1-i0+1

       do j0=1,n,nblk
          j1=min(j0+nblk-1,n)
          nj=j1-j0+1

          if ( j0 > i1 ) then
             cycle
          else if ( j1 <= i0 ) then
#ifdef _DRSDFT_
             call dgemm('T','N',ni,nj,m, dV,a(1,i0),m,b(1,j0),m,zero,c(i0,j0),n)
#else
             call zgemm('C','N',ni,nj,m,zdV,a(1,i0),m,b(1,j0),m,zero,c(i0,j0),n)
#endif
          else
             if ( ni > nblk1 ) then
                allocate( ctmp(ni,ni) ) ; ctmp=0.d0
                nblkh=nblk/2
                call mate_sub(m,ni,nblkh,a(1,i0),b(1,j0),ctmp)
                c(i0:i1,j0:j1)=ctmp(:,:)
                deallocate( ctmp )
             else
                do i=i0,i1
                do j=j0,i
#ifdef _DRSDFT_
                   c(i,j)=sum( a(:,i)*b(:,j) )*dV
#else
                   c(i,j)=sum( conjg(a(:,i))*b(:,j) )*dV
#endif
                end do
                end do
             end if
          end if

       end do ! j0

    end do ! i0

  END SUBROUTINE mate_sub

END MODULE subspace_mate_sl_0_module
