MODULE subspace_mate_eigen_module

  use rgrid_module, only: dV,zdV
  use parallel_module
  use hamiltonian_module
  use wf_module, only: unk,esp,hunk,iflag_hunk &
                      ,ML_0_WF,ML_1_WF,MB_0_WF,MB_1_WF, MB_WF
  use subspace_diag_variables, only: Hsub
#ifdef _EIGEN_
  use eigen_libs
#endif

  implicit none

  PRIVATE
  PUBLIC :: subspace_mate_eigen

#ifdef _DRSDFT_
  real(8),parameter :: zero=0.0d0
  real(8),allocatable :: hpsi(:,:),vtmp2(:,:),wtmp2(:,:)
  character(1),parameter :: TRANSA='T', TRANSB='N'
  integer,parameter :: TYPE_MAIN=MPI_REAL8
#else
  complex(8),parameter :: zero=(0.0d0,0.0d0)
  complex(8),allocatable :: hpsi(:,:),vtmp2(:,:),wtmp2(:,:)
  character(1),parameter :: TRANSA='C', TRANSB='N'
  integer,parameter :: TYPE_MAIN=MPI_COMPLEX16
#endif

  integer :: NBLK1

CONTAINS


  SUBROUTINE subspace_mate_eigen(k,s)

    implicit none
    integer,intent(IN) :: k,s
#ifdef _EIGEN_
    integer :: i,ib1,ib2,j,m,me,mm,ms,MB,ierr,MB_0,MB_1
    integer :: MBLK,MBLKH,ML0,n1,n2,mms,mme,nnn,nns,n,ne,nn,ns
    integer,allocatable :: ir(:),id(:)
    integer :: nnod,inod,x_nnod,x_inod,y_nnod,y_inod,LDR,LDC
    integer :: i_1,j_1,ion,jon,mrnk

    n1    = ML_0_WF
    n2    = ML_1_WF
    ML0   = ML_1_WF-ML_0_WF+1
    mrnk  = id_class(myrank,4)
    MB_0  = MB_0_WF
    MB_1  = MB_1_WF
    MB    = MB_WF
    MBLK  = MB
    NBLK1 = 4

! ---

    call eigen_get_procs( nnod, x_nnod, y_nnod )
    call eigen_get_id( inod, x_inod, y_inod )

! ---

    allocate( id(0:np_band-1),ir(0:np_band-1) )

    id(0:np_band-1) = id_band(0:np_band-1)*ML0
    ir(0:np_band-1) = ir_band(0:np_band-1)*ML0

    call mpi_allgatherv(unk(n1,MB_0,k,s),ir(mrnk),TYPE_MAIN &
                       ,unk(n1,1,k,s),ir,id,TYPE_MAIN,comm_band,ierr)
    if ( iflag_hunk >= 1 ) then
       call mpi_allgatherv(hunk(n1,MB_0,k,s),ir(mrnk),TYPE_MAIN &
                          ,hunk(n1,1,k,s),ir,id,TYPE_MAIN,comm_band,ierr)
    end if

    deallocate( ir,id )

! ---

    allocate( hpsi(n1:n2,MBLK) ) ; hpsi=zero

    do ns=MB_0,MB_1,MBLK

       ne=min(ns+MBLK-1,MB_1)
       nn=ne-ns+1

       if ( iflag_hunk >= 1 ) then

          do ib1=ns,ne
             hpsi(:,ib1-ns+1)=hunk(:,ib1,k,s)
          end do

       else

          do ib1=ns,ne,MB_d
             ib2=min(ib1+MB_d-1,ne)
             call hamiltonian &
                  (k,s,unk(n1,ib1,k,s),hpsi(n1,ib1-ns+1),n1,n2,ib1,ib2)
          end do

       end if

       do ms=ns,ne,MBLK

          me=min(ms+MBLK-1,ne)
          mm=me-ms+1

          allocate( vtmp2(ms:me,ns:ne) ) ; vtmp2=zero
          allocate( wtmp2(ms:me,ns:ne) ) ; wtmp2=zero

          MBLKH = max( MBLK/2, NBLK1 )
          call mate_sub(ms,me,ns,ne,MBLKH,ns,ms,me,k,s)

          call mpi_allreduce &
               (vtmp2,wtmp2,mm*nn,TYPE_MAIN,mpi_sum,comm_grid,ierr)

          do j=ns,ne
          do i=j,me

             ion = eigen_owner_node( i, x_nnod, x_inod )
             jon = eigen_owner_node( j, y_nnod, y_inod )
             if ( ion == x_inod .and. jon == y_inod ) then
                i_1 = eigen_translate_g2l( i, x_nnod, x_inod )
                j_1 = eigen_translate_g2l( j, y_nnod, y_inod )
                Hsub(i_1,j_1) = wtmp2(i,j)
             end if

             if ( i == j ) cycle

             ion = eigen_owner_node( j, x_nnod, x_inod )
             jon = eigen_owner_node( i, y_nnod, y_inod )
             if ( ion == x_inod .and. jon == y_inod ) then
                i_1 = eigen_translate_g2l( j, x_nnod, x_inod )
                j_1 = eigen_translate_g2l( i, y_nnod, y_inod )
                Hsub(i_1,j_1) = wtmp2(i,j)
             end if

          end do ! i
          end do ! j

          deallocate( wtmp2 )
          deallocate( vtmp2 )

       end do ! ms

       do ms=ne+1,MB_1,MBLK

          me=min(ms+MBLK-1,MB_1)
          mm=me-ms+1

          allocate( vtmp2(ms:me,ns:ne) ) ; vtmp2=zero
          allocate( wtmp2(ms:me,ns:ne) ) ; wtmp2=zero

#ifdef _DRSDFT_
          call dgemm(TRANSA,TRANSB,mm,nn,ML0, dV,unk(n1,ms,k,s) &
               ,ML0,hpsi(n1,1),ML0,zero,vtmp2(ms,ns),mm)
#else
          call zgemm(TRANSA,TRANSB,mm,nn,ML0,zdV,unk(n1,ms,k,s) &
               ,ML0,hpsi(n1,1),ML0,zero,vtmp2(ms,ns),mm)
#endif

          call mpi_allreduce &
               (vtmp2,wtmp2,mm*nn,TYPE_MAIN,mpi_sum,comm_grid,ierr)

          do j=ns,ne
          do i=ms,me

             ion = eigen_owner_node( i, x_nnod, x_inod )
             jon = eigen_owner_node( j, y_nnod, y_inod )
             if ( ion == x_inod .and. jon == y_inod ) then
                i_1 = eigen_translate_g2l( i, x_nnod, x_inod )
                j_1 = eigen_translate_g2l( j, y_nnod, y_inod )
                Hsub(i_1,j_1) = wtmp2(i,j)
             end if

             if ( i == j ) cycle

             ion = eigen_owner_node( j, x_nnod, x_inod )
             jon = eigen_owner_node( i, y_nnod, y_inod )
             if ( ion == x_inod .and. jon == y_inod ) then
                i_1 = eigen_translate_g2l( j, x_nnod, x_inod )
                j_1 = eigen_translate_g2l( i, y_nnod, y_inod )
                Hsub(i_1,j_1) = wtmp2(i,j)
             end if

          end do ! i
          end do ! j

          deallocate( wtmp2 )
          deallocate( vtmp2 )

       end do ! ms

    end do ! ns

    deallocate( hpsi )
#endif
  END SUBROUTINE subspace_mate_eigen


  RECURSIVE SUBROUTINE mate_sub(mm1,mm2,nn1,nn2,MBLK,ns0,ld0,ld1,k,s)

    implicit none
    integer,intent(IN) :: mm1,mm2,nn1,nn2,MBLK,ns0,ld0,ld1,k,s
    integer :: n1,n2,ML0,n,ns,ne,nn,m,ms,me,mm,mms,MBLKH,i,ld

    n1  = ML_0_WF
    n2  = ML_1_WF
    ML0 = n2-n1+1
    ld  = ld1-ld0+1

    do ns=nn1,nn2,MBLK
       ne=min(ns+MBLK-1,nn2)
       nn=ne-ns+1

       if ( nn <= 0 ) cycle

       do mms=mm1,mm2,MBLK

          ms=max(ns,mms)
          me=min(mms+MBLK-1,mm2)
          mm=me-ms+1

          if ( mm <= 0 ) cycle

          if ( ms >= ne ) then
#ifdef _DRSDFT_
             call dgemm(TRANSA,TRANSB,mm,nn,ML0, dV,unk(n1,ms,k,s) &
                  ,ML0,hpsi(n1,ns-ns0+1),ML0,zero,vtmp2(ms,ns),ld)
#else
             call zgemm(TRANSA,TRANSB,mm,nn,ML0,zdV,unk(n1,ms,k,s) &
                  ,ML0,hpsi(n1,ns-ns0+1),ML0,zero,vtmp2(ms,ns),ld)
#endif
          else if ( mm <= NBLK1 ) then
             do n=ns,ne
#ifdef _DRSDFT_
                call dgemv(TRANSA,ML0,ne-n+1, dV,unk(n1,n,k,s) &
                     ,ML0,hpsi(n1,n-ns0+1),1,zero,vtmp2(n,n),1)
#else
                call zgemv(TRANSA,ML0,ne-n+1,zdV,unk(n1,n,k,s) &
                     ,ML0,hpsi(n1,n-ns0+1),1,zero,vtmp2(n,n),1)
#endif
             end do

          else

             MBLKH = max( MBLK/2, NBLK1 )
             call mate_sub(ms,me,ns,ne,MBLKH,ns0,ld0,ld1,k,s)

          end if

       end do ! mms

    end do ! ns

    return
  END SUBROUTINE mate_sub


END MODULE subspace_mate_eigen_module
