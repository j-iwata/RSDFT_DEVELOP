MODULE subspace_rotv_eigen_module

#ifdef _EIGEN_
  use eigen_libs
#endif
  use wf_module, only: unk,hunk,iflag_hunk &
                      ,MB_0_WF,MB_1_WF,MB_WF,ML_0_WF,ML_1_WF
  use parallel_module
  use subspace_diag_variables, only: Vsub

  implicit none

  PRIVATE
  PUBLIC :: subspace_rotv_eigen

CONTAINS

  SUBROUTINE subspace_rotv_eigen(k,s)
    implicit none
    integer,intent(IN) :: k,s
#ifdef _EIGEN_
    integer :: ns,ne,nn,ms,me,mm,i,j,ion,jon,il,jl,ierr
    integer :: x_nnod,x_inod,y_nnod,y_inod,inod,nnod
    integer :: MB,MB_0,MB_1,ML_0,ML_1,N_ML,MBLK,NBLK
    real(8),allocatable :: Cmn(:,:),utmp(:,:)
    real(8),parameter :: zero=0.0d0, one=1.0d0

    MB   = MB_WF
    MB_0 = MB_0_WF
    MB_1 = MB_1_WF

    ML_0 = ML_0_WF
    ML_1 = ML_1_WF
    N_ML = ML_1 - ML_0 + 1

    MBLK = MB
    NBLK = MB

    call eigen_get_procs( nnod, x_nnod, y_nnod )
    call eigen_get_id( inod, x_inod, y_inod )

    allocate( utmp(ML_0:ML_1,MB_0:MB_1) ) ; utmp=zero

    do ns=MB_0,MB_1,NBLK
       ne=min(ns+NBLK-1,MB_1)
       nn=ne-ns+1

       do ms=1,MB,MBLK
          me=min(ms+MBLK-1,MB)
          mm=me-ms+1

          allocate( Cmn(ms:me,ns:ne) ) ; Cmn=zero

          do j=ns,ne
             jon = eigen_owner_node( j, y_nnod, y_inod )
             do i=ms,me
                ion = eigen_owner_node( i, x_nnod, x_inod )
                if ( ion == x_inod .and. jon == y_inod ) then
                   il = eigen_translate_g2l( i, x_nnod, x_inod )
                   jl = eigen_translate_g2l( j, y_nnod, y_inod )
                   Cmn(i,j) = Vsub(il,jl)
                end if
             end do ! i
          end do ! j

          call MPI_ALLREDUCE( MPI_IN_PLACE, Cmn, size(Cmn), MPI_REAL8 &
                             ,MPI_SUM, comm_grid, ierr )

          call dgemm('N','N',N_ML,nn,mm,one,unk(ML_0,ms,k,s),N_ML &
                    ,Cmn,mm,one,utmp(ML_0,ns),N_ML)

          deallocate( Cmn )

       end do ! ms

    end do ! ns

    unk(:,MB_0:MB_1,k,s) = utmp(:,:)

    deallocate( utmp )
#endif
  END SUBROUTINE subspace_rotv_eigen

END MODULE subspace_rotv_eigen_module
