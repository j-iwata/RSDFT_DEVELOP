MODULE WFtest
  use InnerProduct
  use wf_module
  use parallel_module, only: COMM_GRID
  use RealComplex, only: TYPE_MAIN
  implicit none

  include 'mpif.h'

  PRIVATE
  PUBLIC :: test_orthnorm_wf

CONTAINS
  SUBROUTINE test_orthnorm_wf(rank)
    implicit none
    integer,intent(IN) :: rank
    integer :: ns,nk,nb1,nb2,ierr
#ifdef _DRSDFT_
    real(8) :: c,d
#else
    complex(8) :: c,d
#endif

    do ns=MS_0_WF,MS_1_WF
      do nk=MK_0_WF,MK_1_WF
        do nb1=MB_0_WF,MB_1_WF
          do nb2=MB_0_WF,nb1
            call get_gSf(unk(ML_0_WF,nb1,nk,ns),unk(ML_0_WF,nb2,nk,ns),ML_0_WF,ML_1_WF,nk,d,0)
            call MPI_ALLREDUCE(d,c,1,TYPE_MAIN,MPI_SUM,COMM_GRID,ierr)
            write(330+rank,'(4I5,4g20.7)') ns,nk,nb1,nb2,d,c
          end do
        end do
      end do
    end do
    return
  END SUBROUTINE test_orthnorm_wf
END MODULE WFtest
