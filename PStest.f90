MODULE PStest
  use VarPSMember
  use VarPSMemberG
  use VarParaPSnonLocG
  use ps_nloc2_variables
  use array_bound_module
  implicit none
  private
  public :: check_all_ps

CONTAINS
  SUBROUTINE check_all_ps(myrank)
    implicit none
    integer,intent(IN) :: myrank
    integer :: i,j,k,s,m,n,l
    integer :: ik,lma,inlop
    write(900,*) '----------------------- nloc'
    write(900,'(A11,I6)') 'Mlma     = ',Mlma
    write(900,'(A11,I6)') 'nzlma    = ',nzlma
    write(900,'(A11,I6)') 'MAXMJJ   = ',MAXMJJ
    write(900,*) '----------------------- qr'
    write(900,'(A11,I6)') 'Mqr      = ',Mqr
    write(900,'(A11,I6)') 'MAXMJJ_Q = ',MAXMJJ_Q

    do ik=MBZ_0,MBZ_1
      do lma=1,nzlma
        do j=1,MJJ_MAP(lma)
          write(910,*) '----- uVk'
          write(901,'(3I5,2g20.7)') ik,lma,j,uVk(j,lma,ik)
        end do
      end do
    end do

    do inlop=1,N_nlop
      write(902,'(I5,g20.7)') inlop,qij(inlop)
    end do
  END SUBROUTINE check_all_ps
END MODULE PStest
