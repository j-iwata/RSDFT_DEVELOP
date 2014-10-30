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
    write(9000+myrank,*) '----------------------- nloc'
    write(9000+myrank,'(A11,I6)') 'Mlma     = ',Mlma
    write(9000+myrank,'(A11,I6)') 'nzlma    = ',nzlma
    write(9000+myrank,'(A11,I6)') 'MAXMJJ   = ',MAXMJJ
    write(9000+myrank,'(A11,I6)') 'MAXMJJ_  = ',maxval(MJJ(1:nzlma))
#ifdef _USPP_
    write(9000+myrank,*) '----------------------- qr'
    write(9000+myrank,'(A11,I6)') 'Mqr      = ',Mqr
    write(9000+myrank,'(A11,I6)') 'MAXMJJ_Q = ',MAXMJJ_Q
#endif

    write(9010+myrank,*) '----- uVk'
    do ik=MBZ_0,MBZ_1
      do lma=1,nzlma
        do j=1,MJJ(lma)
          write(9010+myrank,'(3I5,2g20.7)') ik,lma,j,uVk(j,lma,ik)
        end do
      end do
    end do

#ifdef _USPP_
    do inlop=1,N_nlop
      write(9020+myrank,'(I5,g20.7)') inlop,qij(inlop)
    end do
#endif
  END SUBROUTINE check_all_ps
END MODULE PStest
