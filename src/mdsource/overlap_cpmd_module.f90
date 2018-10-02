module overlap_cpmd_module

  use wf_module, only: unk
  use rgrid_module, only: dV
  use parallel_module
  use cpmd_variables, only: wrk,tau,sig,MBC,MBT,psi_n,psi_v,ir_band_cpmd,id_band_cpmd
  use watch_module
  use calc_overlap_bp_module
  use rsdft_mpi_module

  implicit none

  PRIVATE
  PUBLIC :: overlap2, overlap4, overlap5

CONTAINS


  subroutine overlap2(s,k)
    implicit none
    integer,intent(in) :: s,k
    integer :: i,m1,m2
    real(8) :: ttmp(2),tttt(2,3)

    !call write_border( 1, "overlap2(start)" )
    !call watchb( ttmp, barrier="on" ); tttt=0.0d0

    m1 = id_band_cpmd(myrank_b)+1
    m2 = id_band_cpmd(myrank_b)+ir_band_cpmd(myrank_b)

    !call watchb( ttmp, tttt(:,1), barrier="on" )

    call calc_overlap_bp( MBT, psi_n(:,m1:m2,k,s), psi_n(:,m1:m2,k,s), -dV, sig )

    !call watchb( ttmp, tttt(:,2), barrier="on" )

!$OMP parallel do
    do i=1,MBC
       sig(i,i) = sig(i,i) + 1.0d0
    end do
!$OMP end parallel do

    !call watchb( ttmp, tttt(:,3), barrier="on" )

!    if ( myrank == 0 ) then
!       do i=1,3
!          write(*,'(2x,"time_overlap2(",i1,")",2f10.5)') i,tttt(:,i)
!       end do
!    end if

    !call write_border( 1, "overlap2(end)" )

  end subroutine overlap2


  subroutine overlap4(s,k)
    implicit none
    integer,intent(in) :: s,k
    integer :: i,m1,m2
    real(8) :: ttmp(2),tttt(2,3)

    !call write_border( 1, "overlap4(start)" )
    !call watchb( ttmp, barrier="on" ); tttt=0.0d0

    m1 = id_band_cpmd(myrank_b)+1
    m2 = id_band_cpmd(myrank_b)+ir_band_cpmd(myrank_b)

    !call watchb( ttmp, tttt(:,1), barrier="on" )

    call calc_overlap_bp( MBT, unk(:,m1:m2,k,s), psi_n(:,m1:m2,k,s), -dV, tau )

    !call watchb( ttmp, tttt(:,2), barrier="on" )

!$OMP parallel do
    do i=1,MBC
       tau(i,i) = tau(i,i) + 1.0d0
    end do
!$OMP end parallel do

    !call watchb( ttmp, tttt(:,3), barrier="on" )

!    if ( myrank == 0 ) then
!       do i=1,3
!          write(*,'(2x,"time_overlap4(",i1,")",2f10.5)') i,tttt(:,i)
!       end do
!    end if

    !call write_border( 1, "overlap4(end)" )

  end subroutine overlap4


  subroutine overlap5(s,k)
    implicit none
    integer,intent(IN) :: s,k
    integer :: i,m1,m2
    real(8) :: ttmp(2),tttt(2,2)

    !call write_border( 1, "overlap5(start)" )
    !call watchb( ttmp, barrier="on" ); tttt=0.0d0

    m1    = id_band_cpmd(myrank_b)+1
    m2    = id_band_cpmd(myrank_b)+ir_band_cpmd(myrank_b)

    !call watchb( ttmp, tttt(:,1), barrier="on" )

    call calc_overlap_bp( MBT, unk(:,m1:m2,k,s), psi_v(:,m1:m2,k,s), dV, wrk )

    !call watchb( ttmp, tttt(:,2), barrier="on" )

!    if ( myrank == 0 ) then
!       do i=1,2
!          write(*,'(2x,"time_overlap5(",i1,")",2f10.5)') i, tttt(:,i)
!       end do
!    end if

    !call write_border( 1, "overlap5(end)" )

  end subroutine overlap5


end module overlap_cpmd_module
