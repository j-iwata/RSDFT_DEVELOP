module rotorb_module

  use parallel_module
  use array_bound_module, only: MBZ_0,MBZ_1,MSP_0,MSP_1
  use wf_module, only: unk
  use cpmd_variables, only: gam,scr,wrk,tau,sig,gamn,psi_v,psi_n &
                           ,mstocck,dt,disp_switch,MBC,MBT,MB_0_CPMD,MB_1_CPMD
  use overlap_cpmd_module
  use watch_module
  use rsdft_mpi_module

  implicit none

  PRIVATE
  PUBLIC :: rotorb
  PUBLIC :: rotorb2

  integer,allocatable :: ir_i(:),id_i(:)
  integer :: ls,le,li,n1,n2,ML0,MB0
  real(8),allocatable :: psi_tmp(:,:)

CONTAINS

  subroutine rotorb
    implicit none
    integer,parameter :: maxit=100
    integer :: i,j,k,s,n,it
    real(8),parameter :: eps=1.0d-8
    real(8),parameter :: hf = 0.5d0
    real(8),parameter :: hm =-0.5d0
    real(8),parameter :: on = 1.0d0
    real(8),parameter :: zr = 0.0d0
    real(8) :: error,tmp1,tttt(2,16),ttmp(2)

    !call write_border( 1, "rotorb(start)" )
    !call watchb( ttmp, barrier="on" ); tttt=0.0d0

    if ( .not.allocated(id_i) ) then
       n1  = idisp(myrank)+1
       n2  = idisp(myrank)+ircnt(myrank)
       ML0 = ircnt(myrank)
       MB0 = MB_1_CPMD - MB_0_CPMD + 1
       allocate( id_i(0:nprocs-1) ) ; id_i=0
       allocate( ir_i(0:nprocs-1) ) ; ir_i=0
       ir_i(0:nprocs-1)=MBC/nprocs
       n=MBC-sum(ir_i)
       do i=1,n
          k=mod(i-1,nprocs)
          ir_i(k)=ir_i(k)+1
       end do
       do k=0,nprocs-1
          id_i(k)=sum(ir_i(0:k))-ir_i(k)
       end do
       if ( DISP_SWITCH ) then
          do k=0,nprocs-1
             write(*,*) k,id_i(k)+1,id_i(k)+ir_i(k)
          end do
       end if
       id_i(0:nprocs-1) = id_i(0:nprocs-1)*MBC
       ir_i(0:nprocs-1) = ir_i(0:nprocs-1)*MBC
       allocate( psi_tmp(n1:n2,MB_0_CPMD:MB_1_CPMD) ) ; psi_tmp=0.0d0
       ls = id_i(myrank)/MBC+1
       le = id_i(myrank)/MBC+ir_i(myrank)/MBC
       li = ir_i(myrank)/MBC
    end if

    !call watchb( ttmp, tttt(:,1), barrier="on" )

    do s=MSP_0,MSP_1
    do k=MBZ_0,MBZ_1

       !call watchb( ttmp, barrier="on" )

       MBT=mstocck(k,s)

       call overlap2(s,k) ! sig

       !call watchb( ttmp, tttt(:,2), barrier="on" )

       call overlap4(s,k) ! tau

       !call watchb( ttmp, tttt(:,3), barrier="on" )

!$OMP parallel workshare
       gam(1:MBC,ls:le) = sig(1:MBC,ls:le)*hf      ! (I-A)/2
!$OMP end parallel workshare

       !call watchb( ttmp, tttt(:,4), barrier="on" )

       if ( ls <= le ) then
          call DGEMM &
               ('N','N',MBT,li,MBT,hf,tau,MBC,tau(1,ls),MBC,hf,sig(1,ls),MBC)
       end if

       !call watchb( ttmp, tttt(:,5), barrier="on" )

       do it=1,maxit

          !call watchb( ttmp, barrier="on" )

          do j=ls,le
          do i=1,MBC
             gamn(i,j) = sig(i,j)
          end do
          end do
          do j=ls,le
          do i=1,MBC
             wrk(i,j) = tau(i,j) - gam(i,j)
          end do
          end do

          !call watchb( ttmp, tttt(:,6), barrier="on" )

          call rsdft_allgatherv( wrk(:,ls:le), wrk, ir_i, id_i, MPI_COMM_WORLD )

          !call watchb( ttmp, tttt(:,7), barrier="on" )

          if ( ls <= le ) then
             call DGEMM &
               ('N','N',MBT,li,MBT,hm,wrk,MBC,wrk(1,ls),MBC,on,gamn(1,ls),MBC)
          end if

          !call watchb( ttmp, tttt(:,8), barrier="on" )

          error=0.0d0
          do j=ls,le
          do i=j,MBT
             error = error + ( gamn(i,j)-gam(i,j) )**2
          end do
          end do

          call rsdft_allreduce_sum( error, MPI_COMM_WORLD )
          error=sqrt(error)/MBT
          if ( myrank == 0 .and. error > on ) write(*,*) k,s,error

          !call watchb( ttmp, tttt(:,9), barrier="on" )

          do j=ls,le
          do i=1,MBC
             gam(i,j) = gamn(i,j)
          end do
          end do

          !call watchb( ttmp, tttt(:,10), barrier="on" )

          if ( error <= eps ) exit
          if ( it == maxit .and. myrank == 0 ) then
             write(*,*) "WARNING: iteration was not converged in rotorb"
          end if

       end do ! it

       !call watchb( ttmp, barrier="on" )

       call rsdft_allgatherv( gam(:,ls:le), wrk, ir_i, id_i, MPI_COMM_WORLD )

       !call watchb( ttmp, tttt(:,11), barrier="on" )

       call DGEMM('N','N',ML0,MB0,MBT,on,unk(n1,1,k,s),ML0 &
            ,wrk(1,MB_0_CPMD),MBC,zr,psi_tmp(n1,MB_0_CPMD),ML0)

       !call watchb( ttmp, tttt(:,12), barrier="on" )

!$OMP parallel workshare
       psi_n(n1:n2,MB_0_CPMD:MB_1_CPMD,k,s) = &
       psi_n(n1:n2,MB_0_CPMD:MB_1_CPMD,k,s) + psi_tmp(n1:n2,MB_0_CPMD:MB_1_CPMD)
!$OMP end parallel workshare

       !call watchb( ttmp, tttt(:,13), barrier="on" )

       tmp1=1.0d0/dt
!$OMP parallel workshare
       psi_v(n1:n2,MB_0_CPMD:MB_1_CPMD,k,s) = &
       psi_v(n1:n2,MB_0_CPMD:MB_1_CPMD,k,s) + tmp1*psi_tmp(n1:n2,MB_0_CPMD:MB_1_CPMD)
!$OMP end parallel workshare

       !call watchb( ttmp, tttt(:,14), barrier="on" )

    end do ! k
    end do ! s

    !call watchb( ttmp, barrier="on" )

!$OMP parallel workshare
    unk(:,MB_0_CPMD:MB_1_CPMD,:,:) = psi_n(:,MB_0_CPMD:MB_1_CPMD,:,:)
!$OMP end parallel workshare

    !call watchb( ttmp, tttt(:,15), barrier="on" )

!    if ( myrank == 0 ) then
!       do i=1,15
!          write(*,'("time_rotorb(",i2.2,")",2f10.5)') i,tttt(1:2,i)
!       end do
!    end if

    !call write_border( 1, "rotorb(end)" )

  end subroutine rotorb


!-----------------------------------------------------------------------
!     orthogonalization of wave function velocity
!----------------------------------------------------------------------- 
  subroutine rotorb2

    implicit none
    integer :: k,s,i
    real(8),parameter :: on= 1.0d0
    real(8),parameter :: om=-1.0d0
    real(8) :: ttmp(2),tttt(2,2)

    !call write_border( 1, "rotorb2(start)" )
    !call watchb( ttmp, barrier="on" ); tttt=0.0d0

    do s=MSP_0,MSP_1
    do k=MBZ_0,MBZ_1

       MBT = mstocck(k,s)

       !call watchb( ttmp, barrier="on" )

       call overlap5(s,k)

       !call watchb( ttmp, tttt(:,1), barrier="on" )

       call DGEMM('N','N',ML0,MB0,MBT,om,unk(n1,1,k,s),ML0 &
            ,wrk(1,MB_0_CPMD),MBC,on,psi_v(n1,MB_0_CPMD,k,s),ML0)

       !call watchb( ttmp, tttt(:,2), barrier="on" )

    end do ! k
    end do ! s

!    if ( myrank == 0 ) then
!       do i=1,2
!          write(*,'("time_rotorb2(",i1,")",2f10.5)') i,tttt(:,i)
!       end do
!    end if

    !call write_border( 1, "rotorb2(end)" )

  end subroutine rotorb2


end module rotorb_module
