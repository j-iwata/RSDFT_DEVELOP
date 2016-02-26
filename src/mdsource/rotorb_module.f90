MODULE rotorb_module

  use parallel_module
  use array_bound_module, only: MBZ_0,MBZ_1,MSP_0,MSP_1
  use wf_module, only: unk
  use cpmd_variables, only: gam,scr,wrk,tau,sig,gamn,psi_v,psi_n &
                           ,mstocck,dt,disp_switch,MBC,MBT,MB_0_CPMD,MB_1_CPMD
  use overlap_cpmd_module
  use watch_module

  implicit none

  PRIVATE
  PUBLIC :: rotorb

  integer,allocatable :: ir_i(:),id_i(:)
  real(8),allocatable :: psi_tmp(:,:)

CONTAINS

  SUBROUTINE rotorb
    implicit none
    integer,parameter :: maxit=100
    integer i,j,k,l,m,n,s,ii,it,ij
    integer n1,n2,ML0,ierr,MB0
    integer ls,le,li
    real(8),parameter :: eps=1.0d-8
    real(8),parameter :: hf = 0.5d0
    real(8),parameter :: hm =-0.5d0
    real(8),parameter :: on = 1.0d0
    real(8),parameter :: zr = 0.0d0
    real(8) :: ctime_force(0:10),etime_force(0:10)
    real(8) :: error,error0,tmp1,ct0,ct1,et0,et1
    complex(8),parameter :: zo = (1.0d0,0.0d0)
    complex(8),parameter :: zz = (0.0d0,0.0d0)
    complex(8),allocatable :: zrk(:,:)

    !call watcht(myrank==0,"",0)

    if ( MBC < nprocs ) then
       write(*,*) "MBC<nprocs!!!",MBC,nprocs,myrank
       stop
    end if

    n1    = idisp(myrank)+1
    n2    = idisp(myrank)+ircnt(myrank)
    ML0   = ircnt(myrank)
    MB0   = MB_1_CPMD - MB_0_CPMD + 1

    if ( .not.allocated(id_i) ) then
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
    end if

    ls = id_i(myrank)/MBC+1
    le = id_i(myrank)/MBC+ir_i(myrank)/MBC
    li = ir_i(myrank)/MBC

    !call watcht(myrank==0,"rotorb(00)",1)

    do s=MSP_0,MSP_1
    do k=MBZ_0,MBZ_1

       !call watcht(myrank==0,"",0)

       MBT=mstocck(k,s)

       !call watcht(myrank==0,"rotorb(0)",1)
       !call watch(ct0,et0)

       call overlap4(s,k) ! tau

       !call watch(ct1,et1) ; if ( myrank==0 ) write(*,*) "timet(rotorb(1)=",ct1-ct0,et1-et0

       call overlap2(s,k) ! sig

       !call watch(ct0,et0) ; if ( myrank == 0 ) write(*,*) "timet(rotorb(2)=",ct0-ct1,et0-et1
       !call watcht(myrank==0,"",0)

!$OMP parallel workshare
       gam(1:MBC,ls:le) = sig(1:MBC,ls:le)*hf      ! (I-A)/2
!$OMP end parallel workshare

       !call watcht(myrank==0,"rotorb(3)",1)

       call dgemm('n','n',MBT,li,MBT,hf,tau,MBC,tau(1,ls),MBC,hf,sig(1,ls),MBC)

       !call watcht(myrank==0,"rotorb(5)",1)

       do it=1,maxit

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

          call mpi_allgatherv(wrk(1,ls),ir_i(myrank),mpi_real8 &
               ,wrk,ir_i,id_i,mpi_real8,mpi_comm_world,ierr)

          call dgemm &
               ('n','n',MBT,li,MBT,hm,wrk,MBC,wrk(1,ls),MBC,on,gamn(1,ls),MBC)

          error0=0.0d0
          do j=ls,le
          do i=j,MBT
             error0 = error0 + (gamn(i,j)-gam(i,j))**2
          end do
          end do
          call mpi_allreduce(error0,error,1,MPI_REAL8,mpi_sum,mpi_comm_world,ierr)
          error=sqrt(error)/MBT
          if ( myrank == 0 .and. error > on ) write(*,*) k,s,error

          do j=ls,le
          do i=1,MBC
             gam(i,j) = gamn(i,j)
          end do
          end do

          if ( error <= eps ) exit
          if ( it == maxit .and. myrank == 0 ) then
             write(*,*) "WARNING: iteration was not converged in rotorb"
          end if

       end do ! it

       !call watcht(myrank==0,"rotorb(6)",1)

       call mpi_allgatherv(gam(1,ls),ir_i(myrank),mpi_real8 &
            ,wrk,ir_i,id_i,mpi_real8,mpi_comm_world,ierr)

       !call watcht(myrank==0,"rotorb(7)",1)

       call dgemm('n','n',ML0,MB0,MBT,on,unk(n1,1,k,s),ML0 &
            ,wrk(1,MB_0_CPMD),MBC,zr,psi_tmp(n1,MB_0_CPMD),ML0)

!$OMP parallel workshare
       psi_n(n1:n2,MB_0_CPMD:MB_1_CPMD,k,s) = &
       psi_n(n1:n2,MB_0_CPMD:MB_1_CPMD,k,s) + psi_tmp(n1:n2,MB_0_CPMD:MB_1_CPMD)
!$OMP end parallel workshare

       !call watcht(myrank==0,"rotorb(8)",1)

       tmp1=1.d0/dt
!$OMP parallel workshare
       psi_v(n1:n2,MB_0_CPMD:MB_1_CPMD,k,s) = &
       psi_v(n1:n2,MB_0_CPMD:MB_1_CPMD,k,s) + tmp1*psi_tmp(n1:n2,MB_0_CPMD:MB_1_CPMD)
!$OMP end parallel workshare

       !call watcht(myrank==0,"rotorb(9)",1)

    end do ! k
    end do ! s

    !call watcht(myrank==0,"",0)

!$OMP parallel workshare
    unk(:,MB_0_CPMD:MB_1_CPMD,:,:) = psi_n(:,MB_0_CPMD:MB_1_CPMD,:,:)
!$OMP end parallel workshare

    !call watcht(myrank==0,"rotorb(10)",1)

    return
  END SUBROUTINE rotorb


END MODULE rotorb_module
