!-----------------------------------------------------------------------
subroutine rotorb
  use parallel_module
  use array_bound_module, only: MBZ_0,MBZ_1,MSP_0,MSP_1
  use wf_module, only: unk
  use cpmd_variables, only: gam,scr,wrk,tau,sig,gamn,psi_v,psi_n &
                           ,mstocck,dt,disp_switch,MBC,MBT,MB_0_CPMD,MB_1_CPMD
  use watch_module
  implicit none
  integer,parameter :: maxit=100
  integer i,j,k,l,m,n,s,ii,it,ij
  integer n1,n2,ML0,ierr,MB0
  integer,allocatable :: ir_i(:),id_i(:)
  integer ls,le,li
  real(8),parameter :: eps=1.0d-8
  real(8),parameter :: hf = 0.5d0
  real(8),parameter :: hm =-0.5d0
  real(8),parameter :: on = 1.0d0
  real(8),parameter :: zr = 0.0d0
  real(8) :: ctime_force(0:10),etime_force(0:10)
  real(8) :: error,error0,tmp1
  real(8),allocatable :: psi_tmp(:,:)
  complex(8),parameter :: zo = (1.0d0,0.0d0)
  complex(8),parameter :: zz = (0.0d0,0.0d0)
  complex(8),allocatable :: zrk(:,:)
  logical :: dirprod,lblas

  if ( MBC < nprocs ) then
     write(*,*) "MBC<nprocs!!!",MBC,nprocs,myrank
     stop
  end if

  n1    = idisp(myrank)+1
  n2    = idisp(myrank)+ircnt(myrank)
  ML0   = ircnt(myrank)
  lblas =.true.
  MB0   = MB_1_CPMD-MB_0_CPMD+1

  allocate( id_i(0:nprocs-1),ir_i(0:nprocs-1) )
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

  ls=id_i(myrank)+1
  le=id_i(myrank)+ir_i(myrank)
  li=ir_i(myrank)

  id_i(0:nprocs-1)=id_i(0:nprocs-1)*MBC
  ir_i(0:nprocs-1)=ir_i(0:nprocs-1)*MBC

  allocate( psi_tmp(ML0,MB0) )

  do s=MSP_0,MSP_1
  do k=MBZ_0,MBZ_1

     call watcht(myrank==0,"",0)

     MBT=mstocck(k,s)

     tau(:,:)=zr
     sig(:,:)=zr
     gam(:,:)=zr
     gamn(:,:)=zr
     scr(:,:)=zr
     wrk(:,:)=zr

     call watcht(myrank==0,"rotorb(0)",1)

     call overlap(s,k,4,0,lblas)
     call watcht(myrank==0,"rotorb(1)",1)
     call overlap(s,k,2,0,lblas)
     call watcht(myrank==0,"rotorb(2)",1)

     do i=ls,le
        sig(i,i)=on+sig(i,i) ! I-A
     enddo
     do i=ls,le
        tau(i,i)=on+tau(i,i) ! I-B
     enddo

     gam(:,:)=sig(:,:)*hf    ! (I-A)/2

     wrk(:,:)=zr
     do i=ls,le
     do j=1,MBC
        wrk(j,i)=tau(j,i)
     enddo
     enddo

     call watcht(myrank==0,"rotorb(3)",1)

     call mpi_allgatherv(wrk(1,ls),ir_i(myrank),mpi_real8,wrk &
          ,ir_i,id_i,mpi_real8,mpi_comm_world,ierr)

     call watcht(myrank==0,"rotorb(4)",1)

     call dgemm('n','n',MBT,li,MBT,hf,wrk,MBC,tau(1,ls),MBC,hf,sig(1,ls),MBC)

     call watcht(myrank==0,"rotorb(5)",1)

     do it=1,maxit

        do i=ls,le
        do j=1,MBC
           gamn(j,i)=sig(j,i)
        enddo
        enddo
        do i=ls,le
        do j=1,MBC
           tmp1=tau(j,i)-gam(j,i)
           scr(j,i)=tmp1
           wrk(j,i)=tmp1
        enddo
        enddo

        call mpi_allgatherv(wrk(1,ls),ir_i(myrank),mpi_real8,wrk &
             ,ir_i,id_i,mpi_real8,mpi_comm_world,ierr)

        call dgemm('n','n',MBT,li,MBT,hm,wrk,MBC,scr(1,ls),MBC,on,gamn(1,ls),MBC)

        error0=0.0d0
        do i=ls,le
        do j=i,MBT
           error0=error0+(gamn(j,i)-gam(j,i))**2
        enddo
        enddo
        call mpi_allreduce(error0,error,1,MPI_REAL8,mpi_sum,mpi_comm_world,ierr)
        error=sqrt(error)/MBT
        if ( myrank==0 .and. error > on ) write(*,*) k,s,error

        do i=ls,le
        do j=1,MBC
           gam(j,i)=gamn(j,i)
        enddo
        enddo

        if ( error <= eps ) exit
        if ( it == maxit .and. myrank==0 ) then
           write(*,*) "WARNING: iteration was not converged in rotorb"
        end if

        call watcht(myrank==0,"rotorb(6)",1)

     end do ! it

     call watcht(myrank==0,"",0)

     call mpi_allgatherv(gam(1,ls),ir_i(myrank),mpi_real8,wrk &
          ,ir_i,id_i,mpi_real8,mpi_comm_world,ierr)

     call watcht(myrank==0,"rotorb(7)",1)

!     do ms=MB_0,min(MB_1,MBT),MBLK
!        me=min(ms+MBLK-1,min(MB_1,MBT))
!        mn=me-ms+1
! (org)
!        do i=n1,n2,NBLK2
!           mm=min(NBLK2,n2-i+1)
!           call dgemm('n','t',mm,MBT,mn,on,unk(i,ms,k,s),ML0 &
!                ,wrk(1,ms),MBT,on,psi_n(i,1,k,s),ML0)
!        enddo
! (1)
!        call dgemm('n','t',ML0,MBT,mn,on,unk(n1,ms,k,s),ML0 &
!             ,wrk(1,ms),MBT,on,psi_n(n1,1,k,s),ML0)
! (test)
!        call dgemm('n','n',ML0,mn,MBT,on,unk(n1,1,k,s),ML0 &
!             ,wrk(1,ms),MBC,on,psi_n(n1,ms,k,s),ML0)
!        call dgemm('n','n',ML0,mn,MBT,on,unk(n1,1,k,s),ML0 &
!             ,wrk(1,ms),MBC,on,psi_n(n1,ms,k,s),ML0)
!     enddo
     call dgemm('n','n',ML0,MB0,MBT,on,unk(n1,1,k,s),ML0 &
             ,wrk(1,MB_0_CPMD),MBC,zr,psi_tmp(1,1),ML0)
     psi_n(n1:n2,MB_0_CPMD:MB_1_CPMD,k,s) = psi_n(n1:n2,MB_0_CPMD:MB_1_CPMD,k,s) &
                                          + psi_tmp(1:ML0,1:MB0)
!     call dgemm('n','n',ML0,MB0,MBT,on,unk(n1,1,k,s),ML0 &
!             ,wrk(1,MB_0_CPMD),MBC,on,psi_n(n1,MB_0_CPMD,k,s),ML0)

     call watcht(myrank==0,"rotorb(8)",1)

!     do ms=MB_0,min(MB_1,MBT),MBLK
!        me=min(ms+MBLK-1,min(MB_1,MBT))
!        mn=me-ms+1
! (org)
!        do i=n1,n2,NBLK2
!           mm=min(NBLK2,n2-i+1)
!           call dgemm('n','t',mm,MBT,mn,on/dt,unk(i,ms,k,s),ML0 &
!                ,wrk(1,ms),MBT,on,psi_v(i,1,k,s),ML0)
!        enddo
! (1)
!        call dgemm('n','t',ML0,MBT,mn,on/dt,unk(n1,ms,k,s),ML0 &
!             ,wrk(1,ms),MBT,on,psi_v(n1,1,k,s),ML0)
! (test)
!        call dgemm('n','n',ML0,mn,MBT,on/dt,unk(n1,1,k,s),ML0 &
!             ,wrk(1,ms),MBC,on,psi_v(n1,ms,k,s),ML0)
!     enddo
     tmp1=1.d0/dt
     psi_v(n1:n2,MB_0_CPMD:MB_1_CPMD,k,s) = psi_v(n1:n2,MB_0_CPMD:MB_1_CPMD,k,s) &
                                          + tmp1*psi_tmp(1:ML0,1:MB0)
!     call dgemm('n','n',ML0,MB0,MBT,on/dt,unk(n1,1,k,s),ML0 &
!             ,wrk(1,MB_0_CPMD),MBC,on,psi_v(n1,MB_0_CPMD,k,s),ML0)

     call watcht(myrank==0,"rotorb(9)",1)

  enddo ! k
  enddo ! s

  deallocate( psi_tmp )
  deallocate( ir_i, id_i )

  call watcht(myrank==0,"",0)

  unk(:,MB_0_CPMD:MB_1_CPMD,:,:)=psi_n(:,MB_0_CPMD:MB_1_CPMD,:,:)

  call watcht(myrank==0,"rotorb(10)",1)

  return
end subroutine rotorb
