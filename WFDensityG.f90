MODULE WFDensityG
  use rgrid_module, only: dV
  use electron_module, only: Nelectron
  ! unk(:,n,k,s)
  use wf_module
  use parallel_module, only: COMM_GRID,COMM_BAND,COMM_BZSM,COMM_SPIN,myrank
  use array_bound_module, only: ML_0,ML_1,MB_0,MB_1,MBZ_0,MBZ_1,MSP_0,MSP_1
  use RealComplex, only: RCProduct
    use ps_nloc2_variables, only: Mlma,nzlma,MJJ,JJP,nrlma_xyz,num_2_rank,sbufnl,rbufnl,lma_nsend,sendmap,recvmap,TYPE_MAIN,uVk
    use pseudopot_module, only: pselect
    
    use VarPSMemberG
    ! N_nzqr,nzqr_pair,MJJ_Q(),JJP_Q(),sendmap_Q,recvmap_Q,qr_nsend,sbufnl_Q,rbufnl_Q,num_2_rank_Q,uVunk,uVunk0,QRij
    use VarParaPSnonLocG, only: MJJ_Q,JJP_Q
    use ParaRGridComm, only: do3StepComm

  implicit none
  include 'mpif.h'

  PRIVATE
  PUBLIC :: get_rhonks
  !	PUBLIC :: rho,init_density,normalize_density,calc_density

  ! WARNING_USPP
  real(8),allocatable :: rho(:,:)

CONTAINS

  !---------------------------------------------------------------------------------------
  !	SUBROUTINE init_density
  !		implicit none
  !		real(8) :: c,d

  !		if ( .not. allocated(rho) ) then
  !			allocate( rho(ML_0:ML_1,MSP) )
  !			call random_number( rho )
  !			call normalize_density
  !		end if
  !	END SUBROUTINE init_density

  !---------------------------------------------------------------------------------------

  !	SUBROUTINE normalize_density
  !		implicit none
  !		real(8) :: c,d
  !		integer :: ierr
  !		include 'mpif.h'

  !		c = dV*sum( rho )
  !		call MPI_ALLREDUCE(c,d,1,MPI_REAL8,MPI_SUM,COMM_GRID,ierr )
  !		c = Nelectron/d
  !		rho = c*rho
  !	END SUBROUTINE normalize_density

  !---------------------------------------------------------------------------------------

  !	SUBROUTINE calc_density
  !	! IN:	unk(:.n.k.s),ML_0,ML_1,MSP_0,MSP_1
  !	! OUT:	rho(n1:n2,s)
  !		implicit none
  !		integer :: n,k,s
  !		integer :: n1,n2,n0
  !		real(8),allocatable :: rhonks(:)

  !		n1=ML_0
  !		n2=ML_1
  !		n0 = ML_1 - ML_0 + 1

  !		allocate( rhonks(n1:n2) ) ; rhonks(:)=0.d0

  !		rho(:,:)=0.d0
  !		do s=MSP_0,MSP_1
  !			do k=MBZ_0,MBZ_1
  !				do n=MB_0,MB_1
  !					rhonks(:)=0.d0
  !					call get_rhonks( rhonks,n1,n2,n,k,s )
  !					rho(:,s) = rho(:,s) + occ(n,k,s)*rhonks(:)
  !				end do
  !			end do
  !		end do

  !		call reduce_and_gather

  !		deallocate( rhonks )

  !	END SUBROUTINE calc_density

  !---------------------------------------------------------------------------------------

  SUBROUTINE get_rhonks( rhonks,nn1,nn2,n,k,s )
    ! IN:	nn1,nn2,n,k,s,
    !		unk(nn1:nn2,n,k,s),dV
    !		uVnk(MJJ(lma),lma,k),uVnk0(MJJ(lma),lma,k),
    !		Mlma,nzlma,MJJ(lma),JJP(MJJ(lma),lma),
    !		nrlma_xyz(1:6),num_2_rank(nrlma_xyz(i),i),lma_nsend(num_2_rank(m,i)),sbufnl(i2,irank),rbufnl(i2,jrank),TYPE_MAIN,
    !		sendmap(i1,irank),recvmap(i1,jrank),uVk(:,:,:)
    !		N_nzqr,nzqr_pair(N_nzqr,2),
    !		MJJ_Q(N_nzqr),JJP_Q(MJJ_Q(kk1),kk1),QRij(MJJ_Q(kk1),kk1),
    !		myrank
    ! OUT:	rhonks(:)

    implicit none
    integer,intent(IN) :: nn1,nn2
    integer,intent(IN) :: n,k,s
    real(8),intent(INOUT) :: rhonks(nn1:nn2)
#ifdef _DRSDFT_
    real(8) :: Qrhonks(nn1:nn2)
    real(8) :: p_uVunk1,p_uVunk2
    real(8) :: uVunk1,uVunk2
    real(8),parameter :: zero=0.d0
    real(8) :: uVunk2_z(Mlma)
    real(8) :: tmp
#else
    complex(8) :: Qrhonks(nn1:nn2)
    complex(8) :: p_uVunk1,p_uVunk2
    complex(8) :: uVunk1,uVunk2
    complex(8),parameter :: zero=(0.d0,0.d0)
    complex(8) :: uVunk2_z(Mlma)
    complex(8) :: tmp
#endif
    integer :: kk1,i,j
    integer :: ierr
    integer :: lma,irank,jrank,i1,i2,nreq
    integer :: ib1,ib2,nb,ib,lma1,lma2,m
    integer :: istatus(MPI_STATUS_SIZE,512),ireq(512)

    select case ( pselect )
    case ( 1,2,102 )

!----- term1 -----
      rhonks(nn1:nn2) = abs( unk(nn1:nn2,n,k,s )**2 )
!===== term1 =====

!----- term2 -----
!----- get_uVunk -----
      ib1=n
      ib2=n
      nb = ib2 - ib1 + 1
      allocate( uVunk(nzlma,ib1:ib2)  ) ; uVunk(:,:) =zero
      do ib=ib1,ib2
        do lma=1,nzlma
          do j=1,MJJ(lma)
            i=JJP(j,lma)
#ifdef _DRSDFT_
            uVunk(lma,ib)=uVunk(lma,ib)+uVk(j,lma,k)*unk(i,ib,k,s)
#else
            uVunk(lma,ib)=uVunk(lma,ib)+conjg(uVk(j,lma,k))*unk(i,ib,k,s)
#endif
          end do
          uVunk(lma,ib) = dV*uVunk(lma,ib)
        end do	! lma
      end do	! ib


! 3WayComm
      call do3StepComm(nrlma_xyz,num_2_rank,sendmap,recvmap,lma_nsend,sbufnl,rbufnl,nzlma,ib1,ib2,uVunk)
!===== term2 =====

!----- get Qrhonks -----
      Qrhonks=zero
      do ib=ib1,ib2
        do kk1=1,N_nzqr
          lma1=nzqr_pair(kk1,1)
          lma2=nzqr_pair(kk1,2)
            if ( lma1<lma2 ) stop 'NZQR_PAIR is stange'
            if ( lma1==lma2 ) then
              do j=1,MJJ_Q(kk1)
                i=JJP_Q(j,kk1)
                if ( i<nn1 .or. nn2<i ) then
                  write(*,'(A,6I)') 'STOP *****  myrank,i,nn1,nn2,kk1,j=',myrank,i,nn1,nn2,kk1,j
                  stop
                end if
#ifdef _DRSDFT_
                Qrhonks(i) = Qrhonks(i) + QRij(j,kk1)*uVunk(lma1,ib)*uVunk(lma2,ib)
#else
                Qrhonks(i) = Qrhonks(i) + QRij(j,kk1)*conjg(uVunk(lma1,ib))*uVunk(lma2,ib)
#endif
              end do
            else
              do j=1,MJJ_Q(kk1)
                i=JJP_Q(j,kk1)
#ifdef _DRSDFT_
                Qrhonks(i) = Qrhonks(i) + QRij(j,kk1)*uVunk(lma1,ib)*uVunk(lma2,ib)
                Qrhonks(i) = Qrhonks(i) + QRij(j,kk1)*uVunk(lma2,ib)*uVunk(lma1,ib)
#else
                Qrhonks(i) = Qrhonks(i) + QRij(j,kk1)*conjg(uVunk(lma1,ib))*uVunk(lma2,ib)
                Qrhonks(i) = Qrhonks(i) + QRij(j,kk1)*conjg(uVunk(lma2,ib))*uVunk(lma1,ib)
#endif
              end do
            end if
          end do
        end do	! ib
!===== get Qrhonks =====

!----- total = term1 + term2 -----
        rhonks(nn1:nn2) = rhonks(nn1:nn2) + real( Qrhonks(nn1:nn2) )
!===== total = term1 + term2 =====

        deallocate( uVunk  )

!    case ( 3 )

       !----- term1 -----
!       rhonks(nn1:nn2) = abs( unk(nn1:nn2,n,k,s) )**2
       !===== term1 =====

       !----- term2 -----
!       do i=1,Mlma
!          p_uVunk2=zero
!          uVunk2=zero
!          do nn=nn1,nn2
!             call RCProduct( uVk(nn,i,s),unk(nn,n,k,s),tmp )
!             p_uVunk2 = p_uVunk2 + tmp
!          end do
!          p_uVunk2 = dV*p_uVunk2
!          call MPI_ALLREDUCE( p_uVunk2,uVunk2,1,TYPE_MAIN,MPI_SUM,COMM_GRID,ierr )
!          uVunk2_z(i) = uVunk2
!       end do

!       Qrhonks=zero
!       do kk1=1,N_nzqr
!          i=nzqr_pair(kk1,1)
!          j=nzqr_pair(kk1,2)
!          if ( i==j ) then
!             call RCProduct( uVunk2_z(i),uVunk2_z(j),tmp )
!             Qrhonks = Qrhonks + QRij(nn1:nn2,kk1)*tmp
!          else if ( i>j ) then
!             call RCProduct( uVunk2_z(i),uVunk2_z(j),tmp )
!             Qrhonks = Qrhonks + QRij(nn1:nn2,kk1)*tmp
!             call RCProduct( uVunk2_z(j),uVunk2_z(i),tmp )
!             Qrhonks = Qrhonks + QRij(nn1:nn2,kk1)*tmp
!          else
!             stop 'get_rhonks'
!          end if
!       end do
       !===== term2 =====

       !----- total = term1 + term2 -----
!       rhonks(nn1:nn2) = rhonks(nn1:nn2) + real( Qrhonks(nn1:nn2) )
       !===== total = term1 + term2 =====

    end select

    return
  END SUBROUTINE get_rhonks

  !---------------------------------------------------------------------------------------

  !	SUBROUTINE reduce_and_gather
  !	! IN:	rho(n1:n2,s),ML_0,ML_1,MSP_0,MSP_1,COMM_BAND,COMM_BZSM,COMM_SPIN
  !	! OUT:	rho(n1:n2,s)
  !		use parallel_module, only: id_spin,ir_spin,nprocs_s,myrank_s
  !		implicit none
  !		real(8),allocatable :: w(:)
  !!		integer,allocatable :: ir(:),id(:)
  !		integer :: n,k,s,ierr
  !		integer :: n1,n2,ML0
  !		include 'mpif.h'

  !		n1=ML_0
  !		n2=ML_1
  !		ML0 = n2 - n1 + 1
  !		allocate( w(n1:n2) )
  !		do s=MSP_0,MSP_1
  !			call MPI_ALLREDUCE( rho(n1,s),w(n1),ML0,MPI_REAL8,MPI_SUM,COMM_BAND,ierr )
  !			call MPI_ALLREDUCE( w(n1),rho(n1,s),ML0,MPI_REAL8,MPI_SUM,COMM_BZSM,ierr )
  !		end do

  !	! The following assumes all 'MSP_1-MSP_0+1' are the same
  !	! ?????????
  !		m = m*(MSP_1-MSP_0+1)
  !		call MPI_ALLGATHER( rho(n1,MSP_0),m,MPI_REAL8,rho,m,MPI_REAL8,COMM_SPIN,ierr )
  !	!----- old ver. -----
  !!		allocate( ir(0:nprocs_s-1) ) ; ir=0
  !!		allocate( id(0:nprocs_s-1) ) ; id=0
  !!		id(0:nprocs_s-1) = id_spin(0:nprocs_s-1)*ML0
  !!		ir(0:nprocs_s-1) = ir_spin(0:nprocs_s-1)*ML0
  !!		call MPI_ALLGATHERV( rho(n1,MSP_0),ir(myrank_s),MPI_REAL8,rho,ir,id,MPI_REAL8,COMM_SPIN,ierr )
  !!		deallocate( id )
  !!		deallocate( ir )
  !	!===== old ver. =====

  !		deallocate( w )
  !	END SUBROUTINE reduce_and_gather

END MODULE WFDensityG
