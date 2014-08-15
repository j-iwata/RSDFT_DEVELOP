MODULE InnerProduct
  use watch_module, only: watch
  ! Mlma,nzlma,MMJ(lma),JJP(MJJ(lma),lma),TYPE_MAIN,lma_nsend(irank),sendmap(i1,irank),
  ! sbufnl(i2,irank),rbufnl(1,jrank),uVk(MJJ(lma),lma,k),
  ! nrlma_xyz(1:6),num_2_rank(nrlma_xyz(:),1:6)
  use ps_nloc2_variables, only: nzlma,nrlma_xyz,lma_nsend,num_2_rank,sendmap,recvmap,MJJ,JJP,uVk,sbufnl,rbufnl,TYPE_MAIN,Mlma
  use pseudopot_module, only: pselect
  use rgrid_module, only: dV
  use parallel_module, only: COMM_GRID,myrank
  ! N_nlop,nlop_pair(1:2,:),qij(m)
  use VarPSMemberG, only: N_nlop,nlop_pair,qij,uVunk,uVunk0
  ! unk(nn1,n,k,s)
  use wf_module, only: unk,Sunk
  ! Sunk(nn1,n)
  use RealComplex, only: RCProduct,zero
  
  use ParaRGridComm, only: threeWayComm

  include 'mpif.h'
  PRIVATE
  PUBLIC :: dot_product,get_Sf,get_gf,get_gSf,get_Sunk_Mat


CONTAINS

!---------------------------------------------------------------------------------------	

  SUBROUTINE dot_product( a,b,c,alpha,n,m,nop )
!$ use omp_lib
    implicit none
    integer,intent(IN) :: n,m
    integer :: i
    real(8) :: a(*),b(*),c(*),alpha,nop,tmp

    c(1:m)=0.d0

    tmp=0.d0
!$OMP parallel do reduction(+:tmp)
    do i=1,n
       tmp = tmp + a(i)*b(i)
    end do
!$OMP end parallel do
    c(1) = tmp*alpha
    nop = nop + 2.d0*n +1.d0

    if ( m==2 ) then
       tmp=0.d0
!$OMP parallel do reduction(+:tmp)
       do i=1,n-1,2
          tmp = tmp + a(i)*b(i+1) - a(i+1)*b(i)
       end do
!$OMP end parallel do
       c(m) = tmp*alpha
       nop = nop + 4.d0*(n/2) + 1.d0
    end if
    return
  END SUBROUTINE dot_product

!---------------------------------------------------------------------------------------
  SUBROUTINE get_Sf(fin,nn1,nn2,k,Sf)
    ! IN:	nn1,nn2,k,fin
    !		Mlma,nzlma,pselect,
    !		MMJ(lma),JJP(MJJ(lma),lma)
    !		uVk(MJJ(lma),lma,k),
    !		dV,
    !		nrlma_xyz(1:6),num_2_rank(nrlma_xyz(:),1:6),
    !		lma_nsend(irank),sbufnl(i2,irank),rbufnl(1,jrank),sendmap(i1,irank),
    !		TYPE_MAIN,COMM_GRID
    !		N_nlop,qij(m),nlop_pair(1,N_nlop)
    ! OUT:	Sf(JJP(MMJ(nlop_pair(1,N_nlop)),nlop_pair(1,N_nlop))
!$ use omp_lib
    implicit none

    integer,intent(IN) :: nn1,nn2,k
#ifdef _DRSDFT_
    real(8),intent(IN) :: fin(nn1:nn2)
    real(8),intent(OUT) :: Sf(nn1:nn2)
    real(8) :: p_uVunk2
    real(8) :: uVunk2
    real(8) :: uVunk_z(Mlma)
    real(8) :: tmp
#else
    complex(8),intent(IN) :: fin(nn1:nn2)
    complex(8),intent(OUT) :: Sf(nn1:nn2)
    complex(8) :: p_uVunk2
    complex(8) :: uVunk2
    complex(8) :: uVunk_z(Mlma)
    complex(8) :: tmp
#endif
    integer :: kk1,i,j,ierr,ii
    integer :: ib1,ib2,ib,nb
    integer :: lma,lma1,lma2,m
    integer :: i1,i2,irank,jrank
    integer :: nreq
    integer :: istatus(MPI_STATUS_SIZE,512),ireq(512)
    real(8) :: ctt(0:1),ett(0:1)

! ??????????????????????????????????
    ib1=1
    ib2=1
    nb= ib2 - ib1 + 1

    select case ( pselect )
    case ( 1,2,102 )
#ifdef _SHOWALL_INNER_
write(200+myrank,*) 'get_Sf 1'
#endif
!$OMP parallel do
       do i=nn1,nn2
          Sf(i) = fin(i)
       end do
!$OMP end parallel do
#ifdef _SHOWALL_INNER_
write(200+myrank,*) 'get_Sf 2'
#endif

!----- get uVunk (=<uV|fin>) -----
       allocate( uVunk(nzlma,ib1:ib2) ) ; uVunk=zero
!----- within mygrid -----
#ifndef _OMP_
       do ib=ib1,ib2
          do lma=1,nzlma
             uVunk(lma,ib)=zero
             do j=1,MJJ(lma)
                i=JJP(j,lma)
#ifdef _DRSDFT_
              uVunk(lma,ib)=uVunk(lma,ib)+uVk(j,lma,k)*fin(i)
#else
              uVunk(lma,ib)=uVunk(lma,ib)+conjg(uVk(j,lma,k))*fin(i)
#endif
             end do
             uVunk(lma,ib) = dV*uVunk(lma,ib)
          end do
       end do
#else
!$OMP parallel
!$OMP workshare
       uVunk(:,:)=zero
!$OMP end workshare
       do ib=ib1,ib2
!$OMP do
          do lma=1,nzlma
             do j=1,MJJ(lma)
                i=JJP(j,lma)
#ifdef _DRSDFT_
                uVunk(lma,ib)=uVunk(lma,ib)+uVk(j,lma,k)*fin(i)
#else
                uVunk(lma,ib)=uVunk(lma,ib)+conjg(uVk(j,lma,k))*fin(i)
#endif
             end do
             uVunk(lma,ib) = dV*uVunk(lma,ib)
          end do
!$OMP end do
       end do
!$OMP end parallel
#endif
!===== within mygrid =====
#ifdef _SHOWALL_INNER_
write(200+myrank,*) 'get_Sf 3'
#endif

!----- summation over all grids -----
       call watch( ctt(0),ett(0) )

! 3WayComm
! uVunk
       call threeWayComm( nrlma_xyz,num_2_rank,sendmap,recvmap,lma_nsend,sbufnl,rbufnl,nzlma,ib1,ib2,uVunk )

       call watch( ctt(1),ett(1) )
!===== summation over all grids =====

!----- total = term1 + term2 -----
#ifndef _OMP_
       do ib=ib1,ib2
          do m=1,N_nlop
             lma1 = nlop_pair(1,m)
             lma2 = nlop_pair(2,m)
             do j=1,MJJ(lma1)
                i=JJP(j,lma1)
                Sf(i) = Sf(i) + qij(m)*uVk(j,lma1,k)*uVunk(lma2,ib)
             end do
          end do
       end do
#else
       do ib=ib1,ib2
!$OMP parallel private(lma1,lma2,i)
          do m=1,N_nlop
             lma1 = nlop_pair(1,m)
             lma2 = nlop_pair(2,m)
!$OMP do
             do j=1,MJJ(lma1)
                i=JJP(j,lma1)
                Sf(i) = Sf(i) + qij(m)*uVk(j,lma1,k)*uVunk(lma2,ib)
             end do
!$OMP end do
          end do
!$OMP end parallel
       end do
#endif
!===== total = term1 + term2 =====

       deallocate( uVunk )
#ifdef _SHOWALL_INNER_
write(200+myrank,*) 'get_Sf 4'
#endif

!    case ( 3 )
!----- term1 -----
!!$OMP parallel do
!       do i=nn1,nn2
!          Sf(i) = fin(i)
!       end do
!!$OMP end parallel do
!===== term1 =====

!----- term2 -----
!       do i=1,Mlma
!          uVunk2=zero
!          p_uVunk2=zero
!!$OMP parallel do reduction(+:p_uVunk2)
!          do j=nn1,nn2
!             call RCProduct( uVk(i,j,k),fin(j),tmp )
!             p_uVunk2 = p_uVunk2 + tmp
!          end do
!!$OMP end parallel do
!          p_uVunk2 = dV*p_uVunk2

!          call MPI_ALLREDUCE(p_uVunk2,uVunk2,1,TYPE_MAIN,MPI_SUM,COMM_GRID,ierr)

!          uVunk_z(i) = uVunk2
!       end do

!       do kk1=1,N_nlop
!          i = nlop_pair(1,m)
!          j = nlop_pair(2,m)
!!$OMP parallel do
!          do ii=nn1,nn2
!             Sf(ii) = Sf(ii) + uVk(ii,i,k)*qij(kk1)*uVunk_z(j)
!          end do
!!$OMP end parallel do
!       end do
!===== term2 =====
    end select
#ifdef _SHOWALL_INNER_
write(200+myrank,*) 'get_Sf 5'
#endif

    return
  END SUBROUTINE get_Sf

!---------------------------------------------------------------------------------------

  SUBROUTINE get_gf( gin,fin,nn1,nn2,gf,switch_zp )
    ! IN:	nn1,nn2,gin(nn1:nn2),fin(nn1:nn2),switch_zp,
    !		dV,
    !		TYPE_MAIN,COMM_GRID
    ! OUT:	gf
!$ use omp_lib
    implicit none

    integer,intent(IN) :: nn1,nn2
    integer,intent(IN) :: switch_zp
#ifdef _DRSDFT_
    real(8),intent(IN) :: fin(nn1:nn2),gin(nn1:nn2)
    real(8),intent(OUT) :: gf
    real(8) :: gf_0
    real(8) :: tmp
#else
    complex(8),intent(IN) :: fin(nn1:nn2),gin(nn1:nn2)
    complex(8),intent(OUT) :: gf
    complex(8) :: gf_0
    complex(8) :: tmp
#endif
    integer :: i
    integer :: ierr

    gf_0=zero

!$OMP parallel do reduction(+:gf_0)
#ifdef _DRSDFT_
    do i=nn1,nn2
       gf_0 = gf_0 + gin(i)*fin(i)
    end do
#else
    do i=nn1,nn2
       gf_0 = gf_0 + conjg(gin(i))*fin(i)
    end do
#endif
!$OMP end parallel do
    gf_0 = dV*gf_0

    if ( switch_zp==1 ) then
       call MPI_ALLREDUCE(gf_0,gf,1,TYPE_MAIN,MPI_SUM,COMM_GRID,ierr)
    else
       gf = gf_0
    end if

    return
  END SUBROUTINE get_gf

  !---------------------------------------------------------------------------------------

  SUBROUTINE get_gSf( gin,fin,nn1,nn2,k,gSf,switch_zp )
    ! IN:	nn1,nn2,gin(nn1:nn2),fin(nn1:nn2),k,switch_zp,
    ! OUT:	gSf
!$ use omp_lib
    implicit none

    integer,intent(IN) :: nn1,nn2,k		
    integer,intent(IN) :: switch_zp
#ifdef _DRSDFT_
    real(8),intent(IN) :: gin(nn1:nn2),fin(nn1:nn2)
    real(8),intent(OUT) :: gSf
    real(8) :: Sf(nn1:nn2),gSf_0
#else
    complex(8),intent(IN) :: gin(nn1:nn2),fin(nn1:nn2)
    complex(8),intent(OUT) :: gSf
    complex(8) :: Sf(nn1:nn2),gSf_0
#endif
    integer :: i
    integer :: ierr

!----- |f'> = S|f> = Sf -----
    call get_Sf( fin,nn1,nn2,k,Sf )

!----- <g|f'> = <g|S|f> = gSf
    call get_gf( gin,Sf,nn1,nn2,gSf,switch_zp )

    return
  END SUBROUTINE get_gSf

!---------------------------------------------------------------------------------------

  SUBROUTINE get_Sunk_Mat( nn1,nn2,ns,ne,k,s )
    ! IN:	nn1,nn2,ns,ne,k,s,
    !		unk(nn1,n,k,s)
    ! OUT:	Sunk(nn1,n)
    implicit none

    integer,intent(IN) :: nn1,nn2,ns,ne,k,s
    integer :: n
integer :: i

    do n=ns,ne
      call get_Sf( unk(nn1,n,k,s),nn1,nn2,k,Sunk(nn1,n) )
    end do

    return
  END SUBROUTINE get_Sunk_Mat

END MODULE InnerProduct
