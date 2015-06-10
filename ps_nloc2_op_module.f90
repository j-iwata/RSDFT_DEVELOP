MODULE ps_nloc2_op_module

!$ use omp_lib
  use ps_nloc2_variables
  use watch_module, only: et_hpsi_
  use parallel_module
  use rgrid_module, only: dV,Igrid

  implicit none

  PRIVATE
  PUBLIC :: op_ps_nloc2_hp, init_op_ps_nloc2_hp

  integer,allocatable :: ompMJJ1(:,:), ompMJJ2(:,:)
  integer,allocatable :: omplns1(:,:), omplns2(:,:)
  integer :: ompnzlma1,ompnzlma2,ompmyrank,ompnprocs,ompn1,ompn2

  logical :: init_flag=.false.

!$OMP threadprivate( ompmyrank,ompnprocs,ompnzlma1,ompnzlma2,ompn1,ompn2 )

CONTAINS


  SUBROUTINE init_op_ps_nloc2_hp
    implicit none
    integer :: i,j,k,ompblock,ompblock0,mm

    if ( init_flag ) return
    init_flag = .true.

    mm = Igrid(2,0) - Igrid(1,0) + 1

!$OMP parallel private( i,j,k,ompblock,ompblock0 )

    ompnprocs = 1
!$  ompnprocs = omp_get_num_threads()
    ompmyrank = 0
!$  ompmyrank = omp_get_thread_num()

    ompblock=0
    ompblock0=0
    do i=1,nzlma
       j=mod(i-1,ompnprocs)
       if ( j == ompmyrank ) then
          ompblock = ompblock + 1
       else if ( j < ompmyrank ) then
          ompblock0 = ompblock0 + 1
       end if
    end do

    ompnzlma1 = ompblock0 + 1
    ompnzlma2 = ompnzlma1 + ompblock - 1

! ---

!$OMP master
    allocate( omplns1(0:np_grid-1,0:ompnprocs-1) ) ; omplns1=0
    allocate( omplns2(0:np_grid-1,0:ompnprocs-1) ) ; omplns2=0
    allocate( ompMJJ1(0:nzlma,0:ompnprocs-1)     ) ; ompMJJ1=0
    allocate( ompMJJ2(0:nzlma,0:ompnprocs-1)     ) ; ompMJJ2=0
!$OMP end master
!$OMP barrier

    omplns1(:,ompmyrank)=-1
    omplns2(:,ompmyrank)=-2
    do i=0,np_grid-1
       do j=1,lma_nsend(i)
          k = sendmap(j,i)
          if ( ompnzlma1 <= k .and. k <= ompnzlma2 ) then
             omplns2(i,ompmyrank) = j
          end if
       end do ! j
       do j=lma_nsend(i),1,-1
          k = sendmap(j,i)
          if ( ompnzlma1 <= k .and. k <= ompnzlma2 ) then
             omplns1(i,ompmyrank) = j
          end if
       end do ! j
    end do ! i

! ---

    ompblock=0
    ompblock0=0
    do i=1,mm
       j=mod(i-1,ompnprocs)
       if ( j == ompmyrank ) then
          ompblock = ompblock + 1
       else if ( j < ompmyrank ) then
          ompblock0 = ompblock0 + 1
       end if
    end do

    ompn1 = Igrid(1,0) + ompblock0
    ompn2 = ompn1 + ompblock - 1

! ---

    ompMJJ1(:,ompmyrank)=-1
    ompMJJ2(:,ompmyrank)=-2
    do k=1,nzlma
       do j=1,MJJ(k)
          i=JJP(j,k)
          if ( ompn1 <= i .and. i<= ompn2 ) then
             ompMJJ1(k,ompmyrank) = j
             exit
          end if
       end do
       do j=MJJ(k),1,-1
          i=JJP(j,k)
          if ( ompn1 <= i .and. i<= ompn2 ) then
             ompMJJ2(k,ompmyrank) = j
             exit
          end if
       end do
    end do

!$OMP end parallel

  END SUBROUTINE init_op_ps_nloc2_hp


  SUBROUTINE op_ps_nloc2_hp(k,tpsi,htpsi,n1,n2,ib1,ib2)

    implicit none

    integer,intent(IN) :: k,n1,n2,ib1,ib2
    real(8),intent(IN)  :: tpsi(n1:n2,ib1:ib2)
    real(8),intent(OUT) :: htpsi(n1:n2,ib1:ib2)
    integer :: i,ib,i1,i2,j,jb,lma,m,ML0,n,nb
    integer :: ierr,nreq
    integer :: irank,jrank,istatus(mpi_status_size,512),ireq(512)
    real(8) :: c
    real(8) :: ct0,et0,ct1,et1

    ML0 = n2-n1+1
    nb  = ib2-ib1+1

!$OMP barrier
!$OMP master
    et0=omp_get_wtime()
!$OMP end master

!!$omp master
!    allocate( uVunk(0:nzlma,ib1:ib2),uVunk0(0:nzlma,ib1:ib2) )
!!$omp end master
!!$omp barrier

!$OMP barrier
!$OMP master
    et1=omp_get_wtime()
    et_hpsi_(5)=et_hpsi_(5)+et1-et0
!$OMP end master

    do ib=ib1,ib2
       jb=ib-ib1+1
!$omp do
       do lma=1,nzlma
          uVunk(lma,jb)=zero
          do j=1,MJJ(lma)
             i=JJP(j,lma)
             uVunk(lma,jb)=uVunk(lma,jb)+uVk(j,lma,k)*tpsi(i,ib)
          end do
          uVunk(lma,jb)=iuV(lma)*dV*uVunk(lma,jb)
       end do
!$omp end do nowait
    end do

!$OMP barrier
!$OMP master
    et0=omp_get_wtime()
    et_hpsi_(6)=et_hpsi_(6)+et0-et1
!$OMP end master
         
    do i=1,6

!$omp barrier

       select case(i)
       case(1,3,5)
          j=i+1
          do ib=1,nb
             uVunk0(ompnzlma1:ompnzlma2,ib) &
                  =uVunk(ompnzlma1:ompnzlma2,ib)
          end do
       case(2,4,6)
          j=i-1
       end select

!$omp barrier

       do m=1,nrlma_xyz(i)
          nreq=0
          irank=num_2_rank(m,i)
          jrank=num_2_rank(m,j)
          if( irank>=0 )then
             i2=0
             do ib=1,nb
                i2=omplns1(irank,ompmyrank)-1+lma_nsend(irank)*(ib-1)
                do i1=omplns1(irank,ompmyrank),omplns2(irank,ompmyrank)
                   i2=i2+1
                   sbufnl(i2,irank)=uVunk0(sendmap(i1,irank),ib)
                end do
             end do
!$omp barrier
!$omp master
             nreq=nreq+1
             call mpi_isend(sbufnl(1,irank),lma_nsend(irank)*nb &
                  ,TYPE_MAIN,irank,1,comm_grid,ireq(nreq),ierr)
!$omp end master
          end if
          if( jrank>=0 )then
!$omp master
             nreq=nreq+1
             call mpi_irecv(rbufnl(1,jrank),lma_nsend(jrank)*nb &
                  ,TYPE_MAIN,jrank,1,comm_grid,ireq(nreq),ierr)
!$omp end master
          end if
!$omp master
          call mpi_waitall(nreq,ireq,istatus,ierr)
!$omp end master
!$omp barrier
          if( jrank >= 0 )then
             i2=0
             do ib=1,nb
                i2=omplns1(jrank,ompmyrank)-1 &
                     +lma_nsend(jrank)*(ib-1)
                do i1=omplns1(jrank,ompmyrank), &
                      omplns2(jrank,ompmyrank)
                   i2=i2+1
                   uVunk(recvmap(i1,jrank),ib) &
                  =uVunk(recvmap(i1,jrank),ib)+rbufnl(i2,jrank)
                end do
             end do
          end if

       end do

!$omp barrier

    end do

!$OMP barrier
!$OMP master
    et1=omp_get_wtime()
    et_hpsi_(7)=et_hpsi_(7)+et1-et0
!$OMP end master

!$omp barrier
    do ib=ib1,ib2
       jb=ib-ib1+1
       do lma=1,nzlma
          do j=ompMJJ1(lma,ompmyrank),ompMJJ2(lma,ompmyrank)
             i=JJP(j,lma)
             htpsi(i,ib)=htpsi(i,ib)+uVk(j,lma,k)*uVunk(lma,jb)
          end do
       end do
    end do

!$OMP barrier
!$OMP master
    et0=omp_get_wtime()
    et_hpsi_(8)=et_hpsi_(8)+et0-et1
!$OMP end master

!!$omp barrier
!!$omp master
!    deallocate( uVunk0,uVunk )
!!$omp end master

    return 
      
  END SUBROUTINE op_ps_nloc2_hp


END MODULE ps_nloc2_op_module

