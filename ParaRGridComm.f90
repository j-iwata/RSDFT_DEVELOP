MODULE ParaRGridComm
  use parallel_module, only: node_partition,myrank_g,nprocs_g,COMM_GRID,myrank
  use VarPara
  ! TYPE_MAIN,node_partition,nprocs_g,myrank_g
  implicit none

  PRIVATE
  PUBLIC :: prepThreeWayComm,threeWayComm
  

CONTAINS

!---------------------------------------------------------------------------------------
  SUBROUTINE prepThreeWayComm( nr,NLRankMap,NRxyz,Num2Rank0 )
    implicit none
    
    integer,intent(IN) :: nr
    integer,intent(IN) :: NLRankMap(nr)
    integer,intent(OUT) :: NRxyz(1:6)
    integer,allocatable,intent(OUT) :: Num2Rank0(:,:)
    
    integer :: ir,n,m,i,j,k
    integer :: i1,i2,i3
    integer :: np1,np2,np3
    integer,allocatable :: LLp(:,:)
    integer,allocatable :: itmp(:,:),itmp1(:),itmp2(:),itmp3(:,:),work(:)

write(400+myrank,*) ">>>> prepThreeWayComm"
    np1=node_partition(1)
    np2=node_partition(2)
    np3=node_partition(3)
    
    allocate( LLp(3,0:nprocs_g-1) ) ; LLp=0
    n=-1
    do i3=0,np3-1
        do i2=0,np2-1
            do i1=0,np1-1
                n=n+1
                LLp(1,n)=i1
                LLp(2,n)=i2
                LLp(3,n)=i3
            end do
        end do
    end do
    
    allocate( itmp(3,nr)  ) ; itmp =0
    allocate( itmp1(nr)   ) ; itmp1=0
    allocate( itmp2(nr)   ) ; itmp2=0
    allocate( itmp3(3,nr) ) ; itmp3=0
    allocate( work(nr)    ) ; work =0

    do ir=1,nr
        n=NLRankMap(ir)
        itmp(1,ir)=LLp(1,n)-LLp(1,myrank_g)
        itmp(2,ir)=LLp(2,n)-LLp(2,myrank_g)
        itmp(3,ir)=LLp(3,n)-LLp(3,myrank_g)
    end do
    
    NRxyz(1:6)=0
    
    !----- itmp(1,*)>0 -----
    m=0
    n=0
    do i=1,nr
        if ( itmp(2,i)==0 .and. itmp(3,i)==0 .and. itmp(1,i)>0 ) then
            n=n+1
            work(n)=itmp(1,i)
            itmp2(n)=i
        end if
    end do
    if ( n>0 ) then
        call indexx( n,work,itmp1 )
        do i=1,n
            j=itmp2( itmp1(i) )
            itmp3(:,m+i)=itmp(:,j)
        end do
    end if
    m=m+n
    NRxyz(1)=NRxyz(1)+n
    !===== itmp(1,*)>0 =====
    
    !----- itmp(1,*)<0 -----
    n=0
    do i=1,nr
        if ( itmp(2,i)==0 .and. itmp(3,i)==0 .and. itmp(1,i)<0 ) then
            n=n+1
            work(n)=itmp(1,i)
            itmp2(n)=i
        end if
    end do
    if ( n>0 ) then
        call indexx( n,work,itmp1 )
        do i=1,n
            j=itmp2( itmp1(i) )
            itmp3(:,m+n-i+1)=itmp(:,j)
        end do
    end if
    m=m+n
    NRxyz(2)=NRxyz(2)+n
    !===== itmp(1,*)<0 =====
    
    !----- itmp(2,*)>0 -----
    n=0
    do i=1,nr
        if ( itmp(1,i)==0 .and. itmp(3,i)==0 .and. itmp(2,i)>0 ) then
            n=n+1
            work(n)=itmp(2,i)
            itmp2(n)=i
        end if
    end do
    if ( n>0 ) then
        call indexx( n,work,itmp1 )
        do i=1,n
            j=itmp2( itmp1(i) )
            itmp3(:,m+i)=itmp(:,j)
        end do
    end if
    m=m+n
    NRxyz(3)=NRxyz(3)+n
    !===== itmp(2,*)>0 =====
    
    !----- itmp(2,*)<0 -----
    n=0
    do i=1,nr
        if ( itmp(1,i)==0 .and. itmp(3,i)==0 .and. itmp(2,i)<0 ) then
            n=n+1
            work(n)=itmp(2,i)
            itmp2(n)=i
        end if
    end do
    if ( n>0 ) then
        call indexx( n,work,itmp1 )
        do i=1,n
            j=itmp2( itmp1(i) )
            itmp3(:,m+n-i+1)=itmp(:,j)
        end do
    end if
    m=m+n
    NRxyz(4)=NRxyz(4)+n
    !===== itmp(2,*)<0 =====
    
    !----- itmp(3,*)>0 -----
    n=0
    do i=1,nr
        if ( itmp(1,i)==0 .and. itmp(2,i)==0 .and. itmp(3,i)>0 ) then
            n=n+1
            work(n)=itmp(3,i)
            itmp2(n)=i
        end if
    end do
    if ( n>0 ) then
        call indexx( n,work,itmp1 )
        do i=1,n
            j=itmp2( itmp1(i) )
            itmp3(:,m+i)=itmp(:,j)
        end do
    end if
    m=m+n
    NRxyz(5)=NRxyz(5)+n
    !===== itmp(3,*)>0 =====
    
    !----- itmp(3,*)<0 -----
    n=0
    do i=1,nr
        if ( itmp(1,i)==0 .and. itmp(2,i)==0 .and. itmp(3,i)<0 ) then
            n=n+1
            work(n)=itmp(3,i)
            itmp2(n)=i
        end if
    end do
    if ( n>0 ) then
        call indexx( n,work,itmp1 )
        do i=1,n
            j=itmp2( itmp1(i) )
            itmp3(:,m+n-i+1)=itmp(:,j)
        end do
    end if
    m=m+n
    NRxyz(6)=NRxyz(6)+n
    !===== itmp(3,*)<0 =====
    
    !----- calc. Num2Rank0 -----
    n=maxval( NRxyz )
    if ( allocated(Num2Rank0) ) deallocate( Num2Rank0 )
    allocate( Num2Rank0(n,6) ) ; Num2Rank0=MPI_PROC_NULL
    
    m=0
    do i=1,6
        do j=1,NRxyz(i)
            m=m+1
            i1=itmp3(1,m)+LLp(1,myrank_g)
            i2=itmp3(2,m)+LLp(2,myrank_g)
            i3=itmp3(3,m)+LLp(3,myrank_g)
            k = i1 + i2*np1 + i3*np1*np2
!            k=LLLp(i1,i2,i3)
            Num2Rank0(j,i)=k
        end do
    end do
    !===== calc. Num2Rank0 =====
    
    deallocate( itmp  )
    deallocate( itmp1 )
    deallocate( itmp2 )
    deallocate( itmp3 )
    deallocate( work  )
    deallocate( LLp   )
    
    !----- adjust NRxyz -----
    do i=1,5,2
        n=max( NRxyz(i),NRxyz(i+1) )
        NRxyz(i)=n
        NRxyz(i+1)=n
    end do
    !===== adjust NRxyz =====

write(400+myrank,*) "<<<< prepThreeWayComm"
    return
  END SUBROUTINE prepThreeWayComm

!---------------------------------------------------------------------------------------
  SUBROUTINE threeWayComm( NRxyz,Num2Rank0,SendMap,RecvMap,TarNSend,SbufNL,RbufNL,nz,ib1,ib2,TarIN )
    implicit none
    integer,intent(IN) :: nz,ib1,ib2
    integer,intent(IN) :: NRxyz(1:6)
    integer,allocatable,intent(IN) :: Num2Rank0(:,:)
    integer,allocatable,intent(IN) :: SendMap(:,:),RecvMap(:,:)
    integer,allocatable,intent(INOUT) :: TarNSend(:)
    
#ifdef _DRSDFT_
    real(8),intent(INOUT) :: TarIN(nz,ib1:ib2)
    real(8),allocatable,intent(INOUT) :: SbufNL(:,:),RbufNL(:,:)
    real(8) :: tmp0(nz,ib1:ib2)
#else
    complex(8),intent(INOUT) :: TarIN(nz,ib1:ib2)
    complex(8),allocatable,intent(INOUT) :: SbufNL(:,:),RbufNL(:,:)
    complex(8) :: tmp0(nz,ib1:ib2)
#endif
    integer :: i,j,m,ib,i1,i2
    integer :: irank,jrank
    integer :: nb
    integer :: nreq,istatus(MPI_STATUS_SIZE,512),ireq(512),ierr

#ifdef _SHOWALL_
write(400+myrank,*) ">>>> threeWayComm"
#endif

    nb=ib2-ib1+1

!!$OMP single    
    do i=1,6
        select case ( i )
        case ( 1,3,5 )
            j=i+1
            tmp0(:,:) = TarIN(:,:)
        case ( 2,4,6 )
            j=i-1
        end select
        do m=1,NRxyz(i)
            nreq=0
            irank=Num2Rank0(m,i)
            jrank=Num2Rank0(m,j)
            if ( irank>=0 ) then
                i2=0
                do ib=ib1,ib2
                    do i1=1,TarNSend(irank)
                        i2=i2+1
                        SbufNL(i2,irank)=tmp0(SendMap(i1,irank),ib)
                    end do
                end do
                nreq=nreq+1
                call MPI_ISEND( SbufNL(1,irank),TarNSend(irank)*nb,TYPE_MAIN,irank,1,COMM_GRID,ireq(nreq),ierr )
            end if
            if ( jrank>=0 ) then
                nreq=nreq+1
                call MPI_IRECV( RbufNL(1,jrank),TarNSend(jrank)*nb,TYPE_MAIN,jrank,1,COMM_GRID,ireq(nreq),ierr )
            end if
            call MPI_WAITALL( nreq,ireq,istatus,ierr )
            if ( jrank>=0 ) then
                i2=0
                do ib=ib1,ib2
                    do i1=1,TarNSend(jrank)
                        i2=i2+1
                        TarIN(RecvMap(i1,jrank),ib) = TarIN(RecvMap(i1,jrank),ib) + RbufNL(i2,jrank)
                    end do
                end do
            end if
        end do
    end do
!!$OMP end single

#ifdef _SHOWALL_
write(400+myrank,*) ">>>> threeWayComm"
#endif
    return
  END SUBROUTINE threeWayComm

!---------------------------------------------------------------------------------------

END MODULE ParaRGridComm
