MODULE GScG
  use rgrid_module, only: dV,zdV
  use wf_module, only: unk,Sunk
  use array_bound_module, only: ML_0,ML_1,MB,MB_0
  use parallel_module, only: COMM_GRID,COMM_BAND,ir_band,id_band,id_class,myrank,np_band,myrank_b
  use RealComplex, only: zero,one,TYPE_MAIN,TRANSA,TRANSB
  use InnerProduct
  implicit none
  include 'mpif.h'
  PRIVATE
  PUBLIC :: GramSchmidtG

CONTAINS

!---------------------------------------------------------------------------------------
  ! Gram-Schmidt orthogonalization
  ! ( Takahashi, Block cyclic )
  SUBROUTINE GramSchmidtG(n0,n1,k,s,NBLK,NBLK1)
    implicit none
    integer,intent(IN) :: n0,n1,k,s
    integer,intent(INOUT) :: NBLK,NBLK1
    
    integer :: irank_b,ns,ne,ms,me,n,ML0,ierr
    integer :: nbss,k1,ib,NBAND_BLK,ncycle,mrnk
    integer,allocatable :: ir(:),id(:)

    ! not used???? or in recursive subroutine
    integer :: nn,cc1,cc2
    integer :: nn1,nn2
    integer :: mm1,mm2
write(400+myrank,*) ">>>>>>>> GramSchmidtG"

    ML0 = ML_1 - ML_0 + 1
    mrnk = id_class(myrank,4)

    if ( NBLK  == 0 ) NBLK =MB
    if ( NBLK1 == 0 ) NBLK1=4

    allocate( ir(0:np_band-1),id(0:np_band-1) )

    ir(0:np_band-1)=ir_band(0:np_band-1)*ML0
    id(0:np_band-1)=id_band(0:np_band-1)*ML0

    NBAND_BLK=NBLK
    ncycle=(MB-1)/NBAND_BLK+1

    call MPI_ALLGATHERV( unk(ML_0,MB_0,k,s),ir(mrnk),TYPE_MAIN,unk(ML_0,1,k,s),ir,id,TYPE_MAIN,COMM_BAND,ierr )
    do k1=1,ncycle
       irank_b = mod(k1-1,np_band)
       ns = NBAND_BLK*(k1-1) + 1
       ne = min(ns+NBAND_BLK-1,MB)

       allocate( Sunk(n0:n1,ns:ne) )
write(400+myrank,*) "GramSchmidtG",k1
       call get_Sunk_Mat( n0,n1,ns,ne,k,s )
       if ( id_class(myrank,4)==irank_b ) then
write(400+myrank,*) "GramSchmidtG",k1,"Rec"
          call GramSchmidtGSub(ns,ne,ns,ne,NBLK,k,s,NBLK,NBLK1)
       end if
       n=ML0*(ne-ns+1)
       call MPI_BCAST( unk(ML_0,ns,k,s),n,TYPE_MAIN,irank_b,COMM_BAND,ierr )
       if ( ns <= MB-NBAND_BLK ) then
          do ib=1,(ncycle-1)/np_band+1
             nbss = (ib-1)*np_band + myrank_b + 1
             if ( nbss<=ncycle .and. nbss>= k1+1 ) then
                ms=NBAND_BLK*(nbss-1)+1
                me=min(ms+NBAND_BLK-1,MB)
                if ( ms<=me ) then
write(400+myrank,*) "GramSchmidtG",k1,ib,"Rec"
                   call GramSchmidtGSub(ms,me,ns,ne,NBLK,k,s,NBLK,NBLK1)
                end if
             end if
          end do
       end if
       deallocate(Sunk)
    end do ! k1
    deallocate( id,ir )
write(400+myrank,*) "<<<<<<<< GramSchmidtG"

    return
  END SUBROUTINE GramSchmidtG

!---------------------------------------------------------------------------------------
  ! Gram-Schmidt orthogonalization
  RECURSIVE SUBROUTINE GramSchmidtGSub(mm1,mm2,nn1,nn2,MBLK,k,s,NBLK,NBLK1)
    implicit none
    integer,intent(IN) :: mm1,mm2,nn1,nn2,MBLK,k,s
    integer,intent(IN) :: NBLK,NBLK1
    integer :: n,ns,ne,m,ms,me,m1,m2,mm,nn,MBLKH,ierr
    integer :: ML0,i
#ifdef _DRSDFT_
    real(8),allocatable :: utmp2(:,:),vtmp2(:,:),utmp(:),vtmp(:)
    real(8) :: c,d
#else
    complex(8),allocatable :: utmp2(:,:),vtmp2(:,:),utmp(:),vtmp(:)
    complex(8) :: c,d
#endif

write(400+myrank,*) ">>>>>> GramSchmidtGSub"
    ML0 = ML_1-ML_0+1

    do ms=mm1,mm2,MBLK
       me=min(ms+MBLK-1,mm2)
       mm=me-ms+1
       do ns=nn1,nn2,MBLK
          ne=min(ns+MBLK-1,nn2)
          ne=min(ne,me-1)
          nn=ne-ns+1
          if ( nn<=0 ) cycle
          if ( ms>=ne+1 ) then
write(400+myrank,*) "GramSchmidtGSub IF1"
             allocate( utmp2(ns:ne,ms:me),vtmp2(ns:ne,ms:me) )
#ifdef _DRSDFT_
             call dgemm(TRANSA,TRANSB,nn,mm,ML0, -dV,unk(ML_0,ns,k,s),ML0,Sunk(ML_0,ms),ML0,zero,utmp2,nn)
#else
             call zgemm(TRANSA,TRANSB,nn,mm,ML0,-zdV,unk(ML_0,ns,k,s),ML0,Sunk(ML_0,ms),ML0,zero,utmp2,nn)
#endif

             call MPI_ALLREDUCE( utmp2,vtmp2,nn*mm,TYPE_MAIN,MPI_SUM,COMM_GRID,ierr )

#ifdef _DRSDFT_
             call dgemm(TRANSB,TRANSB,ML0,mm,nn,one,unk(ML_0,ns,k,s),ML0,vtmp2,nn,one,unk(ML_0,ms,k,s),ML0)
#else
             call zgemm(TRANSB,TRANSB,ML0,mm,nn,one,unk(ML_0,ns,k,s),ML0,vtmp2,nn,one,unk(ML_0,ms,k,s),ML0)
#endif

             deallocate( vtmp2,utmp2 )
             if ( ms==ne+1 ) then
write(400+myrank,*) "GramSchmidtGSub IF2"
                call get_gSf(unk(ML_0,ms,k,s),unk(ML_0,ms,k,s),ML_0,ML_1,k,d,0)
                call MPI_ALLREDUCE( d,c,1,mpi_real8,MPI_SUM,COMM_GRID,ierr )
                ! ?????????????????????????????? dV????
                !						c=1.d0/sqrt(c)
                c=1.d0/sqrt(c*dV)
                unk(ML_0:ML_1,ms,k,s) = c*unk(ML_0:ML_1,ms,k,s)
             end if

          else if ( mm<=NBLK1 ) then
             allocate( utmp(NBLK1),vtmp(NBLK1) )
             do m=ms,me
                n = min(m-1,ne)
                if ( n-ns+1>0 ) then
#ifdef _DRSDFT_
                   call dgemv(TRANSA,ML0,n-ns+1, -dV,unk(ML_0,ns,k,s),ML0,Sunk(ML_0,m),1,zero,utmp,1)
#else
                   call zgemv(TRANSA,ML0,n-ns+1,-zdV,unk(ML_0,ns,k,s),ML0,Sunk(ML_0,m),1,zero,utmp,1)
#endif
                   call mpi_allreduce(utmp,vtmp,n-ns+1,TYPE_MAIN,mpi_sum,comm_grid,ierr)

#ifdef _DRSDFT_
                   call dgemv(TRANSB,ML0,n-ns+1,one,unk(ML_0,ns,k,s),ML0,vtmp,1,one,unk(ML_0,m,k,s),1)
#else
                   call zgemv(TRANSB,ML0,n-ns+1,one,unk(ML_0,ns,k,s),ML0,vtmp,1,one,unk(ML_0,m,k,s),1)
#endif
                end if
                if ( m==1 .or. (n==m-1 .and. m/=ns) ) then
                   call get_gSf(unk(ML_0,m,k,s),unk(ML_0,m,k,s),ML_0,ML_1,k,d,0)
                   call mpi_allreduce(d,c,1,mpi_real8,mpi_sum,comm_grid,ierr)
                   !							c=1.d0/sqrt(c)
                   c=1.d0/sqrt(c*dV)
                   unk(ML_0:ML_1,m,k,s)=c*unk(ML_0:ML_1,m,k,s)
                end if
             end do ! m
             deallocate( vtmp,utmp )

          else
             MBLKH=max(MBLK/2,NBLK1)
write(400+myrank,*) "GramSchmidtGSub Rec"
             call GramSchmidtGSub(ms,me,ns,ne,MBLKH,k,s,NBLK,NBLK1)
          end if
       end do ! ns
    end do ! ms
write(400+myrank,*) "<<<<<< GramSchmidtGSub"

    return
  END SUBROUTINE GramSchmidtGSub

END MODULE GScG
