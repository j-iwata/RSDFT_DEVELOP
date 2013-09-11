MODULE subspace_rotv_sl_0_module

  use rgrid_module, only: dV,zdV
  use wf_module, only: unk
  use scalapack_module
  use parallel_module
  use subspace_diag_module
  use array_bound_module, only: ML_0,ML_1,MB_0,MB_1
  use watch_module
  use bcast_module

  implicit none

  PRIVATE
  PUBLIC :: subspace_rotv_sl_0,reset_subspace_rotv_sl_0

#ifdef _DRSDFT_
  real(8),allocatable :: utmp(:,:),utmp2(:,:)
#else
  complex(8),allocatable :: utmp(:,:),utmp2(:,:)
#endif

  integer,allocatable :: array_size(:,:)
  integer,allocatable :: nblock(:)
  integer,allocatable :: mat_indx(:,:,:)
  integer :: comm_sl,myrank_sl,nprocs_sl
  logical :: flag_comm_sl=.true.

CONTAINS

  SUBROUTINE reset_subspace_rotv_sl_0
    if ( allocated(mat_indx)   ) deallocate( mat_indx )
    if ( allocated(nblock)     ) deallocate( nblock )
    if ( allocated(array_size) ) deallocate( array_size )
  END SUBROUTINE reset_subspace_rotv_sl_0

  SUBROUTINE subspace_rotv_sl_0(k,s)
    implicit none
    integer,intent(IN) :: k,s
    integer :: i,i1,i2,ii,n1,n2,i0,j0,j,j1,ns,ne,nn,ms,me,mm,nnn
    integer :: IPCOL,IPROW,iroot1,iroot2,ierr,ML0,ip,MB
    real(8) :: ct0,et0,ct1,et1,ctt(3),ett(3)
    integer :: LLDR0,LLDC0,m,n,iblock
    integer,allocatable :: itmp(:,:)
    integer,allocatable :: idest(:),iorgn(:),ireq(:),istatus(:,:)
    integer :: ndest,norgn,nreq,ig,ib,jk,js

    MB = MB_diag

    if ( NPROW*NPCOL < np_grid*np_band ) then
       if ( .not.allocated(idest) ) then
       allocate( idest(nprocs) ) ; idest=-1
       allocate( iorgn(nprocs) ) ; iorgn=-1
       ndest=0
       norgn=0
       mm=myrank_g+np_grid*myrank_b
       nn=mod(mm,NPROW*NPCOL)
       ii=-1
       do ib=0,np_band-1
       do ig=0,np_grid-1
          ii=ii+1
          do j=0,nprocs-1
             jk=id_class(j,5)
             js=id_class(j,6)
             if ( id_bzsm(jk)+1 <= k .and. k <= id_bzsm(jk)+ir_bzsm(jk) .and. &
                  id_spin(js)+1 <= s .and. s <= id_spin(js)+ir_spin(js) .and. &
                  id_class(j,0) == ig .and. id_class(j,4) == ib ) exit
          end do
          if ( ii < NPROW*NPCOL ) then
             if ( ii == mm ) cycle
             if ( nn == ii ) then
                norgn=norgn+1
                iorgn(norgn)=j
             end if
          else
             if ( mod(ii,NPROW*NPCOL) == mm ) then
                ndest=ndest+1
                idest(ndest)=j
             end if
          end if
       end do ! ig
       end do ! ib
       end if
!       if ( norgn>0 ) then
!          write(*,*) "recv",ndest,norgn,iorgn(1:norgn),myrank,k,s
!       end if
!       if ( ndest>0 ) then
!          write(*,*) "send",ndest,norgn,idest(1:ndest),myrank,k,s
!       end if
!       if ( ndest>0 .and. norgn>0 ) stop "aaaa"
!       call end_mpi_parallel ; stop
       allocate( ireq(nprocs),istatus(MPI_STATUS_SIZE,nprocs) )
       nreq=0
       if ( mm >= NPROW*NPCOL ) then
          Vsub=0.d0
          do i=1,norgn
             nreq=nreq+1
             call mpi_irecv(Vsub,size(Vsub),TYPE_MAIN,iorgn(i),9,mpi_comm_world,ireq(nreq),ierr)
          end do
       else
          do i=1,ndest
             nreq=nreq+1
             call mpi_isend(Vsub,size(Vsub),TYPE_MAIN,idest(i),9,mpi_comm_world,ireq(nreq),ierr)
          end do
       end if
       call mpi_waitall(nreq,ireq,istatus,ierr)
    end if
    if ( .not.allocated(array_size) ) then
       allocate( array_size(2,0:nprocs-1) )
       array_size(1,myrank)=LLD_R
       array_size(2,myrank)=LLD_C
       call mpi_allgather(array_size(1,myrank),2,mpi_integer &
            ,array_size,2,mpi_integer,mpi_comm_world,ierr)
    end if
    if ( .not.allocated(nblock) ) then
       allocate( nblock(0:nprocs-1) ) ; nblock=0
       m=-1
       do i=1,MB,MBSIZE
          m=m+1
          mm=mod(m+NPROW,NPROW)
          n=-1
          do j=1,MB,MBSIZE
             n=n+1
             nn=mod(n+NPCOL,NPCOL)
             nblock(usermap(mm,nn,1))=nblock(usermap(mm,nn,1))+1 
          end do
       end do
    end if
    if ( .not.allocated(mat_indx) ) then
       allocate( mat_indx(8,maxval(nblock),0:nprocs-1) )
       mat_indx(:,:,:)=0
       do ip=0,nprocs-1
          nnn=0
          m=-1
          do i0=1,MB,MBSIZE
             i1=min(i0+MBSIZE-1,MB)
             m=m+1
             mm=mod(m+NPROW,NPROW)
             n=-1
             do j0=1,MB,MBSIZE
                j1=min(j0+MBSIZE-1,MB)
                n=n+1
                nn=mod(n+NPCOL,NPCOL)
                if ( usermap(mm,nn,1) == ip ) then
                   nnn=nnn+1
                   mat_indx(1,nnn,ip)=i0
                   mat_indx(2,nnn,ip)=i1
                   mat_indx(3,nnn,ip)=j0
                   mat_indx(4,nnn,ip)=j1
                end if
             end do
          end do
          if ( nnn /= nblock(ip) ) stop "error!!!"
          mat_indx(5,1,ip)=1
          mat_indx(6,1,ip)=mat_indx(5,1,ip)+mat_indx(2,1,ip)-mat_indx(1,1,ip)
          mat_indx(7,1,ip)=1
          mat_indx(8,1,ip)=mat_indx(7,1,ip)+mat_indx(4,1,ip)-mat_indx(3,1,ip)
          do i=2,nblock(ip)
             if ( mat_indx(1,i,ip) == mat_indx(1,i-1,ip) ) then
                mat_indx(5,i,ip) = mat_indx(5,i-1,ip)
                mat_indx(6,i,ip) = mat_indx(6,i-1,ip)
                mat_indx(7,i,ip) = mat_indx(8,i-1,ip)+1
                mat_indx(8,i,ip) = mat_indx(7,i,ip)+mat_indx(4,i,ip)-mat_indx(3,i,ip)
             else
                mat_indx(5,i,ip) = mat_indx(6,i-1,ip)+1
                mat_indx(6,i,ip) = mat_indx(5,i,ip)+mat_indx(2,i,ip)-mat_indx(1,i,ip)
                mat_indx(7,i,ip) = 1
                mat_indx(8,i,ip) = mat_indx(7,i,ip)+mat_indx(4,i,ip)-mat_indx(3,i,ip)
             end if
          end do
       end do
    end if
    if ( flag_comm_sl ) then
       flag_comm_sl=.false.
       if ( NPROW*NPCOL == np_grid*np_band ) then
          comm_sl=comm_grid
       else
          n =np_grid*np_band
          nn=NPROW*NPCOL
          mm=0
          do i0=0,n-1,nn
             i1=min(i0+nn-1,n-1)
             m=myrank_g+np_grid*myrank_b
             if ( i0 <= m .and. m <= i1 ) mm=i1+np_grid*np_band*myrank_k+np_grid*np_band*np_bzsm*myrank_s
          end do
          call mpi_comm_split(mpi_comm_world,mm,myrank,comm_sl,ierr)
          call mpi_comm_rank(comm_sl,myrank_sl,ierr)
          call mpi_comm_size(comm_sl,nprocs_sl,ierr)
          allocate( itmp(0:NPROW-1,0:NPCOL-1) ) ; itmp=0
          do j=0,NPCOL-1
          do i=0,NPROW-1
             if ( usermap(i,j,1) == myrank ) then
                itmp(i,j) = myrank_sl
             end if
          end do
          end do
          call mpi_allreduce(itmp,usermap(0,0,2),NPROW*NPCOL,MPI_INTEGER,MPI_SUM,comm_sl,ierr)
          deallocate( itmp )
       end if
    end if

    ctt(:)=0.d0
    ett(:)=0.d0

    n1  = ML_0
    n2  = ML_1
    ML0 = ML_1-ML_0+1

    NBLK2 = maxval(ircnt)

    allocate( utmp(NBLK2,MB_0:MB_1) )
    utmp=zero

    do i=1,maxval(ircnt),NBLK2
       i1=n1+i-1
       i2=min(i1+NBLK2-1,n2)
       ii=i2-i1+1

       utmp(:,:)=zero

       do IPCOL=0,NPCOL-1
       do IPROW=0,NPROW-1

          iroot1 = usermap(IPROW,IPCOL,1)
          iroot2 = usermap(IPROW,IPCOL,2)

          LLDR0 = array_size(1,iroot1)
          LLDC0 = array_size(2,iroot1)

          allocate( utmp2(LLDR0,LLDC0) )

          if ( iroot1 == myrank ) utmp2(:,:) = Vsub(:,:)

          call watch(ct0,et0)

!          call mpi_bcast(utmp2,LLDR0*LLDC0,TYPE_MAIN,iroot2,comm_sl,ierr)
#ifdef _DRSDFT_
          call d_rsdft_bcast(utmp2,LLDR0*LLDC0,TYPE_MAIN,iroot2,comm_sl,ierr)
#else
          call z_rsdft_bcast(utmp2,LLDR0*LLDC0,TYPE_MAIN,iroot2,comm_sl,ierr)
#endif
          call watch(ct1,et1) ; ctt(1)=ctt(1)+ct1-ct0 ; ett(1)=ett(1)+et1-et0

          if ( ii>0 ) then
             do iblock=1,nblock(iroot1)
                ms=mat_indx(1,iblock,iroot1)
                mm=mat_indx(2,iblock,iroot1)-mat_indx(1,iblock,iroot1)+1
                ns=mat_indx(3,iblock,iroot1)
                nn=mat_indx(4,iblock,iroot1)-mat_indx(3,iblock,iroot1)+1
                i0=mat_indx(5,iblock,iroot1)
                j0=mat_indx(7,iblock,iroot1)
#ifdef _DRSDFT_
                call dgemm('N','N',ii,nn,mm,one,unk(i1,ms,k,s) &
                     ,ML0,utmp2(i0,j0),LLDR0,one,utmp(1,ns),NBLK2)
#else
                call zgemm('N','N',ii,nn,mm,one,unk(i1,ms,k,s) &
                     ,ML0,utmp2(i0,j0),LLDR0,one,utmp(1,ns),NBLK2)
#endif
             end do ! iblock
          end if

          call watch(ct0,et0) ; ctt(2)=ctt(2)+ct0-ct1 ; ett(2)=ett(2)+et0-et1

          deallocate( utmp2 )

       end do ! IPROW
       end do ! IPCOL

       if ( ii>0 ) then
          unk(i1:i2,MB_0:MB_1,k,s)=utmp(1:ii,MB_0:MB_1)
       end if

    end do ! ii

    deallocate( utmp )

    if ( disp_switch_parallel ) then
       write(*,*) "time(bcast)=",ctt(1),ett(1)
       write(*,*) "time(gemm)=",ctt(2),ett(2)
    end if

  END SUBROUTINE subspace_rotv_sl_0

END MODULE subspace_rotv_sl_0_module
