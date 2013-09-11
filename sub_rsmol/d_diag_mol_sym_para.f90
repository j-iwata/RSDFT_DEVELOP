!--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3--

      SUBROUTINE diag_mol_sym
      use global_variables
      implicit none

      integer :: i,ib1,ib2,j,k,n,n1,n2,ML0,ierr,is,ik,jmax,imax
      integer :: mm0,mm1,mml,nn0,nn1,nnl,j1,j2
      integer :: IPROW,IPCOL,iroot,iroot2,mrnk,nrecv_me,nsend_me
      integer,allocatable :: irecv_me(:,:),isend_me(:,:)
      integer,allocatable :: ir(:),id(:)
      integer :: istatus(MPI_STATUS_SIZE,123),nreq,ireq(123),itag
      integer :: iblk,jblk,nblk_0,nblk_1,nblk_n,mblk_0,mblk_1,mblk_n
      integer :: ns,ne,ms,me,ns_0,ns_1,nns,nne,nnn
      integer :: n3,m,i0,j0,inn,inn1,i1,i2,ii
      complex(8),allocatable :: zwork(:)
      complex(8) :: ztmp
      integer :: mm,nn,nn2
      integer :: mt,nt,mm2,mmm0,mmm1,mmms,nnn0,nnn1,nnns
      integer :: itmp(1),mms,mme,MBLKH,MBLK
      real(8) :: tmp0,tmp
      integer :: irank,s

      real(8) :: ctime0,ctime1,etime0,etime1,mem,memax
      real(8) :: ct(10),et(10),ct0,ct1,ct2,ct3,et0,et1,et2,et3
      real(8) :: s0(37),s1(37,3),nop(10)
      real(8) :: nop1,nop3,nop1t,nop3t

      if ( MB==1 ) return

      call prep_scalapack

      call watch(ctime0,etime0)

      ct(:)=0.d0 ; et(:)=0.d0
      ctime_hpsi=0.d0 ; etime_hpsi=0.d0
      mem=0.d0 ; memax=0.d0
      nop(:)=0.d0

      n1   = idisp(myrank)+1
      n2   = idisp(myrank)+ircnt(myrank)
      ML0  = ircnt(myrank)
      mrnk = id_class(myrank,4)

!
! --- preparation for ScaLAPACK & band-parallel ---
!

      allocate( irecv_me(99,0:8),isend_me(99,0:8) )
      mem=mem+bsintg*99*9*2 ; memax=max(memax,mem)

      allocate( id(0:np_band-1),ir(0:np_band-1) ) ; id=0 ; ir=0
      mem=mem+bsintg*np_band*2 ; memax=max(memax,mem)

      nrecv_me      = 0
      irecv_me(:,:) =-1
      nsend_me      = 0
      isend_me(:,:) =-1

      id(0:np_band-1) = id_band(0:np_band-1)*ML0
      ir(0:np_band-1) = ir_band(0:np_band-1)*ML0

      MBLK = min( NBLK,MBSIZE,NBSIZE )

!
! --- matrix element ---
!
!- allocate ---------------------------------------------
!         allocate( utmp2(LLD_R,LLD_C)   ) ; utmp2=zero
!         allocate( psi_tmp(n1:n2,MBLK) )  ; psi_tmp=zero
!         mem=mem+bdmain*(LLD_R*LLD_C)
!         mem=mem+bdmain*ML0*MBLK
!         memax=max(memax,mem)
!--------------------------------------------------------

!
! --- prep_scalapack ---
!
      NPCOL0 = np_grid
      NPCOL  = np_grid
      NPROW  = 1
      do i=2,np_grid
         j=NPCOL0/i
         if ( i*j==NPCOL0 .and. i<=j ) then
            NPCOL=j
            NPROW=i
         end if
      end do

      allocate( usermap(0:NPROW-1,0:NPCOL-1,2) ) usermap=0

      n=-1
      do is=0,np_2d(6)-1
      do ik=0,np_2d(5)-1
         m=-1
         do ib=0,np_2d(4)-1
            l=-1
            do i3=0,np_2d(3)-1
            do i2=0,np_2d(2)-1
            do i1=0,np_2d(1)-1
               n=n+1
               m=m+1
               l=l+1
               i=mod(l+NPROW,NPROW)
               j=l/NPROW
               if ( id_class(myrank,5)==ik .and. id_class(myrank,6)==is .and. id_class(myrank,4)==ib ) then
                  usermap(i,j,1)=n
                  usermap(i,j,2)=l
               end if
            end do
            end do
            end do
         end do
      end do
      end do

      if ( .not.FLAG_NODE_TEST ) then
         call blacs_get(0,0,ICTXT0)
         ICTXT = ICTXT0
         call blacs_gridmap(ICTXT,usermap(0,0,1),NPROW,NPROW,NPCOL)
         call blacs_gridinfo(ICTXT,NPROW,NPCOL,MYROW,MYCOL)
      end if

      do s=MSP_0,MSP_1
      do k=MBZ_0,MBZ_1

         do ir=Nir_0,Nir_1

            N_excited=iedim(ir)
            nd=irdim(ir)
            fac=dV*nsym/dble(nd)


            i=(N_excited+NPROW-1)/NPROW
            j=(N_excited+NPCOL-1)/NPCOL
            MBSIZE=min(i,j)
            NBSIZE=MBSIZE

            if(DISP_SWITCH)then
               write(*,*) "NPROW,NPCOL,MBSIZE,NBSIZE"
               write(*,*) NPROW,NPCOL,MBSIZE,NBSIZE
            end if

            if ( NBSIZE*NPCOL/=N_excited ) then
               write(*,*) "NBSIZE*NPCOL/=MB!"
               n=max( NBSIZE*NPCOL, N_excited )
               n=min( n, (NBSIZE+1)*NPCOL )
               write(*,*) "recommended value for N_excited =",n
               goto 900
            end if


            j0=0

            do ns=1,N_excited,MBLK
               ne=min(ns+MBLK-1,N_excited)
               nn=ne-ns+1

               do ie=ns,ne
                  do id=1,nd
                     n=irlabel2n(ir,id,ie)
                     vtmp2(n1:n2,id)=unk(n1:n2,n,k,s)
                  end do
                  call hpsi_mol_sym(k,s,vtmp2,wtmp2,n1,n2,nd,ir)
                  do id=1,nd
                     n=(ie-1)*nd+id
                     psi_tmp(n1:n2,n)=wtmp2(n1:n2,id)*wkgrp(n1:n2)
                  end do
               end do

            if ( j0<LLD_C ) i0=0

            if ( j0==0 ) then
               do ms=1,ns-1,MBSIZE
                  me=min(ms+MBSIZE-1,ns-1)
                  mm=me-ms+1
                  IPROW = mod( (ms-1)/MBSIZE, NPROW )
                  iroot = usermap(IPROW,IPCOL,1)
                  if ( iroot==myrank ) i0=i0+mm
               end do
            end if

            do ms=ns,ne,MBLK
               me=min(ms+MBLK-1,ne)
               mm=me-ms+1
               IPROW = mod( (ms-1)/MBSIZE, NPROW )
               iroot = usermap(IPROW,IPCOL,1)
               iroot2= usermap(IPROW,IPCOL,2)
               if ( FLAG_NODE_TEST ) iroot2=myrank_g
               allocate( vtmp2(ms:me,ns:ne),wtmp2(ms:me,ns:ne) ) ; vtmp2=zero ; wtmp2=zero
               mem=mem+bdmain*mm*nn*2 ; memax=max(memax,mem)
               MBLKH=max(MBLK/2,NBLK1)
               call diag_2d_sub(ms,me,ns,ne,MBLKH,ns,ms,me,ct,et,nop,k,s)

               call watch(ct2,et2)

               call mpi_reduce(vtmp2,wtmp2,mm*nn,TYPE_MAIN,mpi_sum,iroot2,comm_grid,ierr)

               call watch(ct3,et3)
               ct(4)=ct(4)+ct3-ct2 ; et(4)=et(4)+et3-et2

               if ( iroot==myrank ) then
                  utmp2(i0+1:i0+mm,j0+1:j0+nn)=wtmp2(ms:me,ns:ne)
                  i0=i0+mm
               end if
               deallocate( vtmp2,wtmp2 ) ; mem=mem-bdmain*mm*nn*2
            end do

            if ( ne+1<=MB_1 ) then
               do ms=ne+1,MB_1,MBLK
                  me=min(ms+MBLK-1,MB_1)
                  mm=me-ms+1
                  IPROW = mod( (ms-1)/MBSIZE, NPROW )
                  iroot = usermap(IPROW,IPCOL,1)
                  iroot2= usermap(IPROW,IPCOL,2)
                  if ( FLAG_NODE_TEST ) iroot2=myrank_g
                  allocate( vtmp2(ms:me,ns:ne),wtmp2(ms:me,ns:ne) ) ; vtmp2=zero ; wtmp2=zero
                  mem=mem+bdmain*mm*nn*2 ; memax=max(memax,mem)
                  call dgemm(TRANSA,TRANSB,mm,nn,ML0,zdV,unk(n1,ms,k,s),ML0,psi_tmp(n1,1),ML0,zero,vtmp2(ms,ns),mm)
                  nop(1)=nop(1)+sop1*mm*nn*ML0

                  call watch(ct2,et2)

                  call mpi_reduce(vtmp2,wtmp2,mm*nn,TYPE_MAIN,mpi_sum,iroot2,comm_grid,ierr)

                  call watch(ct3,et3)
                  ct(4)=ct(4)+ct3-ct2 ; et(4)=et(4)+et3-et2

                  if ( iroot==myrank ) then
                     utmp2(i0+1:i0+mm,j0+1:j0+nn)=wtmp2(ms:me,ns:ne)
                     i0=i0+mm
                  end if
                  deallocate( wtmp2,vtmp2 ) ; mem=mem-bdmain*mm*nn*2
               end do
            end if

            nns = max(ns,MB_1-mat_block(mrnk,2)+1)
            nnn = ne-nns+1
            mme = min(MB,MB_1+mat_block(mrnk,1))

            do ms=MB_1+1,mme,MBLK
               me=min(ms+MBLK-1,mme)
               mm=me-ms+1
               IPROW = mod( (ms-1)/MBSIZE, NPROW )
               iroot = usermap(IPROW,IPCOL,1)
               iroot2= usermap(IPROW,IPCOL,2)
               if ( FLAG_NODE_TEST ) iroot2=myrank_g
               j1=j0
               if ( ns<nns ) then
                  if ( iroot==myrank ) then
                     i=mod( (ns-1)/MBSIZE, NPROW )
                     j=mod( (ms-1)/NBSIZE, NPCOL )
                     n=usermap(i,j,1)
                     if ( n==myrank ) goto 900
                     nrecv_me=nrecv_me+1
                     irecv_me(nrecv_me,0)=n
                     irecv_me(nrecv_me,1)=ms
                     irecv_me(nrecv_me,2)=me
                     irecv_me(nrecv_me,3)=ns
                     irecv_me(nrecv_me,4)=min(ne,nns-1)
                     irecv_me(nrecv_me,5)=i0+1
                     irecv_me(nrecv_me,6)=i0+mm
                     irecv_me(nrecv_me,7)=j1+1
                     irecv_me(nrecv_me,8)=j1+min(ne,nns-1)-ns+1
                     j1=j1+min(ne,nns-1)-ns+1
                  end if
               end if
               if ( nns<=ne ) then
                  allocate( vtmp2(ms:me,nns:ne),wtmp2(ms:me,nns:ne) ) ; vtmp2=zero ; wtmp2=zero
                  mem=mem+bdmain*mm*nnn*2 ; memax=max(memax,mem)
                  call dgemm(TRANSA,TRANSB,mm,nnn,ML0,zdV,unk(n1,ms,k,s),ML0,psi_tmp(n1,nns-ns+1),ML0,zero,vtmp2,mm)
                  nop(1)=nop(1)+sop1*mm*nnn*ML0

                  call watch(ct2,et2)

                  call mpi_reduce(vtmp2,wtmp2,mm*nnn,TYPE_MAIN,mpi_sum,iroot2,comm_grid,ierr)

                  call watch(ct3,et3)
                  ct(4)=ct(4)+ct3-ct2 ; et(4)=et(4)+et3-et2

                  if ( iroot==myrank ) then
                     utmp2(i0+1:i0+mm,j1+1:j1+nnn)=wtmp2(ms:me,nns:ne)
                  end if
                  deallocate( wtmp2,vtmp2 ) ; mem=mem-bdmain*mm*nnn*2
               end if
               if ( iroot==myrank ) then
                  i0=i0+mm
               end if
            end do

            if ( mme+1<=MB ) then
               do ms=mme+1,MB,MBLK
                  me=min(ms+MBLK-1,MB)
                  mm=me-ms+1
                  IPROW = mod( (ms-1)/MBSIZE, NPROW )
                  iroot = usermap(IPROW,IPCOL,1)
                  if ( iroot==myrank ) then
                     i=mod( (ns-1)/MBSIZE, NPROW )
                     j=mod( (ms-1)/NBSIZE, NPCOL )
                     n=usermap(i,j,1)
                     if ( n==myrank ) goto 900
                     nrecv_me=nrecv_me+1
                     irecv_me(nrecv_me,0)=n
                     irecv_me(nrecv_me,1)=ms
                     irecv_me(nrecv_me,2)=me
                     irecv_me(nrecv_me,3)=ns
                     irecv_me(nrecv_me,4)=ne
                     irecv_me(nrecv_me,5)=i0+1
                     irecv_me(nrecv_me,6)=i0+mm
                     irecv_me(nrecv_me,7)=j0+1
                     irecv_me(nrecv_me,8)=j0+nn
                     i0=i0+mm
                  end if
               end do
            end if

            if ( MB+1<=MB_1+mat_block(mrnk,1) ) then
               i0=0
               do mms=MB+1,MB_1+mat_block(mrnk,1),MBLK
                  mme=min(mms+MBLK-1,MB_1+mat_block(mrnk,1))
                  ms=mod(mms+MB-1,MB)+1
                  me=mod(mme+MB-1,MB)+1
                  mm=me-ms+1
                  IPROW = mod( (ms-1)/MBSIZE,NPROW )
                  iroot = usermap(IPROW,IPCOL,1)
                  iroot2= usermap(IPROW,IPCOL,2)
                  if ( FLAG_NODE_TEST ) iroot2=myrank_g
                  allocate( vtmp2(ms:me,nns:ne),wtmp2(ms:me,nns:ne) ) ; vtmp2=zero ; wtmp2=zero
                  mem=mem+bdmain*mm*nnn*2 ; memax=max(memax,mem)
                  call dgemm(TRANSA,TRANSB,mm,nnn,ML0,zdV,unk(n1,ms,k,s),ML0,psi_tmp(n1,nns-ns+1),ML0,zero,vtmp2,mm)
                  nop(1)=nop(1)+sop1*mm*nn*ML0

                  call watch(ct2,et2)

                  call mpi_reduce(vtmp2,wtmp2,mm*nnn,TYPE_MAIN,mpi_sum,iroot2,comm_grid,ierr)

                  call watch(ct3,et3)
                  ct(4)=ct(4)+ct3-ct2 ; et(4)=et(4)+et3-et2

                  if ( iroot==myrank ) then
                     i=mod( (nns-1)/MBSIZE, NPROW )
                     j=mod( (ms-1)/NBSIZE, NPCOL )
                     n=usermap(i,j,1)
                     if ( n==myrank ) goto 900
                     nsend_me=nsend_me+1
                     isend_me(nsend_me,0)=n
                     isend_me(nsend_me,1)=ms
                     isend_me(nsend_me,2)=me
                     isend_me(nsend_me,3)=nns
                     isend_me(nsend_me,4)=ne
                     isend_me(nsend_me,5)=i0+1
                     isend_me(nsend_me,6)=i0+mm
                     isend_me(nsend_me,7)=j0+1
                     isend_me(nsend_me,8)=j0+nnn
                     utmp2(i0+1:i0+mm,j0+1:j0+nnn)=wtmp2(ms:me,nns:ne)
                     i0=i0+mm
                  end if                        
                  deallocate( wtmp2,vtmp2 ) ; mem=mem-bdmain*mm*nnn*2
               end do
            end if

            if ( any(usermap(0:NPROW-1,IPCOL,1)==myrank) ) j0=j0+nn

         end do ! ns

         deallocate( psi_tmp ) ; mem=mem-bdmain*ML0*MBLK

         call watch(ct2,et2)

!
! --- ME-send ---
!
         n=0
         do i=1,nsend_me
            m=(isend_me(i,2)-isend_me(i,1)+1)*(isend_me(i,4)-isend_me(i,3)+1)
            n=max(m,n)
         end do
         allocate( vtmp2(n,nsend_me) )
         mem=mem+bdmain*n*nsend_me ; memax=max(memax,mem)
         n=0
         do i=1,nrecv_me
            m=(irecv_me(i,2)-irecv_me(i,1)+1)*(irecv_me(i,4)-irecv_me(i,3)+1)
            n=max(m,n)
         end do
         allocate( wtmp2(n,nrecv_me) )
         mem=mem+bdmain*n*nrecv_me ; memax=max(memax,mem)

         if ( nsend_me>0 .or. nrecv_me>0 ) then
            nreq=0
            do i=1,nsend_me
               n =isend_me(i,0)
               ms=isend_me(i,1)
               me=isend_me(i,2)
               ns=isend_me(i,3)
               ne=isend_me(i,4)
               i1=isend_me(i,5)
               i2=isend_me(i,6)
               j1=isend_me(i,7)
               j2=isend_me(i,8)
               j=0
               do i0=i1,i2
               do j0=j1,j2
                  j=j+1
                  vtmp2(j,i)=utmp2(i0,j0)
               end do
               end do
               itag=ns+10*ne+100*ms+1000*me
               nreq=nreq+1
               if ( .not.FLAG_NODE_TEST ) then
                  call mpi_isend(vtmp2(1,i),j,TYPE_MAIN,n,itag,mpi_comm_world,ireq(nreq),ierr)
               end if
            end do
            do i=1,nrecv_me
               n =irecv_me(i,0)
               ms=irecv_me(i,1)
               me=irecv_me(i,2)
               ns=irecv_me(i,3)
               ne=irecv_me(i,4)
               i1=irecv_me(i,5)
               i2=irecv_me(i,6)
               j1=irecv_me(i,7)
               j2=irecv_me(i,8)
               j=(me-ms+1)*(ne-ns+1)
               itag=ms+10*me+100*ns+1000*ne
               nreq=nreq+1
               if ( .not.FLAG_NODE_TEST ) then
                  call mpi_irecv(wtmp2(1,i),j,TYPE_MAIN,n,itag,mpi_comm_world,ireq(nreq),ierr)
               end if
            end do

            if ( .not.FLAG_NODE_TEST ) then
               call mpi_waitall(nreq,ireq,istatus,ierr)
            end if

            do i=1,nrecv_me
               i1=irecv_me(i,5)
               i2=irecv_me(i,6)
               j1=irecv_me(i,7)
               j2=irecv_me(i,8)
               if ( TYPE_MAIN==mpi_complex16 ) then
                  j=0
                  do j0=j1,j2
                  do i0=i1,i2
                     j=j+1
                     ztmp=utmp2(i0,j0)
                     ztmp=conjg(ztmp)
                     utmp2(i0,j0)=ztmp
                  end do
                  end do
               else
                  j=0
                  do j0=j1,j2
                  do i0=i1,i2
                     j=j+1
                     utmp2(i0,j0)=wtmp2(j,i)
                  end do
                  end do
               end if
            end do

         end if

         call watch(ct3,et3)
         ct(5)=ct(5)+ct3-ct2 ; et(5)=et(5)+et3-et2

         mem=mem-bdmain*(size(vtmp2)+size(wtmp2))
         deallocate( wtmp2,vtmp2 )

         call watch(ct1,et1)
         ct(6)=ct(6)+ct1-ct0
         et(6)=et(6)+et1-et0

!
! --- Solve eigenvalue problem (SCALAPACK) ---
!

         allocate( vtmp2(LLD_R,LLD_C) ) ; vtmp2=zero
         mem=mem+bdmain*(LLD_R*LLD_C) ; memax=max(memax,mem)

         call ev_solver('L',k,s,mem,memax)

         deallocate( utmp2 )
         mem=mem-bdmain*LLD_R*LLD_C

         call watch(ct0,et0)
         ct(7)=ct(7)+ct0-ct1
         et(7)=et(7)+et0-et1

!
! --- Eigen vectors ---
!

         allocate( psi_tmp(NBLK2,MB_0:MB_1) ) ; psi_tmp=zero
         mem=mem+bdmain*NBLK2*(MB_1-MB_0+1) ; memax=max(mem,memax)

         do i=1,maxval(ircnt),NBLK2
            i1=n1+i-1
            i2=min(i1+NBLK2-1,n2)
            ii=i2-i1+1

            psi_tmp(:,:)=zero

            j0=0
         do ns=MB_0,MB_1,NBSIZE
            ne=min(ns+NBSIZE-1,MB_1)
            nn=ne-ns+1

            IPCOL=mod( (ns-1)/NBSIZE,NPCOL )

            i0=0
         do ms=1,MB,MBSIZE
            me=min(ms+MBSIZE-1,MB)
            mm=me-ms+1

            IPROW=mod( (ms-1)/MBSIZE,NPROW )

            iroot  = usermap(IPROW,IPCOL,1)
            iroot2 = usermap(IPROW,IPCOL,2)

            if ( mm<1 .or. nn<1 ) cycle

            allocate( utmp2(ms:me,ns:ne) )
            mem=mem+bdmain*mm*nn ; memax=max(mem,memax)

            if ( iroot==myrank ) then
               utmp2(ms:me,ns:ne)=vtmp2(i0+1:i0+mm,j0+1:j0+nn)
               i0=i0+mm
            end if

            if ( FLAG_NODE_TEST ) then
               utmp2(:,:)=zero
               do j2=ns,ne
                  do j1=ms,me
                     if ( j1==j2 ) utmp2(j1,j2)=one
                  end do
               end do
            else

               call watch(ct2,et2)

               call mpi_bcast( utmp2(ms,ns),mm*nn,TYPE_MAIN,iroot2,comm_grid,ierr )

               call watch(ct3,et3)
               ct(8)=ct(8)+ct3-ct2 ; et(8)=et(8)+et3-et2

            end if

            if ( ii>0 ) then
               call dgemm(TRANSB,TRANSB,ii,nn,mm,one,unk(i1,ms,k,s),ML0,utmp2(ms,ns),mm,one,psi_tmp(1,ns),NBLK2)
               nop(9)=nop(9)+sop1*nn*mm*ii
            end if

            deallocate( utmp2 ) ; mem=mem-bdmain*mm*nn

         end do ! ms

            if ( any(usermap(0:NPROW-1,IPCOL,1)==myrank) ) j0=j0+nn

         end do ! ns

         if ( ii>0 ) then
            unk(i1:i2,MB_0:MB_1,k,s)=psi_tmp(1:ii,MB_0:MB_1)
         end if

         end do ! ii

         deallocate( psi_tmp )
         deallocate( vtmp2 )
         mem = mem - bdmain*NBLK2*(MB_1-MB_0+1)
         mem = mem - bdmain*LLD_R*LLD_C

         call watch(ct1,et1)
         ct(10)=ct(10)+ct1-ct0
         et(10)=et(10)+et1-et0

      end do ! k
      end do ! s

      deallocate( ir,id )
      mem=mem-bsintg*np_band*2

      mem=mem-bsintg*( size(irecv_me)+size(isend_me) )
      deallocate( isend_me,irecv_me )

!
! --- allgatherv (esp)
!

      allocate( ir(0:np_bzsm-1),id(0:np_bzsm-1) ) ; mem=mem+bsintg*np_bzsm*2 ; memax=max(memax,mem)

      id(0:np_bzsm-1)=id_bzsm(0:np_bzsm-1)*MB
      ir(0:np_bzsm-1)=ir_bzsm(0:np_bzsm-1)*MB
      mrnk           =id_class(myrank,5)

      do s=MSP_0,MSP_1
         call mpi_allgatherv(esp(1,MBZ_0,s),ir(mrnk),mpi_real8,esp(1,1,s),ir,id,mpi_real8,comm_bzsm,ierr)
      end do

      deallocate( id,ir ) ; mem=mem-bsintg*np_bzsm*2

      allocate( ir(0:np_spin-1),id(0:np_spin-1) ) ; mem=mem+bsintg*np_spin*2 ; memax=max(memax,mem)

      id(0:np_bzsm-1)=id_spin(0:np_bzsm-1)*MB*MBZ
      ir(0:np_bzsm-1)=ir_spin(0:np_bzsm-1)*MB*MBZ
      mrnk           =id_class(myrank,5)

      call mpi_allgatherv(esp(1,1,MSP_0),ir(mrnk),mpi_real8,esp(1,1,1),ir,id,mpi_real8,comm_spin,ierr)

      deallocate( id,ir ) ; mem=mem-bsintg*np_spin*2

!

      call watch(ctime1,etime1)

      nop1=sop1*0.5d0*(MB*MB+MB)*ML0
      nop1t=sop1*0.5d0*(MB*MB+MB)*ML
      nop3=sop1*MB*MB*ML0
      nop3t=sop1*MB*MB*ML

      ct(2)=ctime_hpsi
      et(2)=etime_hpsi
      ct(1)=ct(6)-ct(2)-ct(3)-ct(4)-ct(5)
      et(1)=et(6)-et(2)-et(3)-et(4)-et(5)
      ct(9)=ct(10)-ct(8)
      et(9)=et(10)-et(8)

      nop(2)=sum(nop_hpsi(1:8))

      s0(1:10)=ct(1:10)
      s0(11:20)=et(1:10)

      if ( ct(1)>0.d0  ) s0(21)=nop(1)/ct(1)
      if ( ct(2)>0.d0  ) s0(22)=nop(2)/ct(2)
      if ( ct(6)>0.d0  ) s0(23)=(nop(1)+nop(2))/ct(6)
      if ( ct(9)>0.d0  ) s0(24)=nop(9)/ct(9)
      if ( ct(10)>0.d0 ) s0(25)=nop(9)/ct(10)

      s0(26:35)=nop(1:10)
      s0(36)=ctime1-ctime0
      s0(37)=etime1-etime0

      call mpi_allreduce(s0,s1(1,1),37,mpi_real8,mpi_min,mpi_comm_world,ierr)
      call mpi_allreduce(s0,s1(1,2),37,mpi_real8,mpi_max,mpi_comm_world,ierr)
      call mpi_allreduce(s0,s1(1,3),37,mpi_real8,mpi_sum,mpi_comm_world,ierr)

      if (DISP_SWITCH) then
         write(*,'(1x,"TIME(D_DIAG_MOL_SYM_PARA)=",4f10.3)') s1(36,1:2),s1(37,1:2)
         write(*,'(1x," (MATE)      ",6f10.3)') s1(6,1:2),s1(16,1:2),s1(23,1:2)/1.d6
         write(*,'(1x,"  mate       ",6f10.3)') s1(1,1:2),s1(11,1:2),s1(1,1:2)/1.d6
         write(*,'(1x,"  hpsi       ",6f10.3)') s1(2,1:2),s1(12,1:2),s1(2,1:2)/1.d6
         write(*,'(1x,"  allgatherv ",4f10.3)') s1(3,1:2),s1(13,1:2)
         write(*,'(1x,"  reduce     ",4f10.3)') s1(4,1:2),s1(14,1:2)
         write(*,'(1x,"  send/recv  ",4f10.3)') s1(5,1:2),s1(15,1:2)
         write(*,'(1x," (",a6,")    ",4f10.3)') idiag0,s1(7,1:2),s1(17,1:2)
         write(*,'(1x," (ROTV)      ",6f10.3)') s1(10,1:2),s1(20,1:2),s1(25,1:2)/1.d6
         write(*,'(1x,"  mat        ",6f10.3)') s1( 9,1:2),s1(19,1:2),s1(24,1:2)/1.d6
         write(*,'(1x,"  bcast      ",4f10.3)') s1( 8,1:2),s1(18,1:2)
         write(*,*) "MBSIZE,NBSIZE,NBLK2=",MBSIZE,NBSIZE,NBLK2
         write(*,*) "MEM(MB)=",memax*B2MB,mem
         write(*,*) "check # of OP :",s1(26,3)-nop1t,s1(34,3)-nop3t
      end if

      return

 900  call stop_program

      END SUBROUTINE diag_mol_sym

!--------1---------2---------3---------4---------5---------6---------7--

      RECURSIVE SUBROUTINE diag_mol_sym_sub(mm1,mm2,nn1,nn2,MBLK,ns0,ld0,ld1,ct,et,nop,k,s)
      use global_variables
      implicit none
      integer,intent(IN)    :: mm1,mm2,nn1,nn2,MBLK,ns0,ld0,ld1,k,s
      real(8),intent(INOUT) :: ct(*),et(*),nop(*)
      integer :: n1,n2,ML0,n,ns,ne,nn,m,ms,me,mm,mms,MBLKH,i,ld

      n1  = idisp(myrank)+1
      n2  = idisp(myrank)+ircnt(myrank)
      ML0 = ircnt(myrank)
      ld  = ld1-ld0+1

      do ns=nn1,nn2,MBLK
         ne=min(ns+MBLK-1,nn2)
         nn=ne-ns+1
         if ( nn<=0 ) cycle

      do mms=mm1,mm2,MBLK
         ms=max(ns,mms)
         me=min(mms+MBLK-1,mm2)
         mm=me-ms+1
         if ( mm<=0 ) cycle

         if ( ms>=ne ) then

            call dgemm(TRANSA,TRANSB,mm,nn,ML0,zdV,unk(n1,ms,k,s),ML0,psi_tmp(n1,ns-ns0+1),ML0,zero,vtmp2(ms,ns),ld)

            nop(1)=nop(1)+sop1*mm*nn*ML0

         else if ( mm<=NBLK1 ) then

            do n=ns,ne
               call dgemv(TRANSA,ML0,ne-n+1,zdV,unk(n1,n,k,s),ML0,psi_tmp(n1,n-ns0+1),1,zero,vtmp2(n,n),1)
            end do

            nop(1)=nop(1)+sop1*0.5d0*nn*(nn+1.d0)*ML0

         else

            MBLKH=max(MBLK/2,NBLK1)
            call diag_2d_sub(ms,me,ns,ne,MBLKH,ns0,ld0,ld1,ct,et,nop,k,s)

         end if

      end do ! mms
      end do ! ns

      return
      END SUBROUTINE diag_mol_sym_sub
