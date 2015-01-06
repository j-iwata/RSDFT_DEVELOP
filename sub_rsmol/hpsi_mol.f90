!=======================================================================
!========================== Hamiltonian Operation ( for real functions )

      SUBROUTINE hpsi_mol(k,s,tpsi,htpsi,n1,n2,ib1,ib2)
      use global_variables
      implicit none
      integer,intent(IN) :: k,s,n1,n2,ib1,ib2
      real(8),intent(IN)  :: tpsi(n1:n2,ib1:ib2)
      real(8),intent(OUT) :: htpsi(n1:n2,ib1:ib2)
      real(8) :: c,ctime0,ctime1,etime0,etime1
      integer :: i,ib,i1,i2,i3,j,lma,m,ML0,n,nb
      integer :: ierr,nreq
      integer :: irank,jrank,istatus(MPI_STATUS_SIZE,512),ireq(512)

      call watch(ctime0,etime0)

      ML0 = n2-n1+1
      nb  = ib2-ib1+1

!
! --- Kinetic energy ---
!

      do ib=ib1,ib2
         n=ib-ib1+1
         do i=n1,n2
            i1=LL(1,i) ; i2=LL(2,i) ; i3=LL(3,i)
            www(i1,i2,i3,n)=tpsi(i,ib)
         end do
      end do

      call bcset_2d(1,nb,Md,0)

      do ib=ib1,ib2
         do i=n1,n2
            htpsi(i,ib)=coef_lap0*tpsi(i,ib)
         end do
      end do

      do ib=ib1,ib2
         n=ib-ib1+1
         do m=1,Md
            c=coef_lap(1,m)
            do i=n1,n2
               i1=LL(1,i)
               i2=LL(2,i)
               i3=LL(3,i)
               htpsi(i,ib)=htpsi(i,ib)+c*( www(i1+m,i2,i3,n)+www(i1-m,i2,i3,n) &
                                          +www(i1,i2+m,i3,n)+www(i1,i2-m,i3,n) &
                                          +www(i1,i2,i3+m,n)+www(i1,i2,i3-m,n) )
            end do
         end do
      end do

!
! --- local ---
!

      do ib=ib1,ib2
         do i=n1,n2
            htpsi(i,ib)=htpsi(i,ib)+Vloc(i,s)*tpsi(i,ib)
         end do
      end do

      if ( iflag_hf>0 ) then
         do ib=ib1,ib2
            call Fock(n1,n2,ib,tpsi(n1,ib),htpsi(n1,ib))
         end do
      end if

!
! --- Non-local ( Pseudopotential ) ---
!

      select case(pselect)
      case default

      case(1,2)

         allocate( uVunk(nzlma,ib1:ib2),uVunk0(nzlma,ib1:ib2) )

         do ib=ib1,ib2
         do lma=1,nzlma
            uVunk(lma,ib)=zero
            do j=1,MJJ(lma)
               i=JJP(j,lma)
               uVunk(lma,ib)=uVunk(lma,ib)+uVk(j,lma,k)*tpsi(i,ib)
            end do
            uVunk(lma,ib)=iuV(lma)*dV*uVunk(lma,ib)
         end do
         end do

         if ( .not. FLAG_NODE_TEST ) then

         select case(iswitch_eqdiv)
         case default

            do i=1,6
               select case(i)
               case(1,3,5)
                  j=i+1
                  uVunk0(:,:)=uVunk(:,:)
               case(2,4,6)
                  j=i-1
               end select
               do m=1,nrlma_xyz(i)
                  nreq=0
                  irank=num_2_rank(m,i)
                  jrank=num_2_rank(m,j)
                  if( irank>=0 )then
                     i2=0
                     do ib=ib1,ib2
                     do i1=1,lma_nsend(irank)
                        i2=i2+1
                        sbufnl(i2,irank)=uVunk0(sendmap(i1,irank),ib)
                     end do
                     end do
                     nreq=nreq+1
                     call mpi_isend(sbufnl(1,irank),lma_nsend(irank)*nb,TYPE_MAIN,irank,1,comm_grid,ireq(nreq),ierr)
                  end if
                  if( jrank>=0 )then
                     nreq=nreq+1
                     call mpi_irecv(rbufnl(1,jrank),lma_nsend(jrank)*nb,TYPE_MAIN,jrank,1,comm_grid,ireq(nreq),ierr)
                  end if
                  call mpi_waitall(nreq,ireq,istatus,ierr)
                  if ( jrank>=0 ) then
                     i2=0
                     do ib=ib1,ib2
                     do i1=1,lma_nsend(jrank)
                        i2=i2+1
                        uVunk(recvmap(i1,jrank),ib)=uVunk(recvmap(i1,jrank),ib)+rbufnl(i2,jrank)
                     end do
                     end do
                  end if
               end do
            end do

         case(1)

            nreq=0
            do irank=0,nprocs_g-1
               m=lma_nsend(irank)
               if ( irank==myrank_g .or. m<=0 ) cycle
               i2=0
               do ib=ib1,ib2
               do i1=1,m
                  i2=i2+1
                  sbufnl(i2,irank)=uVunk(sendmap(i1,irank),ib)
               end do
               end do
               nreq=nreq+1
               call mpi_isend(sbufnl(1,irank),m*nb,TYPE_MAIN,irank,1,comm_grid,ireq(nreq),ierr)
               nreq=nreq+1
               call mpi_irecv(rbufnl(1,irank),m*nb,TYPE_MAIN,irank,1,comm_grid,ireq(nreq),ierr)
            end do
            if ( nreq>0 ) call mpi_waitall(nreq,ireq,istatus,ierr)
            do irank=0,nprocs_g-1
               m=lma_nsend(irank)
               if ( irank==myrank_g .or. m<=0 ) cycle
               i2=0
               do ib=ib1,ib2
               do i1=1,m
                  i2=i2+1
                  uVunk(recvmap(i1,irank),ib)=uVunk(recvmap(i1,irank),ib)+rbufnl(i2,irank)
               end do
               end do
            end do

         end select ! iswitch_eqdiv

         end if

         do ib=ib1,ib2
         do lma=1,nzlma
            do j=1,MJJ(lma)
               i=JJP(j,lma)
               htpsi(i,ib)=htpsi(i,ib)+uVunk(lma,ib)*uVk(j,lma,k)
            end do
         end do
         end do

         deallocate( uVunk0,uVunk )

      case(3)

         write(*,*) "pselect=",pselect
         write(*,*) "This is only for periodic systems!"
         call stop_program

      end select

      call watch(ctime1,etime1)
      ctime_hpsi=ctime_hpsi+ctime1-ctime0
      etime_hpsi=etime_hpsi+etime1-etime0

      return
      END SUBROUTINE hpsi_mol

!--------1---------2---------3---------4---------5---------6---------7--
!=======================================================================

      SUBROUTINE hpsi_spe_mol(k,s,ib1,ib2,e,iswitch)
      use global_variables
      implicit none
      integer,intent(IN)  :: k,s,ib1,ib2,iswitch
      real(8),intent(OUT) :: e(ib1:ib2)
      integer :: i,ib,i1,i2,i3,j,lma,m,ML0,n,nb,n1,n2
      integer :: ierr,nreq,irank,jrank,istatus(MPI_STATUS_SIZE,512),ireq(512)
      real(8) :: c,mem,ct0,et0,ct1,et1,ct2,et2,ct3,et3

      real(8) :: r1,r2
      integer,allocatable :: itmp(:,:,:)

      n1  = idisp(myrank)+1
      n2  = idisp(myrank)+ircnt(myrank)
      ML0 = n2-n1+1
      nb  = ib2-ib1+1
      mem = 0.d0

      allocate( vtmp2(n1:n2,ib1:ib2) ) ; mem=mem+bdmain*ML0*nb

!
! --- Kinetic energy ---
!

      if ( iswitch==0 .or. iswitch==10 .or. iswitch==1 ) then

      call watch(ct0,et0)

      do ib=ib1,ib2
         n=ib-ib1+1
         do i=n1,n2
            i1=LL(1,i) ; i2=LL(2,i) ; i3=LL(3,i)
            www(i1,i2,i3,n)=unk(i,ib,k,s)
         end do
      end do

      call watch(ct1,et1)

      call bcset_2d(1,nb,Md,0)

      call watch(ct2,et2)

      do ib=ib1,ib2
         do i=n1,n2
            vtmp2(i,ib)=coef_lap0*unk(i,ib,k,s)
         end do
      end do

      do ib=ib1,ib2
         n=ib-ib1+1
         do m=1,Md
            c=coef_lap(1,m)
            do i=n1,n2
               i1=LL(1,i)
               i2=LL(2,i)
               i3=LL(3,i)
               vtmp2(i,ib)=vtmp2(i,ib)+c*( www(i1-m,i2,i3,n)+www(i1+m,i2,i3,n) &
                                          +www(i1,i2-m,i3,n)+www(i1,i2+m,i3,n) &
                                          +www(i1,i2,i3-m,n)+www(i1,i2,i3+m,n) )
            end do
         end do
      end do

      call watch(ct3,et3)

      ct_hpsi(1)=ct_hpsi(1)+ct1-ct0 ; et_hpsi(1)=et_hpsi(1)+et1-et0
      ct_hpsi(2)=ct_hpsi(2)+ct2-ct1 ; et_hpsi(2)=et_hpsi(2)+et2-et1
      ct_hpsi(3)=ct_hpsi(3)+ct3-ct2 ; et_hpsi(3)=et_hpsi(3)+et3-et2
      nop_hpsi(3)=nop_hpsi(3) + sop3*ML0*nb + nb*ML0*Md*(6.d0*sop3+sop3)

      end if

!
! --- single-particle energy 1 ---
!

      if ( iswitch==10 ) then

      call watch(ct0,et0)

      do ib=ib1,ib2
         tsp=tsp+occ(ib,k,s)*sum(vtmp2(n1:n2,ib)*unk(n1:n2,ib,k,s))*dV
      end do

      call watch(ct1,et1)

      ct_hpsi(9)=ct_hpsi(9)+ct1-ct0 ; et_hpsi(9)=et_hpsi(9)+et1-et0
      nop_hpsi(9)=nop_hpsi(9)+sop1*ML0*nb

      end if

!
! --- Local potential ---
!

      if ( iswitch==0 .or. iswitch==10 .or. iswitch==2 ) then

      call watch(ct0,et0)

      do ib=ib1,ib2
         do i=n1,n2
            vtmp2(i,ib)=vtmp2(i,ib)+Vloc(i,s)*unk(i,ib,k,s)
         end do
      end do

      call watch(ct1,et1)
      ct_hpsi(4)=ct_hpsi(4)+ct1-ct0 ; et_hpsi(4)=et_hpsi(4)+et1-et0
      nop_hpsi(4)=nop_hpsi(4)+sop2*ML0*nb

      end if

!
! --- single-particle energy 2 ---
!

      if ( iswitch==10 ) then

      call watch(ct0,et0)

      do ib=ib1,ib2
         elocsp=elocsp+occ(ib,k,s)*sum(abs(unk(n1:n2,ib,k,s))**2*Vloc(n1:n2,s))*dV
      end do

      call watch(ct1,et1)

      ct_hpsi(9)=ct_hpsi(9)+ct1-ct0 ; et_hpsi(9)=et_hpsi(9)+et1-et0
      nop_hpsi(9)=nop_hpsi(9)+sop2*sop2*ML0*nb

      end if

!
! --- Fock ---
!

      if ( iflag_hf>0 ) then
         allocate( utmp2(n1:n2,ib1:ib2) ) ; utmp2=zero
         mem=mem+bdmain*size(utmp2)
         if ( iswitch==0 .or. iswitch==10 .or. iswitch==4 ) then
            do ib=ib1,ib2
               call Fock(n1,n2,ib,unk(n1,ib,k,s),utmp2(n1,ib))
            end do
            do ib=ib1,ib2
               do i=n1,n2
                  vtmp2(i,ib)=vtmp2(i,ib)+utmp2(i,ib)
               end do
            end do
         end if
         if ( iswitch==10 ) then
            do ib=ib1,ib2
               efocksp=efocksp+occ(ib,k,s)*sum(unk(n1:n2,ib,k,s)*utmp2(n1:n2,ib))*dV
            end do
         end if
         mem=mem-bdmain*size(utmp2) ; deallocate(utmp2)
      end if

! 
! --- Non-local pseudopotential ---
!

      select case(pselect)
      case default

      case(1,2)

         if ( iswitch==0 .or. iswitch==10 .or. iswitch==3 ) then

         call watch(ct0,et0)

         allocate( uVunk(nzlma,ib1:ib2),uVunk0(nzlma,ib1:ib2) ) ; mem=mem+bdmain*(nzlma*nb)*2

         do ib=ib1,ib2
         do lma=1,nzlma
            uVunk(lma,ib)=zero
            do j=1,MJJ(lma)
               i=JJP(j,lma)
               uVunk(lma,ib)=uVunk(lma,ib)+uVk(j,lma,k)*unk(i,ib,k,s)
            end do
            uVunk(lma,ib)=iuV(lma)*dV*uVunk(lma,ib)
         end do
         end do

         call watch(ct1,et1)

         if ( .not. FLAG_NODE_TEST ) then

         select case(iswitch_eqdiv)
         case default

            do i=1,6
               select case(i)
               case(1,3,5)
                  j=i+1
                  uVunk0(:,:)=uVunk(:,:)
               case(2,4,6)
                  j=i-1
               end select
               do m=1,nrlma_xyz(i)
                  nreq=0
                  irank=num_2_rank(m,i)
                  jrank=num_2_rank(m,j)
                  if( irank>=0 )then
                     i2=0
                     do ib=ib1,ib2
                     do i1=1,lma_nsend(irank)
                        i2=i2+1
                        sbufnl(i2,irank)=uVunk0(sendmap(i1,irank),ib)
                     end do
                     end do
                     nreq=nreq+1
                     call mpi_isend(sbufnl(1,irank),lma_nsend(irank)*nb,TYPE_MAIN,irank,1,comm_grid,ireq(nreq),ierr)
                  end if
                  if( jrank>=0 )then
                     nreq=nreq+1
                     call mpi_irecv(rbufnl(1,jrank),lma_nsend(jrank)*nb,TYPE_MAIN,jrank,1,comm_grid,ireq(nreq),ierr)
                  end if
                  call mpi_waitall(nreq,ireq,istatus,ierr)
                  if ( jrank>=0 ) then
                     i2=0
                     do ib=ib1,ib2
                     do i1=1,lma_nsend(jrank)
                        i2=i2+1
                        uVunk(recvmap(i1,jrank),ib)=uVunk(recvmap(i1,jrank),ib)+rbufnl(i2,jrank)
                     end do
                     end do
                  end if
               end do
            end do

         case(1)

            nreq=0
            do irank=0,nprocs_g-1
               m=lma_nsend(irank)
               if ( irank==myrank_g .or. m<=0 ) cycle
               i2=0
               do ib=ib1,ib2
               do i1=1,m
                  i2=i2+1
                  sbufnl(i2,irank)=uVunk(sendmap(i1,irank),ib)
               end do
               end do
               nreq=nreq+1
               call mpi_isend(sbufnl(1,irank),m*nb,TYPE_MAIN,irank,1,comm_grid,ireq(nreq),ierr)
               nreq=nreq+1
               call mpi_irecv(rbufnl(1,irank),m*nb,TYPE_MAIN,irank,1,comm_grid,ireq(nreq),ierr)
            end do
            if ( nreq>0 ) call mpi_waitall(nreq,ireq,istatus,ierr)
            do irank=0,nprocs_g-1
               m=lma_nsend(irank)
               if ( irank==myrank_g .or. m<=0 ) cycle
               i2=0
               do ib=ib1,ib2
               do i1=1,m
                  i2=i2+1
                  uVunk(recvmap(i1,irank),ib)=uVunk(recvmap(i1,irank),ib)+rbufnl(i2,irank)
               end do
               end do
            end do

         end select ! iswitch_eqdiv

         end if

         call watch(ct2,et2)

         allocate( utmp2(n1:n2,ib1:ib2) ) ; utmp2=zero
         mem=mem+bdmain*ML0*nb

         do ib=ib1,ib2
            do lma=1,nzlma
               do j=1,MJJ(lma)
                  i=JJP(j,lma)
                  utmp2(i,ib)=utmp2(i,ib)+uVk(j,lma,k)*uVunk(lma,ib)
               end do
            end do
            vtmp2(n1:n2,ib)=vtmp2(n1:n2,ib)+utmp2(n1:n2,ib)
         end do

         call watch(ct3,et3)
 
         ct_hpsi(5)=ct_hpsi(5)+ct1-ct0 ; et_hpsi(5)=et_hpsi(5)+et1-et0
         ct_hpsi(6)=ct_hpsi(6)+ct2-ct1 ; et_hpsi(6)=et_hpsi(6)+et2-et1
         ct_hpsi(7)=ct_hpsi(7)+ct3-ct2 ; et_hpsi(7)=et_hpsi(7)+et3-et2
         nop_hpsi(5)=nop_hpsi(5)+sop1*sum(MJJ(1:nzlma))*nb
         nop_hpsi(7)=nop_hpsi(7)+sop1*sum(MJJ(1:nzlma))*nb+2.d0*ML0*nb

!
! --- single-particle energy 3 ---
!

         if ( iswitch==10 ) then

         call watch(ct0,et0)

         do ib=ib1,ib2
            enlsp=enlsp+occ(ib,k,s)*sum(unk(n1:n2,ib,k,s)*utmp2(n1:n2,ib))*dV
         end do

         call watch(ct1,et1)

         ct_hpsi(9)=ct_hpsi(9)+ct1-ct0 ; et_hpsi(9)=et_hpsi(9)+et1-et0
         nop_hpsi(9)=nop_hpsi(9)+sop1*ML0*nb

         end if

         mem=mem-bdmain*(size(utmp2)+size(uVunk0)+size(uVunk))
         deallocate( utmp2,uVunk0,uVunk )

         end if ! iswitch

      case(3)

         write(*,*) "pselect=",pselect
         write(*,*) "This is only for periodic systems!"
         call stop_program

      end select ! pselect

!
! --- single-particle energy 4 ---
!

      if ( iswitch==10 ) then

      call watch(ct0,et0)

      do ib=ib1,ib2
         e(ib)=sum(unk(n1:n2,ib,k,s)*vtmp2(n1:n2,ib))*dV
      end do

      call watch(ct1,et1)

      ct_hpsi(9)=ct_hpsi(9)+ct1-ct0 ; et_hpsi(9)=et_hpsi(9)+et1-et0
      nop_hpsi(9)=nop_hpsi(9)+sop1*ML0*nb

      end if

      mem=mem-bdmain*size(vtmp2) ; deallocate( vtmp2 )

      return
      END SUBROUTINE hpsi_spe_mol
