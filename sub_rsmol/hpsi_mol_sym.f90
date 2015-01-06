!=======================================================================
!========================== Hamiltonian Operation ( for real functions )

      SUBROUTINE hpsi_mol_sym(k,s,tpsi,htpsi,n1,n2,nd,ir)
      use global_variables
      implicit none
      integer,intent(IN) :: k,s,n1,n2,nd,ir
      real(8),intent(IN)  :: tpsi(n1:n2,nd)
      real(8),intent(OUT) :: htpsi(n1:n2,nd)
      real(8) :: c,ctime0,ctime1,etime0,etime1
      integer :: i,id,i1,i2,i3,j,ii,lma,m,ML0,n,jd,isym
      integer :: ierr,nreq
      integer :: irank,jrank,istatus(MPI_STATUS_SIZE,512),ireq(512)
      real(8),allocatable :: w(:,:,:),w0(:,:,:),tphi(:,:)

      call watch(ctime0,etime0)

      allocate( tphi(ML_irreducible+MK_irreducible,nd) ) ; tphi=0.d0
      allocate( uVunk(nzlma,1),uVunk0(nzlma,1) )
!      allocate( w(-ML1-Md:ML1+Md,-ML2-Md:ML2+Md,-ML3-Md:ML3+Md) )
!      allocate( w0(-ML1-Md:ML1+Md,-ML2-Md:ML2+Md,-ML3-Md:ML3+Md) )

      ML0 = n2-n1+1

      do id=1,nd
         call mpi_allgatherv(tpsi(n1,id),ML0,mpi_real8,tphi(1,id) &
              ,ir_grid,id_grid,mpi_real8,comm_grid,ierr)
      end do

      do id=1,nd

!
! --- Kinetic energy ---
!

         do i=n1,n2
            htpsi(i,id)=coef_lap0*tpsi(i,id)
         end do

         goto 10
         w0(:,:,:)=zero
         do i=n1,n2
            do isym=1,nsym
               i1=rga(1,1,isym)*LL(1,i)+rga(1,2,isym)*LL(2,i)+rga(1,3,isym)*LL(3,i)
               i2=rga(2,1,isym)*LL(1,i)+rga(2,2,isym)*LL(2,i)+rga(2,3,isym)*LL(3,i)
               i3=rga(3,1,isym)*LL(1,i)+rga(3,2,isym)*LL(2,i)+rga(3,3,isym)*LL(3,i)
               ii=LLL(i1,i2,i3)
               c=0.d0
               do jd=1,nd
                  c=c+Rir(id,jd,isym,ir)*tpsi(i,jd)
               end do
               if ( w0(i1,i2,i3)==0.d0 ) w0(i1,i2,i3)=c
!               w0(i1,i2,i3)=w0(i1,i2,i3)+c*wkgrp(i)
            end do
         end do
         m=size(w0)
         call mpi_allreduce(w0,w,m,mpi_real8,mpi_sum,comm_grid,ierr)
         do m=1,Md
            c=coef_lap(1,m)
            do i=n1,n2
               i1=LL(1,i)
               i2=LL(2,i)
               i3=LL(3,i)
               htpsi(i,id)=htpsi(i,id)+c*( w(i1+m,i2,i3)+w(i1-m,i2,i3) &
                                          +w(i1,i2+m,i3)+w(i1,i2-m,i3) &
                                          +w(i1,i2,i3+m)+w(i1,i2,i3-m) )
            end do
         end do
 10      continue

         do m=1,Md
            do jd=1,nd
               c=coef_lap(1,m)
               do i=n1,n2
                  i1=LL(1,i)
                  i2=LL(2,i)
                  i3=LL(3,i)
                  j=LLL(i1-m,i2,i3)
                  isym=LLR(i1-m,i2,i3)
                  htpsi(i,id)=htpsi(i,id)+c*Rir(id,jd,isym,ir)*tphi(j,jd)
                  j=LLL(i1+m,i2,i3)
                  isym=LLR(i1+m,i2,i3)
                  htpsi(i,id)=htpsi(i,id)+c*Rir(id,jd,isym,ir)*tphi(j,jd)
                  j=LLL(i1,i2-m,i3)
                  isym=LLR(i1,i2-m,i3)
                  htpsi(i,id)=htpsi(i,id)+c*Rir(id,jd,isym,ir)*tphi(j,jd)
                  j=LLL(i1,i2+m,i3)
                  isym=LLR(i1,i2+m,i3)
                  htpsi(i,id)=htpsi(i,id)+c*Rir(id,jd,isym,ir)*tphi(j,jd)
                  j=LLL(i1,i2,i3-m)
                  isym=LLR(i1,i2,i3-m)
                  htpsi(i,id)=htpsi(i,id)+c*Rir(id,jd,isym,ir)*tphi(j,jd)
                  j=LLL(i1,i2,i3+m)
                  isym=LLR(i1,i2,i3+m)
                  htpsi(i,id)=htpsi(i,id)+c*Rir(id,jd,isym,ir)*tphi(j,jd)
               end do
            end do
         end do

!
! --- local ---
!

         do i=n1,n2
            htpsi(i,id)=htpsi(i,id)+Vloc(i,s)*tpsi(i,id)
         end do

!
! --- Non-local ( Pseudopotential ) ---
!

         select case(pselect)
         case default

         case(1,2)

            goto 20
            do lma=1,nzlma
               uVunk0(lma,1)=zero
               do j=1,MJJ(lma)
                  i1=JJP3(1,j,lma)
                  i2=JJP3(2,j,lma)
                  i3=JJP3(3,j,lma)
                  uVunk0(lma,1)=uVunk0(lma,1)+uVk(j,lma,k)*w(i1,i2,i3)
               end do
               uVunk0(lma,1)=iuV(lma)*dV*uVunk0(lma,1)
               uVunk(lma,1)=uVunk0(lma,1)
            end do
 20         continue

            do lma=1,nzlma
               uVunk0(lma,1)=zero
               do jd=1,nd
                  do j=1,MJJ(lma)
                     i1=JJP3(1,j,lma)
                     i2=JJP3(2,j,lma)
                     i3=JJP3(3,j,lma)
                     i=LLL(i1,i2,i3)
                     isym=LLR(i1,i2,i3)
                     uVunk0(lma,1)=uVunk0(lma,1)+uVk(j,lma,k)*Rir(id,jd,isym,ir)*tphi(i,jd)
                  end do
               end do
               uVunk0(lma,1)=iuV(lma)*dV*uVunk0(lma,1)
               uVunk(lma,1)=uVunk0(lma,1)
            end do

            if ( .not. FLAG_NODE_TEST ) then

            nreq=0
            do irank=0,nprocs_g-1
               m=lma_nsend(irank)
               if ( irank==myrank_g .or. m<=0 ) cycle
               do i1=1,m
                  sbufnl(i1,irank)=uVunk0(sendmap(i1,irank),1)
               end do
               nreq=nreq+1
               call mpi_isend(sbufnl(1,irank),m,TYPE_MAIN,irank,1,comm_grid,ireq(nreq),ierr)
               nreq=nreq+1
               call mpi_irecv(rbufnl(1,irank),m,TYPE_MAIN,irank,1,comm_grid,ireq(nreq),ierr)
            end do
            if ( nreq>0 ) call mpi_waitall(nreq,ireq,istatus,ierr)
            do irank=0,nprocs_g-1
               m=lma_nsend(irank)
               if ( irank==myrank_g .or. m<=0 ) cycle
               do i1=1,m
                  uVunk(recvmap(i1,irank),1)=uVunk(recvmap(i1,irank),1)+rbufnl(i1,irank)
               end do
            end do

            end if

            goto 30
            do lma=1,nzlma
               w0(:,:,:)=0.d0
               do j=1,MJJ(lma)
                  i1=JJP3(1,j,lma)
                  i2=JJP3(2,j,lma)
                  i3=JJP3(3,j,lma)
                  w0(i1,i2,i3)=uVunk(lma,1)*uVk(j,lma,k)
               end do
               do i=n1,n2
                  htpsi(i,id)=htpsi(i,id)+w0(LL(1,i),LL(2,i),LL(3,i))
               end do
            end do
 30         continue

            do lma=1,nzlma
               do j=1,MJJ0(lma)
                  i=JJP(j,lma)
                  htpsi(i,id)=htpsi(i,id)+uVunk(lma,1)*uVk(j,lma,k)
               end do
            end do

         case(3)

            write(*,*) "pselect=",pselect
            write(*,*) "This is only for periodic systems!"
            call stop_program

         end select ! pselect

      end do ! id

!      deallocate( w0 )
!      deallocate( w )
      deallocate( uVunk0,uVunk )
      deallocate( tphi )

      call watch(ctime1,etime1)
      ctime_hpsi=ctime_hpsi+ctime1-ctime0
      etime_hpsi=etime_hpsi+etime1-etime0

      return

 900  call stop_program

      END SUBROUTINE hpsi_mol_sym

!=======================================================================
!========================== Hamiltonian Operation ( for real functions )

      SUBROUTINE hpsi_spe_mol_sym(k,s,tpsi,n1,n2,nd,ir,sb)
      use global_variables
      implicit none
      integer,intent(IN) :: k,s,n1,n2,nd,ir
      real(8),intent(IN)  :: tpsi(n1:n2,nd)
      real(8),intent(OUT) :: sb(4,nd)
      real(8) :: c,ctime0,ctime1,etime0,etime1,fac,nop
      integer :: i,id,i1,i2,i3,j,ii,lma,m,ML0,n,jd,isym
      integer :: ierr,nreq
      integer :: irank,jrank,istatus(MPI_STATUS_SIZE,512),ireq(512)
      real(8),allocatable :: tphi(:,:),work(:,:),htpsi(:,:),sa(:)

      call watch(ctime0,etime0)

      allocate( tphi(ML_irreducible+MK_irreducible,nd) ) ; tphi=0.d0
      allocate( uVunk(nzlma,1),uVunk0(nzlma,1) )
      allocate( work(n1:n2,nd) )
      allocate( htpsi(n1:n2,nd) )
      allocate( sa(nd) )

      ML0 = n2-n1+1
!      fac = dV*nsym/dble(nd)
      fac = dV/dble(nd)

      do id=1,nd
         call mpi_allgatherv(tpsi(n1,id),ML0,mpi_real8,tphi(1,id) &
              ,ir_grid,id_grid,mpi_real8,comm_grid,ierr)
      end do

      do id=1,nd

!
! --- Kinetic energy ---
!

         do i=n1,n2
            htpsi(i,id)=coef_lap0*tpsi(i,id)
         end do

         do m=1,Md
            do jd=1,nd
               c=coef_lap(1,m)
               do i=n1,n2
                  i1=LL(1,i)
                  i2=LL(2,i)
                  i3=LL(3,i)
                  j=LLL(i1-m,i2,i3)
                  isym=LLR(i1-m,i2,i3)
                  htpsi(i,id)=htpsi(i,id)+c*Rir(id,jd,isym,ir)*tphi(j,jd)
                  j=LLL(i1+m,i2,i3)
                  isym=LLR(i1+m,i2,i3)
                  htpsi(i,id)=htpsi(i,id)+c*Rir(id,jd,isym,ir)*tphi(j,jd)
                  j=LLL(i1,i2-m,i3)
                  isym=LLR(i1,i2-m,i3)
                  htpsi(i,id)=htpsi(i,id)+c*Rir(id,jd,isym,ir)*tphi(j,jd)
                  j=LLL(i1,i2+m,i3)
                  isym=LLR(i1,i2+m,i3)
                  htpsi(i,id)=htpsi(i,id)+c*Rir(id,jd,isym,ir)*tphi(j,jd)
                  j=LLL(i1,i2,i3-m)
                  isym=LLR(i1,i2,i3-m)
                  htpsi(i,id)=htpsi(i,id)+c*Rir(id,jd,isym,ir)*tphi(j,jd)
                  j=LLL(i1,i2,i3+m)
                  isym=LLR(i1,i2,i3+m)
                  htpsi(i,id)=htpsi(i,id)+c*Rir(id,jd,isym,ir)*tphi(j,jd)
               end do
            end do
         end do

      end do ! id

!      call dot_product_mol_sym(tpsi,htpsi,sa,wkgrp,fac,nd,ML0,nop)
      call dot_product_mol_sym(tpsi,htpsi,sa,Nstar,fac,nd,ML0,nop)
      sb(1,1:nd)=sa(1:nd)

!
! --- local ---
!
      do id=1,nd
         work(n1:n2,id)=Vloc(n1:n2,s)*tpsi(n1:n2,id)
         htpsi(n1:n2,id)=htpsi(n1:n2,id)+work(n1:n2,id)
      end do

!      call dot_product_mol_sym(tpsi,work,sa,wkgrp,fac,nd,ML0,nop)
      call dot_product_mol_sym(tpsi,work,sa,Nstar,fac,nd,ML0,nop)
      sb(2,1:nd)=sa(1:nd)

!
! --- Non-local ( Pseudopotential ) ---
!
      work(n1:n2,1:nd)=0.d0

      do id=1,nd

         select case(pselect)
         case default

         case(1,2)

            do lma=1,nzlma
               uVunk0(lma,1)=zero
               do jd=1,nd
                  do j=1,MJJ(lma)
                     i1=JJP3(1,j,lma)
                     i2=JJP3(2,j,lma)
                     i3=JJP3(3,j,lma)
                     i=LLL(i1,i2,i3)
                     isym=LLR(i1,i2,i3)
                     uVunk0(lma,1)=uVunk0(lma,1)+uVk(j,lma,k)*Rir(id,jd,isym,ir)*tphi(i,jd)
                  end do
               end do
               uVunk0(lma,1)=iuV(lma)*dV*uVunk0(lma,1)
               uVunk(lma,1)=uVunk0(lma,1)
            end do

            if ( .not. FLAG_NODE_TEST ) then

            nreq=0
            do irank=0,nprocs_g-1
               m=lma_nsend(irank)
               if ( irank==myrank_g .or. m<=0 ) cycle
               do i1=1,m
                  sbufnl(i1,irank)=uVunk0(sendmap(i1,irank),1)
               end do
               nreq=nreq+1
               call mpi_isend(sbufnl(1,irank),m,TYPE_MAIN,irank,1,comm_grid,ireq(nreq),ierr)
               nreq=nreq+1
               call mpi_irecv(rbufnl(1,irank),m,TYPE_MAIN,irank,1,comm_grid,ireq(nreq),ierr)
            end do
            if ( nreq>0 ) call mpi_waitall(nreq,ireq,istatus,ierr)
            do irank=0,nprocs_g-1
               m=lma_nsend(irank)
               if ( irank==myrank_g .or. m<=0 ) cycle
               do i1=1,m
                  uVunk(recvmap(i1,irank),1)=uVunk(recvmap(i1,irank),1)+rbufnl(i1,irank)
               end do
            end do

            end if

            do lma=1,nzlma
               do j=1,MJJ0(lma)
                  i=JJP(j,lma)
                  work(i,id)=work(i,id)+uVunk(lma,1)*uVk(j,lma,k)
               end do
            end do

         case(3)

            write(*,*) "pselect=",pselect
            write(*,*) "This is only for periodic systems!"
            call stop_program

         end select ! pselect

         htpsi(n1:n2,id)=htpsi(n1:n2,id)+work(n1:n2,id)

      end do ! id

!      call dot_product_mol_sym(tpsi,work,sa,wkgrp,fac,nd,ML0,nop)
      call dot_product_mol_sym(tpsi,work,sa,Nstar,fac,nd,ML0,nop)
      sb(3,1:nd)=sa(1:nd)
!      call dot_product_mol_sym(tpsi,htpsi,sa,wkgrp,fac,nd,ML0,nop)
      call dot_product_mol_sym(tpsi,htpsi,sa,Nstar,fac,nd,ML0,nop)
      sb(4,1:nd)=sa(1:nd)

      deallocate( sa )
      deallocate( htpsi )
      deallocate( work )
      deallocate( uVunk0,uVunk )
      deallocate( tphi )

      call watch(ctime1,etime1)
      ctime_hpsi=ctime_hpsi+ctime1-ctime0
      etime_hpsi=etime_hpsi+etime1-etime0

      return

 900  call stop_program

      END SUBROUTINE hpsi_spe_mol_sym
