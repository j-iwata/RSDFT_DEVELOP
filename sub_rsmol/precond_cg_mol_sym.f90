!--------1---------2---------3---------4---------5---------6---------7--

      SUBROUTINE precond_cg_mol_sym(E,k,s,nn,ML0,ir,Npc)
      use global_variables
      implicit none
      integer,intent(IN)    :: k,s,nn,ML0,ir
      integer,intent(INOUT) :: Npc
      real(8),intent(IN)    :: E(nn)
      real(8),parameter :: ep=1.d-24
      real(8) :: rr0(nn),rr1(nn),sb(nn),pAp(nn),a(nn),E0(nn),b,fac
      real(8) :: mem,memax
      integer :: n,n1,n2,mloop,ipc,ierr,iloop,mm

      n1    = idisp(myrank)+1
      n2    = idisp(myrank)+ircnt(myrank)
      mloop = NBLK3
      ipc   = 2 !NBLK4
      mm    = ML0  ; if (TYPE_MAIN==mpi_complex16) mm=mm*2
      E0(:) = 0.d0
      mem   = 0.d0
      memax = 0.d0
!      fac   = dble(nsym)/dble(nn)
      fac   = 1.d0/dble(nn)

      if ( mloop<=0 ) return
      if ( ipc<1 .or. 3<ipc ) return

      allocate( rk_pc(ML0,nn),pk_pc(ML0,nn) )
      allocate( gtmp2(ML0,nn),ftmp2(ML0,nn) )
      mem=mem+bdmain*ML0*nn*4 ; memax=max(memax,mem)

      gtmp2(1:ML0,1:nn)=gk(n1:n2,1:nn)
      select case(ipc)
      case(1) ; call precond_cg_mat1_mol_sym(E0,k,s,ML0,nn,ir)
      case(2) ; call precond_cg_mat2_mol_sym(E0,k,s,ML0,nn,ir)
      case(3) ; call precond_cg_mat3_mol_sym(E0,k,s,ML0,nn,ir)
      end select

      do n=1,nn
         rk_pc(1:ML0,n)=gk(n1:n2,n)-ftmp2(1:ML0,n)
         pk_pc(1:ML0,n)=rk_pc(1:ML0,n)
      end do
      nop_pc=nop_pc+sop3*ML0*nn

!      do n=1,nn
!         sb(n)=sum(abs(rk_pc(1:ML0,n))**2)
!      end do
!      call dot_product_mol_sym(rk_pc,rk_pc,sb,wkgrp,fac,nn,ML0,nop_pc)
      call dot_product_mol_sym(rk_pc,rk_pc,sb,Nstar,fac,nn,ML0,nop_pc)
      nop_pc=nop_pc+sop2*ML0*nn

      call mpi_allreduce(sb,rr0,nn,mpi_real8,mpi_sum,comm_grid,ierr)

      if ( all(rr0(1:nn)<ep) ) then
         deallocate( ftmp2,gtmp2,pk_pc,rk_pc ) ; mem=mem-bdmain*ML0*nn*4
         return
      end if

      Npc = Npc + nn

      do iloop=1,mloop !+1

         gtmp2(1:ML0,1:nn)=pk_pc(1:ML0,1:nn)
         select case(ipc)
         case(1) ; call precond_cg_mat1_mol_sym(E0,k,s,ML0,nn,ir)
         case(2) ; call precond_cg_mat2_mol_sym(E0,k,s,ML0,nn,ir)
         case(3) ; call precond_cg_mat3_mol_sym(E0,k,s,ML0,nn,ir)
         end select

!         do n=1,nn
!            call dot_product(pk_pc(1,n),ftmp2(1,n),sb(n),1.d0,mm,1,nop_pc)
!         end do
!         call dot_product_mol_sym(pk_pc,ftmp2,sb,wkgrp,fac,nn,ML0,nop_pc)
         call dot_product_mol_sym(pk_pc,ftmp2,sb,Nstar,fac,nn,ML0,nop_pc)
         call mpi_allreduce(sb,pAp,nn,mpi_real8,mpi_sum,comm_grid,ierr)

         do n=1,nn
            a(n)=rr0(n)/pAp(n)
            rk_pc(1:ML0,n)=rk_pc(1:ML0,n)-a(n)*ftmp2(1:ML0,n)
         end do
         nop_pc=nop_pc+(sop3*ML0+sop3*ML0)*nn

!         do n=1,nn
!            sb(n)=sum(abs(rk_pc(1:ML0,n))**2)
!         end do
!         call dot_product_mol_sym(rk_pc,rk_pc,sb,wkgrp,fac,nn,ML0,nop_pc)
         call dot_product_mol_sym(rk_pc,rk_pc,sb,Nstar,fac,nn,ML0,nop_pc)
         nop_pc=nop_pc+sop2*ML0*nn

         call mpi_allreduce(sb,rr1,nn,mpi_real8,mpi_sum,comm_grid,ierr)

         if ( iloop==mloop+1 ) then
            exit
         end if

         do n=1,nn
            Pgk(n1:n2,n)=Pgk(n1:n2,n)+a(n)*pk_pc(1:ML0,n)
            b=rr1(n)/rr0(n)
            rr0(n)=rr1(n)
            pk_pc(1:ML0,n)=rk_pc(1:ML0,n)+b*pk_pc(1:ML0,n)
         end do
         nop_pc=nop_pc+(sop3*ML0+sop3*ML0+sop3*ML0+sop3*ML0)*nn

      end do ! iloop

      deallocate( ftmp2,gtmp2,pk_pc,rk_pc ) ; mem=mem-bdmain*ML0*nn*4

      return

      END SUBROUTINE precond_cg_mol_sym

!--------1---------2---------3---------4---------5---------6---------7--

      SUBROUTINE precond_cg_mat1_mol_sym(E,k,s,mm,nn,ir)
      use global_variables
      implicit none
      integer,intent(IN) :: k,s,mm,nn,ir
      real(8) :: E(nn)
      integer :: n,n1,n2

      n1 = idisp(myrank)+1
      n2 = idisp(myrank)+ircnt(myrank)

      call hpsi_mol_sym(k,s,gtmp2,ftmp2,n1,n2,nn,ir)

      if ( any(E(1:nn)/=0.d0) ) then
         do n=1,nn
            ftmp2(1:mm,n)=ftmp2(1:mm,n)-E(n)*gtmp2(1:mm,n)
         end do
      end if

      return
      END SUBROUTINE precond_cg_mat1_mol_sym

!--------1---------2---------3---------4---------5---------6---------7--

      SUBROUTINE precond_cg_mat2_mol_sym(E,k,s,mm,nn,ir)
      use global_variables
      integer,intent(IN) :: k,s,mm,nn,ir
      real(8) :: E(nn)
      real(8) :: c,d
      real(8),allocatable :: work(:,:)
      integer :: n,n1,n2,i,i1,i2,i3,ii,j,id,jd,isym,ierr

      n1 = idisp(myrank)+1
      n2 = idisp(myrank)+ircnt(myrank)

      allocate( work(ML_irreducible+MK_irreducible,nn) ) ; work=0.d0
      do id=1,nn
         call mpi_allgatherv(gtmp2(1,id),mm,mpi_real8,work(1,id) &
              ,ir_grid,id_grid,mpi_real8,comm_grid,ierr)
      end do

      do n=1,nn
         d=coef_lap0-E(n)
         do i=1,mm
!            ftmp2(i,n)=(Vloc(i+n1-1,s)+d)*gtmp2(i,n)
            ftmp2(i,n)=d*gtmp2(i,n)
         end do
      end do
      nop_pc=nop_pc+sop3*mm*nn

      do id=1,nn
         do jd=1,nn
            c=coef_lap(1,1)
            do i=n1,n2
               ii=i-n1+1
               i1=LL(1,i)
               i2=LL(2,i)
               i3=LL(3,i)
               j=LLL(i1-1,i2,i3)
               isym=LLR(i1-1,i2,i3)
               ftmp2(ii,id)=ftmp2(ii,id)+c*Rir(id,jd,isym,ir)*work(j,jd)
               j=LLL(i1+1,i2,i3)
               isym=LLR(i1+1,i2,i3)
               ftmp2(ii,id)=ftmp2(ii,id)+c*Rir(id,jd,isym,ir)*work(j,jd)
               j=LLL(i1,i2-1,i3)
               isym=LLR(i1,i2-1,i3)
               ftmp2(ii,id)=ftmp2(ii,id)+c*Rir(id,jd,isym,ir)*work(j,jd)
               j=LLL(i1,i2+1,i3)
               isym=LLR(i1,i2+1,i3)
               ftmp2(ii,id)=ftmp2(ii,id)+c*Rir(id,jd,isym,ir)*work(j,jd)
               j=LLL(i1,i2,i3-1)
               isym=LLR(i1,i2,i3-1)
               ftmp2(ii,id)=ftmp2(ii,id)+c*Rir(id,jd,isym,ir)*work(j,jd)
               j=LLL(i1,i2,i3+1)
               isym=LLR(i1,i2,i3+1)
               ftmp2(ii,id)=ftmp2(ii,id)+c*Rir(id,jd,isym,ir)*work(j,jd)
            end do ! i
         end do ! jd
      end do ! id
      nop_pc=nop_pc+(3.d0*sop3*mm+3.d0*sop3*mm+3.d0*sop3*mm)*nn

      deallocate( work )

      return
      END SUBROUTINE precond_cg_mat2_mol_sym

!--------1---------2---------3---------4---------5---------6---------7--

      SUBROUTINE precond_cg_mat3_mol_sym(E,k,s,mm,nn,ir)
      use global_variables
      integer,intent(IN) :: k,s,mm,nn,ir
      real(8) :: E(nn)
      real(8) :: c,d
      real(8),allocatable :: work(:,:)
      integer :: n,n1,n2,i,i1,i2,i3,ii,j,m,id,jd,isym,ierr

      n1 = idisp(myrank)+1
      n2 = idisp(myrank)+ircnt(myrank)

      allocate( work(ML_irreducible+MK_irreducible,nn) ) ; work=0.d0
      do id=1,nn
         call mpi_allgatherv(gtmp2(1,id),mm,mpi_real8,work(1,id) &
              ,ir_grid,id_grid,mpi_real8,comm_grid,ierr)
      end do

      do n=1,nn
         d=coef_lap0-E(n)
         do i=1,mm
!            ftmp2(i,n)=(Vloc(i+n1-1,s)+d)*gtmp2(i,n)
            ftmp2(i,n)=d*gtmp2(i,n)
         end do
      end do
      nop_pc=nop_pc+sop3*mm*nn

      do id=1,nn
         do m=1,Md
            do jd=1,nn
               c=coef_lap(1,m)
               do i=n1,n2
                  ii=i-n1+1
                  i1=LL(1,i)
                  i2=LL(2,i)
                  i3=LL(3,i)
                  j=LLL(i1-m,i2,i3)
                  isym=LLR(i1-m,i2,i3)
                  ftmp2(ii,id)=ftmp2(ii,id)+c*Rir(id,jd,isym,ir)*work(j,jd)
                  j=LLL(i1+m,i2,i3)
                  isym=LLR(i1+m,i2,i3)
                  ftmp2(ii,id)=ftmp2(ii,id)+c*Rir(id,jd,isym,ir)*work(j,jd)
                  j=LLL(i1,i2-m,i3)
                  isym=LLR(i1,i2-m,i3)
                  ftmp2(ii,id)=ftmp2(ii,id)+c*Rir(id,jd,isym,ir)*work(j,jd)
                  j=LLL(i1,i2+m,i3)
                  isym=LLR(i1,i2+m,i3)
                  ftmp2(ii,id)=ftmp2(ii,id)+c*Rir(id,jd,isym,ir)*work(j,jd)
                  j=LLL(i1,i2,i3-m)
                  isym=LLR(i1,i2,i3-m)
                  ftmp2(ii,id)=ftmp2(ii,id)+c*Rir(id,jd,isym,ir)*work(j,jd)
                  j=LLL(i1,i2,i3+m)
                  isym=LLR(i1,i2,i3+m)
                  ftmp2(ii,id)=ftmp2(ii,id)+c*Rir(id,jd,isym,ir)*work(j,jd)
               end do ! i
            end do ! jd
         end do ! m
      end do ! id
      nop_pc=nop_pc+(3.d0*sop3*mm+3.d0*sop3*mm+3.d0*sop3*mm)*nn

      deallocate( work )

      return
      END SUBROUTINE precond_cg_mat3_mol_sym
