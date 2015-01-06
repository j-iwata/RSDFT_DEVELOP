!--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3--
!
      SUBROUTINE init_parallel_mol_sym
      use global_variables
      implicit none

      integer :: i,j,ix,iy,iz,i1,i2,i3,i4,m,n,ierr,isym,LL_tmp(3),irank,n1,n2,mark_KK
      integer :: Mx,My,Mz
      real(8),parameter :: eps=1.d-10
      real(8) :: mem,memax,ctime0,ctime1,etime0,etime1,Rc2,r2,sum0,sum1,H
      real(8),allocatable :: work(:)
      integer,allocatable :: iwork(:),itmp(:),itmp2(:,:)

      if ( isymmetry/=1 ) return

      call watch(ctime0,etime0)
      if (DISP_SWITCH) then
         write(*,'(a60," init_parallel_mol_sym")') repeat("-",60)
      end if

      mem   = 0.d0
      memax = 0.d0

      H  = H1
      Mx = ML1+Md
      My = ML2+Md
      Mz = ML3+Md

! pinfo_grid_2d(7:8,0:np_grid-1)

      pinfo_grid_2d(8,0:np_grid-1) = ML_irreducible/np_grid
      if ( any(pinfo_grid_2d(8,:)==0) ) then
         write(*,*) "Too much CPUs!"
         write(*,*) "ML_irreducible,np_grid=",ML_irreducible,np_grid
         goto 900
      end if
      m=ML_irreducible-sum(pinfo_grid_2d(8,0:np_grid-1))
      do n=1,m
         irank=mod(n-1,np_grid)
         pinfo_grid_2d(8,irank)=pinfo_grid_2d(8,irank)+1
      end do
      m=sum(pinfo_grid_2d(8,:))
      if ( m/=ML_irreducible ) then
         write(*,*) "# of grid points is incosistent"
         goto 900
      end if
      do n=0,np_grid-1
         pinfo_grid_2d(7,n)=sum(pinfo_grid_2d(8,0:n))-pinfo_grid_2d(8,n)
      end do

      allocate( LL(3,ML_irreducible) ) ; LL=0

      do iz=-Mz,Mz
      do iy=-My,My
      do ix=-Mx,Mx
         i=LLL(ix,iy,iz)
         if ( i>0 ) then
            if ( any(LL(1:3,i)/=0) ) cycle
            LL(1,i)=ix
            LL(2,i)=iy
            LL(3,i)=iz
         end if
      end do
      end do
      end do

      goto 1000
      n=0
      j=0
      m=0
      do iz=0,ML3
      do iy=0,ML2
      do ix=0,ML1
         m=m+1
         i=LLL(ix,iy,iz)
         if ( i>0 ) then
            if ( any(LL(1:3,i)/=0) ) cycle
            j=j+1
            LL(1,i)=ix
            LL(2,i)=iy
            LL(3,i)=iz
         end if
      end do
      end do
      end do
      if ( n==0 .and. j==ML_irreducible ) n=m
      do iz=-1,-ML3,-1
      do iy=0,ML2
      do ix=0,ML1
         m=m+1
         i=LLL(ix,iy,iz)
         if ( i>0 ) then
            if ( any(LL(1:3,i)/=0) ) cycle
            j=j+1
            LL(1,i)=ix
            LL(2,i)=iy
            LL(3,i)=iz
         end if
      end do
      end do
      end do
      if ( n==0 .and. j==ML_irreducible ) n=m
      do iz=0,ML3
      do iy=-1,-ML2,-1
      do ix=0,ML1
         m=m+1
         i=LLL(ix,iy,iz)
         if ( i>0 ) then
            if ( any(LL(1:3,i)/=0) ) cycle
            j=j+1
            LL(1,i)=ix
            LL(2,i)=iy
            LL(3,i)=iz
         end if
      end do
      end do
      end do
      if ( n==0 .and. j==ML_irreducible ) n=m
      do iz=0,ML3
      do iy=0,ML2
      do ix=-1,-ML1,-1
         m=m+1
         i=LLL(ix,iy,iz)
         if ( i>0 ) then
            if ( any(LL(1:3,i)/=0) ) cycle
            j=j+1
            LL(1,i)=ix
            LL(2,i)=iy
            LL(3,i)=iz
         end if
      end do
      end do
      end do
      if ( n==0 .and. j==ML_irreducible ) n=m
      do iz=-1,-ML3,-1
      do iy=-1,-ML2,-1
      do ix=0,ML1
         m=m+1
         i=LLL(ix,iy,iz)
         if ( i>0 ) then
            if ( any(LL(1:3,i)/=0) ) cycle
            j=j+1
            LL(1,i)=ix
            LL(2,i)=iy
            LL(3,i)=iz
         end if
      end do
      end do
      end do
      if ( n==0 .and. j==ML_irreducible ) n=m
      do iz=-1,-ML3,-1
      do iy=0,ML2
      do ix=-1,-ML1,-1
         m=m+1
         i=LLL(ix,iy,iz)
         if ( i>0 ) then
            if ( any(LL(1:3,i)/=0) ) cycle
            j=j+1
            LL(1,i)=ix
            LL(2,i)=iy
            LL(3,i)=iz
         end if
      end do
      end do
      end do
      if ( n==0 .and. j==ML_irreducible ) n=m
      do iz=0,ML3
      do iy=-1,-ML2,-1
      do ix=-1,-ML1,-1
         m=m+1
         i=LLL(ix,iy,iz)
         if ( i>0 ) then
            if ( any(LL(1:3,i)/=0) ) cycle
            j=j+1
            LL(1,i)=ix
            LL(2,i)=iy
            LL(3,i)=iz
         end if
      end do
      end do
      end do
      if ( n==0 .and. j==ML_irreducible ) n=m
      do iz=-1,-ML3,-1
      do iy=-1,-ML2,-1
      do ix=-1,-ML1,-1
         m=m+1
         i=LLL(ix,iy,iz)
         if ( i>0 ) then
            if ( any(LL(1:3,i)/=0) ) cycle
            j=j+1
            LL(1,i)=ix
            LL(2,i)=iy
            LL(3,i)=iz
         end if
      end do
      end do
      end do
      if ( n==0 .and. j==ML_irreducible ) n=m
      if ( DISP_SWITCH ) then
         write(*,*) "j,ML_irreducible,n=",j,ML_irreducible,n
         write(*,*) "m,(2*ML1+1)*(2*ML2+1)*(2*ML3+1) =",m,(2*ML1+1)*(2*ML2+1)*(2*ML3+1)
      end if
 1000 continue


!- allocate -
      allocate( iwork(ML_irreducible) )
      iwork=0 ; mem=mem+bsintg*size(iwork) ; memax=max(mem,memax)
      allocate( work(ML_irreducible) )
      work=0.d0 ; mem=mem+bdreal*size(work) ; memax=max(mem,memax)
!-
      do i=1,ML_irreducible
         work(i)=H*H*(LL(1,i)*LL(1,i)+LL(2,i)*LL(2,i)+LL(3,i)*LL(3,i))
      end do
      call indexx(ML_irreducible,work,iwork)
!- deallocate & allocate -
      mem=mem-bdreal*size(work) ; deallocate( work )
      allocate( itmp(ML_irreducible) )
      itmp=-1 ; mem=mem+bsintg*size(itmp) ; memax=max(mem,memax)
      allocate( itmp2(3,ML_irreducible) )
      itmp2=0 ; mem=mem+bsintg*size(itmp2) ; memax=max(mem,memax)
!-
      do i=1,ML_irreducible
         itmp(iwork(i))=mod(i-1,np_grid)
      end do
      j=0
      do n=0,np_grid-1
         do i=1,ML_irreducible
            if ( itmp(iwork(i))==n ) then
               j=j+1
               itmp2(1:3,j)=LL(1:3,iwork(i))
            end if
         end do
      end do
      LL(1:3,1:ML_irreducible)=itmp2(1:3,1:ML_irreducible)
!      if(myrank==0)then
!         do n=0,np_grid-1
!            write(*,*) n,count(itmp==n),count(itmp>=0)
!         end do
!         write(*,*)
!         do i=1,ML_irreducible
!            write(*,*) i,H*H*sum(LL(1:3,i)**2),mod(i-1,np_grid)
!         end do
!      end if
!- deallocate -
      mem=mem-bsintg*size(itmp2) ; deallocate( itmp2 )
      mem=mem-bsintg*size(itmp)  ; deallocate( itmp )
      mem=mem-bsintg*size(iwork) ; deallocate( iwork )
!-


      mark_KK=minval(LLL)*2
      do i=1,ML_irreducible
         i1=LL(1,i)
         i2=LL(2,i)
         i3=LL(3,i)
         do m=-Md,Md
            if ( m==0 ) cycle
            if ( LLL(i1+m,i2,i3)<0 ) LLL(i1+m,i2,i3)=mark_KK
            if ( LLL(i1,i2+m,i3)<0 ) LLL(i1,i2+m,i3)=mark_KK
            if ( LLL(i1,i2,i3+m)<0 ) LLL(i1,i2,i3+m)=mark_KK
         end do
      end do
      MK_irreducible=count(LLL==mark_KK)
      if ( DISP_SWITCH ) then
         write(*,*) "# of outer grid points reduced by symmetry",MK_irreducible,MK,minval(LLL)
      end if

! idisp2,ircnt2,ir_grid2,id_grid2

      call gv_alloc("idisp2")

      ir_grid2(0:np_grid-1) = MK_irreducible/np_grid
      if ( any(ir_grid2==0) ) then
         write(*,*) "Too much CPUs!"
         write(*,*) "MK_irreducible,np_grid=",MK_irreducible,np_grid
         goto 900
      end if
      m=MK_irreducible-sum(ir_grid2(0:np_grid-1))
      do n=1,m
         irank=mod(n-1,np_grid)
         ir_grid2(irank)=ir_grid2(irank)+1
      end do
      m=sum(ir_grid2)
      if ( m/=MK_irreducible ) then
         write(*,*) "# of grid points is incosistent"
         goto 900
      end if
      do n=0,np_grid-1
         id_grid2(n)=sum(ir_grid2(0:n))-ir_grid2(n)
      end do

      do n=0,nprocs-1
         idisp2(n)=id_grid2(id_class(n,0))
         ircnt2(n)=ir_grid2(id_class(n,0))
      end do

      allocate( KK(3,ML_irreducible+1:ML_irreducible+MK_irreducible) ) ; KK=0

      j=ML_irreducible
      do i=1,ML_irreducible
         i1=LL(1,i)
         i2=LL(2,i)
         i3=LL(3,i)
         do m=-Md,Md
            if ( m==0 ) cycle
            if ( LLL(i1+m,i2,i3)==mark_KK ) then
               j=j+1
               KK(1,j)=i1+m               
               KK(2,j)=i2
               KK(3,j)=i3
               LLL(i1+m,i2,i3)=0
            end if
            if ( LLL(i1,i2+m,i3)==mark_KK ) then
               j=j+1
               KK(1,j)=i1               
               KK(2,j)=i2+m
               KK(3,j)=i3
               LLL(i1,i2+m,i3)=0
            end if
            if ( LLL(i1,i2,i3+m)==mark_KK ) then
               j=j+1
               KK(1,j)=i1               
               KK(2,j)=i2
               KK(3,j)=i3+m
               LLL(i1,i2,i3+m)=0
            end if
         end do
      end do
      j=j-ML_irreducible
      if ( j/=MK_irreducible ) then
         write(*,*) "j,MK_irreducible=",j,MK_irreducible
         goto 900
      end if

!
! ---
!

      n1=pinfo_grid_2d(7,id_class(myrank,0))+1
      n2=n1+pinfo_grid_2d(8,id_class(myrank,0))-1

      allocate( Nstar(n1:n2) ) ; Nstar=0
      allocate( wkgrp(n1:n2) ) ; wkgrp=0.d0

      LLL(:,:,:)=0
      do i=n1,n2
         do isym=1,nsym
            LL_tmp(1:3)=matmul( rga(1:3,1:3,isym),LL(1:3,i) )
            LLL( LL_tmp(1),LL_tmp(2),LL_tmp(3) )=i
         end do
         Nstar(i)=count(LLL==i)
         wkgrp(i)=dble(Nstar(i))/dble(nsym)
      end do

      if (DISP_SWITCH) then
         write(*,*) "Nstar(min,max)  =",minval(Nstar),maxval(Nstar)
         write(*,*) "wkgrp(min,max)  =",minval(wkgrp),maxval(wkgrp)
         write(*,*) "1/wkgrp(min,max)=",1.d0/maxval(wkgrp),1.d0/minval(wkgrp)
      end if

      m=maxval(Nstar)
      allocate( LL_star(1:3,m,n1:n2) ) ; LL_star=0

      LLL(:,:,:)=0
      do i=n1,n2
         j=0
         do isym=1,nsym
            LL_tmp(1:3)=matmul( rga(1:3,1:3,isym),LL(1:3,i) )
            if ( LLL( LL_tmp(1),LL_tmp(2),LL_tmp(3) )==0 ) then
               j=j+1
               LL_star(1:3,j,i)=LL_tmp(1:3)
               LLL( LL_tmp(1),LL_tmp(2),LL_tmp(3) )=i
            end if
         end do
         if ( j/=Nstar(i) ) goto 900
      end do

!
! ---
!

      allocate( LLR(-Mx:Mx,-My:My,-Mz:Mz) ) ; LLR=0


      LLL(:,:,:)=0
      do i=1,ML_irreducible
         LLL( LL(1,i),LL(2,i),LL(3,i) )=i
         do isym=1,nsym
            LL_tmp(1:3)=matmul( rga(1:3,1:3,isym),LL(1:3,i) )
            LLL(LL_tmp(1),LL_tmp(2),LL_tmp(3))=i
            LLR(LL_tmp(1),LL_tmp(2),LL_tmp(3))=isym
         end do
      end do
      if ( DISP_SWITCH ) then
         write(*,*) count(LLR/=0),maxval(LLR)
      end if
      do i=ML_irreducible+1,ML_irreducible+MK_irreducible
         LLL( KK(1,i),KK(2,i),KK(3,i) )=i
         do isym=1,nsym
            LL_tmp(1:3)=matmul( rga(1:3,1:3,isym),KK(1:3,i) )
            LLL(LL_tmp(1),LL_tmp(2),LL_tmp(3))=i
            LLR(LL_tmp(1),LL_tmp(2),LL_tmp(3))=isym
         end do
      end do
      if ( DISP_SWITCH ) then
         write(*,*) count(LLL/=0),count(1<=LLL.and.LLL<=ML_irreducible),count(ML_irreducible<LLL)
         write(*,*) count(LLR/=0),maxval(LLR)
      end if

!
! --- integral weight ---
!

      allocate( weight(n1:n2) ) ; weight=0.d0

      do i=n1,n2
         weight(i)=count(LLL==i)
      end do

      sum0=sum( weight(n1:n2) )
      call mpi_allreduce(sum0,sum1,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
      sum1=sum1/dble(np_band*np_bzsm*np_spin)

      if (DISP_SWITCH) then
         write(*,*) "sum(weight)=",sum1
      end if

      sum1=sum1-count(0<LLL.and.LLL<=ML_irreducible)
      if ( abs(sum1)>1.d-10 ) then
         write(*,*) "weight is invalid"
         goto 900
      end if

      call watch(ctime1,etime1)
      if (DISP_SWITCH) then
         write(*,*) "TIME(INIT_PARALLEL_MOL_SYM)=",ctime1-ctime0,etime1-etime0
         write(*,*) "MEM(MB)=",mem,memax*B2MB
      end if

      return

 900  call stop_program1("init_parallel_mol_sym",1)

      END SUBROUTINE init_parallel_mol_sym
