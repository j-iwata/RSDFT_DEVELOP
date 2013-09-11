!--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3--
!
! Grid generation and numbering
!
      SUBROUTINE init_mesh_sym
      use global_variables
      implicit none

      integer :: i,j,ix,iy,iz,i1,i2,i3,i4,m,n,ierr,isym,LL_tmp(3),irank,n1,n2,mark_KK
      integer :: Mx,My,Mz,max_read,ie,ir,k,id,jd,itmp(2),jr
      integer,allocatable :: chi(:,:)
      real(8),parameter :: eps=1.d-10
      real(8) :: mem,memax,ctime0,ctime1,etime0,etime1,Rc2,r2,sum0,sum1,H
      real(8),allocatable :: tmp(:,:),tmp3(:,:,:)
      character(8) :: label

      if ( isymmetry/=1 ) return

      call watch(ctime0,etime0)
      if (DISP_SWITCH) then
         write(*,'(a60," init_mesh_sym")') repeat("-",60)
      end if

      mem   = 0.d0
      memax = 0.d0

      H  = H1
      Mx = ML1+Md
      My = ML2+Md
      Mz = ML3+Md

!
! --- Symmetry ---
!

      if ( SYStype/=1 ) then
         write(*,*) "SYStype/=1 is not available yet"
         call stop_program
      end if

!- Read -

      if ( myrank==0 ) then

         write(*,*) "Read symmetry data from ",file_symdat
         open(95,file=file_symdat,status='old')
         read(95,*) nsym
         write(*,*) "nsym=",nsym

!- allocate ------------------------------------
         allocate( Rsym(3,3,nsym) ) ; Rsym=0.d0
!-----------------------------------------------

         do isym=1,nsym
            read(95,*) Rsym(1:3,1,isym)
            read(95,*) Rsym(1:3,2,isym)
            read(95,*) Rsym(1:3,3,isym)
         end do
         close(95)

      end if

!- Broad cast -

      call mpi_bcast(nsym,1,mpi_integer,0,mpi_comm_world,ierr)
      if ( myrank/=0 ) then
         allocate( Rsym(3,3,nsym) ) ; Rsym=0.d0
      end if
      call mpi_bcast(Rsym,9*nsym,mpi_real8,0,mpi_comm_world,ierr)

!- Read -

      if ( myrank==0 ) then
         open(94,file='sym_irmat_Td',status='old')
         read(94,*)
         read(94,*)
         read(94,*)
         read(94,*) nsym,Nir,max_irdim
!- allocate ------------------------------------------------------
         allocate( irdim(Nir) ) ; irdim=0
         allocate( Rir(max_irdim,max_irdim,nsym,Nir) ) ; Rir=0.d0
!-----------------------------------------------------------------
         do i=1,Nir
            read(94,*) i1,n
            irdim(i)=n
            do isym=1,nsym
               do j=1,n
                  read(94,*) Rir(j,1:n,isym,i)
               end do
            end do
         end do
      end if

!- Broad cast -

      call mpi_bcast(Nir,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(max_irdim,1,mpi_integer,0,mpi_comm_world,ierr)
!- allocate -------------------------------------------------------
      if ( myrank/=0 ) then
         allocate( irdim(Nir) ) ; irdim=0
         allocate( Rir(max_irdim,max_irdim,nsym,Nir) ) ; Rir=0.d0
      end if
!------------------------------------------------------------------
      n=size(Rir)
      call mpi_bcast(irdim,Nir,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(Rir,n,mpi_real8,0,mpi_comm_world,ierr)

!- state labeling -

!- allocate ------------------------------
      allocate( iedim(Nir) ) ; iedim=0
      allocate( irindx(Nir) ) ; irindx=0
!-----------------------------------------

      if ( myrank==0 ) then
         max_read=max(2*MI,1000)
         do i=1,max_read
            read(unit,'(a8)') label ; write(*,*) "label=",label
            if ( label=='# SYMMOL' ) exit
            if ( i>=max_read ) stop "Label '# SYMMOL' was not found."
         end do
         read(unit,*) iedim(1:Nir)
         write(*,*) "# of excited states for each irreducible rep."
         write(*,*) "Nir=",Nir
         write(*,'(1x,"iedim(1:Nir)=",10i6)') iedim(1:Nir)
      end if

      call mpi_bcast(iedim,Nir,mpi_integer,0,mpi_comm_world,ierr,ierr)

      n=sum( irdim(1:Nir)*iedim(1:Nir) )
      if (DISP_SWITCH) then
         write(*,*) "Nir=",Nir
         write(*,*) "sum(irdim*iedim) (total # of states)=",n
         write(*,*) "MB=",MB
      end if

      if ( np_2d(4)>0 ) then
         if ( myrank==0 ) then
            read(unit,*) irindx(1:Nir)
            write(*,'(1x,"irindx(1:Nir)=",9i3)') irindx(1:Nir)
         end if
         call mpi_bcast(irindx,Nir,mpi_integer,0,mpi_comm_world,ierr)
!- sorting -
         allocate( tmp3(max_irdim,max_irdim,nsym) )
         do ir=1,Nir
            jr=irindx(ir)
            itmp(1)=irdim(ir)
            itmp(2)=iedim(ir)
            tmp3(:,:,:)=Rir(:,:,:,ir)
            irdim(ir)=irdim(jr)
            iedim(ir)=iedim(jr)
            Rir(:,:,:,ir)=Rir(:,:,:,jr)
            irdim(jr)=itmp(1)
            iedim(jr)=itmp(2)
            Rir(:,:,:,jr)=tmp3(:,:,:)
            do jr=1,Nir
               if ( irindx(jr)==ir ) then
                  irindx(jr)=irindx(ir)
                  exit
               end if
            end do
            irindx(ir)=ir
         end do
         deallocate( tmp3 )
      end if

!- character index table -

      allocate( chi(Nir,nsym) ) ; chi=0
      do isym=1,nsym
         do ir=1,Nir
            sum0=0.d0
            do i=1,irdim(ir)
               sum0=sum0+Rir(i,i,isym,ir)
            end do
            chi(ir,isym)=nint(sum0)
         end do
      end do
      if (DISP_SWITCH) then
         do ir=1,Nir
            write(*,'(1x,i3,3x,24i3)') ir,chi(ir,1:nsym)
         end do
      end if
      deallocate( chi )

!- allocate --------------------------------------------------------
      n=max(n,MB)
      N_excited=maxval( iedim(1:Nir) )
      allocate( irlabel(3,n) ) ; irlabel=0
      allocate( irlabel2n(Nir,max_irdim,N_excited) ) ; irlabel2n=0
!-------------------------------------------------------------------

      n=0
      do ir=1,Nir
         N_excited=iedim(ir)
         do id=1,irdim(ir)
            do ie=1,N_excited
               n=n+1
               if ( n>MB ) cycle
               irlabel(1,n)=ir
               irlabel(2,n)=id
               irlabel(3,n)=ie
               irlabel2n(ir,id,ie)=n
            end do
         end do
      end do

      if ( n/=MB ) then
         write(*,*) "MB should be ",n
         goto 900
      end if
      if (DISP_SWITCH) then
         write(*,*) "n,MB=",n,MB
         write(*,'(1x,a6,3a6)') "i","ir","id","ie"
         do i=1,n
            if ( all(irlabel(1:3,i)==0) ) cycle
            write(*,'(1x,i6,3i6)') i,irlabel(1:3,i)
         end do
      end if

! check(D-matrix orthonomalization relations)                                     

      if (DISP_SWITCH) write(*,*) "--- check D-matrix ortho-normalization ---"
      n=sum( irdim(1:Nir)**2 )
      allocate( tmp(nsym,n) ) ; tmp=0.d0
      mem=mem+bdreal*size(tmp) ; memax=max(mem,memax)
      k=0
      do ir=1,Nir
         do jd=1,irdim(ir)
         do id=1,irdim(ir)
            k=k+1
            do isym=1,nsym
               tmp(isym,k)=Rir(id,jd,isym,ir)
            end do
         end do
         end do
      end do
      if ( k/=n ) stop "k/=n"
      do j=1,n
      do i=1,n
         sum0=0.d0
         do isym=1,nsym
            sum0=sum0+tmp(isym,i)*tmp(isym,j)
         end do
         if ( abs(sum0)>1.d-10 ) then
            if ( DISP_SWITCH ) write(*,*) i,j,sum0
         end if
      end do
      end do
      mem=mem-bdreal*size(tmp) ; deallocate( tmp )

!- allocate -------------------------------------------
      allocate( rga(3,3,nsym) ) ; rga=0
      rga(1:3,1:3,1:nsym)=nint( Rsym(1:3,1:3,1:nsym) )
      allocate( LLL(-Mx:Mx,-My:My,-Mz:Mz) ) ; LLL=0
!------------------------------------------------------

      Rc2 = Rsize**2

      i1=0
      i3=0
      do iz=-Mz,Mz
      do iy=-My,My
      do ix=-Mx,Mx
         r2=H*H*(ix*ix+iy*iy+iz*iz)
         if ( r2>=Rc2+eps ) cycle
         i2=0
         do isym=1,nsym
            LL_tmp(1)=rga(1,1,isym)*ix+rga(1,2,isym)*iy+rga(1,3,isym)*iz
            LL_tmp(2)=rga(2,1,isym)*ix+rga(2,2,isym)*iy+rga(2,3,isym)*iz
            LL_tmp(3)=rga(3,1,isym)*ix+rga(3,2,isym)*iy+rga(3,3,isym)*iz
            j=LLL(LL_tmp(1),LL_tmp(2),LL_tmp(3))
            if ( j==0 ) then
               if ( i2==0 ) then
                  i1=i1+1
                  i2=1
               end if
               LLL(LL_tmp(1),LL_tmp(2),LL_tmp(3))=i1
            end if
         end do
         do m=-Md,Md
            if ( m==0 ) cycle
            r2=H*H*((ix+m)**2+iy*iy+iz*iz)
            if ( r2>=Rc2+eps ) then
               i4=0
               do isym=1,nsym
                  LL_tmp(1)=rga(1,1,isym)*(ix+m)+rga(1,2,isym)*iy+rga(1,3,isym)*iz
                  LL_tmp(2)=rga(2,1,isym)*(ix+m)+rga(2,2,isym)*iy+rga(2,3,isym)*iz
                  LL_tmp(3)=rga(3,1,isym)*(ix+m)+rga(3,2,isym)*iy+rga(3,3,isym)*iz
                  j=LLL( LL_tmp(1),LL_tmp(2),LL_tmp(3) )
                  if ( j==0 ) then
                     if ( i4==0 ) then
                        i3=i3-1
                        i4=1
                     end if
                     LLL( LL_tmp(1),LL_tmp(2),LL_tmp(3) )=i3
                  end if
               end do
            end if
            r2=H*H*(ix*ix+(iy+m)**2+iz*iz)
            if ( r2>=Rc2+eps ) then
               i4=0
               do isym=1,nsym
                  LL_tmp(1)=rga(1,1,isym)*ix+rga(1,2,isym)*(iy+m)+rga(1,3,isym)*iz
                  LL_tmp(2)=rga(2,1,isym)*ix+rga(2,2,isym)*(iy+m)+rga(2,3,isym)*iz
                  LL_tmp(3)=rga(3,1,isym)*ix+rga(3,2,isym)*(iy+m)+rga(3,3,isym)*iz
                  j=LLL( LL_tmp(1),LL_tmp(2),LL_tmp(3) )
                  if ( j==0 ) then
                     if ( i4==0 ) then
                        i3=i3-1
                        i4=1
                     end if
                     LLL( LL_tmp(1),LL_tmp(2),LL_tmp(3) )=i3
                  end if
               end do
            end if
            r2=H*H*(ix*ix+iy*iy+(iz+m)**2)
            if ( r2>=Rc2+eps ) then
               i4=0
               do isym=1,nsym
                  LL_tmp(1)=rga(1,1,isym)*ix+rga(1,2,isym)*iy+rga(1,3,isym)*(iz+m)
                  LL_tmp(2)=rga(2,1,isym)*ix+rga(2,2,isym)*iy+rga(2,3,isym)*(iz+m)
                  LL_tmp(3)=rga(3,1,isym)*ix+rga(3,2,isym)*iy+rga(3,3,isym)*(iz+m)
                  j=LLL( LL_tmp(1),LL_tmp(2),LL_tmp(3) )
                  if ( j==0 ) then
                     if ( i4==0 ) then
                        i3=i3-1
                        i4=1
                     end if
                     LLL( LL_tmp(1),LL_tmp(2),LL_tmp(3) )=i3
                  end if
               end do
            end if
         end do ! m
      end do
      end do
      end do
      ML_irreducible=i1
!      MK_irreducible=abs(i3)
      ML=count(LLL>0)
      MK=count(LLL<0)
      if ( DISP_SWITCH ) then
         write(*,*) "total # of inner and outer grid points",ML,MK
         write(*,*) "# of inner grid points reduced by symmetry",ML_irreducible,maxval(LLL)
!         write(*,*) "# of outer grid points reduced by symmetry",MK_irreducible,minval(LLL)
      end if

      call watch(ctime1,etime1)
      if (DISP_SWITCH) then
         write(*,*) "TIME(INIT_MESH_SYM)=",ctime1-ctime0,etime1-etime0
         write(*,*) "MEM(MB)=",mem,memax*B2MB
      end if

      return

 900  call stop_program1("init_mesh_sym",1)

      END SUBROUTINE init_mesh_sym
