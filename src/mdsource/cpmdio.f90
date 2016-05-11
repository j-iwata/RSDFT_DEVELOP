!--------1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3--

      SUBROUTINE write_data2(flag)
!      use global_variables
      use cpmd_variables
      implicit none
      logical,intent(IN) :: flag
      integer :: i,j,k,n,n1,n2,ierr,i1,i2,i3,j1,j2,j3,isym,ir,ie,id,i0,jd
      integer :: a1,a2,a3,b1,b2,b3,ML0,irank,ispin,lt(3)
      integer :: istatus(MPI_STATUS_SIZE,123)
      integer,allocatable :: irc(:),ids(:),irtmp(:),itmp3(:,:,:),idtmp(:)
      integer,save :: icount=0
      real(8) :: c,fs,ct0,ct1,et0,et1,mem,memax
      real(8),allocatable :: rtmp(:)
      logical :: flag_related

      if ( DISP_SWITCH ) then
         write(*,*)
         write(*,'(a60," write_data2")') repeat("-",60)
      end if

      call bwatch(ct0,et0)

      n1  = idisp(myrank)+1
      n2  = idisp(myrank)+ircnt(myrank)
      ML0 = ircnt(myrank)

      mem   = 0.d0
      memax = 0.d0

      if ( MBwr1<1 .or. MBwr2<MBwr1 .or. MB<MBwr1 ) MBwr1=1
      if ( MBwr2<1 .or. MBwr2<MBwr1 .or. MB<MBwr2 ) MBwr2=MB

      if(.not.allocated(LL2)) then
         allocate( LL2(3,ML) )
!         mem=mem+bsintg*ML*3 ; memax=max(mem,memax)
      endif
      LL2=0

      if ( SYStype/=0 .and. isymmetry==1 ) then

         if(.not.allocated(irtmp)) then
            allocate( irtmp(0:np_grid-1) )
!            mem=mem+bsintg*size(irtmp) ; memax=max(mem,memax)
         endif
         if(.not.allocated(idtmp)) then
            allocate( idtmp(0:np_grid-1) )
!            mem=mem+bsintg*size(idtmp) ; memax=max(mem,memax)
         endif
         if(.not.allocated(itmp3)) then
            allocate( itmp3(-ML1:ML1,-ML2:ML2,-ML3:ML3) )
!            mem=mem+bsintg*size(itmp3) ; memax=max(mem,memax)
         endif
         irtmp=0 ; idtmp=0 ; itmp3=0

         i0=0
         itmp3(:,:,:)=0
         do i=n1,n2
            do isym=1,nsym
               lt(1:3)=matmul( rga(1:3,1:3,isym),LL(1:3,i) )
               j=itmp3(lt(1),lt(2),lt(3))
               if ( j==0 ) then
                  i0=i0+1
                  itmp3(lt(1),lt(2),lt(3))=i0
               end if
            end do
         end do
         j=count(itmp3/=0)
         if ( j/=i0 ) goto 900
         irtmp(myrank_g)=i0
         call rsdft_allgather( irtmp(myrank_g:myrank_g),irtmp,comm_grid )
!         call mpi_allgather(irtmp(myrank_g),1,mpi_integer,irtmp,1,mpi_integer,comm_grid,ierr)
         do n=0,np_grid-1
            idtmp(n)=sum(irtmp(0:n))-irtmp(n)
         end do

         itmp3(:,:,:)=0
         i0=idtmp(myrank_g)
         do i=n1,n2
            do isym=1,nsym
               lt(1:3)=matmul( rga(1:3,1:3,isym),LL(1:3,i) )
               j=itmp3(lt(1),lt(2),lt(3))
               if ( j==0 ) then
                  i0=i0+1
                  itmp3(lt(1),lt(2),lt(3))=i0
                  LL2(1:3,i0)=lt(1:3)
               end if
            end do
         end do

         i1=idtmp(myrank_g)+1
         i2=idtmp(myrank_g)+irtmp(myrank_g)

         irtmp(0:np_grid-1)=irtmp(0:np_grid-1)*3
         idtmp(0:np_grid-1)=idtmp(0:np_grid-1)*3

         call rsdft_allgatherv( LL2(:,i1:i2),LL2,irtmp,idtmp,comm_grid )
!         call mpi_allgatherv(LL2(1,i1),irtmp(myrank_g),mpi_integer,LL2,irtmp,idtmp,mpi_integer,comm_grid,ierr)

         irtmp(0:np_grid-1)=irtmp(0:np_grid-1)/3
         idtmp(0:np_grid-1)=idtmp(0:np_grid-1)/3

         i0=0
         do i=1,ML
            if ( LL2(1,i)/=0 .or. LL2(2,i)/=0 .or. LL2(3,i)/=0 ) then
               i0=i0+1
            end if
         end do
         if ( i0+1/=ML ) goto 900

      else

         call Make_GridMap_1(LL2(1,n1),n1,n2)
         if(.not.allocated(irc)) then
            allocate( irc(0:np_grid-1) ) !; mem=mem+bsintg*np_grid ; memax=max(mem,memax)
         endif
         if(.not.allocated(ids)) then
            allocate( ids(0:np_grid-1) ) !; mem=mem+bsintg*np_grid ; memax=max(mem,memax)
         endif
         irc(0:np_grid-1)=3*ir_grid(0:np_grid-1)
         ids(0:np_grid-1)=3*id_grid(0:np_grid-1)
         call rsdft_algatherv( LL2(:,n1:n2),LL2,irc,idc,comm_grid )
!         call mpi_allgatherv(LL2(1,n1),irc(myrank_g),mpi_integer,LL2,irc,ids,mpi_integer,comm_grid,ierr)
!         mem=mem-bsintg*np_grid-1
         deallocate( ids,irc )

      end if

!
! --- Wave function ---
!

      if(.not.allocated(utmp)) then
         allocate( utmp(ML) )
!         mem=mem+bdmain*ML ; memax=max(mem,memax)
      endif
      utmp=zero

      if ( IO_CTRL==0 ) then
         if ( myrank==0 ) then
            open(1,file="cpmdwf.dat0",form="unformatted")
            write(1) ML,ML1,ML2,ML3
            write(1) MB,MBwr1,MBwr2
            write(1) LL2(:,:)
            write(1) occ(:,:,:)
         end if
      else
         open(1,file="cpmdwf.dat0",form="unformatted")
         write(1) ML,ML1,ML2,ML3
         write(1) MB,MBwr1,MBwr2
         write(1) LL2(:,:)
         write(1) occ(:,:,:)
      end if

      do ispin=1,nspin
      do k=1,MBZ
      do n=MBwr1,MBwr2

         do irank=0,nprocs-1
            i1=id_band(id_class(irank,4))+1
            j1=id_band(id_class(irank,4))+ir_band(id_class(irank,4))
            i2=id_bzsm(id_class(irank,5))+1
            j2=id_bzsm(id_class(irank,5))+ir_bzsm(id_class(irank,5))
            i3=id_spin(id_class(irank,6))+1
            j3=id_spin(id_class(irank,6))+ir_spin(id_class(irank,6))
            if ( id_grid(id_class(irank,0))==0 .and. i1<=n .and. n<=j1 .and. &
&                i2<=k .and. k<=j2 .and. i3<=ispin .and. ispin<=j3 ) exit
         end do

         if ( irank>=nprocs ) then
            write(*,*) "ERROR(write_data)",myrank
!            call stop_program
            stop
         end if

         flag_related = .false.
         if ( MBZ_0<=k .and. k<=MBZ_1 .and. MB_0 <=n .and. n<=MB_1 .and. MSP_0<=ispin .and. ispin<=MSP_1 ) then
            flag_related=.true.
         end if

         if ( IO_ctrl==0 ) then

            if ( SYStype/=0 .and. isymmetry==1 ) then

               if ( flag_related ) then
                  ir=irlabel(1,n)
                  id=irlabel(2,n)
                  ie=irlabel(3,n)
                  i0=idtmp(myrank_g)
                  itmp3(:,:,:)=0
                  do i=n1,n2
                     do isym=1,nsym
                        lt(1:3)=matmul( rga(1:3,1:3,isym),LL(1:3,i) )
                        j=itmp3(lt(1),lt(2),lt(3))
                        if ( j==0 ) then
                           c=0.d0
                           do jd=1,irdim(ir)
                              c=c+Rir(id,jd,isym,ir)*unk(i,irlabel2n(ir,jd,ie),k,ispin)
                           end do
                           i0=i0+1
                           utmp(i0)=c
                           itmp3(lt(1),lt(2),lt(3))=i0
                        end if
                     end do
                  end do
                  i0=idtmp(myrank_g)+1
                  call mpi_gatherv(utmp(i0),irtmp(myrank_g),TYPE_MAIN,utmp,irtmp,idtmp,TYPE_MAIN,0,comm_grid,ierr)
               end if
                  
            else

               if ( flag_related ) then
                  call mpi_gatherv(unk(n1,n,k,ispin),ML0,TYPE_MAIN,utmp,ir_grid,id_grid,TYPE_MAIN,0,comm_grid,ierr)
               end if

            end if

            call mpi_barrier(mpi_comm_world,ierr)

            if ( irank/=0 ) then
               if ( irank==myrank ) then
                  call mpi_send(utmp,ML,TYPE_MAIN,0,0,mpi_comm_world,ierr)
               end if
               if ( myrank==0 ) then
                  call mpi_recv(utmp,ML,TYPE_MAIN,irank,0,mpi_comm_world,istatus,ierr)
               end if
            end if
            if ( myrank==0 ) then
               write(1) utmp(:)
            end if

         else

            if ( flag_related ) then
               write(1) unk(n1:n2,n,k,ispin)
            end if

         end if

      end do ! n
      end do ! k
      end do ! ispin

      do ispin=1,nspin
      do k=1,MBZ
      do n=MBwr1,MBwr2

         do irank=0,nprocs-1
            i1=id_band(id_class(irank,4))+1
            j1=id_band(id_class(irank,4))+ir_band(id_class(irank,4))
            i2=id_bzsm(id_class(irank,5))+1
            j2=id_bzsm(id_class(irank,5))+ir_bzsm(id_class(irank,5))
            i3=id_spin(id_class(irank,6))+1
            j3=id_spin(id_class(irank,6))+ir_spin(id_class(irank,6))
            if ( id_grid(id_class(irank,0))==0 .and. i1<=n .and. n<=j1 .and. &
&                i2<=k .and. k<=j2 .and. i3<=ispin .and. ispin<=j3 ) exit
         end do

         if ( irank>=nprocs ) then
            write(*,*) "ERROR(write_data)",myrank
!            call stop_program
            stop
         end if

         flag_related = .false.
         if ( MBZ_0<=k .and. k<=MBZ_1 .and. MB_0 <=n .and. n<=MB_1 .and. MSP_0<=ispin .and. ispin<=MSP_1 ) then
            flag_related=.true.
         end if

         if ( IO_ctrl==0 ) then

            if ( SYStype/=0 .and. isymmetry==1 ) then

               if ( flag_related ) then
                  ir=irlabel(1,n)
                  id=irlabel(2,n)
                  ie=irlabel(3,n)
                  i0=idtmp(myrank_g)
                  itmp3(:,:,:)=0
                  do i=n1,n2
                     do isym=1,nsym
                        lt(1:3)=matmul( rga(1:3,1:3,isym),LL(1:3,i) )
                        j=itmp3(lt(1),lt(2),lt(3))
                        if ( j==0 ) then
                           c=0.d0
                           do jd=1,irdim(ir)
                              c=c+Rir(id,jd,isym,ir)*psi_v(i,irlabel2n(ir,jd,ie),k,ispin)
                           end do
                           i0=i0+1
                           utmp(i0)=c
                           itmp3(lt(1),lt(2),lt(3))=i0
                        end if
                     end do
                  end do
                  i0=idtmp(myrank_g)+1
                  call mpi_gatherv(utmp(i0),irtmp(myrank_g),TYPE_MAIN,utmp,irtmp,idtmp,TYPE_MAIN,0,comm_grid,ierr)
               end if
                  
            else

               if ( flag_related ) then
                  call mpi_gatherv(psi_v(n1,n,k,ispin),ML0,TYPE_MAIN,utmp,ir_grid,id_grid,TYPE_MAIN,0,comm_grid,ierr)
               end if

            end if

            call mpi_barrier(mpi_comm_world,ierr)

            if ( irank/=0 ) then
               if ( irank==myrank ) then
                  call mpi_send(utmp,ML,TYPE_MAIN,0,0,mpi_comm_world,ierr)
               end if
               if ( myrank==0 ) then
                  call mpi_recv(utmp,ML,TYPE_MAIN,irank,0,mpi_comm_world,istatus,ierr)
               end if
            end if
            if ( myrank==0 ) then
               write(1) utmp(:)
            end if

         else

            if ( flag_related ) then
               write(1) psi_v(n1:n2,n,k,ispin)
            end if

         end if

      end do ! n
      end do ! k
      end do ! ispin
!!

      if ( IO_ctrl==0 ) then

         if ( myrank==0 ) then
            close(1)
         end if

      else

         close(1)

      end if

!      fs = bsintg*( 7.d0 + 3.d0*ML ) + bdreal*( MB*MBZ*nspin ) + bdmain*( (MBwr2-MBwr1+1)*MBZ*nspin*ML )

      if (DISP_SWITCH) then
         write(*,*) "write to ",FILE_WF
         write(*,*) "MBwr1,MBwr2 =",MBwr1,MBwr2
!         write(*,*) " File size (MB) =",fs*B2MB
      end if

      deallocate( utmp ) !; mem=mem-bdmain*ML

      if ( SYStype/=0 .and. isymmetry==1 ) then
!         mem=mem-bsintg*size(itmp3) ;
         deallocate( itmp3 )
!         mem=mem-bsintg*size(idtmp) ;
         deallocate( idtmp )
!         mem=mem-bsintg*size(irtmp) ;
         deallocate( irtmp )
      end if

      deallocate( LL2 ) !; mem=mem-bsintg*ML*3

      icount=0

!      Max_mem_inst = max( Max_mem+memax,Max_mem_inst )

      call bwatch(ct1,et1)
      if (DISP_SWITCH) then
         write(*,*) "TIME(WRITE_DATA2)=",ct1-ct0,et1-et0
!         write(*,*) "MEM(MB)",memax*B2MB,Max_mem_inst*B2MB,mem
         write(*,*)
      end if

      return

 900  stop "stop@write_data2"

      end subroutine write_data2

!--------1---------2---------3---------4---------5---------6---------7--

      SUBROUTINE read_data2(flag)
!      use global_variables
      use cpmd_variables
      implicit none
      logical,intent(IN) :: flag
      integer :: k,n,i,j,ML_tmp,MB_tmp,MB1_tmp,MB2_tmp,n1,n2,ML0,irank
      integer :: ML1_tmp,ML2_tmp,ML3_tmp,ierr,i1,i2,i3,j1,j2,j3,ispin,i0
      integer :: itmp(7),a1,a2,a3,b1,b2,b3,istatus(MPI_STATUS_SIZE),lt(3),isym
      integer,allocatable :: LL_tmp(:,:),ir(:),id(:),itmp3(:,:,:),idtmp(:),irtmp(:)
      real(8) :: fs,mem,memax,ct0,et0,ct1,et1
      real(8),allocatable :: rtmp(:),rtmp3(:,:,:)
      logical :: flag_related

      call bwatch(ct0,et0)

      n1    = idisp(myrank)+1
      n2    = idisp(myrank)+ircnt(myrank)
      ML0   = ircnt(myrank)
      mem   = 0.d0
      memax = 0.d0

      if (DISP_SWITCH) then
         write(*,'(a60," read_data2")') repeat("-",60)
      end if

      if(.not.allocated(LL2)) then
         allocate( LL2(3,ML)    )
!         mem=mem+bsintg*ML*3 ; memax=max(memax,mem)
      endif
      if(.not.allocated(LL_tmp)) then
         allocate( LL_tmp(3,ML) )
!         mem=mem+bsintg*ML*3 ; memax=max(memax,mem)
      endif
      LL2=0 ; LL_tmp=0

      if ( SYStype/=0 .and. isymmetry==1 ) then

         if(.not.allocated(irtmp)) then
            allocate( irtmp(0:np_grid-1) )
!            mem=mem+bsintg*size(irtmp) ; memax=max(mem,memax)
         endif
         if(.not.allocated(idtmp)) then
            allocate( idtmp(0:np_grid-1) )
!            mem=mem+bsintg*size(idtmp) ; memax=max(mem,memax)
         endif
         if(.not.allocated(itmp3)) then
            allocate( itmp3(-ML1:ML1,-ML2:ML2,-ML3:ML3) )
!            mem=mem+bsintg*size(itmp3) ; memax=max(mem,memax)
         endif
         irtmp=0 ; idtmp=0 ; itmp3=0

         i0=0
         itmp3(:,:,:)=0
         do i=n1,n2
            do isym=1,nsym
               lt(1:3)=matmul( rga(1:3,1:3,isym),LL(1:3,i) )
               j=itmp3(lt(1),lt(2),lt(3))
               if ( j==0 ) then
                  i0=i0+1
                  itmp3(lt(1),lt(2),lt(3))=i0
               end if
            end do
         end do
         j=count(itmp3/=0)
         if ( j/=i0 ) goto 900
         irtmp(myrank_g)=i0
         call rsdft_allgather( irtmp(myrank_g:myrank_g),irtmp,comm_grid )
!         call mpi_allgather(irtmp(myrank_g),1,mpi_integer,irtmp,1,mpi_integer,comm_grid,ierr)
         do n=0,np_grid-1
            idtmp(n)=sum(irtmp(0:n))-irtmp(n)
         end do

         itmp3(:,:,:)=0
         i0=idtmp(myrank_g)
         do i=n1,n2
            do isym=1,nsym
               lt(1:3)=matmul( rga(1:3,1:3,isym),LL(1:3,i) )
               j=itmp3(lt(1),lt(2),lt(3))
               if ( j==0 ) then
                  i0=i0+1
                  itmp3(lt(1),lt(2),lt(3))=i0
                  LL2(1:3,i0)=lt(1:3)
               end if
            end do
         end do

         i1=idtmp(myrank_g)+1
         i2=idtmp(myrank_g)+irtmp(myrank_g)

         irtmp(0:np_grid-1)=irtmp(0:np_grid-1)*3
         idtmp(0:np_grid-1)=idtmp(0:np_grid-1)*3

         call rsdft_allgatherv( LL2(:,i1:i2),LL2,irtmp,idtmp,comm_grid )
!         call mpi_allgatherv(LL2(1,i1),irtmp(myrank_g),mpi_integer &
!                            ,LL2,irtmp,idtmp,mpi_integer,comm_grid,ierr)

         irtmp(0:np_grid-1)=irtmp(0:np_grid-1)/3
         idtmp(0:np_grid-1)=idtmp(0:np_grid-1)/3

         i0=0
         do i=1,ML
            if ( LL2(1,i)/=0 .or. LL2(2,i)/=0 .or. LL2(3,i)/=0 ) then
               i0=i0+1
            end if
         end do
         if ( i0+1/=ML ) goto 900

      else

         call Make_GridMap_1(LL2(1,n1),n1,n2)

         if(.not.allocated(ir)) then
            allocate( ir(0:np_grid-1) )
!            mem=mem+bsintg*np_grid ; memax=max(mem,memax)
         endif
         if(.not.allocated(id)) then
            allocate( id(0:np_grid-1) )
!            mem=mem+bsintg*np_grid ; memax=max(mem,memax)
         endif

         ir(0:np_grid-1)=3*ir_grid(0:np_grid-1)
         id(0:np_grid-1)=3*id_grid(0:np_grid-1)
         call rsdft_allgatherv( LL2(:,n1:n2),ir,id,LL2,comm_grid )
!         call mpi_allgatherv(LL2(1,n1),ir(myrank_g),mpi_integer &
!              ,LL2,ir,id,mpi_integer,comm_grid,ierr)

!         mem=mem-bsintg*np_grid-1 ;
         deallocate( id,ir )

      end if

!
! --- Read WF ---
!

      if ( IO_ctrl==0 ) then
         if ( myrank==0 ) then
            open(3,file="cpmdwf.dat",form='unformatted')
            read(3) ML_tmp,ML1_tmp,ML2_tmp,ML3_tmp
            read(3) MB_tmp,MB1_tmp,MB2_tmp
            itmp(1)=ML_tmp
            itmp(2)=ML1_tmp
            itmp(3)=ML2_tmp
            itmp(4)=ML3_tmp
            itmp(5)=MB_tmp
            itmp(6)=MB1_tmp
            itmp(7)=MB2_tmp
         end if
         call mpi_bcast(itmp,7,mpi_integer,0,mpi_comm_world,ierr)
         ML_tmp=itmp(1)
         ML1_tmp=itmp(2)
         ML2_tmp=itmp(3)
         ML3_tmp=itmp(4)
         MB_tmp=itmp(5)
         MB1_tmp=itmp(6)
         MB2_tmp=itmp(7)
      else
         open(3,file="cpmdwf.dat",form='unformatted')
         read(3) ML_tmp,ML1_tmp,ML2_tmp,ML3_tmp
         read(3) MB_tmp,MB1_tmp,MB2_tmp
      end if

      if (DISP_SWITCH) then
         write(*,*) "ML,ML_tmp  =",ML,ML_tmp
         write(*,*) "MB,MB_tmp=",MB,MB_tmp
      end if
      if ( ML_tmp/=ML ) then
         write(*,*) "ML,ML_tmp  =",ML,ML_tmp
         goto 900
      end if
      if ( ML1_tmp/=ML1 .or. ML2_tmp/=ML2 .or. ML3_tmp/=ML3 ) then
         if ( DISP_SWITCH ) then
            write(*,*) "ML1_tmp,ML2_tmp,ML3_tmp =",ML1_tmp,ML2_tmp,ML3_tmp
            write(*,*) "ML1,ML2,ML3 =",ML1,ML2,ML3
         end if
         goto 900
      end if

      if ( 1/=MB1_tmp .or. MB2_tmp/=MB ) then
         if (DISP_SWITCH) then
            write(*,*) "******** WARNING! ********"
            write(*,*) "MB1_tmp,MB2_tmp=",MB1_tmp,MB2_tmp
            write(*,*) "stop!"
         end if
         goto 900
      end if

      if ( IO_ctrl==0 ) then
         if ( myrank==0 ) then
            read(3) LL_tmp(:,:)
         end if
         call mpi_bcast(LL_tmp,3*ML,mpi_integer,0,mpi_comm_world,ierr)
      else
         read(3) LL_tmp(:,:)
      end if

      if(DISP_SWITCH)then
         write(*,*) "minval(LL_tmp),maxval(LL_tmp)"
         write(*,*) minval(LL_tmp(1,1:ML) ),maxval( LL_tmp(1,1:ML))
         write(*,*) minval(LL_tmp(2,1:ML) ),maxval( LL_tmp(2,1:ML))
         write(*,*) minval(LL_tmp(3,1:ML) ),maxval( LL_tmp(3,1:ML))
      endif

      i=sum(abs(LL_tmp(:,:)-LL2(:,:)))
      if ( i/=0 ) then
         if (DISP_SWITCH) then
            write(*,*) "LL and LL_tmp is different"
         end if
         if ( IO_ctrl /= 0 ) goto 900
      end if

      if ( IO_ctrl==0 ) then
         if ( myrank==0 ) then
            read(3) occ(:,:,:)
         end if
         call mpi_bcast(occ,MB*MBZ*nspin,mpi_real8,0,mpi_comm_world,ierr)
      else
         read(3) occ(:,:,:)
      end if

      if ( myrank==0 ) then
         write(*,*) "sum(occ)=",sum(occ)
      end if

      select case(SYStype)
      case default
         if(.not.allocated(utmp3)) then
            allocate( utmp3(0:ML1-1,0:ML2-1,0:ML3-1) )
!            mem=mem+bdmain*ML
         endif
         utmp3=zero
      case(1,2)
         if(.not.allocated(utmp3)) then
            allocate( utmp3(-ML1:ML1,-ML2:ML2,-ML3:ML3) )
!            mem=mem+bdmain*(2*ML1+1)*(2*ML2+1)*(2*ML2+1)
         endif
         utmp3=zero
      end select

      if(.not.allocated(utmp)) then
         allocate( utmp(ML) )
!         mem=mem+bdmain*ML ; memax=max(mem,memax)
      endif
      utmp=zero

      do ispin=1,nspin
      do k=1,MBZ
      do n=MB1_tmp,min(MB2_tmp,MB)

         do irank=0,nprocs-1
            i1=id_band(id_class(irank,4))+1
            j1=id_band(id_class(irank,4))+ir_band(id_class(irank,4))
            i2=id_bzsm(id_class(irank,5))+1
            j2=id_bzsm(id_class(irank,5))+ir_bzsm(id_class(irank,5))
            i3=id_spin(id_class(irank,6))+1
            j3=id_spin(id_class(irank,6))+ir_spin(id_class(irank,6))
            if ( id_grid(id_class(irank,0))==0 .and. i1<=n .and. n<=j1 .and. &
&                i2<=k .and. k<=j2 .and. i3<=ispin .and. ispin<=j3 ) exit
         end do
         if ( irank>=nprocs ) then
            write(*,*) "ERROR(read_data)",myrank
!            call stop_program
            stop
         end if
            
         flag_related = .false.
         if ( MBZ_0<=k .and. k<=MBZ_1 .and. MB_0 <=n .and. n<=MB_1 .and. &
              MSP_0<=ispin .and. ispin<=MSP_1 ) then
            flag_related=.true.
         end if

         if ( IO_ctrl==0 ) then

            if ( SYStype/=0 .and. isymmetry==1 ) then

               if ( myrank==0 ) then
                  read(3) utmp
                  do i=1,ML
                     i1=LL_tmp(1,i) ; i2=LL_tmp(2,i) ; i3=LL_tmp(3,i)
                     utmp3(i1,i2,i3)=utmp(i)
                  end do
                  do i=1,ML_irreducible
                     i1=LL(1,i) ; i2=LL(2,i) ; i3=LL(3,i)
                     utmp(i)=utmp3(i1,i2,i3)
                  end do
               end if

               call mpi_barrier(mpi_comm_world,ierr)

               if ( irank/=0 ) then
                  if ( myrank==0 ) then
                     call mpi_send(utmp,ML_irreducible,TYPE_MAIN,irank,0,mpi_comm_world,ierr)
                  end if
                  if ( myrank==irank ) then
                     call mpi_recv(utmp,ML_irreducible,TYPE_MAIN,0,0,mpi_comm_world,istatus,ierr)
                  end if
               end if

               call mpi_barrier(mpi_comm_world,ierr)

               if ( flag_related ) then
                  call mpi_scatterv(utmp,ir_grid,id_grid,TYPE_MAIN,unk(n1,n,k,ispin),ML0,TYPE_MAIN,0,comm_grid,ierr)
               end if

            else !( SYStype==0 )

               if ( myrank==0 ) then
                  read(3) utmp
                  do i=1,ML
                     i1=LL_tmp(1,i) ; i2=LL_tmp(2,i) ; i3=LL_tmp(3,i)
                     utmp3(i1,i2,i3)=utmp(i)
                  end do
                  do i=1,ML
                     i1=LL2(1,i) ; i2=LL2(2,i) ; i3=LL2(3,i)
                     utmp(i)=utmp3(i1,i2,i3)
                  end do
               end if

               call mpi_barrier(mpi_comm_world,ierr)

               if ( irank/=0 ) then
                  if ( myrank==0 ) then
                     call mpi_send(utmp,ML,TYPE_MAIN,irank,0,mpi_comm_world,ierr)
                  end if
                  if ( myrank==irank ) then
                     call mpi_recv(utmp,ML,TYPE_MAIN,0,0,mpi_comm_world,istatus,ierr)
                  end if
               end if

               call mpi_barrier(mpi_comm_world,ierr)

               if ( flag_related ) then
                  call mpi_scatterv(utmp,ir_grid,id_grid,TYPE_MAIN &
                       ,unk(n1,n,k,ispin),ML0,TYPE_MAIN,0,comm_grid,ierr)
               end if

            end if

         else

            if ( flag_related ) then
               read(3) unk(n1:n2,n,k,ispin)
            end if

         end if

      end do ! n
      end do ! k
      end do ! ispin

      do ispin=1,nspin
      do k=1,MBZ
      do n=MB1_tmp,min(MB2_tmp,MB)

         do irank=0,nprocs-1
            i1=id_band(id_class(irank,4))+1
            j1=id_band(id_class(irank,4))+ir_band(id_class(irank,4))
            i2=id_bzsm(id_class(irank,5))+1
            j2=id_bzsm(id_class(irank,5))+ir_bzsm(id_class(irank,5))
            i3=id_spin(id_class(irank,6))+1
            j3=id_spin(id_class(irank,6))+ir_spin(id_class(irank,6))
            if ( id_grid(id_class(irank,0))==0 .and. i1<=n .and. n<=j1 .and. &
&                i2<=k .and. k<=j2 .and. i3<=ispin .and. ispin<=j3 ) exit
         end do
         if ( irank>=nprocs ) then
            write(*,*) "ERROR(read_data)",myrank
!            call stop_program
            stop
         end if
            
         flag_related = .false.
         if ( MBZ_0<=k .and. k<=MBZ_1 .and. MB_0 <=n .and. n<=MB_1 .and. &
              MSP_0<=ispin .and. ispin<=MSP_1 ) then
            flag_related=.true.
         end if

         if ( IO_ctrl==0 ) then

            if ( SYStype/=0 .and. isymmetry==1 ) then

               if ( myrank==0 ) then
                  read(3) utmp
                  do i=1,ML
                     i1=LL_tmp(1,i) ; i2=LL_tmp(2,i) ; i3=LL_tmp(3,i)
                     utmp3(i1,i2,i3)=utmp(i)
                  end do
                  do i=1,ML_irreducible
                     i1=LL(1,i) ; i2=LL(2,i) ; i3=LL(3,i)
                     utmp(i)=utmp3(i1,i2,i3)
                  end do
               end if

               call mpi_barrier(mpi_comm_world,ierr)

               if ( irank/=0 ) then
                  if ( myrank==0 ) then
                     call mpi_send(utmp,ML_irreducible,TYPE_MAIN,irank,0,mpi_comm_world,ierr)
                  end if
                  if ( myrank==irank ) then
                     call mpi_recv(utmp,ML_irreducible,TYPE_MAIN,0,0,mpi_comm_world,istatus,ierr)
                  end if
               end if

               call mpi_barrier(mpi_comm_world,ierr)

               if ( flag_related ) then
                  call mpi_scatterv(utmp,ir_grid,id_grid,TYPE_MAIN,psi_v(n1,n,k,ispin),ML0,TYPE_MAIN,0,comm_grid,ierr)
               end if

            else !( SYStype==0 )

               if ( myrank==0 ) then
                  read(3) utmp
                  do i=1,ML
                     i1=LL_tmp(1,i) ; i2=LL_tmp(2,i) ; i3=LL_tmp(3,i)
                     utmp3(i1,i2,i3)=utmp(i)
                  end do
                  do i=1,ML
                     i1=LL2(1,i) ; i2=LL2(2,i) ; i3=LL2(3,i)
                     utmp(i)=utmp3(i1,i2,i3)
                  end do
               end if

               call mpi_barrier(mpi_comm_world,ierr)

               if ( irank/=0 ) then
                  if ( myrank==0 ) then
                     call mpi_send(utmp,ML,TYPE_MAIN,irank,0,mpi_comm_world,ierr)
                  end if
                  if ( myrank==irank ) then
                     call mpi_recv(utmp,ML,TYPE_MAIN,0,0,mpi_comm_world,istatus,ierr)
                  end if
               end if

               call mpi_barrier(mpi_comm_world,ierr)

               if ( flag_related ) then
                  call mpi_scatterv(utmp,ir_grid,id_grid,TYPE_MAIN &
                       ,psi_v(n1,n,k,ispin),ML0,TYPE_MAIN,0,comm_grid,ierr)
               end if

            end if

         else

            if ( flag_related ) then
               read(3) psi_v(n1:n2,n,k,ispin)
            end if

         end if

      end do ! n
      end do ! k
      end do ! ispin
!!

      if ( IO_ctrl==0 ) then
         if ( myrank==0 ) then
            close(3)
         end if
      else
         close(3)
      end if

!      mem=mem-bdmain*size(utmp)-bdmain*size(utmp3) ;
      deallocate( utmp,utmp3 )

      if (DISP_SWITCH) then
         write(*,*) "read from ",FILE_WF2
      end if

!      mem=mem-bsintg*size(LL_tmp) ;
      deallocate( LL_tmp )
!      mem=mem-bsintg*size(LL2)    ;
      deallocate( LL2 )

!      Max_mem_inst = max( Max_mem+memax, Max_mem_inst )

      if ( SYStype/=0 .and. isymmetry==1 ) then
!         mem=mem-bsintg*size(irtmp) ; memax=max(mem,memax)
!         mem=mem-bsintg*size(idtmp) ; memax=max(mem,memax)
!         mem=mem-bsintg*size(itmp3) ; memax=max(mem,memax)
         deallocate(irtmp,idtmp,itmp3)
      endif


      call bwatch(ct1,et1)
      if (DISP_SWITCH) then
         write(*,*) "TIME(READ_DATA2)=",ct1-ct0,et1-et0
!         write(*,*) "MEM(MB)",memax*B2MB,Max_mem_inst*B2MB,mem
      end if

      return

 900  stop "stop@read_data2"

      END SUBROUTINE read_data2
