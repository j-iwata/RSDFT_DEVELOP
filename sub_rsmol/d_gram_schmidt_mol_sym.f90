!--------1---------2---------3---------4---------5---------6---------7--
! Gram-Schmidt orthogonalization
! ( Takahashi, Block cyclic )
!
      SUBROUTINE Gram_Schmidt_Mol_Sym(iswitch_gs_baka)
      use global_variables
      implicit none

      integer :: iswitch_gs_baka
      integer :: k,s
      integer :: n1,n2,ML0,nn1,nn2,irank_b,ns,ne,ms,me
      integer :: mm1,mm2,n,ierr
      real(8) :: memax,mem,nop_tot,nop_max,nop_min,nop_0
      real(8) :: nop(13),flops(9),ct(13),et(13),nop1,nop2,check(9)
      real(8) :: ctime0,ctime1,etime0,etime1,ctt,ett
      real(8) :: ct0,ct1,ct2,et0,et1,et2
      real(8) :: tmp(41),tmp1(41),tmp2(41),tmp3(41)
      real(8) :: flops_0_max,flops_0_min,flops_1_max,flops_1_min
      real(8) :: ctt_max,ctt_min,ett_max,ett_min,flops_0,flops_1
      integer,allocatable :: ir(:),id(:)
      integer :: nbss,k1,ib,NBAND_BLK,ncycle,mrnk
      integer :: ir,id

      call watch(ctime0,etime0)
      ct=0.d0 ; et=0.d0 ; nop=0.d0 ; flops=0.d0

      n1    = idisp(myrank)+1
      n2    = idisp(myrank)+ircnt(myrank)
      ML0   = n2-n1+1
      mrnk  = id_class(myrank,4)
      mem   = 0.d0
      memax = 0.d0

      if ( np_band/=1 ) then
         write(*,*) "Band-parallel calc is temporary unavailable"
      end if

      allocate( ir(0:np_band-1),id(0:np_band-1) )
      mem=mem+bsintg*np_band*2 ; memax=max(mem,memax)
      ir(0:np_band-1)=ir_band(0:np_band-1)*ML0
      id(0:np_band-1)=id_band(0:np_band-1)*ML0
!
!      NBAND_BLK=NBLK
      NBAND_BLK=MB
      ncycle=(MB-1)/NBAND_BLK+1
!

      do s=MSP_0,MSP_1
      do k=MBZ_0,MBZ_1

         call watch1(ct(10),et(10),0)

         call mpi_allgatherv(unk(n1,MB_0,k,s),ir(mrnk),TYPE_MAIN,unk(n1,1,k,s),ir,id,TYPE_MAIN,comm_band,ierr)

         call watch1(ct(10),et(10),1)
         nop(10)=nop(10)+ir(mrnk)
         nop(11)=nop(11)+1.d0

         do ir=1,Nir
         do id=1,irdim(ir)

         do k1=1,ncycle

            irank_b=mod(k1-1,np_band)

!            ns=NBAND_BLK*(k1-1)+1
!            ne=min(ns+NBAND_BLK-1,MB)
            ns=irlabel2n(ir,id,1)
            ne=irlabel2n(ir,id,Nexcited)

            if ( id_class(myrank,4)==irank_b ) then

               call Gram_Schmidt_Mol_Sub(ns,ne,ns,ne,NBLK,k,s,ct,et,nop,mem,memax)

            end if

            call watch1(ct(12),et(12),0)

            n=ML0*(ne-ns+1)
            call mpi_bcast(unk(n1,ns,k,s),n,TYPE_MAIN,irank_b,comm_band,ierr)

            call watch1(ct(12),et(12),1)
            nop(12)=nop(12)+n
            nop(13)=nop(13)+1.d0

            if ( ns <= MB-NBAND_BLK ) then

               do ib=1,(ncycle-1)/np_band+1

                  nbss=(ib-1)*np_band+myrank_b+1

                  if ( nbss<=ncycle .and. nbss>= k1+1 ) then

                     ms=NBAND_BLK*(nbss-1)+1
                     me=min(ms+NBAND_BLK-1,MB)

                     if ( ms<=me ) then
                        call Gram_Schmidt_Mol_Sub(ms,me,ns,ne,NBLK,k,s,ct,et,nop,mem,memax)
                     end if

                  end if

               enddo

            end if

         end do ! k1

      end do ! k
      end do ! s

      mem=mem-bsintg*(size(ir)+size(id)) ; deallocate( id,ir )

      call watch(ctime1,etime1)
      ctt=ctime1-ctime0
      ett=etime1-etime0

!
! --- performance data ---
!
      nop_0=nop(1)+nop(4)+nop(7)

      nop1=( sop1*MB*(MB-1)*ML )*MBZ*nspin
      nop2=( sop2*MB*ML+sop3*MB*ML )*MBZ*nspin

!      flops_0 = nop_0/ctt/1.d6
!      flops_1 = (nop1+nop2)/ctt/1.d9

      ct(3)=ct(1)-ct(2)
      ct(6)=ct(4)-ct(5)
      ct(9)=ct(7)-ct(8)

!      flops(1) = nop(1)/ct(1)/1.d6
!      flops(2) = nop(1)/ct(3)/1.d6
!      flops(4) = nop(4)/ct(4)/1.d6
!      flops(5) = nop(4)/ct(6)/1.d6
!      flops(7) = nop(7)/ct(7)/1.d6
!      flops(8) = nop(7)/ct(9)/1.d6

      nop(2)=nop(2)*bdmain*B2MB
      nop(5)=nop(5)*bdmain*B2MB
      nop(8)=nop(8)*bdreal*B2MB
      nop(10)=nop(10)*bdmain*B2MB
      nop(12)=nop(12)*bdmain*B2MB

      tmp(1)=ctt
      tmp(2)=ett
      tmp(3)=nop_0
      tmp(4)=0.d0 !flops_0
      tmp(5)=0.d0 !flops_1
      tmp(6:13)=flops(1:8)
      tmp(14:26)=nop(1:13)
      tmp(27:39)=ct(1:13)
      tmp(40)=mem
      tmp(41)=memax
      call mpi_allreduce(tmp,tmp1,41,mpi_real8,mpi_sum,mpi_comm_world,ierr)
      call mpi_allreduce(tmp,tmp2,41,mpi_real8,mpi_max,mpi_comm_world,ierr)
      call mpi_allreduce(tmp,tmp3,41,mpi_real8,mpi_min,mpi_comm_world,ierr)

      ctt_max=tmp2(1)
      ctt_min=tmp3(1)
      ett_max=tmp2(2)
      ett_min=tmp3(2)
      nop_tot=tmp1(3)
      flops_0_max=tmp2(4)
      flops_0_min=tmp3(4)
      flops_1_max=tmp2(5)
      flops_1_min=tmp3(5)

      if (DISP_SWITCH) then
         write(*,'(1x,"TIME(D_GRAM_SCHMIDT_MOL_SYM)=",4f10.3)') ctt_max,ctt_min,ett_max,ett_min
!         write(*,'(1x," FLOPS   =",4f12.3)') flops_0_max,flops_0_min,flops_1_max,flops_1_min
         write(*,'(1x," # of OP =",2g20.8)') nop_tot,nop1+nop2
         write(*,'(1x," MM = ",6f10.3)') tmp2(29),tmp3(29),tmp2(6),tmp3(6),tmp2(7),tmp3(7)
         write(*,'(1x," -allreduce=    ",4f10.3)') tmp2(28),tmp3(28),tmp2(15),tmp2(16)
         write(*,'(1x," MV = ",6f10.3)') tmp2(32),tmp3(32),tmp2(9),tmp3(9),tmp2(10),tmp3(10)
         write(*,'(1x," -allreduce=    ",4f10.3)') tmp2(31),tmp3(31),tmp2(18),tmp2(19)
         write(*,'(1x," NOR= ",6f10.3)') tmp2(35),tmp3(35),tmp2(12),tmp3(12),tmp2(13),tmp3(13)
         write(*,'(1x," -allreduce=    ",4f10.3)') tmp2(34),tmp3(34),tmp2(21),tmp2(22)
         write(*,'(1x," %(band)allgatherv= ",4f10.3)') tmp2(36),tmp3(36),tmp2(23),tmp2(24)
         write(*,'(1x," %(band)bcast     = ",4f10.3)') tmp2(38),tmp3(38),tmp2(25),tmp2(26)
         write(*,'(1x," (mem)(MB)=",2f12.3,g24.15)') tmp2(41)*B2MB,tmp3(41)*B2MB,tmp1(40)
         write(*,*) "NBLK,NBLK1,MB,MBZ=",NBLK,NBLK1,MB,MBZ
      end if

      return

 900  call stop_program

      END SUBROUTINE Gram_Schmidt_Mol_Sym

!--------1---------2---------3---------4---------5---------6---------7--
! Gram-Schmidt orthogonalization
!
      RECURSIVE SUBROUTINE Gram_Schmidt_Mol_Sub(mm1,mm2,nn1,nn2,MBLK,k,s,ct,et,nop,mem,memax)
      use global_variables
      implicit none
      integer,intent(IN) :: mm1,mm2,nn1,nn2,MBLK,k,s
      real(8),intent(INOUT) :: ct(13),et(13),nop(13),mem,memax
      integer :: n,ns,ne,m,ms,me,m1,m2,mm,nn,MBLKH,ierr
      integer :: n1,n2,ML0,i
      real(8) :: c,d
      real(8) :: ct0,ct1,ct2,ct3
      real(8) :: et0,et1,et2,et3

      n1  = idisp(myrank)+1
      n2  = idisp(myrank)+ircnt(myrank)
      ML0 = n2-n1+1

      do ms=mm1,mm2,MBLK

         me=min(ms+MBLK-1,mm2)
         mm=me-ms+1

      do ns=nn1,nn2,MBLK

         ne=min(ns+MBLK-1,nn2)
         ne=min(ne,me-1)
         nn=ne-ns+1

         if ( nn<=0 ) cycle

         if ( ms>=ne+1 ) then

            allocate( utmp2(ns:ne,ms:me),vtmp2(ns:ne,ms:me) )
            mem=mem+bdmain*nn*mm*2 ; memax=max(mem,memax)

            call watch(ct0,et0)

            call dgemm(TRANSA,TRANSB,nn,mm,ML0,-zdV,unk(n1,ns,k,s),ML0,unk(n1,ms,k,s),ML0,zero,utmp2,nn)

            call watch(ct2,et2)

            call mpi_allreduce(utmp2,vtmp2,nn*mm,TYPE_MAIN,mpi_sum,comm_grid,ierr)

            call watch(ct3,et3)

            call dgemm(TRANSB,TRANSB,ML0,mm,nn,one,unk(n1,ns,k,s),ML0,vtmp2,nn,one,unk(n1,ms,k,s),ML0)

            call watch(ct1,et1)

            deallocate( vtmp2,utmp2 ) ; mem=mem-bdmain*nn*mm*2

            nop(1)=nop(1)+sop1*mm*nn*ML0+sop1*mm*nn*ML0
            nop(2)=nop(2)+nn*mm
            nop(3)=nop(3)+1.d0
            ct(1)=ct(1)+ct1-ct0
            et(1)=et(1)+et1-et0
            ct(2)=ct(2)+ct3-ct2
            et(2)=et(2)+et3-et2

            if ( ms==ne+1 ) then

               call watch(ct0,et0)

               d=sum(abs(unk(n1:n2,ms,k,s))**2)*dV

               call watch(ct2,et2)

               call mpi_allreduce(d,c,1,mpi_real8,mpi_sum,comm_grid,ierr)

               call watch(ct3,et3)

               c=1.d0/sqrt(c)
               unk(n1:n2,ms,k,s)=c*unk(n1:n2,ms,k,s)

               call watch(ct1,et1)

               nop(7)=nop(7)+sop2*ML0+sop3*ML0
               nop(8)=nop(8)+1.d0
               nop(9)=nop(9)+1.d0
               ct(7)=ct(7)+ct1-ct0
               et(7)=et(7)+et1-et0
               ct(8)=ct(8)+ct3-ct2
               et(8)=et(8)+et3-et2

            end if

         else if ( mm<=NBLK1 ) then

            allocate( utmp(NBLK1),vtmp(NBLK1) )
            mem=mem+bdmain*NBLK1*2 ; memax=max(memax,mem)

            do m=ms,me

               n=min(m-1,ne)

               if ( n-ns+1>0 ) then

                  call watch(ct0,et0)

                  call dgemv(TRANSA,ML0,n-ns+1,-zdV,unk(n1,ns,k,s),ML0,unk(n1,m,k,s),1,zero,utmp,1)

                  call watch(ct2,et2)

                  call mpi_allreduce(utmp,vtmp,n-ns+1,TYPE_MAIN,mpi_sum,comm_grid,ierr)

                  call watch(ct3,et3)

                  call dgemv(TRANSB,ML0,n-ns+1,one,unk(n1,ns,k,s),ML0,vtmp,1,one,unk(n1,m,k,s),1)

                  call watch(ct1,et1)

                  nop(4)=nop(4)+sop1*(n-ns+1)*ML0+sop1*(n-ns+1)*ML0
                  nop(5)=nop(5)+n-ns+1
                  nop(6)=nop(6)+1.d0
                  ct(4)=ct(4)+ct1-ct0
                  et(4)=et(4)+et1-et0
                  ct(5)=ct(5)+ct3-ct2
                  et(5)=et(5)+et3-et2

               end if

               if ( m==1 .or. (n==m-1 .and. m/=ns) ) then

                  call watch(ct0,et0)

                  d=sum(abs(unk(n1:n2,m,k,s))**2)*dV

                  call watch(ct2,et2)

                  call mpi_allreduce(d,c,1,mpi_real8,mpi_sum,comm_grid,ierr)

                  call watch(ct3,et3)

                  c=1.d0/sqrt(c)
                  unk(n1:n2,m,k,s)=c*unk(n1:n2,m,k,s)

                  call watch(ct1,et1)

                  nop(7)=nop(7)+sop2*ML0+sop3*ML0
                  nop(8)=nop(8)+1.d0
                  nop(9)=nop(9)+1.d0
                  ct(7)=ct(7)+ct1-ct0
                  et(7)=et(7)+et1-et0
                  ct(8)=ct(8)+ct3-ct2
                  et(8)=et(8)+et3-et2

               end if

            end do ! m

            deallocate( vtmp,utmp ) ; mem=mem-bdmain*NBLK1*2

         else

            MBLKH=max(MBLK/2,NBLK1)
            call Gram_Schmidt_Mol_Sub(ms,me,ns,ne,MBLKH,k,s,ct,et,nop,mem,memax)

         end if

      end do ! ns

      end do ! ms

      return

      END SUBROUTINE Gram_Schmidt_Mol_Sub
